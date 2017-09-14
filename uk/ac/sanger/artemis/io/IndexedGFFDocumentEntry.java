/* IndexedGFFDocumentEntry.java
 *
 * created: 2012
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2012  Genome Research Limited
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or(at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 **/
package uk.ac.sanger.artemis.io;

import java.io.File;
import java.io.IOException;
import java.io.Writer;
import java.util.Comparator;
import java.util.Date;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Set;
import java.util.Vector;

import net.sf.samtools.util.BlockCompressedInputStream;

import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.components.FeatureDisplay;
import uk.ac.sanger.artemis.components.genebuilder.GeneUtils;
import uk.ac.sanger.artemis.components.variant.FeatureContigPredicate;
import uk.ac.sanger.artemis.components.variant.TabixReader;
import uk.ac.sanger.artemis.util.CacheHashMap;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.Document;
import uk.ac.sanger.artemis.util.FileDocument;
import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.util.ReadOnlyException;
import uk.ac.sanger.artemis.util.StringVector;


public class IndexedGFFDocumentEntry implements DocumentEntry
{
   private TabixReader reader;
   private String sequenceNames[];
   private LinkedHashMap<String, IndexContig> contigHash;
   private String name;
   
   private String contig;
   private boolean combinedReference = false;
   
   private Document document;
   private EntryInformation entryInfo;
   private EntryGroup entryGroup;
   private int featureCount = -1;
   
   private boolean isGTF = false;
   // cache used by getFeatureAtIndex() and indexOf()
   private CacheHashMap gffCache = new CacheHashMap(150,5);
   
   public static org.apache.log4j.Logger logger4j = 
       org.apache.log4j.Logger.getLogger(IndexedGFFDocumentEntry.class);
   
  /**
   *  Create a new IndexedGFFDocumentEntry object associated with the given
   *  Document.
   *  @param document This is the file that we will read from.  This is also
   *    used for saving the entry back to the file it came from and to give
   *    the new object a name.
   *  @param listener The object that will listen for ReadEvents.
   *  @exception IOException thrown if there is a problem reading the entry -
   *    most likely ReadFormatException.
   **/
  public IndexedGFFDocumentEntry(final Document document) 
  {
    this.document = document;
    entryInfo = new GFFEntryInformation();

    try
    {
      final File gffFile = ((FileDocument)getDocument()).getFile();
      setName(gffFile.getName());
      
      reader = new TabixReader(gffFile.getAbsolutePath());
      sequenceNames = reader.getSeqNames();

      final BlockCompressedInputStream in = 
           (BlockCompressedInputStream) getDocument().getInputStream();
      String ln;
      contigHash = new LinkedHashMap<String, IndexContig>(sequenceNames.length);
      int offset = 0;
      int cnt = 0;
      while( (ln = in.readLine()) != null && ln.startsWith("#"))
      {
        // ##sequence-region seqid start end
        if(ln.startsWith("##sequence-region "))
        {
          logger4j.debug(ln);
          final String parts[] = ln.split(" ");

          try
          {
            contigHash.put(parts[1], new IndexContig(parts[1], 1, Integer.parseInt(parts[3]), offset));
            offset+=Integer.parseInt(parts[3]);
            cnt++;
          }
          catch(Exception ae)
          {
            contigHash.clear();
            cnt = 0;
            break;
          }
        }
      }
      in.close();
      
      // no GFF header found
      if(cnt < 1)
      {
        logger4j.debug("No GFF header found for "+gffFile.getAbsolutePath());
        for(int i=0; i<sequenceNames.length; i++)
          contigHash.put(sequenceNames[i], new IndexContig(sequenceNames[i], 1, Integer.MAX_VALUE, 0));
      }
    }
    catch (IOException e)
    {
      e.printStackTrace();
    }
  }
  
  /**
   * Used when editing a subsequence and features.
   * @param constraint  base range to edit
   * @param entry       new entry to add features to
   */
  public void truncate(final Range constraint, final uk.ac.sanger.artemis.Entry entry)
  {
    final FeatureVector features = getFeaturesInRange(constraint);
    try
    {
      for(int i=0; i<features.size(); i++)
      {
        final GFFStreamFeature f = (GFFStreamFeature)features.get(i);
        f.setLocation(f.getLocation().truncate(constraint));
        f.setReadOnlyFeature(false);
        entry.getEMBLEntry().forcedAdd(f);
      }
    }
    catch (ReadOnlyException e)
    {
      e.printStackTrace();
    }
    catch (OutOfRangeException e)
    {
      e.printStackTrace();
    }
  }

  /**
   *  Return a vector containing the references of the Feature objects within
   *  the given range.
   *  @param range Return features that overlap this range - ie the start of
   *    the feature is less than or equal to the end of the range and the end
   *    of the feature is greater than or equal to the start of the range.
   *  @return The features of this feature table the are within
   *    the given range.  The returned object is a copy - changes will not
   *    effect the FeatureTable object itself.
   **/
  public FeatureVector getFeaturesInRange(Range range) 
  {
    if(contig == null)
      initContig();
    
    final FeatureVector featuresInRange = new FeatureVector();
    final List<IndexContig> contigs = getContigsInRange(range);

    for(IndexContig c: contigs)
    {
      try
      {
        getFeaturesInRange(c, range, featuresInRange);
      }
      catch(IOException ioe)
      {
        ioe.printStackTrace();
      }
    }

    if(featuresInRange.size() > 0 && GFFStreamFeature.isGTF((Feature)featuresInRange.get(0)))
    {
      isGTF = true;
      // GTF
      try
      {
        mergeGtfFeatures(featuresInRange, "CDS");
        mergeGtfFeatures(featuresInRange, "exon");
      }
      catch (ReadOnlyException e)
      {
        e.printStackTrace();
      }
      
    }
    else 
    {
      // GFF
      combineGeneFeatures(featuresInRange);
    }

    //combineGeneFeatures(featuresInRange);
    return featuresInRange;
  }
  
  private void getFeaturesInRange(IndexContig c, Range range, FeatureVector features) throws NumberFormatException, IOException
  {
    int start = getCoordInContigCoords(range.getStart(), c);
    int end = getCoordInContigCoords(range.getEnd(), c);
    
    if(isGTF)
    {
      // for GTF grab a larger range so Artemis knows about any
      // connecting exons outside the view
      start -= 500000;
      if(start < 1)
        start = 1;
      end += 500000;
      try
      {
        range = new Range(start, end);
      }
      catch (OutOfRangeException e){}
      if(end < start)
        return;
    }
    
    String r = c.chr+":"+start+"-"+end;
    TabixReader.Iterator tabixIterator = null;
    try
    {
      tabixIterator = reader.query(r);
    }
    catch(NullPointerException npe){}
     
    if(tabixIterator == null)
      return;
    
    FeatureVector featuresInRange = new FeatureVector();
    int pos[] = iterate(c, range.getStart(), range.getEnd(), tabixIterator, featuresInRange);

    if(pos[0] < range.getStart() || pos[1] > range.getEnd())
    {
      start = getCoordInContigCoords(pos[0], c);
      end = getCoordInContigCoords(pos[1], c);
      r = c.chr+":"+start+"-"+end;

      tabixIterator = reader.query(r);
      if(tabixIterator == null)
        return;
      featuresInRange.clear();

      iterate(c, pos[0], pos[1], tabixIterator, featuresInRange);
    }
    features.addAll(featuresInRange);
  }
  
  private int[] iterate(final IndexContig c, 
                        int min, int max,
                        final TabixReader.Iterator tabixIterator, 
                        final FeatureVector features) throws NumberFormatException, ReadFormatException, IOException
  {
    String ln;
    while( (ln = tabixIterator.next()) != null )
    {
      StringVector parts = StringVector.getStrings(ln, "\t", true);
      ln = getGffInArtemisCoordinates(ln, parts, c);
      parts = StringVector.getStrings(ln, "\t", true);
      
      int sbeg = Integer.parseInt(parts.elementAt(3).trim());
      int send = Integer.parseInt(parts.elementAt(4).trim());

      if( (sbeg < min && send < min) || (sbeg > max && send > max) )
        continue;

      GFFStreamFeature gff = new GFFStreamFeature(ln);
      gff.setReadOnlyFeature(true);
      features.add(gff);

      if( parts.elementAt(2).equals("gene") )
      {
        if(sbeg < min)
          min = sbeg;
        if(send > max)
          max = send;
      }
    }
    return new int[]{min, max};
  }
  
  /**
   * Return the sequences that lie within a given range
   * @param range
   * @return
   */
  private List<IndexContig> getContigsInRange(Range range)
  {
    final List<IndexContig> list = new Vector<IndexContig>();
    if(!combinedReference)
    {
      if(contig != null)
      {
        if(contigHash.get(contig) == null)
          System.err.println(contig+" not found in "+this.getName());
        else
          list.add(contigHash.get(contig));
      }
      else
        list.add(contigHash.get(sequenceNames[0]));
      return list;
    }
    
    for (String key : contigHash.keySet())
    {
      IndexContig contig = contigHash.get(key);
      
      if( (range.getStart() >= contig.getOffsetStart() && range.getStart() <= contig.getOffsetEnd()) ||
          (range.getEnd()   >= contig.getOffsetStart() && range.getEnd()   <= contig.getOffsetEnd()) ||
          (contig.getOffsetStart() >= range.getStart() && contig.getOffsetStart() <= range.getEnd()) ||
          (contig.getOffsetEnd() >= range.getStart()   && contig.getOffsetEnd() <= range.getEnd()) )
      {
        list.add(contig);
      }
    }
    return list;
  }
  
  /**
   * Get the list of contigs in the feature display.
   * @return
   */
  private List<IndexContig> getListOfContigs()
  {
    List<IndexContig> contigs = new Vector<IndexContig>();
    for (String key : contigHash.keySet())
    {
      IndexContig c = contigHash.get(key);
      if(combinedReference || c.chr.equals(contig))
        contigs.add(c);
    }
    return contigs;
  }
  
  /**
   * Get the features start coordinate.
   * @param gffParts
   * @param c
   * @return
   */
  private int getStartInArtemisCoords(final StringVector gffParts, final IndexContig c)
  {
    int sbeg = Integer.parseInt(((String)gffParts.elementAt(3)).trim());
    if(combinedReference)
      sbeg += c.getOffsetStart() - 1;
    return sbeg;
  }
  
  /**
   * Get the features start coordinate.
   * @param gffParts
   * @param c
   * @return
   */
  private int getEndInArtemisCoords(final StringVector gffParts, final IndexContig c)
  {
    int send = Integer.parseInt(((String)gffParts.elementAt(4)).trim());
    if(combinedReference)
      send += c.getOffsetStart() - 1;
    return send;
  }
  
  /**
   * Get coordinate on the contig.
   * @param start
   * @param c
   * @return
   */
  private int getCoordInContigCoords(int coord, final IndexContig c)
  {
    if(combinedReference)
      coord+=-c.getOffsetStart()+1;
    if(coord<1)
      coord = 1;
    return coord;
  }
  
  /**
   * Get the GFF line for this feature, adjusting the coordinates if contigs
   * are concatenated.
   * @param ln
   * @param gffParts
   * @param c
   * @return
   */
  private String getGffInArtemisCoordinates(String gffLine, final StringVector gffParts, final IndexContig c)
  {
    if(combinedReference)
    {
      int sbeg = Integer.parseInt(((String)gffParts.elementAt(3)).trim());
      int send = Integer.parseInt(((String)gffParts.elementAt(4)).trim());
      
      sbeg += c.getOffsetStart() - 1;
      send += c.getOffsetStart() - 1;
      final StringBuffer newLn = new StringBuffer();
      for(int i=0; i<gffParts.size(); i++)
      {
        if(i==3)
          newLn.append(sbeg);
        else if(i==4)
          newLn.append(send);
        else
          newLn.append((String)gffParts.elementAt(i));
        newLn.append("\t");
      }
      gffLine = newLn.toString();
    }
    return gffLine;
  }

  private boolean isTranscript(Key key)
  {
    if(key.getKeyString().indexOf("RNA") > -1 || 
       key.getKeyString().indexOf("transcript") > -1)
      return true;
    if(GeneUtils.isNonCodingTranscripts(key))
      return true;
    return false;
  }
  
  private void combineGeneFeatures(FeatureVector original_features)
  {
    Feature this_feature;
    HashMap<String, ChadoCanonicalGene> chado_gene = new HashMap<String, ChadoCanonicalGene>();
    try
    {
      // find the genes
      for(int i = 0 ; i < original_features.size() ; ++i) 
      {
        this_feature = original_features.featureAt(i);
        final String key = this_feature.getKey().getKeyString();
        if(this_feature instanceof GFFStreamFeature &&
           (GeneUtils.isHiddenFeature(key) ||
            GeneUtils.isObsolete((GFFStreamFeature)this_feature)))
          ((GFFStreamFeature)this_feature).setVisible(false);
        
        if(key.equals("gene") || key.equals("pseudogene"))
        {
          final Qualifier idQualifier = this_feature.getQualifierByName("ID");
          if(idQualifier != null)
          {
            String id = (String)this_feature.getQualifierByName("ID").getValues().get(0);
            ChadoCanonicalGene gene = new ChadoCanonicalGene();
            gene.setGene(this_feature);
            chado_gene.put(id, gene);
            ((GFFStreamFeature)this_feature).setChadoGene(gene);
          }
        }
      }

      // find the transcripts
      HashMap<String, ChadoCanonicalGene> transcripts_lookup = new HashMap<String, ChadoCanonicalGene>();
      for(int i = 0 ; i < original_features.size() ; ++i) 
      {
        this_feature = original_features.featureAt(i);
        // transcript 
        Qualifier parent_qualifier = this_feature.getQualifierByName("Parent");
        if(parent_qualifier == null || !isTranscript(this_feature.getKey()))
          continue;

        StringVector parents = parent_qualifier.getValues();
        for(int j=0; j<parents.size(); j++)
        {
          String parent = (String)parents.get(j);
          if(chado_gene.containsKey(parent))
          {
            ChadoCanonicalGene gene = (ChadoCanonicalGene)chado_gene.get(parent);

            // store the transcript ID with its ChadoCanonicalGene object
            try
            {
              transcripts_lookup.put((String)this_feature.getQualifierByName("ID").getValues().get(0),
                                   gene);
              ((GFFStreamFeature)this_feature).setChadoGene(gene);
              gene.addTranscript(this_feature);
            }
            catch(NullPointerException npe)
            {
              System.err.println(gene.getGeneUniqueName()+" "+this_feature.getKey().toString()+" "+this_feature.getLocation());
            }
            continue;
          }
        }
      }

      // find exons & protein
      String key;
      for(int i = 0 ; i < original_features.size() ; ++i) 
      {
        this_feature = original_features.featureAt(i);
        // exons
        key = this_feature.getKey().getKeyString();

        final Qualifier parent_qualifier  = this_feature.getQualifierByName("Parent");
        final Qualifier derives_qualifier = this_feature.getQualifierByName("Derives_from");
        if(parent_qualifier == null && derives_qualifier == null)
          continue;    
          
        final Qualifier featureRelationship = 
          this_feature.getQualifierByName("feature_relationship_rank");
        // compare this features parent_id's to transcript id's in the 
        // chado gene hash to decide if it is part of it
        final StringVector parent_id;
        
        if(parent_qualifier != null)
          parent_id = parent_qualifier.getValues();
        else
          parent_id = derives_qualifier.getValues();
        
        for(int j=0; j<parent_id.size(); j++)
        {
          final String parent = (String)parent_id.get(j);
         
          if(transcripts_lookup.containsKey(parent))
          {
            final ChadoCanonicalGene gene = (ChadoCanonicalGene)transcripts_lookup.get(parent);
            ((GFFStreamFeature)this_feature).setChadoGene(gene);
            
            if(parent_qualifier == null)
              gene.addProtein(parent, this_feature);
            else if(key.equals("three_prime_UTR"))
              gene.add3PrimeUtr(parent, this_feature);
            else if(key.equals("five_prime_UTR"))
              gene.add5PrimeUtr(parent, this_feature);
            else if(key.equals(DatabaseDocument.EXONMODEL) || key.equals("exon") || 
                    featureRelationship != null ||
                    key.equals("pseudogenic_exon"))
              gene.addSplicedFeatures(parent, this_feature);
            else
              gene.addOtherFeatures(parent, this_feature);
          }
        } 
      }
  
      // now join exons
      Iterator<String> enum_genes = chado_gene.keySet().iterator();
      while(enum_genes.hasNext())
      {
        ChadoCanonicalGene gene = chado_gene.get(enum_genes.next());
        combineChadoExons(gene, original_features);
      } 

    }
    catch(InvalidRelationException e)
    {
      e.printStackTrace();
    }
  }
  
  /**
   *  Combine the features (which are exons) and delete the orignals from this
   *  Entry.  The key of this hash will be the group name and the value is a
   *  FeatureVector containing the feature that are in that group.  Groups
   *  that have more than one member will be combined.
   **/
  private void combineChadoExons(ChadoCanonicalGene gene, FeatureVector features) 
  {
    final List<Feature> transcripts = gene.getTranscripts();
    gene.correctSpliceSiteAssignments();
    
    for(int i=0; i<transcripts.size(); i++)
    {
      GFFStreamFeature transcript = (GFFStreamFeature)transcripts.get(i);
      String transcript_id = (String)(transcript.getQualifierByName("ID").getValues().get(0));
      Set<String> splicedSiteTypes = gene.getSpliceTypes(transcript_id);
      if(splicedSiteTypes == null)
        continue;

      Iterator<String> it = splicedSiteTypes.iterator();
      Vector<Feature> new_set = new Vector<Feature>();
      while(it.hasNext())
      {
        String type = (String)it.next();
        List<Feature> splicedSites = gene.getSpliceSitesOfTranscript(transcript_id, type);
        if(splicedSites == null)
          continue;
      
        mergeFeatures(splicedSites, new_set, 
                      (String)(transcript.getQualifierByName("ID").getValues().get(0)));
        features.removeAll(splicedSites);
      }
      
      for(int j=0; j<new_set.size(); j++)
      {
        features.add(new_set.get(j));
        if(j == 0)
          gene.addSplicedFeatures(transcript_id, new_set.get(j), true );
        else
          gene.addSplicedFeatures(transcript_id, new_set.get(j));
      }
    }   
  }
  
  private void mergeFeatures(final List<Feature> gffFeatures, 
                             final List<Feature> new_set,
                             final String transcript_id)
  {
    final Hashtable<String, Range> id_range_store = new Hashtable<String, Range>();
    final RangeVector new_range_vector = new RangeVector();
    QualifierVector qualifier_vector = new QualifierVector();

    for (int j = 0; j < gffFeatures.size(); j++)
    {
      final GFFStreamFeature this_feature = (GFFStreamFeature) gffFeatures.get(j);
      final Location this_feature_location = this_feature.getLocation();

      if (this_feature_location.getRanges().size() > 1)
      {
        System.err.println("error - new location should have "
            + "exactly one range " + transcript_id + " "
            + this_feature.getKey().toString() + " "
            + this_feature_location.toStringShort());
        return;
      }

      final Range new_range = (Range) this_feature_location.getRanges().elementAt(0);

      Qualifier id_qualifier = this_feature.getQualifierByName("ID");
      if (id_qualifier != null)
      {
        String id = (String) (id_qualifier.getValues()).elementAt(0);
        id_range_store.put(id, new_range);
      }
      else
        logger4j.warn("NO ID FOUND FOR FEATURE AT: "
            + this_feature.getLocation().toString());

      if (this_feature_location.isComplement())
        new_range_vector.insertElementAt(new_range, 0);
      else
        new_range_vector.add(new_range);
      qualifier_vector.addAll(this_feature.getQualifiers());
    }

    final GFFStreamFeature first_old_feature = (GFFStreamFeature) gffFeatures.get(0);

    final Location new_location = new Location(new_range_vector,
        first_old_feature.getLocation().isComplement());

    qualifier_vector = mergeQualifiers(qualifier_vector, first_old_feature.getLocation().isComplement());

    final GFFStreamFeature new_feature = new GFFStreamFeature(
        first_old_feature.getKey(), new_location, qualifier_vector);

    if (first_old_feature.getChadoGene() != null)
      new_feature.setChadoGene(first_old_feature.getChadoGene());

    new_feature.setSegmentRangeStore(id_range_store);
    new_feature.setGffSource(first_old_feature.getGffSource());
    new_feature.setGffSeqName(first_old_feature.getGffSeqName());
    new_feature.setReadOnlyFeature(first_old_feature.isReadOnly());
    
    // set the ID
    String ID;
    try
    {
      ID = new_feature.getSegmentID(new_feature.getLocation().getRanges());
    }
    catch (NullPointerException npe)
    {
      if (new_feature.getQualifierByName("Parent") != null)
        ID = ((String) new_feature.getQualifierByName("Parent").getValues()
            .get(0))
            + ":"
            + new_feature.getKey().getKeyString()
            + ":"
            + new_feature.getLocation().getFirstBase();
      else
        ID = new_feature.getKey().getKeyString();
    }
    final Qualifier id_qualifier = new_feature.getQualifierByName("ID");
    id_qualifier.removeValue((String) (id_qualifier.getValues()).elementAt(0));
    id_qualifier.addValue(ID);

    // set visibility
    if (GeneUtils.isHiddenFeature(new_feature.getKey().getKeyString())
        || GeneUtils.isObsolete(new_feature))
      new_feature.setVisible(false);

    try
    {
      new_feature.setLocation(new_location);
      final Qualifier gene_qualifier = new_feature.getQualifierByName("gene");

      if (gene_qualifier != null
          && gene_qualifier.getValues().size() > 0
          && ((String) (gene_qualifier.getValues()).elementAt(0))
              .startsWith("Phat"))
      {
        // special case to handle incorrect output of the Phat gene
        // prediction tool
        new_feature.removeQualifierByName("codon_start");
      }
      
      new_set.add(new_feature);
    }
    catch (ReadOnlyException e)
    {
      throw new Error("internal error - unexpected exception: " + e);
    }
    catch (OutOfRangeException e)
    {
      throw new Error("internal error - unexpected exception: " + e);
    }
    catch (EntryInformationException e)
    {
      throw new Error("internal error - unexpected exception: " + e);
    }
  }

  private QualifierVector mergeQualifiers(final QualifierVector qualifier_vector,
                                          final boolean complement)
  {
    QualifierVector merge_qualifier_vector = new QualifierVector();
    boolean seen = false;

    for (int i = 0; i < qualifier_vector.size(); ++i)
    {
      Qualifier qual = (Qualifier) qualifier_vector.elementAt(i);
      if (qual.getName().equals("codon_start"))
      {
        if (!complement && !seen)
        {
          merge_qualifier_vector.addElement(qual);
          seen = true;
        }
        else if (complement)
          merge_qualifier_vector.setQualifier(qual);
      }
      else if (qual.getName().equals("Alias"))
      {
        final Qualifier id_qualifier = merge_qualifier_vector.getQualifierByName("Alias");
        if (id_qualifier == null)
          merge_qualifier_vector.addElement(qual);
        else
        {
          String id1 = (String) (id_qualifier.getValues()).elementAt(0);
          String id2 = (String) (qual.getValues()).elementAt(0);
          id_qualifier.removeValue(id1);
          id_qualifier.addValue(id1 + "," + id2);
        }
      }
      else if (!qual.getName().equals("ID")
          && !qual.getName().equals("feature_id"))
        merge_qualifier_vector.setQualifier(qual);
    }
    return merge_qualifier_vector;
  }

  /**
   * Merge function for GTF features
   * @param original_features
   * @param keyStr
   * @throws ReadOnlyException
   */
  private void mergeGtfFeatures(FeatureVector original_features, String keyStr) throws ReadOnlyException
  {
    Hashtable<String, Vector<GFFStreamFeature>> group = new Hashtable<String, Vector<GFFStreamFeature>>();
    for(int i=0; i<original_features.size(); i++)
    {
      GFFStreamFeature feature = (GFFStreamFeature)original_features.get(i);
      if(!feature.getKey().getKeyString().equals(keyStr))
        continue;
      String transcriptId = 
          ((String) feature.getQualifierByName("transcript_id").getValues().get(0)).replaceAll("'", "");
      if(group.containsKey(transcriptId))
        group.get(transcriptId).add(feature);
      else
      {
        Vector<GFFStreamFeature> this_group = new Vector<GFFStreamFeature>();
        this_group.add(feature);
        group.put(transcriptId, this_group);
      }
    }
    
    Enumeration<String> enumGroup = group.keys();
    while(enumGroup.hasMoreElements())
    {
      String transcriptId = enumGroup.nextElement();
      Vector<GFFStreamFeature> this_group = group.get(transcriptId);
      QualifierVector qualifier_vector = new QualifierVector();
      final RangeVector new_range_vector = new RangeVector();
      
      for(GFFStreamFeature this_feature: this_group)
      {
        qualifier_vector.addAll(this_feature.getQualifiers());
        
        final Range new_range = (Range) this_feature.getLocation().getRanges().elementAt(0);
        if(this_feature.getLocation().isComplement())
          new_range_vector.insertElementAt(this_feature.getLocation().getTotalRange(), 0);
        else
          new_range_vector.add(new_range);
        
        original_features.remove(this_feature);
      }
      final GFFStreamFeature old_feature = (GFFStreamFeature)this_group.get(0);

      final Location new_location = new Location(new_range_vector,
          old_feature.getLocation().isComplement());
      
      qualifier_vector = mergeQualifiers(qualifier_vector, new_location.isComplement());
      if(qualifier_vector.getQualifierByName("gene_id") != null)
        qualifier_vector.addQualifierValues(new Qualifier("ID",
            keyStr+":"+qualifier_vector.getQualifierByName("gene_id").getValues().get(0)));
      
      final GFFStreamFeature new_feature = new GFFStreamFeature(old_feature
          .getKey(), new_location, qualifier_vector);
      original_features.add(new_feature);
    }
  }

  public boolean hasUnsavedChanges()
  {
    return false;
  }

  public boolean isReadOnly()
  {
    return true;
  }

  public String getHeaderText()
  {
    return null;
  }

  public boolean setHeaderText(String new_header) throws IOException
  {
    return true;
  }

  public String getName()
  {
    return name;
  }

  public boolean setName(String name)
  {
    this.name = name;
    return true;
  }

  public Feature createFeature(Key key, Location location,
      QualifierVector qualifiers) throws EntryInformationException,
      ReadOnlyException, OutOfRangeException
  {
    // not for read only entry
    return null;
  }

  public int getFeatureCount()
  {
    if(contig == null)
      initContig();
 
    if(featureCount > -1)
      return featureCount;

    featureCount = 0;
    List<IndexContig> contigs = getListOfContigs();
    for(IndexContig c: contigs)
    {
      int nfeatures = 0;
      final String r = c.chr+":"+1+"-"+Integer.MAX_VALUE;

      TabixReader.Iterator tabixIterator = reader.query(r);
      if(tabixIterator == null)
        continue;

      try
      {
        while( tabixIterator.next() != null )
        {
          featureCount++;
          nfeatures++;
        }
        c.nfeatures = nfeatures;
      }
      catch(IOException ioe){}      
    }

    return featureCount;
  }

  public Feature add(Feature feature) throws EntryInformationException,
      ReadOnlyException
  {
    return null;
  }

  public Feature forcedAdd(Feature feature) throws ReadOnlyException
  {
    return null;
  }

  public boolean remove(Feature feature) throws ReadOnlyException
  {
    return false;
  }


  public Feature getFeatureAtIndex(int idx)
  {    
    Object cachedGFF = gffCache.get(idx);
    if(cachedGFF != null)
      return (GFFStreamFeature)cachedGFF;

    int cnt = 0;
    int start = 1;

    final List<IndexContig> contigs = getListOfContigs();

    for(IndexContig c: contigs)
    {
      int nfeatures = c.nfeatures;
      if(idx > cnt+nfeatures)
      {
        cnt+=nfeatures;
        continue;
      }
      String r = c.chr+":"+start+"-"+Integer.MAX_VALUE;

      TabixReader.Iterator tabixIterator = reader.query(r);
      if(tabixIterator == null)
        return null;
      try
      {
        String ln;
        while( (ln = tabixIterator.next()) != null )
        {
          if(idx == cnt++)
          {
            StringVector parts = StringVector.getStrings(ln, "\t", true);
            final GFFStreamFeature gff = new GFFStreamFeature(
                getGffInArtemisCoordinates(ln, parts, c));

            gffCache.put(idx, gff);
            
            // see if the following line is cached and if not cache the
            // next block of lines - this speeds up the generation of the
            // feature list
            if(gffCache.get(idx+1) == null)
            {
              cnt = 1;
              while(cnt < 32 && (ln = tabixIterator.next()) != null) 
              {
                parts = StringVector.getStrings(ln, "\t", true);
                gffCache.put(idx+cnt, new GFFStreamFeature(
                  getGffInArtemisCoordinates(ln, parts, c)));
                cnt++;
              }
            }
            return gff;
          }
        }
      }
      catch(IOException ioe){}
    }

    return null;
  }

  public int indexOf(Feature feature)
  {
    if(gffCache.containsValue(feature))
    {
      // retrieve from GFF cache
      for (Object key : gffCache.keySet())
      {
        Feature f = (Feature)gffCache.get(key);
        if(f.equals(feature))
          return (Integer)key;
      }
    }

    final List<IndexContig> contigs = getListOfContigs();
    int cnt = 0;
    
    final String keyStr = feature.getKey().getKeyString();
    final int sbeg1 = feature.getFirstBase();
    final int send1 = feature.getLastBase();

    for(IndexContig c: contigs)
    {
      if(combinedReference && sbeg1 > c.getOffsetEnd() && send1 > c.getOffsetStart())
      {
        cnt+=c.nfeatures;
        continue;
      }
      
      String r = c.chr+":"+1+"-"+Integer.MAX_VALUE;
      TabixReader.Iterator tabixIterator = reader.query(r);
      if(tabixIterator == null)
        continue;
      try
      {
        String ln;
        while( (ln = tabixIterator.next()) != null )
        { 
          final StringVector parts = StringVector.getStrings(ln, "\t", true);
          int sbeg2 = getStartInArtemisCoords(parts, c);
          int send2 = getEndInArtemisCoords(parts, c);

          if(sbeg1 == sbeg2 && parts.get(2).equals(keyStr))
          {
            if(send1 == send2 || feature.getLocation().getRanges().size() > 1)
            {
              if(gffCache.get(cnt) == null)
              {
                // add to cache
                final GFFStreamFeature gff = new GFFStreamFeature(
                  getGffInArtemisCoordinates(ln, parts, c));
                gffCache.put(cnt, gff);
              }
              return cnt;
            }
          }
          cnt++;
        }
      }
      catch(IOException ioe){}
    }
    return -1;
  }

  public boolean contains(Feature feature)
  {
    return (indexOf(feature)>-1);
  }

  public FeatureEnumeration features()
  {
    return new IndexGFFFeatureEnumeration();
  }

  public FeatureVector getAllFeatures()
  {
    return new FeatureVector(){
      private static final long serialVersionUID = 1L;

      public int size() 
      {
        return getFeatureCount();
      }
      
      public Feature featureAt(int index) 
      {
        return getFeatureAtIndex(index);
      }
    };
  }

  public Sequence getSequence()
  {
    return null;
  }


  public void dispose()
  {
  }

  public void save() throws IOException
  {
    save(getDocument());
  }

  public void save(Document document) throws IOException
  {
    try
    {
      final Writer out = document.getWriter();
      writeToStream(out);
      out.close();
    }
    catch(NullPointerException npe)
    {
      return;
    }
  }

  public void writeToStream(Writer writer) throws IOException
  {
  }

  public void setDirtyFlag()
  {
  }

  public Date getLastChangeTime()
  {
    return null;
  }

  public Document getDocument()
  {
    return document;
  }
  
  public EntryInformation getEntryInformation()
  {
    return entryInfo;
  }
  
  /**
   * Test if the tabix (.tbi) index is present.
   * @param f
   * @return
   */
  public static boolean isIndexed(File f)
  {
    File index = new File(f.getAbsolutePath() + ".tbi");
    return index.exists();
  }

  public void updateReference(String contig, boolean combinedReference)
  {
    this.contig = contig;
    this.combinedReference = combinedReference;
    featureCount = -1;
    
    gffCache.clear();
  }

  public void setEntryGroup(EntryGroup entryGroup)
  {
    this.entryGroup = entryGroup;
  }
  
  private void initContig()
  {
    Entry entry = entryGroup.getSequenceEntry().getEMBLEntry();
    if(entry.getSequence() instanceof IndexFastaStream)
      updateReference(((IndexFastaStream)entry.getSequence()).getContig(), false);
    else
    {
      int len = 0;
      int off = 0;
      for (String key : contigHash.keySet())
      {
        IndexContig contig = contigHash.get(key);
        int clen = contig.getOffsetEnd();
        if(clen > len)
          len = clen;
        off = contig.offset;
      }

      if(contigHash.size() > 1 && (len == entry.getSequence().length() || off == 0))
      {
        // check the order of the contigs/chromosomes
        checkOffset();
        updateReference(sequenceNames[0], true);
      }
      else
        updateReference(sequenceNames[0], false);
    }
  }
  
  /**
   * For concatenated sequences check the offset is the same as that in
   * the header of the GFF.
   */
  private void checkOffset()
  {
    if(entryGroup.getSequenceEntry().getFeatureCount() == sequenceNames.length)
    {
      final List<IndexContig> list = new Vector<IndexContig>();
      uk.ac.sanger.artemis.FeatureVector features = entryGroup.getSequenceEntry().getAllFeatures();
      
      for (String key : contigHash.keySet())
      {
        IndexContig contig = contigHash.get(key);
        list.add(contig);
        FeatureContigPredicate predicate = new FeatureContigPredicate(contig.chr);
        for(int j=0; j<features.size(); j++)
        {
          if(predicate.testPredicate(features.elementAt(j)))
          {
            // correct offset
            if(contig.getOffsetStart() != features.elementAt(j).getFirstBase() ||
               contig.offset == 0)
            {
              contig.offset = features.elementAt(j).getFirstBase()-1;
              
              // this needs to be set when the GFF header is missing
              contig.end = features.elementAt(j).getLastBase()-contig.offset;
            }
            
            break;
          }
        }
      }
      // sort the list by the contig offset
/*      Collections.sort(list, new ContigCompare());
      contigHash = new LinkedHashMap<String, IndexContig>(sequenceNames.length);
      for(IndexContig c: list)
        contigHash.put(c.chr, c);*/
    }  
  }
  
  /**
   *  Return true if the FeatureVector contains the given Feature.
   **/
  public static boolean contains(final uk.ac.sanger.artemis.Feature f, final uk.ac.sanger.artemis.FeatureVector fs)
  {
    final String id = f.getIDString();
    final String keyStr = f.getKey().toString();
    final String pId = FeatureDisplay.getParentQualifier(f);
    final String loc = f.getLocation().toStringShort();

    for(int i=0; i<fs.size(); i++)
      if(contains(fs.elementAt(i), id, keyStr, pId, loc))
        return true;
    return false;
  }
  
  private static boolean contains(final uk.ac.sanger.artemis.Feature f, 
                                  final String id, 
                                  final String keyStr, 
                                  final String pId,
                                  final String loc)
  {
    if(keyStr.equals(f.getKey().getKeyString()))
    {
      final String thisParentId = FeatureDisplay.getParentQualifier(f);
      if( f.getIDString().equals(id) &&
          (pId != null && thisParentId.equals(pId)) )
        return true;
      else if(id.indexOf("{")>-1)
      {
        int ind = f.getIDString().lastIndexOf(":");
        if( ind > -1 && id.startsWith(f.getIDString().substring(0, ind)) && 
           (pId == null || thisParentId.equals(pId)) )
          return true;
      }
      else if(f.getIDString().indexOf("{")>-1)
      {
        int ind = id.lastIndexOf(":");
        if( ind > -1 && f.getIDString().startsWith(id.substring(0, ind)) && 
           (pId == null || thisParentId.equals(pId)) )
          return true;
      }
      else if(id != null && id.equals(f.getIDString()))
      {
        if(loc.equals(f.getLocation().toStringShort()))
          return true;
      }
    }
    return false;
  }
  

  class IndexContig
  {
    private String chr;
    private int start;
    private int end;
    private int offset;
    private int nfeatures = 0;

    IndexContig(String chr, int s, int e, int off)
    {
      this.chr = chr;
      this.start = s;
      this.end = e;
      this.offset = off;
    }
    
    public boolean equals(Object obj)
    {
      IndexContig c = (IndexContig)obj;
      if(chr.equals(c.chr) && start == c.start && end == c.end)
        return true;
      return false;
    }
    
    private int getOffsetStart()
    {
      return start+offset;
    }
    
    private int getOffsetEnd()
    {
      return end+offset;
    }
  }

  class ContigCompare implements Comparator<IndexContig>
  {
    public int compare(IndexContig c1, IndexContig c2)
    {
      if(c1.offset < c2.offset)
        return -1;
      else if(c2.offset > c1.offset)
        return 1;
      return 0;
    }
  }
  
  class IndexGFFFeatureEnumeration implements FeatureEnumeration
  { 
    private FeatureVector features;
    private int idx = 0;
    private int contigIdx = 0;
    private List<IndexContig> contigs;
    
    public boolean hasMoreFeatures()
    {
      if(features == null)
      {
        if(entryGroup == null)
          return false;

        if(contigs == null)
        {
          try
          {
            contigs =  getContigsInRange(
                new Range(1, entryGroup.getSequenceLength()));
          }
          catch (OutOfRangeException e){}
        }
        getFeaturesInContig();
      }

      if(idx < features.size())
        return true;
      else
      {
        idx = 0;
        contigIdx++;
        if(contigIdx < contigs.size())
        {
          getFeaturesInContig();
          if(idx < features.size())
            return true;
        }
      }
      return false;
    }
    
    public Feature nextFeature() throws NoSuchElementException
    {
      idx++;
      return features.elementAt(idx-1);
    }
    
    private void getFeaturesInContig()
    {
      try
      {
        features = getFeaturesInRange(
            new Range(contigs.get(contigIdx).getOffsetStart(), 
                      contigs.get(contigIdx).getOffsetEnd()));
      }
      catch (OutOfRangeException e){}
    }
  }
}
