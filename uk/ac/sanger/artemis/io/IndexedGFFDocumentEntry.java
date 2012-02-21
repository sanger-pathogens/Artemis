package uk.ac.sanger.artemis.io;

import java.io.File;
import java.io.IOException;
import java.io.Writer;
import java.util.Date;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.Vector;

import net.sf.samtools.util.BlockCompressedInputStream;

import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.components.genebuilder.GeneUtils;
import uk.ac.sanger.artemis.components.variant.TabixReader;
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
   private Hashtable<String, Contig> contigHash;
   private String name;
   
   private String contig;
   private boolean combinedReference = false;
   
   private Document document;
   private EntryInformation entryInfo;
   private EntryGroup entryGroup;
   private int featureCount = -1;
   
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
      StringBuilder hdr = new StringBuilder();
      contigHash = new Hashtable<String, Contig>(sequenceNames.length);
      int offset = 0;
      int cnt = 0;
      while( (ln = in.readLine()) != null && ln.startsWith("#"))
      {
        hdr.append(ln);
        
        if(ln.startsWith("##sequence-region "))
        {
          logger4j.debug(ln);
          final String parts[] = ln.split(" ");

          sequenceNames[cnt++] = parts[1];
          contigHash.put(parts[1], new Contig(parts[1], offset+1, offset+=Integer.parseInt(parts[3])));
        }
      }
      in.close();
    }
    catch (IOException e)
    {
      e.printStackTrace();
    }
  }
 
  
  private List<Contig> getChr(Range range)
  {
    final List<Contig> list = new Vector<Contig>();
    if(!combinedReference)
    {
      if(contig != null)
        list.add(contigHash.get(contig));
      else
        list.add(contigHash.get(sequenceNames[0]));
      return list;
    }
    
    Enumeration<String> seqs = contigHash.keys();
    while(seqs.hasMoreElements())
    {
      Contig contig = contigHash.get(seqs.nextElement());
      if( (range.getStart() >= contig.start && range.getStart() <= contig.end) ||
          (range.getEnd()   >= contig.start && range.getEnd()   <= contig.end))
      {
        //System.out.println(chr+" getChr() "+r.toString()+" RANGE "+range.toString());
        list.add(contig);
      }
    }
    return list;
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
    final List<Contig> contigs = getChr(range);
    for(Contig c: contigs)
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

    combineGeneFeatures(featuresInRange);
    return featuresInRange;
  }
  
  private void getFeaturesInRange(Contig c, Range range, FeatureVector features) throws NumberFormatException, IOException
  {
    int start = range.getStart()-c.start+1;
    if(start<1)
      start = 1;
    int end = range.getEnd(); 
    if(combinedReference)
      end-=c.start+1;
    String r = c.chr+":"+start+"-"+end;

    TabixReader.Iterator tabixIterator = reader.query(r);
    if(tabixIterator == null)
      return;
    int pos[] = iterate(c, range.getStart(), range.getEnd(), tabixIterator, features);

    if(pos[0] < start || pos[1] > end)
    {
      r = c.chr+":"+pos[0]+"-"+pos[1];
      tabixIterator = reader.query(r);
      if(tabixIterator == null)
        return;
      features.clear();

      iterate(c, pos[0], pos[1], tabixIterator, features);
    }
  }
  
  private int[] iterate(final Contig c, 
                        int min,
                        int max,
                        final TabixReader.Iterator tabixIterator, 
                        final FeatureVector features) throws NumberFormatException, ReadFormatException, IOException
  {
    String ln;

    while( (ln = tabixIterator.next()) != null )
    {
      final StringVector parts = StringVector.getStrings(ln, "\t", true);
      int sbeg = Integer.parseInt(((String)parts.elementAt(3)).trim());
      int send = Integer.parseInt(((String)parts.elementAt(4)).trim());

      if( (sbeg < min && send < min) || (sbeg > max && send > max) )
        continue;
      
      if(combinedReference)
      {
        sbeg += c.start - 1;
        send += c.start - 1;
        final StringBuffer newLn = new StringBuffer();
        for(int i=0; i<parts.size(); i++)
        {
          if(i==3)
            newLn.append(sbeg);
          else if(i==4)
            newLn.append(send);
          else
            newLn.append((String)parts.elementAt(i));
          newLn.append("\t");
        }
        ln = newLn.toString();
      }
      
      features.add(new GFFStreamFeature(ln));

      if( ((String)parts.elementAt(2)).equals("gene") )
      {
        if(sbeg < min)
          min = sbeg;
        if(send > max)
          max = send;
      }
    }
    return new int[]{min, max};
  }
  
  private boolean isTranscript(Key key)
  {
    if(key.getKeyString().equals("mRNA") || key.getKeyString().equals("transcript"))
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
    // TODO Auto-generated method stub
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
    // TODO Auto-generated method stub
    return null;
  }

  public int getFeatureCount()
  {
    if(contig == null)
      initContig();
 
    if(featureCount > -1)
      return featureCount;
    
    final String r = contig+":"+1+"-"+Integer.MAX_VALUE;
    TabixReader.Iterator tabixIterator = reader.query(r);
    if(tabixIterator == null)
      return 0;

    try
    {
      while( (tabixIterator.next()) != null )
        featureCount++;
    }
    catch(IOException ioe){}

    return featureCount;
  }

  public Feature add(Feature feature) throws EntryInformationException,
      ReadOnlyException
  {
    // TODO Auto-generated method stub
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

  public Feature getFeatureAtIndex(int i)
  {
    String r = contig+":"+1+"-"+Integer.MAX_VALUE;
    TabixReader.Iterator tabixIterator = reader.query(r);
    if(tabixIterator == null)
      return null;
    try
    {
      int cnt = 0;
      String ln;
      while( (ln = tabixIterator.next()) != null )
      {
        if(i == cnt++)
          return new GFFStreamFeature(ln.toString());
      }
    }
    catch(IOException ioe){}
    return null;
  }

  public int indexOf(Feature feature)
  {
    return 0;
  }

  public boolean contains(Feature feature)
  {
    // TODO Auto-generated method stub
    return false;
  }

  public FeatureEnumeration features()
  {
    System.out.println("HERE features()");
    return new FeatureEnumeration() {
      public boolean hasMoreFeatures() {
        return false;
      }

      public Feature nextFeature() {
        return null;
      }
    };
  }

  public FeatureVector getAllFeatures()
  {
    System.out.println("HERE getAllFeatures()");
    return null;
  }

  public Sequence getSequence()
  {
    return null;
  }


  public void dispose()
  {
    // TODO Auto-generated method stub 
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
    // TODO Auto-generated method stub
    
  }

  public void setDirtyFlag()
  {
    // TODO Auto-generated method stub
  }

  public Date getLastChangeTime()
  {
    // TODO Auto-generated method stub
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
  
  public static boolean isIndexed(File f)
  {
    File index = new File(f.getAbsolutePath() + ".tbi");
    return index.exists();
  }
  
  class Contig
  {
    protected String chr;
    protected int start;
    protected int end;
    protected int offset;

    Contig(String chr, int s, int e)
    {
      this.chr = chr;
      this.start = s;
      this.end = e;
      this.offset = s;
    }
    
    public boolean equals(Object obj)
    {
      Contig c = (Contig)obj;
      if(chr.equals(c.chr) && start == c.start && end == c.end)
        return true;
      return false;
    }
  }

  public void updateReference(String contig, boolean combinedReference)
  {
    // TODO Auto-generated method stub
    this.contig = contig;
    this.combinedReference = combinedReference;
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
      Enumeration<Contig> contigs = contigHash.elements();

      while(contigs.hasMoreElements())
      {
        Contig contig = contigs.nextElement();
        int clen = contig.end;
        if(clen > len)
          len = clen;
      }

      if(contigHash.size() > 1 && len == entry.getSequence().length())
      {
        updateReference(sequenceNames[0], true);
      }
      else
        updateReference(sequenceNames[0], false);
    }
  }
}