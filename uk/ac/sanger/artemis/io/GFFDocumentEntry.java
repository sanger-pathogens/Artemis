/* GFFDocumentEntry.java
 *
 * created: Tue Sep 14 1999
 *
 * This file is part of Artemis
 *
 * Copyright (C) 1999-2005  Genome Research Limited
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/GFFDocumentEntry.java,v 1.69 2009-09-03 13:47:31 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.chado.FeatureLocLazyQualifierValue;
import uk.ac.sanger.artemis.components.Splash;
import uk.ac.sanger.artemis.components.filetree.LocalAndRemoteFileManager;
import uk.ac.sanger.artemis.components.genebuilder.GeneUtils;
import uk.ac.sanger.artemis.util.*;

import java.io.IOException;
import java.util.Hashtable;
import java.util.Enumeration;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.Vector;
import java.sql.Timestamp;

import org.gmod.schema.cv.CvTerm;
import org.gmod.schema.sequence.FeatureLoc;

/**
 *  A DocumentEntry that can read an GFF entry from a Document.
 *
 *  @author Kim Rutherford
 *  @version $Id: GFFDocumentEntry.java,v 1.69 2009-09-03 13:47:31 tjc Exp $
 **/

public class GFFDocumentEntry extends SimpleDocumentEntry
    implements DocumentEntry 
{
  private boolean finished_constructor = false;
  private boolean isReadOnly = false;
  
  /**
   *  Create a new GFFDocumentEntry object associated with the given
   *  Document.
   *  @param document This is the file that we will read from.  This is also
   *    used for saving the entry back to the file it came from and to give
   *    the new object a name.
   *  @param listener The object that will listen for ReadEvents.
   *  @exception IOException thrown if there is a problem reading the entry -
   *    most likely ReadFormatException.
   **/
  GFFDocumentEntry(final Document document, final ReadListener listener)
      throws IOException, EntryInformationException 
  {
    super(new GFFEntryInformation(), document, listener);
    super.in_constructor = true;
    // join the separate exons into one feature (if appropriate)
    final FeatureVector original_features = getAllFeatures();
    if(original_features.size() > 0 && GFFStreamFeature.isGTF((Feature)original_features.get(0)))
    {
      // GTF
      mergeGtfFeatures(original_features, "CDS");
      mergeGtfFeatures(original_features, "exon");
    }
    else 
    {
      // GFF
      combineGeneFeatures(original_features);
    }
    super.in_constructor = false;
    finished_constructor = true;
  }

  /**
   *  Create a new GFFDocumentEntry that will be a copy of the given Entry and
   *  has no Document associated with it.  The new GFFDocumentEntry cannot be
   *  saved to a file with save() unless save(Document) has been called
   *  first.  Some qualifier and location information will be lost.
   *  @param force If true then invalid qualifiers and any features with
   *    invalid keys in the new Entry will be quietly thrown away.  "Invalid"
   *    means that the key/qualifier is not allowed to occur in an Entry of
   *    this type (probably determined by the EntryInformation object of this
   *    Entry).  If false an EntryInformationException will be thrown for
   *    invalid keys or qualifiers.
   **/
  public GFFDocumentEntry(final Entry new_entry, final boolean force)
      throws EntryInformationException 
  {
    super(new GFFEntryInformation(), new_entry, force);
    finished_constructor = true;
  }

  /**
   *  Create a new empty GFFDocumentEntry object that has no Document
   *  associated with it.  The new GFFDocumentEntry cannot be saved to a
   *  file with save() unless save(Document) has been called first.  The
   *  save(Document) method will assign a Document.
   **/
  public GFFDocumentEntry(final EntryInformation entry_information) 
  {
    super(new GFFEntryInformation());
    finished_constructor = true;
  }

  /**
   *  Returns true if and only if this entry is read only.  For now this
   *  always returns true - GFFDocumentEntry objects can't be changed.
   **/
  /**
   *  Returns true if and only if this entry is read only.  For now this
   *  always returns true - BlastDocumentEntry objects can't be changed.
   **/
  public boolean isReadOnly() 
  {
    return isReadOnly;
  }
  
  public void setReadOnly(final boolean isReadOnly)
  {
    this.isReadOnly = isReadOnly;
  }

  /**
   *  If the given feature can be added directly to this Entry, then return
   *  it, otherwise create and return a new feature of the appropriate type.
   *  @param copy if true then always new a new copy of the Feature.
   **/
  protected Object makeNativeFeature(final Feature feature,
                                     final boolean copy) 
  {
    if(!copy && feature instanceof GFFStreamFeature) 
      return (GFFStreamFeature)feature;
    else 
    {
      if(PublicDBDocumentEntry.IGNORE_OBSOLETE_FEATURES)
      {
        Qualifier isObsoleteQualifier = 
          feature.getQualifiers().getQualifierByName("isObsolete");
        if(isObsoleteQualifier != null)
        {
          String value = (String)isObsoleteQualifier.getValues().get(0);
          if(Boolean.parseBoolean(value))
            return null;
        }
      }
      return new GFFStreamFeature(feature);
    }
  }

  /**
   *  If the given Sequence can be added directly to this Entry, then return a
   *  copy of it, otherwise create and return a new feature of the appropriate
   *  type for this Entry.
   **/
  protected StreamSequence makeNativeSequence(final Sequence sequence)
  {
    return new FastaStreamSequence(sequence);
  }

  private void combineGeneFeatures(FeatureVector original_features)
  {
    Feature this_feature;
    Hashtable chado_gene = new Hashtable();
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
      Hashtable transcripts_lookup = new Hashtable();
      for(int i = 0 ; i < original_features.size() ; ++i) 
      {
        this_feature = original_features.featureAt(i);
        
        // transcript 
        Qualifier parent_qualifier = this_feature.getQualifierByName("Parent");
        
        if(parent_qualifier == null)
          continue;

        StringVector parents = parent_qualifier.getValues();
        for(int j=0; j<parents.size(); j++)
        {
          String parent = (String)parents.get(j);
          
          if(chado_gene.containsKey(parent))
          {
            // store transcript
            ChadoCanonicalGene gene = (ChadoCanonicalGene)chado_gene.get(parent);
            gene.addTranscript(this_feature);
            ((GFFStreamFeature)this_feature).setChadoGene(gene);
            
            // store the transcript ID with its ChadoCanonicalGene object
            
            if(this_feature.getQualifierByName("ID") != null)
              transcripts_lookup.put((String)this_feature.getQualifierByName("ID").getValues().get(0),
                                   gene);
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
            {
              //((GFFStreamFeature)this_feature).setVisible(false);
              gene.addProtein(parent, this_feature);
            }
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
  
      //
      // 
      if(getDocument() instanceof DatabaseDocument)
      {
        DatabaseDocument doc = (DatabaseDocument)getDocument();
        loadFeatureLocLazyData(original_features); 
        if(doc.isLazyFeatureLoad())
        {
          // using lazy loading - add the lazy chado feature to GFFStreamFeature
          final Hashtable idFeatureStore = doc.getIdFeatureStore();
          for(int i = 0 ; i < original_features.size() ; ++i) 
          {
            this_feature = original_features.featureAt(i);
            String featureId = (String) this_feature.getQualifierByName("feature_id").getValues().get(0);
            org.gmod.schema.sequence.Feature chadoLazyFeature = 
              (org.gmod.schema.sequence.Feature)idFeatureStore.get(featureId);
            ((GFFStreamFeature)this_feature).setChadoLazyFeature(chadoLazyFeature);
          }
          idFeatureStore.clear();
        }
      }

      // now join exons
      //combineFeatures();
      Enumeration enum_genes = chado_gene.elements();
      while(enum_genes.hasMoreElements())
      {
        ChadoCanonicalGene gene = (ChadoCanonicalGene)enum_genes.nextElement();
        combineChadoExons(gene);

        // inferring CDS and UTRs
        if(DatabaseDocument.CHADO_INFER_CDS)
        {
          final Vector transcripts = (Vector)gene.getTranscripts();
          gene.correctSpliceSiteAssignments();
          
          for(int i=0; i<transcripts.size(); i++)
          {
            GFFStreamFeature transcript = (GFFStreamFeature)transcripts.get(i);
            String transcript_id = null;
            transcript_id = GeneUtils.getUniqueName(transcript);
            
            List exons = gene.getSpliceSitesOfTranscript(transcript_id, "exon");
            if(exons == null)
              continue;
            
            Iterator it = exons.iterator();
            while(it.hasNext())
            {
              final GFFStreamFeature exonFeature = (GFFStreamFeature)it.next();

              QualifierVector qualifiers = new QualifierVector();
              qualifiers.add(new Qualifier("ID", transcript_id+":CDS"));
              qualifiers.add(new Qualifier("Parent", transcript_id));
              
              DatabaseInferredFeature cdsFeature = new DatabaseInferredFeature(
                  Key.CDS, exonFeature.getLocation(), qualifiers, gene);
              
              try
              {
                gene.addSplicedFeatures(transcript_id, cdsFeature);
                forcedAdd(cdsFeature);
              }
              catch (ReadOnlyException e)
              {
                // TODO Auto-generated catch block
                e.printStackTrace();
              }
              
            }
          }
        }
      } 

    }
    catch(InvalidRelationException e)
    {
      e.printStackTrace();
    }
  }
  
  /**
   * Get 'similarity', polypeptide_domain qualifiers
   * @param fv
   * @throws InvalidRelationException
   */
  private void loadFeatureLocLazyData(final FeatureVector fv)
          throws InvalidRelationException
  {
    final DatabaseDocument doc = (DatabaseDocument)getDocument();
    List matches;
 
    if(fv.size() < 30 && fv.size() > 0)  // if just a few features to look up e.g. for gene editor
    {
      List featureIds = new Vector(fv.size());
      for(int i=0;i<fv.size(); i++)
      {
        Qualifier featureIdQualifier = fv.featureAt(i).getQualifierByName("feature_id");
        featureIds.add( (String)featureIdQualifier.getValues().get(0) );
      }
      matches = doc.getSimilarityMatches(featureIds);
    }
    else
      matches = doc.getSimilarityMatches(null);
    
    if(matches == null || matches.size() < 1)
      return;
    final Hashtable temp_lookup_hash = new Hashtable(matches.size()/2);
    String f_id;
    for(int i=0; i<fv.size(); i++)
    {
      Feature f = (Feature)fv.elementAt(i);
      Qualifier qualifier = ((Feature)f).getQualifierByName("feature_id");
      
      if(qualifier != null)
      {
        f_id = (String)qualifier.getValues().get(0);
        temp_lookup_hash.put(f_id, f);
      }
    }
    
    // bulk load featureLocs and create lookup hash 
    final Hashtable hashFeatureLocs = getFeatureLocsHash(doc, matches);
    if(hashFeatureLocs == null)
      return;
    
    final Hashtable cvTermCache = new Hashtable();
    final Feature f = (Feature)fv.elementAt(0);
    for(int i=0; i<matches.size(); i++)
    {
      org.gmod.schema.sequence.Feature matchFeature =
           (org.gmod.schema.sequence.Feature)matches.get(i);
      
      final CvTerm cvTerm = getCvTermFromCache(cvTermCache, matchFeature, f);
      final String qualifierName;
      if(cvTerm.getName().indexOf("match") > -1)
        qualifierName = "similarity";
      else
        qualifierName = cvTerm.getName();
      

      final List featureLocs = 
        (List) hashFeatureLocs.get(new Integer(matchFeature.getFeatureId()));
      if(featureLocs == null)
        continue;
      
      matchFeature.setFeatureLocsForFeatureId(featureLocs);
      //java.util.Collection featureLocs = matchFeature.getFeatureLocsForFeatureId();
      java.util.Iterator it = featureLocs.iterator();
      while(it.hasNext())
      {
        final FeatureLoc featureLoc = (FeatureLoc)it.next();
        final Feature queryFeature = 
          (Feature)temp_lookup_hash.get(Integer.toString(featureLoc.getSrcFeatureId()));
        
        if(queryFeature != null)
        {
          Qualifier qualifier = queryFeature.getQualifierByName(qualifierName);
          final FeatureLocLazyQualifierValue sim = 
            new FeatureLocLazyQualifierValue(matchFeature, featureLoc.getSrcFeatureId());
          
          if(qualifier == null || !(qualifier instanceof QualifierLazyLoading))
            qualifier = new QualifierLazyLoading(qualifierName, sim);
          else
            ((QualifierLazyLoading)qualifier).addValue(sim);

          try
          {
            queryFeature.setQualifier(qualifier);
          }
          catch(ReadOnlyException e)
          {
            e.printStackTrace();
          }
          catch(EntryInformationException e)
          {
            e.printStackTrace();
          }
          
          if(qualifierName.equals("polypeptide_domain") && 
             LocalAndRemoteFileManager.domainLoad.isSelected()) 
            addDomain(queryFeature, featureLoc, matchFeature);

          break;
        }
      }
    }
    cvTermCache.clear();
    temp_lookup_hash.clear(); 
  }
  
  /**
   * Retrieve a CvTerm from a Hashtable with keys equal to the cvterm_id
   * and values of the corresponding  CvTerm. If the term is not in the cache 
   * then look it up in the main DatabaseDocument cache.
   * @param cvTermCache
   * @param matchFeature
   * @param f
   * @return
   */
  private CvTerm getCvTermFromCache(final Hashtable cvTermCache,
                                    final org.gmod.schema.sequence.Feature matchFeature,
                                    final Feature f)
  {
    final Integer cvTermId = new Integer(matchFeature.getCvTerm().getCvTermId());
    final CvTerm cvTerm;
    
    if(cvTermCache.containsKey(cvTermId))
      cvTerm = (CvTerm) cvTermCache.get(cvTermId);
    else
    {
      cvTerm = DatabaseDocument.getCvTermByCvTermId(
          matchFeature.getCvTerm().getCvTermId(), f);
      cvTermCache.put(cvTermId, cvTerm);
    }
    
    matchFeature.setCvTerm(cvTerm);
    return cvTerm;
  }
  
  /**
   * Add domain features as read-only features
   * @param queryFeature
   * @param featureLoc
   * @param matchFeature
   */
  private void addDomain(final Feature queryFeature, 
                         final FeatureLoc featureLoc,
                         org.gmod.schema.sequence.Feature matchFeature)
  {
    try
    {
      int start = queryFeature.getLocation().getFirstBase();
      Location location = null;
      
      ChadoCanonicalGene chadoGene = ((GFFStreamFeature)queryFeature).getChadoGene();
      if(chadoGene != null)
        location = chadoGene.getNucLocation(queryFeature, featureLoc);
      else if(queryFeature.getLocation().isComplement())
        location = new Location("complement("+
            (start+(featureLoc.getFmin()*3)+1)+".."+(start+(featureLoc.getFmax()*3))+")");
      if(location == null)
        location = new Location(
            (start+(featureLoc.getFmin()*3)+1)+".."+(start+(featureLoc.getFmax()*3)));
   
      final GFFStreamFeature newFeature = new GFFStreamFeature(
          new Key("polypeptide_domain"), 
          location, null);
      newFeature.setReadOnlyFeature(true);
      newFeature.setChadoLazyFeature(matchFeature);
      add(newFeature);
    }
    catch (Exception e)
    {
      e.printStackTrace();
    }
  }
  

  /**
   * Bulk load match features featureLoc's and create a Hashtable with the
   * feature_id's as keys and a list of the corresponding featureLocs as values.
   * @param doc
   * @param matches
   * @return the hashtable; null if no featureLocs are found
   */
  private Hashtable getFeatureLocsHash(final DatabaseDocument doc, final List matches)
  {
    final List matchFeatureIds = new Vector(matches.size());
    for(int i=0; i< matches.size(); i++)
    {
      String matchFeatureId = Integer.toString( 
          ((org.gmod.schema.sequence.Feature)matches.get(i)).getFeatureId() );
      matchFeatureIds.add( matchFeatureId );
    }
    
    final List allFeatureLocs = doc.getFeatureLocsByListOfIds(matchFeatureIds);
    if(allFeatureLocs == null)
      return null;
    
    Hashtable hashFeatureLocs = new Hashtable();
    for(int i=0;i<allFeatureLocs.size(); i++)
    {
      FeatureLoc featureLoc = (FeatureLoc)allFeatureLocs.get(i);
      Integer featureId = new Integer(featureLoc.getFeatureByFeatureId().getFeatureId());
      List list;
      if(hashFeatureLocs.containsKey(featureId))
        list = (List) hashFeatureLocs.get(featureId);
      else
        list = new Vector();
      list.add(featureLoc);
      hashFeatureLocs.put(featureId, list);
    }
    return hashFeatureLocs;
  }

  /**
   *  Combine the features (which are exons) and delete the orignals from this
   *  Entry.  The key of this hash will be the group name and the value is a
   *  FeatureVector containing the feature that are in that group.  Groups
   *  that have more than one member will be combined.
   **/
  private void combineChadoExons(ChadoCanonicalGene gene) 
  {
    final List<Feature> transcripts = gene.getTranscripts();
    gene.correctSpliceSiteAssignments();
    
    for(int i=0; i<transcripts.size(); i++)
    {
      GFFStreamFeature transcript = (GFFStreamFeature)transcripts.get(i);
      String transcript_id = null;
      
      if(transcript.getQualifierByName("ID") == null)
        continue;
      transcript_id = (String)(transcript.getQualifierByName("ID").getValues().get(0));
      
      Set<String> splicedSiteTypes = gene.getSpliceTypes(transcript_id);
      if(splicedSiteTypes == null)
        continue;
      
      Iterator<String> it = splicedSiteTypes.iterator();
      Vector<Feature> new_set = new Vector<Feature>();
      while(it.hasNext())
      {
        String type = it.next();
        List<Feature> splicedSites = gene.getSpliceSitesOfTranscript(transcript_id, type);
           
        if(splicedSites == null)
          continue;
      
        mergeFeatures(splicedSites, new_set, 
                      (String)(transcript.getQualifierByName("ID").getValues().get(0)),
                      transcript.getLocation().isComplement());
      }
      
      for(int j=0; j<new_set.size(); j++)
      {
        if(j == 0)
          gene.addSplicedFeatures(transcript_id, 
                new_set.get(j), true );
        else
          gene.addSplicedFeatures(transcript_id, 
                new_set.get(j));
      }
      
    }   
  }
  
  
  private void mergeFeatures(final List<Feature> gffFeatures,
                             final List<Feature> new_set, 
                             final String transcript_id,
                             final boolean isComplement)
  {
    final Hashtable<String, Integer> feature_relationship_rank_store = new Hashtable<String, Integer>();
    final Hashtable<String, Range> id_range_store = new Hashtable<String, Range>();
    final RangeVector new_range_vector = new RangeVector();
    QualifierVector qualifiers = new QualifierVector();
    Timestamp lasttimemodified = null;

    final Qualifier codon_start = getCodonStart(gffFeatures, isComplement);
    for(int j = 0; j < gffFeatures.size(); j++)
    {
      final GFFStreamFeature this_feature = (GFFStreamFeature)gffFeatures.get(j);
      
      Integer rank;
      Qualifier rankQualifier = this_feature
          .getQualifierByName("feature_relationship_rank");
      if(rankQualifier == null)
        rank = new Integer(0);
      else
      {
        rank = new Integer((String) (rankQualifier.getValues().get(0)));
        this_feature.getQualifiers().removeQualifierByName("feature_relationship_rank");
      }
      
      // use the most current lastmodified datestamp
      if(this_feature.getLastModified() != null
          && (lasttimemodified == null || this_feature.getLastModified()
              .compareTo(lasttimemodified) > 0))
        lasttimemodified = this_feature.getLastModified();

      final Location this_feature_location = this_feature.getLocation();

      if(this_feature_location.getRanges().size() > 1)
      {
        String id= "";
        try
        {
          id = (String)this_feature.getQualifierByName("ID").getValues().get(0);
        }
        catch(Exception e){}
        throw new Error("internal error - new location should have "
            + "exactly one range (there may be non-unique ID's):\n"+
            this_feature_location.toStringShort()+"\n"+id);
      }

      final Range new_range = (Range) this_feature_location.getRanges()
          .elementAt(0);

      Qualifier id_qualifier = this_feature.getQualifierByName("ID");
      if(id_qualifier != null)
      {
        String id = (String) (id_qualifier.getValues()).elementAt(0);
        id_range_store.put(id, new_range);
        feature_relationship_rank_store.put(id, rank);
      }
      else
        Splash.logger4j.warn("NO ID FOUND FOR FEATURE AT: "+
                     this_feature.getLocation().toString());
      
      if(this_feature_location.isComplement())
        new_range_vector.insertElementAt(new_range, 0);
      else
        new_range_vector.add(new_range);

      removeInternal(this_feature);
      qualifiers.addAll(this_feature.getQualifiers());
    }

    final GFFStreamFeature first_old_feature = (GFFStreamFeature)gffFeatures.get(0);

    final Location new_location = new Location(new_range_vector,
        first_old_feature.getLocation().isComplement());

    if(codon_start != null)
    {
      QualifierVector tmp_qualifier_vector = new QualifierVector();

      for(Qualifier q: qualifiers)
        if(!q.getName().equals("codon_start"))
          tmp_qualifier_vector.addElement(q);
      qualifiers = tmp_qualifier_vector;
      qualifiers.setQualifier(codon_start);
    }

    qualifiers = mergeQualifiers(qualifiers);

    final GFFStreamFeature new_feature = new GFFStreamFeature(first_old_feature
        .getKey(), new_location, qualifiers);

    if(lasttimemodified != null)
      new_feature.setLastModified(lasttimemodified);
    
    if(first_old_feature.getChadoGene() != null)
      new_feature.setChadoGene(first_old_feature.getChadoGene());
    
    new_feature.setSegmentRangeStore(id_range_store);
    new_feature
        .setFeature_relationship_rank_store(feature_relationship_rank_store);
    new_feature.setGffSource(first_old_feature.getGffSource());
    new_feature.setGffSeqName(first_old_feature.getGffSeqName());
    
//  set the ID
    String ID;
    try
    {
      ID = new_feature.getSegmentID(new_feature.getLocation().getRanges());
    }
    catch(NullPointerException npe)
    {
      if(new_feature.getQualifierByName("Parent") != null)
        ID = ((String)new_feature.getQualifierByName("Parent").getValues().get(0)) +
            ":"+new_feature.getKey().getKeyString();
      else
        ID = new_feature.getKey().getKeyString();
    }
    final Qualifier id_qualifier = new_feature.getQualifierByName("ID");
    id_qualifier.removeValue((String)(id_qualifier.getValues()).elementAt(0));
    id_qualifier.addValue(ID);
    
    
    // set visibility
    if(GeneUtils.isHiddenFeature(new_feature.getKey().getKeyString()) ||
       GeneUtils.isObsolete(new_feature))
      new_feature.setVisible(false);
    
    try
    {
      new_feature.setLocation(new_location);

      final Qualifier gene_qualifier = new_feature.getQualifierByName("gene");

      if(gene_qualifier != null
          && gene_qualifier.getValues().size() > 0
          && ((String) (gene_qualifier.getValues()).elementAt(0))
              .startsWith("Phat"))
      {
        // special case to handle incorrect output of the Phat gene
        // prediction tool
        new_feature.removeQualifierByName("codon_start");
      }

      forcedAdd(new_feature);
      //gene.addExon(transcript_id, new_feature, true );
      new_set.add(new_feature);
    }
    catch(ReadOnlyException e)
    {
      throw new Error("internal error - unexpected exception: " + e);
    }
    catch(OutOfRangeException e)
    {
      throw new Error("internal error - unexpected exception: " + e);
    }
    catch(EntryInformationException e)
    {
      throw new Error("internal error - unexpected exception: " + e);
    }
  }

  private QualifierVector mergeQualifiers(QualifierVector qualifier_vector)
  {
    QualifierVector merge_qualifier_vector = new QualifierVector();
    boolean seen = false;

    for(int i = 0 ; i < qualifier_vector.size() ; ++i)
    {
      Qualifier qual = (Qualifier)qualifier_vector.elementAt(i);
 
      if(qual.getName().equals("codon_start"))
      {
        if(!seen)
        {
          merge_qualifier_vector.addElement(qual);
          seen = true;
        }
      }
      else if(qual.getName().equals("Alias"))
      { 
        final Qualifier id_qualifier = 
          merge_qualifier_vector.getQualifierByName("Alias");

        if(id_qualifier == null)
          merge_qualifier_vector.addElement(qual);
        else
        {
          String id1 = (String)(id_qualifier.getValues()).elementAt(0);
          String id2 = (String)(qual.getValues()).elementAt(0);
          id_qualifier.removeValue(id1);
          id_qualifier.addValue(id1+","+id2);
        }
      }
      else if(!qual.getName().equals("ID") &&
              !qual.getName().equals("feature_id"))
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
        removeInternal(this_feature);
        qualifier_vector.addAll(this_feature.getQualifiers());
        
        final Range new_range = (Range) this_feature.getLocation().getRanges().elementAt(0);
        if(this_feature.getLocation().isComplement())
          new_range_vector.insertElementAt(this_feature.getLocation().getTotalRange(), 0);
        else
          new_range_vector.add(new_range);
      }
      final GFFStreamFeature old_feature = (GFFStreamFeature)this_group.get(0);

      final Location new_location = new Location(new_range_vector,
          old_feature.getLocation().isComplement());
      
      qualifier_vector = mergeQualifiers(qualifier_vector);
      if(qualifier_vector.getQualifierByName("gene_id") != null)
        qualifier_vector.addQualifierValues(new Qualifier("ID",
            keyStr+":"+qualifier_vector.getQualifierByName("gene_id").getValues().get(0)));
      
      final GFFStreamFeature new_feature = new GFFStreamFeature(old_feature
          .getKey(), new_location, qualifier_vector);
      forcedAdd(new_feature);
    }
  }
  
  /**
   * Get the phase/codon_start for the first feature segment
   * @param gffFeatures
   * @param isComplement
   * @return
   */
  private Qualifier getCodonStart(final List<Feature> gffFeatures, final boolean isComplement)
  {
    int fstart = (isComplement ? 0 : Integer.MAX_VALUE);
    Feature firstFeature = null;
    for(Feature f: gffFeatures)
    {
      final GFFStreamFeature this_feature = (GFFStreamFeature)f;
      if(isComplement && this_feature.getFirstBase() > fstart)
      {
        firstFeature = this_feature;
        fstart = this_feature.getFirstBase();
      }
      else if(!isComplement && this_feature.getFirstBase() < fstart)
      {
        firstFeature = this_feature;
        fstart = this_feature.getFirstBase();
      }
    }
    
    if(firstFeature == null)
      return null;
    try
    {
      Qualifier codon_start = firstFeature.getQualifierByName("codon_start");
      if(codon_start != null)
        return codon_start.copy();
    }
    catch (InvalidRelationException e){}
    return null;
  }

  /**
   * Adjust feature coordinates to match the contig positions when loaded
   * with a multiple fasta. 
   * @param sequenceEntry   sequence entry
   */
  public void adjustCoordinates(uk.ac.sanger.artemis.Entry sequenceEntry)
  {
    final Entry entry;
    if(sequenceEntry != null)
      entry = sequenceEntry.getEMBLEntry();
    else
      entry = this;
    if(entry instanceof SimpleDocumentEntry)
    {
      // adjust feature coordinates to match contig positions
      final Hashtable<String, Range> contig_ranges = ((SimpleDocumentEntry)entry).contig_ranges;
      if(contig_ranges != null)
      {
        final FeatureVector gff_regions = getAllFeatures();
        final Enumeration<Feature> gff_features  = gff_regions.elements();

        while(gff_features.hasMoreElements())
        {
          final Feature f = gff_features.nextElement();
          if( !(f instanceof GFFStreamFeature) )
            continue;
          
          final String name = ((GFFStreamFeature)f).getGffSeqName();
          if(name == null)
            continue;
          if(contig_ranges.containsKey(name))
          {
            
            try
            {
              final Range new_range = contig_ranges.get(name);
              final RangeVector new_ranges = new RangeVector();
              final RangeVector ranges = f.getLocation().getRanges();
              for(int i = 0 ; i<ranges.size () ; ++i) 
              {
                final Range r = (Range)ranges.elementAt(i);

                new_ranges.add(new Range(r.getStart()+new_range.getStart()-1, 
                                         r.getEnd()+new_range.getStart()-1));
              }

              Location l = new Location(new_ranges, f.getLocation().isComplement());
              f.setLocation(l);
              ((uk.ac.sanger.artemis.Feature)f.getUserData()).setLocation(l);
            }
            catch(OutOfRangeException e)
            {
              throw new Error("internal error - unexpected exception: " + e);
            }
            catch (ReadOnlyException e)
            {
              e.printStackTrace();
            }
          }
        }
        // store so these can be used when writing out
        GFFStreamFeature.contig_ranges = contig_ranges;
      }
    }
  }

}