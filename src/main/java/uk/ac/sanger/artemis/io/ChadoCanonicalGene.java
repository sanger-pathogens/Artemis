/* ChadoCanonicalGene.java
 *
 * created: 2006
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2006 Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/ChadoCanonicalGene.java,v 1.34 2009-08-11 08:59:46 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.components.genebuilder.GeneUtils;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.ReadOnlyException;
import uk.ac.sanger.artemis.util.StringVector;

import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.Vector;
import java.util.Hashtable;
import java.util.Enumeration;
import java.util.List;
import java.util.Set;
import java.util.HashSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.gmod.schema.sequence.FeatureLoc;

/**
 *  Used by GFFStreamFeature to represent the chado canonical gene.
 *  Contains gene, transcript, exons and proteins.
 **/
public class ChadoCanonicalGene  
{
  private Feature gene;
  
  // part_of gene
  private List<Feature> transcripts = new Vector<Feature>();
  
  // part_of transcripts
  private Hashtable<String, List<Feature>> splicedFeatures =
    new Hashtable<String, List<Feature>>();
  
  // derives_from transript
  private Hashtable<String, Feature> proteins = new Hashtable<String, Feature>();

  // utr features
  private Hashtable<String,  List<Feature>> three_prime_utr =
    new Hashtable<String,  List<Feature>>();
  private Hashtable<String,  List<Feature>> five_prime_utr  = 
    new Hashtable<String,  List<Feature>>();
  
  // other child features of transcript
  private Hashtable<String,  List<Feature>> other_features = 
    new Hashtable<String,  List<Feature>>();
  
  // srcfeature
  private int srcfeature_id;
 
  // srcfeature length
  private int seqlen;
  
  
  /**
   * Get the gene feaure object.
   * @return
   */
  public Feature getGene()
  {
    return gene;
  } 
 
  /**
   * Set the gene feature object.
   * @param gene
   */
  public void setGene(Feature gene)
  {
    this.gene = gene;
  }
 
  public String getGeneUniqueName()
  {
    try
    {
      return getQualifier(getGene(), "ID");
    }
    catch(InvalidRelationException e)
    {
      return null;
    } 
  }
  
  /**
   * Add a transcript to the model
   * @param transcript
   */
  public void addTranscript(Feature transcript)
  {
    transcripts.add(transcript);
  }
  
  /**
   * Delete a transcript and child features.
   * @param transcript_name
   */
  public void deleteTranscript(String transcript_name)
  {
    for(int i=0; i<transcripts.size(); i++)
    {
      try
      {
        Feature transcript = (Feature)transcripts.get(i);
        
        if( transcript_name.equals(getQualifier(transcript, "ID")) )
        {
          transcripts.remove(transcript);
          splicedFeatures.remove(transcript_name);
          three_prime_utr.remove(transcript_name);
          five_prime_utr.remove(transcript_name);
          other_features.remove(transcript_name);
          proteins.remove(transcript_name);
        }
      }
      catch(InvalidRelationException e)
      {
        e.printStackTrace();
      }
    }  
  }
  
  /**
   * This should be called if the uniqueName of a gene model 
   * feature is changed.
   * @param oldName
   * @param newName
   * @param children
   */
  public void updateUniqueName(final String oldName, 
                               final String newName,
                               final Set<Feature> children)
  {
    updateNames(splicedFeatures,oldName,newName);
    updateNames(proteins,oldName,newName);
    updateNames(three_prime_utr,oldName,newName);
    updateNames(five_prime_utr,oldName,newName);
    updateNames(other_features,oldName,newName);
    
    if(children != null)
      GeneUtils.fixParentQualifier(oldName, newName, children);
  }
  
  /**
   * Utility for changing the key used in a Hashtable
   * @param hash
   * @param oldName
   * @param newName
   */
  private static void updateNames(final Hashtable hash,
                                  final String oldName, 
                                  final String newName)
  {
    Object features = hash.get(oldName);
    if(features != null)
    {
      hash.remove(oldName);
      hash.put(newName, features);
    }
  }
  
  /**
   * Delete features.
   * @param embl_feature
   */
  public void deleteFeature(final Feature embl_feature)
  {
    try
    {
      final String name = getQualifier(embl_feature, "ID");
      Object feature = getSplicedFeatures(name);

      if(feature != null)
      {
        String transcript_name = getQualifier((Feature) feature, "Parent");
        splicedFeatures.remove(transcript_name);
        return;
      }
      
      final Enumeration<String> enum_protein = proteins.keys();
      while(enum_protein.hasMoreElements())
      {
        final String transcriptName = (String)enum_protein.nextElement();
        Feature protein = (Feature)proteins.get(transcriptName);
        if(getQualifier(protein, "ID").equals(name))
        {
          proteins.remove(transcriptName);
          return;
        }
      }
      
      feature = getFeatureFromHash(name, three_prime_utr);
      if(feature != null)
      {
        String transcript_name = getQualifier((Feature) feature, "Parent");
        List<Feature> utr = get3UtrOfTranscript(transcript_name);
        utr.remove(feature);
        return;
      }
      
      feature = getFeatureFromHash(name, five_prime_utr);
      if(feature != null)
      {
        String transcript_name = getQualifier((Feature) feature, "Parent");
        List<Feature> utr = get5UtrOfTranscript(transcript_name);
        utr.remove(feature);
        return;
      }
      
      feature = getFeatureFromHash(name, other_features);
      if(feature != null)
      {
        String transcript_name = getQualifier((Feature) feature, "Parent");
        List<Feature> others = getOtherFeaturesOfTranscript(transcript_name);
        others.remove(feature);
        return;
      }
      
      deleteTranscript(name);
    }
    catch(InvalidRelationException e1)
    {
      e1.printStackTrace();
    }
  }
  
  /**
   * Get all child members of a feature
   * @param embl_feature
   * @return
   */
  public Set<Feature> getChildren(Feature embl_feature)
  {
    Set<Feature> children = new HashSet<Feature>();
    try
    {
      String name = getQualifier(embl_feature, "ID");
      
      String gene_name = getQualifier(getGene(), "ID");
      if(name.equals(gene_name))
      {
        List<Feature> transcripts = getTranscripts();
        for(int i=0; i<transcripts.size(); i++)
        {
          Feature transcript = transcripts.get(i);
          children.add(transcript);
          children.addAll( getChildren(transcript) );
        }
        return children;
      }
      
      searchForChildren(splicedFeatures, name, children);
      searchForChildren(three_prime_utr, name, children);
      searchForChildren(five_prime_utr, name, children);
      searchForChildren(other_features, name, children);
      
      // protein
      Enumeration<Feature> pep_enum = proteins.elements();
      while(pep_enum.hasMoreElements())
      {
        Feature child = pep_enum.nextElement();
        String parent = getQualifier(child, "Derives_from");
        if(parent != null && parent.equals(name))
          children.add(child);
      }
      return children;
    }
    catch(InvalidRelationException e)
    {
      e.printStackTrace();
    }
    return null;
  }
  
  /**
   * Search in a <code>Hashtable</code> for child Features with a
   * matching parent ID. Child features are added to the <code>Set</code>
   * that is passed into this method.
   * @param hash        Hashtable to search for children in
   * @param parent_id   uniquname to look for      
   * @param children    collection to add child features to
   * @throws InvalidRelationException
   */
  private void searchForChildren(Hashtable<String, List<Feature>> hash,
                                 String parent_id,
                                 Set<Feature> children)
               throws InvalidRelationException
  {
    Enumeration<List<Feature>> feature_enum = hash.elements();
    String parent;
    
    while(feature_enum.hasMoreElements())
    {
      List<Feature> child_list = feature_enum.nextElement();
      
      for(int i=0; i<child_list.size(); i++)
      {       
        Feature child = child_list.get(i);
        //if(children.contains(child))
        //  continue;
        
        parent = getQualifier(child, "Parent");
        if(parent != null && parent.equals(parent_id))
          children.add(child);
        else 
        {
          parent = getQualifier(child, "Derives_from");
          if(parent != null && parent.equals(parent_id))
            children.add(child);
        }
      } 
    }
  }
  
  /**
   * Add exon feature to the chado gene model.
   * @param transcript_name
   * @param exon
   * @param reset
   * @throws InvalidRelationException
   */
  public void addSplicedFeatures(final String transcript_name, 
                      final Feature exon, boolean reset) 
  {
    if(reset)
      splicedFeatures.remove(transcript_name);
    addSplicedFeatures(transcript_name, exon);
  }
  
  /**
   * Add exon feature to the chado gene model.
   * @param transcript_name
   * @param v_spliced
   * @throws InvalidRelationException
   */
  public void addSplicedFeatures(final String transcript_name, 
                                 final Feature spliced)
  {   
    final List<Feature> v_spliced;
    if(splicedFeatures.containsKey(transcript_name))
      v_spliced = (Vector<Feature>)splicedFeatures.get(transcript_name);
    else
      v_spliced = new Vector<Feature>();
    
    v_spliced.add(spliced);
    splicedFeatures.put(transcript_name, v_spliced);
  }
  
  public void correctSpliceSiteAssignments()
  {
    Enumeration<String> enumSplicedFeatures = splicedFeatures.keys();
    while(enumSplicedFeatures.hasMoreElements())
    {
      String transcriptId = enumSplicedFeatures.nextElement();
      Vector<Feature> v_spliced = (Vector<Feature>)splicedFeatures.get(transcriptId);
      Set<String> splicedTypes = getSpliceTypes(transcriptId);
      Iterator<String> it = splicedTypes.iterator();
      while(it.hasNext())
      {
        String type = it.next();
        if(!type.equals(DatabaseDocument.EXONMODEL) &&
           !type.equals("pseudogenic_exon") &&
           !type.equals("exon"))
        {
          List<Feature> splicedFeatures = getSpliceSitesOfTranscript(transcriptId, type);
          if(splicedFeatures.size() == 1)
          {
            Feature f = (Feature)splicedFeatures.get(0);
            addOtherFeatures(transcriptId, f);
            v_spliced.remove(f);
            try
            {
              f.removeQualifierByName("feature_relationship_rank");
            }
            catch(ReadOnlyException e){}
            catch(EntryInformationException e){}
          }
        }
      }
      splicedFeatures.put(transcriptId, v_spliced);
    }
  }
  
  /**
   * Add protein feature to the chado gene model.
   * @param transcript_name
   * @param protein
   * @throws InvalidRelationException
   */
  public void addProtein(final String transcript_name,
                         final Feature protein) 
  {   
    proteins.put(transcript_name, protein);
  }
  
  /**
   * Add 3'UTR to chado gene model.
   * @param transcript_name
   * @param utr
   * @throws InvalidRelationException
   */
  public void add3PrimeUtr(final String transcript_name, 
                           final Feature utr) 
  {  
    final List<Feature> utr_list;
    if(three_prime_utr.containsKey(transcript_name))
      utr_list = three_prime_utr.get(transcript_name);
    else
      utr_list = new Vector<Feature>();
    
    utr_list.add(utr);
    three_prime_utr.put(transcript_name, utr_list);
  }
  
  /**
   * Add 5'UTR to chado gene model.
   * @param transcript_name
   * @param utr
   * @throws InvalidRelationException
   */
  public void add5PrimeUtr(final String transcript_name, 
                           final Feature utr)
  {
    final List<Feature> utr_list;
    if(five_prime_utr.containsKey(transcript_name))
      utr_list = (Vector<Feature>)five_prime_utr.get(transcript_name);
    else
      utr_list = new Vector<Feature>();
    
    utr_list.add(utr);
    five_prime_utr.put(transcript_name, utr_list);
  }
  
  /**
   * Add other child features of a transcript to the chado
   * gene model.
   * @param transcript_name
   * @param other_feature
   */
  public void addOtherFeatures(final String transcript_name, 
                               final Feature other_feature)
  {
    final List<Feature> v_other_features;
    if(other_features.containsKey(transcript_name))
      v_other_features = (Vector<Feature>)other_features.get(transcript_name);
    else
      v_other_features = new Vector<Feature>();
    v_other_features.add(other_feature);
    other_features.put(transcript_name, v_other_features);
  }
  
  /**
   * Check if this gene model contains a transcript with an ID equal to
   * any of the names in the <code>StringVector</code>. If it does find 
   * it returns the transcript feature, otherwise it returns null.
   * @param names
   * @return
   */
  public Feature containsTranscript(final StringVector names)
  {
	final int numTranscripts = transcripts.size();
    for(int i=0; i < numTranscripts; i++)
    {
      try
      {
        Feature transcript = (Feature)transcripts.get(i);
        
        if( names.contains(getQualifier(transcript, "ID")) )
          return transcript;
      }
      catch(InvalidRelationException e)
      {
        e.printStackTrace();
      }
    }
    return null;
  }
 
  
  public List<Feature> getSpliceSitesOfTranscript(final String transcript_name,
                                         final String type)
  {
    if(splicedFeatures.containsKey(transcript_name))
    {
      List<Feature> splicedFeaturesOfTranscript = splicedFeatures.get(transcript_name);
      List<Feature> results = new Vector<Feature>();
      
      int numSplicedFeatures = splicedFeaturesOfTranscript.size();
      for(int i=0; i < numSplicedFeatures; i++)
      {
        Feature feature = (Feature)splicedFeaturesOfTranscript.get(i);
        if(feature.getKey().getKeyString().equals(type))
          results.add(feature);
      }
      return results;
    }
 
    return null;   
  }
  
  /**
   * Get a list of the feature keys of the types that are splice sites
   * @param transcript_name
   * @return
   */
  public Set<String> getSpliceTypes(final String transcript_name)
  {
    if(splicedFeatures.containsKey(transcript_name))
    {
      List<Feature> splicedFeaturesOfTranscript = splicedFeatures.get(transcript_name);
      Set<String> splicedTypes = new HashSet<String>();
      
      int numSplicedfeatures = splicedFeaturesOfTranscript.size();
      for(int i=0; i < numSplicedfeatures; i++)
      {
        Feature feature = (Feature)splicedFeaturesOfTranscript.get(i);
        splicedTypes.add( feature.getKey().getKeyString() );
      }
      return splicedTypes;
    }
 
    return null; 
  }

  /**
   * Return the exons of a given transcript as a <code>List</code>.
   * @param transcript_name
   * @return
   */
  public List<Feature> getSplicedFeaturesOfTranscript(final String transcript_name)
  {
    if(splicedFeatures.containsKey(transcript_name))
    {
      return splicedFeatures.get(transcript_name);
    }
 
    return null;
  }
  
  /**
   * Return the transcript from the name of a constituent feature
   * @param constituent feature name
   * @return transcript 
   */
  public Feature getTranscriptFeatureFromName(final String name)
  {
    String transcriptName = getTranscriptFromName(name);
    if(transcriptName == null)
      return null;
    
    try
    {
      for (int i = 0; i < transcripts.size(); i++)
      {
        Feature feature = (Feature) transcripts.get(i);
        if (getQualifier(feature, "ID").equals(transcriptName))
          return feature;
      }
    }
    catch (InvalidRelationException ire){}
    return null;
  }
  
  /**
   * Return the transcript from the name of a constituent feature
   * @param constituent feature name
   * @return transcript name
   */
  public String getTranscriptFromName(final String name)
  {
    //  check transcript
    StringVector sv = new StringVector();
    sv.add(name);    
    Feature feature = containsTranscript(sv);
    
    if(feature != null)
      return name;

    // check exons
    List<String> transcriptNames = getTranscriptNames();
    feature = getSplicedFeatures(name);
    
    if(feature != null)
    {
      for(int i=0; i<transcriptNames.size(); i++)
      {
        String transcriptName = (String)transcriptNames.get(i);
        List<Feature> splicedSegments = getSplicedFeaturesOfTranscript(transcriptName);
        
        if(splicedSegments != null)
        {
          for(int j=0; j<splicedSegments.size(); j++)
          {
            Feature segment = splicedSegments.get(j);
            try
            {
              String segmentName = (String)segment.getQualifierByName("ID").getValues().get(0);
              if(name.equals(segmentName))
                return transcriptName;
            }
            catch(InvalidRelationException e)
            {
              // TODO Auto-generated catch block
              e.printStackTrace();
            }
          }
        }
      }
    }
    
    feature = getProtein(name);
    
    if(feature != null)
    {
      for(int i=0; i<transcriptNames.size(); i++)
      {
        String transcriptName = (String)transcriptNames.get(i);
        Feature protein = getProteinOfTranscript(transcriptName);
        try
        {
          String proteinsName = (String)protein.getQualifierByName("ID").getValues().get(0);
          if(name.equals(proteinsName))
            return transcriptName;
        }
        catch(InvalidRelationException e)
        {
          // TODO Auto-generated catch block
          e.printStackTrace();
        }
      }
    }
    
    // search children of all transcripts
    List<Feature> transcripts = getTranscripts();
    for(int i=0;i<transcripts.size(); i++)
    {
      Feature transcript = transcripts.get(i);
      Set<Feature> children = getChildren(transcript);
      Iterator<Feature> it = children.iterator();
      while(it.hasNext())
      {
        Feature f = it.next();
        if(name.equals(GeneUtils.getUniqueName(f)))
          return GeneUtils.getUniqueName(transcript);
      }
    }
 
    return null;
  }
  
  /**
   * Return the protein feature of a transcipt.
   * @param transcript_name
   * @return
   */
  public Feature getProteinOfTranscript(final String transcript_name)
  {
    if(proteins.containsKey(transcript_name))
      return (Feature)proteins.get(transcript_name);;
 
    return null;
  }

  /**
   * Return the 3'UTR features of a transcriot as a <code>List</code>.
   * @param transcript_name
   * @return
   */
  public List<Feature> get3UtrOfTranscript(final String transcript_name)
  {
    if(three_prime_utr.containsKey(transcript_name))
      return (List<Feature>)three_prime_utr.get(transcript_name);
 
    return null;
  }
  
  /**
   * Return the 5'UTR features of a transcriot as a <code>List</code>.
   * @param transcript_name
   * @return
   */
  public List<Feature> get5UtrOfTranscript(final String transcript_name)
  {
    if(five_prime_utr.containsKey(transcript_name))
      return (List<Feature>)five_prime_utr.get(transcript_name);
 
    return null;
  }
  
  /**
   * Utility to determine if this is the first or only UTR, so that
   * partial qualifiers can be added to the correct UTR feature.
   * @param utrName
   * @param isFwd
   * @return
   */
  public boolean isFirstUtr(final String utrName, final boolean isFwd)
  {
    try
    {
      Feature this5Utr = getFeatureFromHash(utrName, five_prime_utr);
      if (this5Utr != null)
      {
        String transcript_name = getQualifier(this5Utr, "Parent");
        List<Feature> utrs = get5UtrOfTranscript(transcript_name);
        if (utrs.size() == 1)
          return true;

        for (Feature utr : utrs)
        {
          if (isFwd && utr.getFirstBase() < this5Utr.getFirstBase())
            return false;
          else if (!isFwd && utr.getLastBase() > this5Utr.getLastBase())
            return false;
        }
        return true;
      }
      
      
      Feature this3Utr = getFeatureFromHash(utrName, three_prime_utr);
      if (this3Utr != null)
      {
        String transcript_name = getQualifier(this3Utr, "Parent");
        List<Feature> utrs = get3UtrOfTranscript(transcript_name);
        if (utrs.size() == 1)
          return true;

        for (Feature utr : utrs)
        {
          if (!isFwd && utr.getFirstBase() < this3Utr.getFirstBase())
            return false;
          else if (isFwd && utr.getLastBase() > this3Utr.getLastBase())
            return false;
        }
        return true;
      }
    }
    catch(InvalidRelationException ire){}
    return false;
  }
  
  /**
   * Return the other child features of a transcriot as a <code>List</code>.
   * @param transcript_name
   * @return
   */
  public List<Feature> getOtherFeaturesOfTranscript(final String transcript_name)
  {
    if(other_features.containsKey(transcript_name))
      return other_features.get(transcript_name);
    return null;
  }
  
  /**
   * Get a list of trancripts.
   * @return
   */
  public List<Feature> getTranscripts()
  {
    return transcripts;
  }
  
  
  /**
   * Get a list of trancripts.
   * @return
   */
  private List<String> getTranscriptNames()
  {
    List<String> names = new Vector<String>();
    for(int i=0; i<transcripts.size(); i++)
    {
      Feature f = (Feature)transcripts.get(i);
      try
      {
        names.add( (String)f.getQualifierByName("ID").getValues().get(0) );
      }
      catch(InvalidRelationException e)
      {
        // TODO Auto-generated catch block
        e.printStackTrace();
      }
      
    }
    
    return names;
  }
  
  /**
   * Test if a name is already used in this gene model
   * @param name
   * @return
   */
  private boolean isUniqueName(final String name)
  {
    if(isTranscript(name))
      return false;
    if(isSplicedFeatures(name))
      return false;

    try
    {
      if(getFeatureFromHash(name, three_prime_utr) != null)
        return false;
      if(getFeatureFromHash(name, five_prime_utr) != null)
        return false;
      if(getFeatureFromHash(name, other_features) != null)
        return false;
      
      final Enumeration<Feature> enum_pp = proteins.elements();
      while(enum_pp.hasMoreElements())
      {
        final Feature pp = enum_pp.nextElement();
        if( getQualifier(pp, "ID").equals(name) )
          return false;
      }  
      
      if( getQualifier(getGene(), "ID").equals(name) )
        return false;
    }
    catch(InvalidRelationException e)
    {
      e.printStackTrace();
    }
    
    return true;
  }
  
  /**
   * Test if the name is a transcript in this gene model.
   * @param feature_id
   * @return true if a transcript
   */
  public boolean isTranscript(final String feature_id)
  {
    try
    {
      for(int i=0; i<transcripts.size(); i++)
      {
        if(feature_id.equals(getQualifier((Feature)transcripts.get(i), "ID")))
          return true;
      }
    }
    catch(InvalidRelationException e)
    {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    
    return false;
  }
  
  /**
   * Test if this is an exon of transcript.
   * @param feature_id    exon feature
   * @param transcript_id transcript feature
   * @return
   */
  private boolean isSplicedFeatures(final String feature_id)
  {
    List<Feature> splicedFeatures = new Vector<Feature>();
    List<Feature> transcripts = getTranscripts();
    
    try
    {
      for(int i = 0; i < transcripts.size(); i++)
      {
        Feature transcript = (Feature) transcripts.get(i);
        String transcript_id = getQualifier(transcript, "ID");
        List<Feature> splicedSites = getSplicedFeaturesOfTranscript(transcript_id);
        if(splicedSites != null)
          splicedFeatures.addAll(splicedSites);
      }

      if(splicedFeatures == null)
        return false;

      for(int i=0; i<splicedFeatures.size(); i++)
      {
        GFFStreamFeature feature = (GFFStreamFeature)splicedFeatures.get(i);
        RangeVector rv = feature.getLocation().getRanges();
        for(int j=0; j<rv.size(); j++)
        {
          String this_feature_id = feature.getSegmentID((Range)rv.get(j));
          if(feature_id.equals(this_feature_id))
            return true;
        }
      }
    }
    catch(InvalidRelationException e)
    {
      e.printStackTrace();
    }
    
    return false;
  }
  
  /**
   * Method to automatically generate ID's for transcripts
   * @param transcript_key
   * @return
   */
  public String autoGenerateTanscriptName(String transcript_key)
  {
    try
    {
      String name = getQualifier(getGene(), "ID");
      int auto = 1;
      while( isTranscript( name + "." + auto ) &&
             auto < 50)
        auto++;
      return name + "." + auto;
    }
    catch(InvalidRelationException e)
    {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    return null;
  }
  
  /**
   * Generate new names for exon features for this gene model
   * @param transcript_id
   * @return
   */
  public String autoGenerateSplicedFeatureName(final String transcript_id)
  {
    try
    {
      int index = transcript_id.lastIndexOf('.');
      if(index == -1)
        index = transcript_id.lastIndexOf(':');
      int transcript_number = -1;
      String name = (String)getGene().getQualifierByName("ID").getValues().get(0);
      
      if(index > -1)
      {
        try
        {
          transcript_number = Integer.parseInt(transcript_id.substring(index+1));
        }
        catch(NumberFormatException nfe)
        {
          transcript_number = -1;
        }
      }
      
      if(transcript_number < 1)
      {
        for(transcript_number = 0; transcript_number <= transcripts.size(); 
            transcript_number++)
        {
          Feature transcript = (Feature) transcripts.get(transcript_number);
          if(transcript_id.equals(getQualifier(transcript, "ID")))
            break;
        }
      }
      if(transcript_number == 0)
        name = name + ":exon:";
      else
        name = name + "." + transcript_number + ":exon:";
      
      int auto = 1;
      while( isSplicedFeatures(name + auto) && auto < 50)
        auto++;
      return name + auto;
    }
    catch(InvalidRelationException e)
    {
      e.printStackTrace();
    }
    return null;
  }
  
  
  /**
   * Generate new names for peptide features for this gene model
   * @param transcript_id
   * @return
   */
  public String autoGeneratePepName(final String transcript_id)
  {
    try
    {
      int index = transcript_id.lastIndexOf('.');
      if(index == -1)
        index = transcript_id.lastIndexOf(':');
      int transcript_number = -1;
      
      if(index > -1)
      {
        try
        {
          transcript_number = Integer.parseInt(transcript_id.substring(index+1));
        }
        catch(NumberFormatException nfe)
        {
          transcript_number = -1;
        }
      }
      
      if(transcript_number < 1)
      {
        for(transcript_number = 1; transcript_number <= transcripts.size(); 
            transcript_number++)
        {
          Feature transcript = (Feature) transcripts.get(transcript_number - 1);
          if(transcript_id.equals(getQualifier(transcript, "ID")))
            break;
        }
      }
      
      String name = (String)getGene().getQualifierByName("ID").getValues().get(0);

      if(isUniqueName(name+ "." + transcript_number + ":pep"))
        return name+ "." + transcript_number + ":pep";
      else
        return name + "." + transcript_number + "a:pep";
    }
    catch(InvalidRelationException e)
    {
      e.printStackTrace();
    }
    return null;
  }
  
  /**
   * Generate new names for generic region features for this gene model
   * @param transcript_id
   * @return
   */
  public String autoGenerateFeatureName(final String transcript_id, 
                                        final String keyName)
  {
    String featureName = "";
    try
    {
      featureName = 
        (String)getGene().getQualifierByName("ID").getValues().get(0);
    }
    catch(InvalidRelationException e){}
    
    final Pattern pattern = Pattern.compile("\\d+$");
    final Matcher matcher = pattern.matcher(transcript_id);
    if(matcher.find())
      featureName = featureName+"."+matcher.group()+":"+keyName;
    else
      featureName = featureName+":"+keyName;
    
    if(!isUniqueName(featureName))
    {
      int num = 1;
      while(!isUniqueName(featureName + ":" + num) && num < 100)
        num++;
      featureName = featureName + ":" + num;
    }
    
    return featureName;
  }
  
  /**
   * Search for the feature with a particular uniquename
   * @param name  uniquename
   * @return
   */
  public Object getFeatureFromId(final String name)
  {
    Object feature = null;
    
    // check gene
    try
    {
      final String uniquename = getQualifier(gene, "ID");    
      
      if(uniquename.equals(name))
        return gene;
    }
    catch(InvalidRelationException e)
    {
      e.printStackTrace();
    }
    
    // check transcript
    StringVector sv = new StringVector();
    sv.add(name);
    
    feature = containsTranscript(sv);
    
    if(feature != null)
      return feature;

    // check exons
    feature = getSplicedFeatures(name);
    
    if(feature != null)
      return feature;
    
    feature = getProtein(name);
    
    if(feature != null)
      return feature;
    
    try
    {
      feature = getFeatureFromHash(name, three_prime_utr);
      if(feature != null)
        return feature;
      
      feature = getFeatureFromHash(name, five_prime_utr);
      if(feature != null)
        return feature;
      
      feature = getFeatureFromHash(name, other_features);
    }
    catch(InvalidRelationException e)
    {
      e.printStackTrace();
    }
   
    return feature;
  }
  
  /**
   * Routine to look for a exon with a particular 
   * uniquename
   * @param name
   * @return
   */
  private Feature getSplicedFeatures(final String name)
  {
    Enumeration<List<Feature>> enum_exons = splicedFeatures.elements();
    try
    {
      while(enum_exons.hasMoreElements())
      {
        List<Feature> exons = enum_exons.nextElement();
        
        for(int i=0; i<exons.size(); i++)
        {
          String uniquename = getQualifier((Feature)exons.get(i), "ID");
          
          if(uniquename.equals(name))
            return (Feature)exons.get(i);
        }
      }
    }
    catch(InvalidRelationException e)
    {
      e.printStackTrace();
    }
    return null;
  }
  
  private Feature getProtein(final String id)
  {
    Enumeration<Feature> enum_proteins = proteins.elements();
    try
    {
      while(enum_proteins.hasMoreElements())
      {
        Feature protein = enum_proteins.nextElement();
        if(getQualifier(protein, "ID").equals(id))
          return protein;
      }
    }
    catch(InvalidRelationException e)
    {
      e.printStackTrace();
    }
    return null;
  }

  /**
   * Search for a feature uniquename 
   * @param id
   * @param UTR
   * @return
   * @throws InvalidRelationException
   */
  private Feature getFeatureFromHash
                         (final String id,
                         final Hashtable<String, List<Feature>> UTR) 
          throws InvalidRelationException
  {
    Enumeration<List<Feature>> enum_utr = UTR.elements();

    while(enum_utr.hasMoreElements())
    {
      List<Feature> utrs = enum_utr.nextElement();
      
      for(int i=0; i<utrs.size(); i++)
      {
        Feature utr = utrs.get(i);
        if(getQualifier(utr, "ID").equals(id))
          return utr;
      }
    }

    return null;
  }
  
  /**
   * Utility for get feature ID and Parent qualifiers.
   * @param feature
   * @param name
   * @return
   * @throws InvalidRelationException
   */
  public String getQualifier(final Feature feature,
                              final String name) 
          throws InvalidRelationException
  {
    Qualifier qualifier = feature.getQualifierByName(name);
    if(qualifier == null)
      return null;
    
    return (String)(qualifier.getValues().get(0));
  }
  
  /**
   * Get the srcfeature residue length
   * @return
   */
  public int getSeqlen()
  {
    return seqlen;
  }

  /**
   * Set the srcfeature residue length
   * @param seqlen
   */
  public void setSeqlen(int seqlen)
  {
    this.seqlen = seqlen;
  }

  public int getSrcfeature_id()
  {
    return srcfeature_id;
  }

  public void setSrcfeature_id(int srcfeature_id)
  {
    this.srcfeature_id = srcfeature_id;
  }

  public Hashtable<String, List<Feature>> getSplicedFeatures()
  {
    return splicedFeatures;
  }
  
  /**
   * Get the nucleotide location for a featureloc in amino acid
   * coordinates.
   * @param proteinFeature
   * @param featureLocToProtein
   * @return
   * @throws LocationParseException
   */
  public Location getNucLocation(final Feature proteinFeature, 
                                 final FeatureLoc featureLocToProtein) 
         throws LocationParseException
  {
    String transcriptName = getTranscriptFromName(
                         GeneUtils.getUniqueName(proteinFeature));
    List<Feature> spliced = getSplicedFeaturesOfTranscript(transcriptName);
    if(spliced == null)
      return null;
    
    RangeVector ranges = new RangeVector();
    for(int i=0; i<spliced.size(); i++)
    {
      Feature f = spliced.get(i);
      if(f.getKey().getKeyString().equals(DatabaseDocument.EXONMODEL))
        ranges.addAll(f.getLocation().getRanges());
    }
    
    int start = proteinFeature.getLocation().getFirstBase();
    int fmin = start+(featureLocToProtein.getFmin()*3)+1;
    int fmax = start+(featureLocToProtein.getFmax()*3);

    int len = proteinFeature.getEntry().getSequence().length();
    if(fmax > len)
      fmax = len;

    if(ranges.size()>1)
    {
      Collections.sort(ranges, new RangeComparator());

      for(int i=0;i<ranges.size()-1; i++)
      {
        Range range1 = (Range) ranges.get(i);
        Range range2 = (Range) ranges.get(i+1);
        if(fmin > range1.getEnd())
          fmin += range2.getStart()-range1.getEnd();
        if(fmax > range1.getEnd())
          fmax += range2.getStart()-range1.getEnd();
      }
    }

    Location location;
    if(proteinFeature.getLocation().isComplement())
      location = new Location("complement("+fmin+".."+fmax+")");
    else
      location = new Location(fmin+".."+fmax); 
    return location;
  }
  
  
  class RangeComparator implements Comparator<Range>
  {
    public int compare(Range o1, Range o2)
    {
      int start1 = o1.getStart();
      int start2 = o2.getStart();
      return start1-start2;
    }
  }
}