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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/ChadoCanonicalGene.java,v 1.16 2006-08-17 13:20:57 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.StringVector;
import java.util.Vector;
import java.util.Hashtable;
import java.util.Enumeration;
import java.util.List;
import java.util.Set;
import java.util.HashSet;

/**
 *  Used by GFFStreamFeature to represent the chado canonical gene.
 *  Contains gene, transcript, exons and proteins.
 **/
public class ChadoCanonicalGene  
{
  private Feature gene;
  
  // part_of gene
  private List transcripts = new Vector();
  
  // part_of transcripts
  private Hashtable exons = new Hashtable();
  
  // derives_from transript
  private Hashtable proteins = new Hashtable();

  // utr features
  private Hashtable three_prime_utr = new Hashtable();
  private Hashtable five_prime_utr  = new Hashtable();
  
  // other child features of transcript
  private Hashtable other_features = new Hashtable();
  
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
          exons.remove(transcript_name);
          three_prime_utr.remove(transcript_name);
          five_prime_utr.remove(transcript_name);
          other_features.remove(transcript_name);
        }
      }
      catch(InvalidRelationException e)
      {
        e.printStackTrace();
      }
    }  
  }
  
  /**
   * Delete features.
   * @param embl_feature
   */
  public void deleteFeature(Feature embl_feature)
  {
    try
    {
      String name = getQualifier(embl_feature, "ID");
      Object feature = getExon(name);

      if(feature != null)
      {
        String transcript_name = getQualifier((Feature) feature, "Parent");
        exons.remove(transcript_name);
        return;
      }
      
      feature = getFeatureFromHash(name, three_prime_utr);
      if(feature != null)
      {
        String transcript_name = getQualifier((Feature) feature, "Parent");
        List utr = get3UtrOfTranscript(transcript_name);
        utr.remove(feature);
        return;
      }
      
      feature = getFeatureFromHash(name, five_prime_utr);
      if(feature != null)
      {
        String transcript_name = getQualifier((Feature) feature, "Parent");
        List utr = get5UtrOfTranscript(transcript_name);
        utr.remove(feature);
        return;
      }
      
      feature = getFeatureFromHash(name, other_features);
      if(feature != null)
      {
        String transcript_name = getQualifier((Feature) feature, "Parent");
        List others = getOtherFeaturesOfTranscript(transcript_name);
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
  public Set getChildren(Feature embl_feature)
  {
    Set children = new HashSet();
    try
    {
      String name = getQualifier(embl_feature, "ID");
      
      String gene_name = getQualifier(getGene(), "ID");
      if(name.equals(gene_name))
      {
        List transcripts = getTranscripts();
        for(int i=0; i<transcripts.size(); i++)
        {
          Feature transcript = (Feature)transcripts.get(i);
          children.add(transcript);
          children.addAll( getChildren(transcript) );
        }
        return children;
      }
      
      searchForChildren(exons, name, children);
      searchForChildren(three_prime_utr, name, children);
      searchForChildren(five_prime_utr, name, children);
      searchForChildren(other_features, name, children);
      
      // protein
      Enumeration pep_enum = proteins.elements();
      while(pep_enum.hasMoreElements())
      {
        Feature child = (Feature)pep_enum.nextElement();
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
  private void searchForChildren(Hashtable hash, String parent_id,
                                 Set children)
               throws InvalidRelationException
  {
    Enumeration feature_enum = hash.elements();
    String parent;
    
    while(feature_enum.hasMoreElements())
    {
      List child_list = (List)feature_enum.nextElement();
      
      for(int i=0; i<child_list.size(); i++)
      {       
        Feature child = (Feature)child_list.get(i);
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
  public void addExon(final String transcript_name, 
                      final Feature exon, boolean reset) 
         throws InvalidRelationException
  {
    if(reset)
      exons = new Hashtable();
    addExon(transcript_name, exon);
  }
  
  /**
   * Add exon feature to the chado gene model.
   * @param transcript_name
   * @param exon
   * @throws InvalidRelationException
   */
  public void addExon(final String transcript_name, 
                      final Feature exon) 
         throws InvalidRelationException
  {   
    final List v_exons;
    if(exons.containsKey(transcript_name))
      v_exons = (Vector)exons.get(transcript_name);
    else
      v_exons = new Vector();
    
    v_exons.add(exon);
    exons.put(transcript_name, v_exons);
  }
  
  /**
   * Add protein feature to the chado gene model.
   * @param transcript_name
   * @param protein
   * @throws InvalidRelationException
   */
  public void addProtein(final String transcript_name,
                         final Feature protein) 
         throws InvalidRelationException
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
         throws InvalidRelationException
  {  
    final List utr_list;
    if(three_prime_utr.containsKey(transcript_name))
      utr_list = (Vector)three_prime_utr.get(transcript_name);
    else
      utr_list = new Vector();
    
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
         throws InvalidRelationException
  {
    final List utr_list;
    if(five_prime_utr.containsKey(transcript_name))
      utr_list = (Vector)five_prime_utr.get(transcript_name);
    else
      utr_list = new Vector();
    
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
    final List v_other_features;
    if(other_features.containsKey(transcript_name))
      v_other_features = (Vector)other_features.get(transcript_name);
    else
      v_other_features = new Vector();
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
    for(int i=0; i<transcripts.size(); i++)
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

  /**
   * Return the exons of a given transcript as a <code>List</code>.
   * @param transcript_name
   * @return
   */
  public List getExonsOfTranscript(final String transcript_name)
  {
    if(exons.containsKey(transcript_name))
      return (List)exons.get(transcript_name);
 
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
  public List get3UtrOfTranscript(final String transcript_name)
  {
    if(three_prime_utr.containsKey(transcript_name))
      return (List)three_prime_utr.get(transcript_name);
 
    return null;
  }
  
  /**
   * Return the 5'UTR features of a transcriot as a <code>List</code>.
   * @param transcript_name
   * @return
   */
  public List get5UtrOfTranscript(final String transcript_name)
  {
    if(five_prime_utr.containsKey(transcript_name))
      return (List)five_prime_utr.get(transcript_name);
 
    return null;
  }
  
  /**
   * Return the other child features of a transcriot as a <code>List</code>.
   * @param transcript_name
   * @return
   */
  public List getOtherFeaturesOfTranscript(final String transcript_name)
  {
    if(other_features.containsKey(transcript_name))
      return (List)other_features.get(transcript_name);
    return null;
  }
  
  /**
   * Get a list of trancripts.
   * @return
   */
  public List getTranscripts()
  {
    return transcripts;
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
  private boolean isExon(final String feature_id, 
                         final String transcript_id)
  {
    List exons = getExonsOfTranscript(transcript_id);
    if(exons == null)
      return false;
    try
    {
      for(int i=0; i<exons.size(); i++)
      {
        if(feature_id.equals(getQualifier((Feature)exons.get(i), "ID")))
          return true;
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
      while( isTranscript( name + ":" + transcript_key + ":" + auto ) &&
             auto < 50)
        auto++;
      return name + ":" + transcript_key + ":" + auto;
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
  public String autoGenerateExonName(final String transcript_id)
  {
    try
    {
      //List exons = getExonsOfTranscript(transcript_id);
      int transcript_number = 1;
      for(transcript_number=1; transcript_number<=transcripts.size(); 
          transcript_number++)
      {
        Feature transcript = (Feature)transcripts.get(transcript_number-1);
        if(transcript_id.equals( getQualifier(transcript, "ID") ))
          break;
      }
      
      String name = (String)getGene().getQualifierByName("ID").getValues().get(0);
      int auto = 1;
      while( isExon( name + ":" + transcript_number + ":exon:" + auto, transcript_id ) && auto < 50)
        auto++;
      return name + ":" + transcript_number + ":exon:" + auto;
    }
    catch(InvalidRelationException e)
    {
      e.printStackTrace();
    }
    return null;
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
    feature = getExon(name);
    
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
  private Feature getExon(final String name)
  {
    Enumeration enum_exons = exons.elements();
    try
    {
      while(enum_exons.hasMoreElements())
      {
        Vector exons = (Vector)enum_exons.nextElement();
        
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
    Enumeration enum_proteins = proteins.elements();
    try
    {
      while(enum_proteins.hasMoreElements())
      {
        Feature protein = (Feature)enum_proteins.nextElement();
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
                         final Hashtable UTR) 
          throws InvalidRelationException
  {
    Enumeration enum_utr = UTR.elements();

    while(enum_utr.hasMoreElements())
    {
      List utrs = (List)enum_utr.nextElement();
      
      for(int i=0; i<utrs.size(); i++)
      {
        Feature utr = (Feature)utrs.get(i);
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
  private String getQualifier(final Feature feature,
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
}