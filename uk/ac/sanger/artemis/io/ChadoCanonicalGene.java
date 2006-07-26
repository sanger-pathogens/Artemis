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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/ChadoCanonicalGene.java,v 1.8 2006-07-26 13:50:16 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.chado.ChadoFeature;
import uk.ac.sanger.artemis.util.StringVector;
import java.util.Vector;
import java.util.Hashtable;
import java.util.Enumeration;
import java.util.List;

/**
 *  Used by GFFStreamFeature to represent the chado canonical gene.
 *  Contains gene, transcript, exons and proteins.
 **/
public class ChadoCanonicalGene  
{
  private Feature gene;
  
  // part_of gene
  private List transcripts = new Vector();
  
  // part_of transcrips
  private Hashtable exons = new Hashtable();
  
  // derives_from transript
  private Hashtable proteins = new Hashtable();

  private Hashtable three_prime_utr = new Hashtable();
  private Hashtable five_prime_utr  = new Hashtable();
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
  
  public void addTranscript(Object transcript)
  {
    transcripts.add(transcript);
  }
  
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
        }
      }
      catch(InvalidRelationException e)
      {
        e.printStackTrace();
      }
    }  
  }
  
  public void deleteFeature(Feature embl_feature)
  {
    try
    {
      String name = getQualifier(embl_feature, "ID");
      Object feature = getExon(name);

      if(feature != null)
      {
        String trenascript_name = getQualifier((Feature) feature, "Parent");
        exons.remove(trenascript_name);

        return;
      }
      deleteTranscript(name);
    }
    catch(InvalidRelationException e1)
    {
      e1.printStackTrace();
    }
  }
  
  public void addExon(String transcript_name, 
                      Object exon, boolean reset) 
         throws InvalidRelationException
  {
    if(reset)
      exons = new Hashtable();
    addExon(transcript_name, exon);
  }
  
  public void addExon(String transcript_name, Object exon) 
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
  
  public void addProtein(String transcript_name, Object protein) 
         throws InvalidRelationException
  {   
    proteins.put(transcript_name, protein);
  }
  
  public void add3PrimeUtr(String transcript_name, Object protein) 
         throws InvalidRelationException
  {   
    three_prime_utr.put(transcript_name, protein);
  }
  
  public void add5PrimeUtr(String transcript_name, Object protein) 
         throws InvalidRelationException
  {
    five_prime_utr.put(transcript_name, protein);
  }
  
  public Object containsTranscript(final StringVector names)
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

  public List getExonsOfTranscript(final String transcript_name)
  {
    if(exons.containsKey(transcript_name))
      return (List)exons.get(transcript_name);
 
    return null;
  }
  
  public Object getProteinOfTranscript(final String transcript_name)
  {
    if(proteins.containsKey(transcript_name))
      return proteins.get(transcript_name);;
 
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
  
  public boolean isTranscript(final String name)
  {
    try
    {
      for(int i=0; i<transcripts.size(); i++)
      {
        if(name.equals(getQualifier((Feature)transcripts.get(i), "ID")))
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
      feature = getUTR(name, three_prime_utr);
      if(feature != null)
        return feature;
      
      feature = getUTR(name, five_prime_utr);
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
  private Object getExon(final String name)
  {
    Enumeration enum_exons = exons.elements();
    try
    {
      while(enum_exons.hasMoreElements())
      {
        Vector exons = (Vector)enum_exons.nextElement();
        
        for(int i=0; i<exons.size(); i++)
        {
          String uniquename;
          
          if(exons.get(i) instanceof ChadoFeature)
            uniquename = ((ChadoFeature)exons.get(i)).getUniquename();
          else
            uniquename = getQualifier((Feature)exons.get(i), "ID");
          
          if(uniquename.equals(name))
            return exons.get(i);
        }
      }
    }
    catch(InvalidRelationException e)
    {
      e.printStackTrace();
    }
    return null;
  }
  
  private Object getProtein(final String id)
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

  private Object getUTR(final String id,
                        final Hashtable UTR) 
          throws InvalidRelationException
  {
    Enumeration enum_utr = UTR.elements();

    while(enum_utr.hasMoreElements())
    {
      Feature utr = (Feature)enum_utr.nextElement();
      if(getQualifier(utr, "ID").equals(id))
        return utr;
    }

    return null;
  }
  
  
  private String getQualifier(Feature feature, String name) 
          throws InvalidRelationException
  {
    return (String)(feature.getQualifierByName(name).getValues().get(0));
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