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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/ChadoCanonicalGene.java,v 1.3 2006-07-04 11:06:33 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.chado.ChadoFeature;
import uk.ac.sanger.artemis.util.StringVector;
import java.util.Vector;
import java.util.Hashtable;
import java.util.Enumeration;
import java.util.List;

/**
 *  Used by GFFStreamFaeture to represent the chado canonical gene.
 **/
public class ChadoCanonicalGene  
{
  private Object gene;
  private String gene_id;
  
  // part_of gene
  private List transcripts = new Vector();
  
  // part_of transcrips
  private Hashtable exons = new Hashtable();
  
  // derives_from transript
  private Hashtable proteins = new Hashtable();

  // srcfeature length
  private int seqlen;
  
  public void addObject(Object obj, String parent_id, String type)
               throws InvalidRelationException
  {
    if(type.equalsIgnoreCase("polypeptide"))
      addProtein(parent_id, obj);
    else if(type.equalsIgnoreCase("exon"))
      addExon(parent_id, obj);
  }
  
  public Object getGene()
  {
    return gene;
  }

  public void setGene(Object gene)
  {
    this.gene = gene;
  }
  
  public void addTranscript(Object transcript)
  {
    transcripts.add(transcript);
  }
  
  public void addExon(String transcript_id, Object exon, boolean reset) 
         throws InvalidRelationException
  {
    exons = new Hashtable();
    addExon(transcript_id, exon);
  }
  
  public void addExon(String transcript_id, Object exon) 
         throws InvalidRelationException
  {   
    final List v_exons;
    if(exons.containsKey(transcript_id))
      v_exons = (Vector)exons.get(transcript_id);
    else
      v_exons = new Vector();
    
    v_exons.add(exon);
    exons.put(transcript_id, v_exons);
  }
  
  public void addProtein(String transcript_id, Object protein) 
         throws InvalidRelationException
  {   
    proteins.put(transcript_id, protein);
  }
  
  public Object containsTranscript(final StringVector ids)
  {
    for(int i=0; i<transcripts.size(); i++)
    {
      try
      {
        Feature transcript = (Feature)transcripts.get(i);
        
        if( ids.contains((String)transcript.getQualifierByName("ID").getValues().get(0)) )
          return transcript;
      }
      catch(InvalidRelationException e)
      {
        e.printStackTrace();
      }
    }
    return null;
  }

  public List getExonsOfTranscript(final String transcript_id)
  {
    if(exons.containsKey(transcript_id))
      return (List)exons.get(transcript_id);;
 
    return null;
  }
  
  public Object getProteinOfTranscript(final String transcript_id)
  {
    if(proteins.containsKey(transcript_id))
      return proteins.get(transcript_id);;
 
    return null;
  }
  
  public Hashtable getExons()
  {
    return exons;
  }

  public List getTranscripts()
  {
    return transcripts;
  }
  
  public Object getFeatureFromId(final String id)
  {
    Object feature = null;
    
    // check gene
    try
    {
      if(((String)(((Feature)gene).getQualifierByName("ID").getValues().get(0))).equals(id))
        return gene;
    }
    catch(InvalidRelationException e)
    {
      e.printStackTrace();
    }
    
    // check transcript
    StringVector sv = new StringVector();
    sv.add(id);
    
    feature = containsTranscript(sv);
    
    if(feature != null)
      return feature;

    // check exons
    feature = getExon(id);
    
    if(feature != null)
      return feature;
    
    return getProtein(id);
  }
  
  private Object getExon(final String id)
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
            uniquename = (String)((Feature)exons.get(i)).getQualifierByName("ID").getValues().get(0);
          
          if(uniquename.equals(id))
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
  
  private Feature getProtein(final String id)
  {
    Enumeration enum_proteins = proteins.elements();
    try
    {
      while(enum_proteins.hasMoreElements())
      {
        Feature protein = (Feature)enum_proteins.nextElement();
        if(((String)(protein.getQualifierByName("ID").getValues().get(0))).equals(id))
          return protein;
      }
    }
    catch(InvalidRelationException e)
    {
      e.printStackTrace();
    }
    return null;
    
  }

  public int getSeqlen()
  {
    return seqlen;
  }

  public void setSeqlen(int seqlen)
  {
    this.seqlen = seqlen;
  }
}