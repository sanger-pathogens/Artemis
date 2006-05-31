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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/ChadoCanonicalGene.java,v 1.2 2006-05-31 15:41:17 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.StringVector;
import java.util.Vector;
import java.util.Hashtable;
import java.util.Enumeration;

/**
 *  Used by GFFStreamFaeture to represent the chado canonical gene.
 **/
public class ChadoCanonicalGene  
{
  private Feature gene;
  
  // part_of gene
  private Vector transcripts = new Vector();
  
  // part_of transcrips
  private Hashtable exons = new Hashtable();
  
  // derives_from transript
  private Hashtable proteins = new Hashtable();

  public Feature getGene()
  {
    return gene;
  }

  public void setGene(Feature gene)
  {
    this.gene = gene;
  }
  
  public void addTranscript(Feature transcript)
  {
    transcripts.add(transcript);
  }
  
  public void addExon(Feature transcript, Feature exon, boolean reset) 
         throws InvalidRelationException
  {
    exons = new Hashtable();
    addExon(transcript, exon);
  }
  
  public void addExon(Feature transcript, Feature exon) 
         throws InvalidRelationException
  {
    final String transcript_id = 
         (String)transcript.getQualifierByName("ID").getValues().get(0);
    
    final Vector v_exons;
    if(exons.containsKey(transcript_id))
      v_exons = (Vector)exons.get(transcript_id);
    else
      v_exons = new Vector();
    
    v_exons.add(exon);
    exons.put(transcript_id, v_exons);
  }
  
  public void addProtein(Feature transcript, Feature protein) 
         throws InvalidRelationException
  {
    final String transcript_id = 
      (String)transcript.getQualifierByName("ID").getValues().get(0);
    
    proteins.put(transcript_id, protein);
  }
  
  public Feature containsTranscript(final StringVector ids)
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

  public Vector getExonsOfTranscript(final String transcript_id)
  {
    if(exons.containsKey(transcript_id))
      return (Vector)exons.get(transcript_id);;
 
    return null;
  }
  
  public Feature getProteinOfTranscript(final String transcript_id)
  {
    if(proteins.containsKey(transcript_id))
      return (Feature)proteins.get(transcript_id);;
 
    return null;
  }
  
  public Hashtable getExons()
  {
    return exons;
  }

  public Vector getTranscripts()
  {
    return transcripts;
  }
  
  public Feature getFeatureFromId(final String id)
  {
    Feature feature = null;
    
    // check gene
    try
    {
      if(((String)(gene.getQualifierByName("ID").getValues().get(0))).equals(id))
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
  
  private Feature getExon(final String id)
  {
    Enumeration enum_exons = exons.elements();
    try
    {
      while(enum_exons.hasMoreElements())
      {
        Vector exons = (Vector)enum_exons.nextElement();
        
        for(int i=0; i<exons.size(); i++)
        {
          Feature exon = (Feature)exons.get(i);
          if(((String)(exon.getQualifierByName("ID").getValues().get(0))).equals(id))
            return exon;
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
}