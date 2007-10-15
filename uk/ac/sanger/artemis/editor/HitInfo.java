/*
 *
 * created: Wed Aug 3 2004
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2000  Genome Research Limited
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
 */

package uk.ac.sanger.artemis.editor;

import java.util.StringTokenizer;
import java.util.NoSuchElementException;
import java.util.Vector;
import java.util.Hashtable;

/**
*
* Hit information is contained in this object.
*
*/
public class HitInfo
{
  /** Hit in database */
  private String db     = null;
  /** Hit ID */
  private String id     = null;
  /** Hit accession number */
  private String acc    = null;
  /** Organism source */
  private String org    = null;
  /** Gene name */
  private String geneName = null;
  /** Description     */
  private String desc   = null;
  /** Sequence length */
  private String length = null;

  /** */
  //private String opt    = null;
  //private String zscore = null;
  private String evalue = null;
  private String header = null;
  /** EMBL ID from linking to EMBL in SRS */
  private String emblID = null;
  private String score  = null;
  /** percentage identity */
  private String identity  = null;
  /** ungapped identity   */
  private String ungapped  = null;
  /** amino acid overlap  */
  private String aaOverlap = null;
  /** query range in the alignment   */
  private String queryRange   = null;
  /** subject range in the alignment */
  private String subjectRange = null;
  /** start position of this hit in the results output */
  private int startPosition = 0;
  /** end position of this hit in the results output   */
  private int endPosition   = 0;
  /** collection of GO terms */
  private Vector go;
  /** EC_number */
  private String ec_number = null;

  /** query hit range */
  private Vector queryPosition = new Vector();

  /** hastable for go terms */
  private Hashtable go_hash;

  public HitInfo(String header, String format)
  {
    header = header.trim();
    this.header = header;
    if(format.equals("fasta"))
      setFastaHitInfo(header);
    else if(format.equals("blastp"))
      setBLASTPInfo(header);
  }


  /**
  *
  * Extract id, accession numbers etc. from the
  * header of BLASTP results.
  * @param header	results header line from BLASTP.
  *
  */
  protected void setBLASTPInfo(String header)
  {
    int ind1 = header.indexOf(" ");
    if(ind1 > -1)
    {
      id = header.substring(0,ind1);

      int indDot;
      
      if( (indDot = id.indexOf(".")) > 5)
      {
        //
        // uniprot fasta headers
        // >pac.sv id desc
        //
        String version = id.substring(indDot+1);
        try
        {
          int version_no = Integer.parseInt(version);
          // version number looks like uniprot so strip this
          if(version_no < 50)
            id = id.substring(0,indDot);
        }
        catch(NumberFormatException nfe){}

        acc = id;
      }
    }
    else
      return;

    int ind1a = id.indexOf(":");
    if(ind1a > -1)
    {
      db = id.substring(0,ind1a);

      if(db.equals("UNIPROT"))
        db = "UniProt";

      id = id.substring(ind1a+1);
    }

    int ind2 = header.indexOf(" ",ind1+1);
    if(ind2 > -1)
    {
      if(acc == null)
        acc = header.substring(ind1+1,ind2);
//    else
//      id = header.substring(ind1+1,ind2);
    }
    else
      return;

    if(acc.startsWith("|")) // genedb
      acc = id;

    ind1 = header.lastIndexOf(" ");
    setEValue(header.substring(ind1).trim());
  
    ind2  = header.substring(0,ind1).trim().lastIndexOf(" ");
    if(ind2>-1)
      score = header.substring(ind2,ind1).trim();
  }

  /**
  *
  * Extract id, accession numbers etc. from the
  * header of FASTA results.
  * @param header       results header line from FASTA.
  *
  */
  protected void setFastaHitInfo(String header)
  {
    int ind1 = header.indexOf(" ");
    if(ind1 > -1)
    {
      id = header.substring(0,ind1);
      int indDot;
      if( (indDot = id.indexOf(".")) > 5)
      {
        //
        // uniprot fasta headers
        // >pac.sv id desc
        //
        String version = id.substring(indDot+1);
        try
        {
          int version_no = Integer.parseInt(version);
          // version number looks like uniprot so strip this
          if(version_no < 50)
            id = id.substring(0,indDot);
        }
        catch(NumberFormatException nfe){}
        acc = id;
      }
    }
    else
      return;

    int ind2 = header.indexOf(" ",ind1+1);
    if(ind2 > -1)
    {
     if(acc == null)
        acc = header.substring(ind1+1,ind2);
//   else
//      id = header.substring(ind1+1,ind2);
    }
    else
      return;

    ind1 = ind2;
    ind2 = header.lastIndexOf("(");

    if(ind2 > -1)
      desc = "";
//    desc = header.substring(ind1,ind2).trim();
    else
      return;

    ind1 = ind2+1;
    ind2 = header.indexOf(")",ind1);
    if(ind2 > -1)
      length = header.substring(ind1,ind2).trim();
    else
      return;

    
    StringTokenizer tok = new StringTokenizer(header.substring(ind2+1));
    try
    {
      tok.nextToken();  // opt
      tok.nextToken();  // zscore
      setEValue(tok.nextToken());
    }
    catch(NoSuchElementException exp){}
    
  }

  
  /**
  *
  * Set EC_number
  *
  */
  protected void setEC_number(String ec_number)
  {
    this.ec_number = ec_number;
  }

 
  /**
  *
  * Get EC_number
  *
  */
  protected String getEC_number()
  {
    return ec_number;
  }


  /**
  *
  * Set the start position for the query sequence in the alignment.
  * @param startQuery 	startQuery
  *
  */
  protected void setQueryPosition(int startQuery, int endQuery)
  {
    queryPosition.add(new Integer(startQuery));
    queryPosition.add(new Integer(endQuery));
  }


  /**
  *
  * Get the start position for the query sequence in the alignment.
  * @param startQuery   startQuery
  *
  */
  protected Vector getQueryPosition()
  {
    return queryPosition;
  }

  /**
  *
  * Set the start position of the alignment in the
  * results.
  * @param startPosition	start position.
  *
  */
  protected void setStartPosition(int startPosition)
  {
    this.startPosition = startPosition;
  }

  /**
  *
  * Set the end position of the alignment in the
  * results.
  * @param endPosition		end position.
  *
  */
  protected void setEndPosition(int endPosition)
  {
    this.endPosition = endPosition;
  }


  /**
  *
  * Set the sequence length.
  * @param length	sequence length.
  *
  */
  protected void setLength(String length)
  {
    this.length = length;
  }

  /**
  *
  * Get the sequence length.
  * @return sequence length.
  *
  */
  protected String getLength()
  {
    return length;
  }

  /**
  *
  * Get the results header.
  * @return results header.
  *
  */
  protected String getHeader()
  {
    return header; 
  }

  /**
  *
  * Get the alignment start position.
  * @return start position.
  *
  */
  protected int getStartPosition()
  {
   return startPosition;
  }

  /**
  *
  * Get the alignment end position.
  * @return end position.
  *
  */
  protected int getEndPosition()
  {
   return endPosition;
  }

  /**
  *
  * Get the database.
  * @return database.
  *
  */
  protected String getDB()
  {
    return db;
  }


  /**
  *
  * Get the sequence ID.
  * @return sequence ID.
  *
  */
  protected String getID()
  {
    return id;
  }
  
  /**
  *
  * Get the sequence accession number.
  * @return sequence accession number.
  *
  */
  protected String getAcc()
  {
    return acc;
  }

  /**
  *
  * Get the organism source (OS).
  * @return organism source.
  *
  */
  protected String getOrganism()
  {
    return org;
  }

  /**
  *
  * Set the organism source (OS).
  * @param org 	organism source.
  *
  */
  protected void setOrganism(final String orgStr)
  {
    if(org == null || org.equals(""))
      this.org = orgStr;
    else 
      this.org = this.org + " " + orgStr; 
  }

  /**
  *
  * Set the description (DE).
  * @param desc	sequence description.
  *
  */
  protected void setDescription(String desc)
  {
    this.desc = desc;
  }


  /**
  *
  * Append to the description.
  * @param description to append.
  *
  */
  protected void appendDescription(String s)
  {
    if(desc == null || desc.equals(""))
      desc = new String(s);
    else 
      desc = desc + " " + s;

    desc = desc.trim();

    // check for EC number
    int ind1 = s.indexOf("(EC ");
    if(ind1 > -1)
    {
      int ind2 = s.indexOf(")",ind1);
      if(ind2 > -1)
        setEC_number(s.substring(ind1+4,ind2).trim());
    }
  }


  /**
  *
  * Get the description (DE).
  * @return sequence description.
  *
  */
  protected String getDescription()
  {
    return desc;
  }

  /**
  *
  * Set the EMBL linked entry (obtained from database 
  * linking in SRS).
  * @param emblID	EMBL ID.
  *
  */
  protected void setEMBL(String emblID)
  {
    if(emblID.endsWith(";"))
      emblID = emblID.substring(0, emblID.length()-1);
    this.emblID = emblID;
  }

  /**
  *
  * Get the EMBL linked entry (obtained from database 
  * linking in SRS).
  * @return EMBL ID.
  *
  */
  protected String getEMBL()
  {
    return emblID;
  }

  /**
  *
  * Set the sequence alignment score.
  * @param score 	alignment score.
  *
  */
  protected void setScore(String score)
  {
    this.score = score;
  }

  /**
  *
  * Get the sequence alignment score.
  * @return alignment score.
  *
  */
  protected String getScore()
  {
    return score;
  }

  /**
  *
  * Set the percentage identity found in the alignment.
  * @param identity	percentage identity.
  *
  */
  protected void setIdentity(String identity)
  {
    this.identity = identity;
  }

  /**
  *
  * Get the percentage identity found in the alignment.
  * @return percentage identity.
  *
  */
  protected String getIdentity()
  {
    return identity;
  }

  /**
  *
  * Set the percentage ungapped identity found in the alignment.
  * @param ungapped     percentage ungapped identity.
  *
  */
  protected void setUngapped(String ungapped)
  {
    this.ungapped = ungapped;
  }

  /**
  *
  * Get the percentage ungapped identity found in the alignment.
  * @return percentage ungapped identity.
  *
  */
  protected String getUngapped()
  {
    return ungapped;
  }

  /**
  *
  * Set the e-value for this hit.
  * @param evalue  e-value.
  *
  */
  protected void setEValue(String evalue)
  {
    if(evalue.startsWith("e"))
      evalue = "0"+evalue;
    this.evalue = evalue;
  }

  /**
  *
  * Get the e-value for this hit.
  * @return e-value.
  *
  */
  protected String getEValue()
  {
    return evalue;
  }

  /**
  *
  * Set the number of amino acids that overlap in this
  * alignment.
  * @param aaOverlap  number of aa-overlaps.
  *
  */
  protected void setOverlap(String aaOverlap)
  {
    this.aaOverlap = aaOverlap;
  }

  /**
  *
  * Get the number of amino acids that overlap in this
  * alignment.
  * @return number of aa-overlaps.
  *
  */
  protected String getOverlap()
  {
    return aaOverlap;
  }

  /**
  *
  * Set the query range in this alignment.
  * @param queryRange 	query range.
  *
  */
  protected void setQueryRange(String queryRange)
  {
    this.queryRange = queryRange;
  }

  /**
  *
  * Get the query range in this alignment.
  * @return query range.
  *
  */
  protected String getQueryRange()
  {
    return queryRange;
  }

  /**
  *
  * Set the subject range in this alignment.
  * @param subjectRange   subject range.
  *
  */
  protected void setSubjectRange(String subjectRange)
  {
    this.subjectRange = subjectRange;
  }

  /**
  *
  * Get the subject range in this alignment.
  * @return subject range.
  *
  */
  protected String getSubjectRange()
  {
    return subjectRange;
  }

  /**
  *
  * Set the gene name (for ortholog).
  * @param geneName	gene name.
  *
  */
  protected void setGeneName(String geneName)
  {
    this.geneName = geneName;
  }

  /**
  *
  * Get the gene name (for ortholog).
  * @return gene name.
  *
  */
  protected String getGeneName()
  {
    return geneName;
  }


  /**
  *
  * Set the GO term
  *
  */
  protected void setGO(String go_terms)
  {
    go_terms = go_terms.trim();
    StringTokenizer tok = new StringTokenizer(go_terms);
    while(tok.hasMoreElements())
    {
      if(go == null)
        go = new Vector();
      go.add(tok.nextToken());
    }
  }

  
  /**
  *
  * Get GO collection.
  *
  */
  protected Vector getGO()
  {
    return go;
  }


  /**
  *
  * Set a GO association.
  *
  */
  protected void setGoAssociation(String go_id, String assoc)
  {
    if(go_hash == null)
      go_hash = new Hashtable();
    go_hash.put(go_id,assoc);
  }

  /**
  *
  * Get a GO association.
  *
  */
  protected String getGoAssociation(String go_id)
  {
    if(go_hash != null && go_hash.containsKey(go_id))
      return (String)go_hash.get(go_id);
    return null;
  }

}

