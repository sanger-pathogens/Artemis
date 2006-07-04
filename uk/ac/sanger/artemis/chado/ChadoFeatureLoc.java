/* 
 *
 * created: 2006
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2006  Genome Research Limited
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
 */

package uk.ac.sanger.artemis.chado;

/**
 * Representation of featureloc table.
 */
public class ChadoFeatureLoc
{
  /** id of the parent feature */
  private int srcfeature_id = 0;
  /** start position */
  private int fmin;
  /** end position */
  private int fmax;
  /** +1 or -1 depending if on the forward or reverse */
  private int strand;
  /** phase */
  private int phase;
  /** Alternative residues, differing from feature.residues e.g. SNP.*/
  private String residue_info;
  /** used to manifest redundant, derivable extra locations for a feature*/
  private int locgroup;
  /** Used when a feature has >1 location */
  private int rank;

  /**
   * Get the strand. The orientation/directionality of the location. 
   * Should be 0, -1 or +1. Use 0 if transpliced.
   * @return the strand
   */
  public int getStrand()
  {
    return strand;
  }

  /**
   * Set the strand. The orientation/directionality of the location. 
   * Should be 0, -1 or +1. Use 0 if transpliced.
   * @param strand the strand
   */
  public void setStrand(int strand)
  {
    this.strand = strand;
  }

  /**
   * Get the fmin value. The leftmost/minimal boundary in the linear range 
   * represented by the featureloc. To convert this to the leftmost position 
   * in a base-oriented system (e.g. GFF) add 1 to fmin.
   * @return the fmin value
   */
  public int getFmin()
  {
    return fmin;
  }

  /**
   * Set the fmin value. The leftmost/minimal boundary in the linear range 
   * represented by the featureloc. To convert this to the leftmost position 
   * in a base-oriented system (e.g. GFF) add 1 to fmin.
   * @param fmin the fmin value
   *
   */
  public void setFmin(int fmin)
  {
    this.fmin = fmin;
  }

  /**
   * Get the fmax value. The rightmost/maximal boundary in the linear range 
   * represented by the featureloc. No conversion is required to go from fmax 
   * to the rightmost coordinate in a base-oriented system that counts from 1 
   * (eg GFF).
   * @return the fmax value
   */
  public int getFmax()
  {
    return fmax;
  }

  /**
   * Set the fmax value. The rightmost/maximal boundary in the linear range 
   * represented by the featureloc. No conversion is required to go from fmax 
   * to the rightmost coordinate in a base-oriented system that counts from 1 
   * (eg GFF).
   * @param fmax the fmax value
   */
  public void setFmax(int fmax)
  {
    this.fmax = fmax;
  }

  /**
   * The source feature which the location in featureloc is relative to.
   * @return    the source feature
   */
 public int getSrcfeature_id()
 {
   return srcfeature_id;
 }

  /**
   * The source feature which the location in featureloc is relative to.
   * @param srcfeature_id the source feature
   */
  public void setSrcfeature_id(int srcfeature_id)
  {
    this.srcfeature_id = srcfeature_id;
  }

  /**
   * The phase of translation wrt srcfeature_id. Values are 0,1,2. 
   * @return  the phase
   */
  public int getPhase()
  {
    return phase;  
  }

  /**
   * The phase of translation wrt srcfeature_id. Values are 0,1,2. 
   * @param phase the phase
   */
  public void setPhase(int phase)
  {
    this.phase = phase;
  }

  /**
   * This is used to manifest redundant, derivable extra locations for a feature. 
   * The default locgroup=0 is used for the DIRECT location of a feature. 
   * MOST CHADO USERS MAY NEVER USE featurelocs WITH logroup>0.
   * For example, the position of an exon on a BAC and in global chromosome coordinates. 
   * This column is used to differentiate these groupings of locations.
   * @return
   */
  public int getLocgroup()
  {
    return locgroup;
  }

  /**
   * Set the locgroup.
   * @param locgroup
   */
  public void setLocgroup(int locgroup)
  {
    this.locgroup = locgroup;
  }

  /**
   * Used when a feature has >1 location, otherwise the default rank 0 
   * is used. Some features (eg blast hits and HSPs) have two locations - 
   * one on the query and one on the subject.
   * @return
   */
  public int getRank()
  {
    return rank;
  }

  /**
   * Used when a feature has >1 location, otherwise the default rank 0 
   * is used. Some features (eg blast hits and HSPs) have two locations - 
   * one on the query and one on the subject.
   * @param rank
   */
  public void setRank(int rank)
  {
    this.rank = rank;
  }

  public String getResidue_info()
  {
    return residue_info;
  }

  public void setResidue_info(String residue_info)
  {
    this.residue_info = residue_info;
  }

}