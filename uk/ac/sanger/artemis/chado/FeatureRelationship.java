/* ChadoFeatureRelationship
 *
 * created: 2006
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2005  Genome Research Limited
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
 * Representation of feature_relationship table.
 */
public class FeatureRelationship
{
  
  public Feature subject;
  /** id of the parent feature */
  public Feature object;
  /** cv term of feature relationship (e.g. part_of) */
  public Cvterm cvterm;
  /** Additional notes/comments */
  public String value;
  /** ordering of subject features */
  public int rank;
  
  /**
   * Get the cv term for the relationship type between 
   * subject and object (e.g. part_of).
   * @return cvterm
   */
  public Cvterm getCvterm()
  {
    return cvterm;
  }
  
  /**
   * Set the cv term for the relationship type between 
   * subject and object (e.g. part_of).
   * @param cvterm
   */
  public void setCvterm(Cvterm cvterm)
  {
    this.cvterm = cvterm;
  }
  
  /**
   * Get the parent feature_id.
   * @return the parent feature_id
   */
  public Feature getObject()
  {
    return object;
  }
  
  /**
   * Set the parent feature_id.
   * @param object_id  the parent feature_id
   */
  public void setObject(Feature object)
  {
    this.object = object;
  }
  
  /**
   * The ordering of subject features with respect to the object
   * feature may be important (for example, exon ordering on a transcript
   *  - not always derivable if you take trans spliced genes into consideration). 
   *  rank is used to order these; starts from zero.
   * @return rank the rank
   */
  public int getRank()
  {
    return rank;
  }
  
  /**
   * The ordering of subject features with respect to the object
   * feature may be important (for example, exon ordering on a transcript
   *  - not always derivable if you take trans spliced genes into consideration). 
   *  rank is used to order these; starts from zero.
   * @param rank
   */
  public void setRank(int rank)
  {
    this.rank = rank;
  }
  
  public Feature getSubject()
  {
    return subject;
  }
  
  public void setSubject(Feature subject)
  {
    this.subject = subject;
  }
  
  /**
   * Additional notes/comments
   * @return
   */
  public String getValue()
  {
    return value;
  }
  
  /**
   * Additional notes/comments
   * @param value
   */
  public void setValue(String value)
  {
    this.value = value;
  }
  
}