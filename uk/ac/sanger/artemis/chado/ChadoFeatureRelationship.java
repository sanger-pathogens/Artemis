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

public class ChadoFeatureRelationship
{
  
  public int subject_id;
  /** id of the parent feature */
  public int object_id;
  /** cv term of feature relationship (e.g. part_of) */
  public Cvterm cvterm;
  /** Additional notes/comments */
  public String value;
  /** ordering of subject features */
  public int rank;
  
  public Cvterm getCvterm()
  {
    return cvterm;
  }
  
  public void setCvterm(Cvterm cvterm)
  {
    this.cvterm = cvterm;
  }
  
  /**
   * Get the parent feature_id.
   * @return the parent feature_id
   */
  public int getObject_id()
  {
    return object_id;
  }
  
  /**
   * Set the parent feature_id.
   * @param object_id  the parent feature_id
   */
  public void setObject_id(int object_id)
  {
    this.object_id = object_id;
  }
  
  public int getRank()
  {
    return rank;
  }
  
  public void setRank(int rank)
  {
    this.rank = rank;
  }
  
  public int getSubject_id()
  {
    return subject_id;
  }
  
  public void setSubject_id(int subject_id)
  {
    this.subject_id = subject_id;
  }
  
  public String getValue()
  {
    return value;
  }
  
  public void setValue(String value)
  {
    this.value = value;
  }
  
}