/* ChadoCv.java
 *
 * created: 2005
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
 * Representation of cv table.
 */
public class ChadoCv
{
  private int cvId;
  private String name;
  private String definition;
  
  public int getCvId()
  {
    return cvId;
  }
  
  public void setCvId(int cvId)
  {
    this.cvId = cvId;
  }
  
  /**
   * A description of the criteria for membership of this ontology.
   * @return
   */
  public String getDefinition()
  {
    return definition;
  }
  
  /**
   * A description of the criteria for membership of this ontology.
   * @param definition
   */
  public void setDefinition(String definition)
  {
    this.definition = definition;
  }
  
  /**
   * The name of the ontology. 
   * @return
   */
  public String getName()
  {
    return name;
  }
  
  /**
   * The name of the ontology.
   * @param name
   */
  public void setName(String name)
  {
    this.name = name;
  }
}