/* CvTerm.java
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
 * Representation of cvterm table.
 */
public class CvTerm
{
  private long cvTermId;
  private String name;
  private Cv cv;
  private DbXRef dbXRef;
  private String definition;
  private int isObsolete;
  private int isRelationshiptype;
  
  public long getCvTermId()
  {
    return cvTermId;
  }

  public void setCvTermId(long cvTermId)
  {
    this.cvTermId = cvTermId;
  }

  /**
   * Primary dbxref - The unique global OBO identifier for this cvterm. 
   * Note that a cvterm may have multiple secondary dbxrefs.
   * @return
   */
  public DbXRef getDbXRef()
  {
    return dbXRef;
  }

  /**
   * Primary dbxref - The unique global OBO identifier for this cvterm. 
   * Note that a cvterm may have multiple secondary dbxrefs.
   * @param dbxref
   */
  public void setDbXRef(DbXRef dbXRef)
  {
    this.dbXRef = dbXRef;
  }

  public Cv getCv()
  {
    return cv;
  }

  public void setCv(Cv cv)
  {
    this.cv = cv;
  }

  /**
   * A human-readable text definition.
   * @return
   */
  public String getDefinition()
  {
    return definition;
  }

  /**
   * A human-readable text definition.
   * @param definition
   */
  public void setDefinition(String definition)
  {
    this.definition = definition;
  }

  public int getIsObsolete()
  {
    return isObsolete;
  }

  public void setIsObsolete(int isObsolete)
  {
    this.isObsolete = isObsolete;
  }

  /**
   * Use this flag to indicate whether this cvterm is an actual 
   * term/concept or a relationship type.
   * @return
   */
  public int getIsRelationshiptype()
  {
    return isRelationshiptype;
  }

  /**
   * Use this flag to indicate whether this cvterm is an actual 
   * term/concept or a relationship type.
   * @param isRelationshiptype
   */
  public void setIsRelationshiptype(int isRelationshiptype)
  {
    this.isRelationshiptype = isRelationshiptype;
  }

  /**
   * A concise human-readable name describing the meaning of the cvterm
   * @return
   */
  public String getName()
  {
    return name;
  }

  /**
   * A concise human-readable name describing the meaning of the cvterm
   * @param name
   */
  public void setName(String name)
  {
    this.name = name;
  }

}
