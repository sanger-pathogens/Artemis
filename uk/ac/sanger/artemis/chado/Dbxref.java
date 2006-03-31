/* Dbxref.java
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
 *
 * Java Object used by ORM to represent the chado dbxref table.
 *
 */
public class Dbxref
{
  private String name;
  private String accession;
  private int feature_id;

  public Dbxref()
  {
  }
  
  public String getName()
  {
    return name; 
  }
  
  public void setName(final String name)
  {
    this.name = name;
  }

  public String getAccession()
  {
    return accession;
  }

  public void setAccession(String accession)
  {
    this.accession = accession;
  }

  public int getFeature_id()
  {
    return feature_id;
  }

  public void setFeature_id(int feature_id)
  {
    this.feature_id = feature_id;
  }

}

