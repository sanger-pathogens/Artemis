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
 * Representation of the chado feature_dbxref, dbxref and db tables.
 */
public class ChadoDbxref
{
  private String name;
  private String accession;
  private int feature_id;
  /** database id */ 
  private int db_id;
  /** database cross reference id */
  private int dbxref_id;
  private boolean current = true;
  private String schema;

  public ChadoDbxref()
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

  public boolean isCurrent()
  {
    return current;
  }

  public void setCurrent(boolean current)
  {
    this.current = current;
  }

  public String getSchema()
  {
    return schema;
  }

  public void setSchema(String schema)
  {
    this.schema = schema;
  }

  public int getDb_id()
  {
    return db_id;
  }

  public void setDb_id(int db_id)
  {
    this.db_id = db_id;
  }

  public int getDbxref_id()
  {
    return dbxref_id;
  }

  public void setDbxref_id(int dbxref_id)
  {
    this.dbxref_id = dbxref_id;
  }
  
  
}

