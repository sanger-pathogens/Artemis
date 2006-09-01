/* ChadoDbxref.java
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
 * Representation of the chado dbxref table.
 */
public class DbXRef
{ 
  private int dbXRefId;
  private String accession;
  private String version;
  private String description;
  private Db db;
  
  public String getAccession()
  {
    return accession;
  }

  public void setAccession(String accession)
  {
    this.accession = accession;
  }

  public Db getDb()
  {
    return db;
  }

  public void setDb(Db db)
  {
    this.db = db;
  }

  public int getDbXRefId()
  {
    return dbXRefId;
  }

  public void setDbXRefId(int dbXRefId)
  {
    this.dbXRefId = dbXRefId;
  }

  public String getDescription()
  {
    return description;
  }

  public void setDescription(String description)
  {
    this.description = description;
  }

  public String getVersion()
  {
    return version;
  }

  public void setVersion(String version)
  {
    this.version = version;
  }  
}