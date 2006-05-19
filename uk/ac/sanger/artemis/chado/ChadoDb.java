/* ChadoDb.java
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
 * Representation of the chado db table.
 */
public class ChadoDb
{ 
  private int db_id;
  private String name;
  private String description;
  private String urlprefix;
  private String url;
  
  public String getDescription()
  {
    return description;
  }
  
  public void setDescription(String description)
  {
    this.description = description;
  }
  
  public int getDb_id()
  {
    return db_id;
  }

  public void setDb_id(int db_id)
  {
    this.db_id = db_id;
  }
  
  public String getName()
  {
    return name;
  }
  
  public void setName(String name)
  {
    this.name = name;
  }
  
  public String getUrl()
  {
    return url;
  }
  
  public void setUrl(String url)
  {
    this.url = url;
  }
  
  public String getUrlprefix()
  {
    return urlprefix;
  }
  
  public void setUrlprefix(String urlprefix)
  {
    this.urlprefix = urlprefix;
  }

}