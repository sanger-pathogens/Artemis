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
 * Alias GFF tags - values are stored in the synonym table and
 * linked to feature via feature_synonym table.
 */
public class Alias
{
  private Integer feature_id;
  private String name;
  private String schema;
  private String uniquename;
  private String cvterm_name;
  
  public String getSchema()
  {
    return schema;
  }

  public void setSchema(String schema)
  {
    this.schema = schema;
  }

  public String getUniquename()
  {
    return uniquename;
  }

  public void setUniquename(String uniquename)
  {
    this.uniquename = uniquename;
  }

  public Integer getFeature_id()
  {
    return feature_id;
  }
  
  public void setFeature_id(Integer feature_id)
  {
    this.feature_id = feature_id;
  }
  
  public String getName()
  {
    return name;
  }
  
  public void setName(String name)
  {
    this.name = name;
  }

  public String getCvterm_name()
  {
    return cvterm_name;
  }

  public void setCvterm_name(String cvterm_name)
  {
    this.cvterm_name = cvterm_name;
  }

}
