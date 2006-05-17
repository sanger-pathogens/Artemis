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
public class ChadoSynonym
{
  private String schema;
  private String uniquename;
  
  // feature_synonym
  private Integer feature_id;
  private Integer synonym_id;
  private Integer pub_id;
  private boolean is_current;
  private boolean is_internal;
 
  // synonym
  private String name;
  private ChadoCvterm cvterm;
  private String synonym_sgml;
  
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

  public Integer getSynonym_id()
  {
    return synonym_id;
  }

  public void setSynonym_id(Integer synonym_id)
  {
    this.synonym_id = synonym_id;
  }
  
  public ChadoCvterm getCvterm()
  {
    return cvterm;
  }

  public void setCvterm(ChadoCvterm cvterm)
  {
    this.cvterm = cvterm;
  }

  public boolean isIs_current()
  {
    return is_current;
  }

  public void setIs_current(boolean is_current)
  {
    this.is_current = is_current;
  }

  public boolean isIs_internal()
  {
    return is_internal;
  }

  public void setIs_internal(boolean is_internal)
  {
    this.is_internal = is_internal;
  }

  public Integer getPub_id()
  {
    return pub_id;
  }

  public void setPub_id(Integer pub_id)
  {
    this.pub_id = pub_id;
  }
  
  public String getSynonym_sgml()
  {
    return synonym_sgml;
  }

  public void setSynonym_sgml(String synonym_sgml)
  {
    this.synonym_sgml = synonym_sgml;
  }

}
