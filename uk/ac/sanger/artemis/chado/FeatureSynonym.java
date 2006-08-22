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
public class FeatureSynonym
{
  private String schema;
  private String uniquename;
  
  // feature_synonym
  private Integer feature_id;
  private Integer synonym_id;
  private Integer pub_id;
  private boolean current;
  private boolean internal;
  
  private Synonym synonym;
  
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

  public Integer getSynonym_id()
  {
    return synonym_id;
  }

  public void setSynonym_id(Integer synonym_id)
  {
    this.synonym_id = synonym_id;
  }

  public Integer getPub_id()
  {
    return pub_id;
  }

  public void setPub_id(Integer pub_id)
  {
    this.pub_id = pub_id;
  }
  
  public boolean isCurrent()
  {
    return current;
  }

  public void setCurrent(boolean current)
  {
    this.current = current;
  }

  public boolean isInternal()
  {
    return internal;
  }

  public void setInternal(boolean internal)
  {
    this.internal = internal;
  }

  public Synonym getSynonym()
  {
    return synonym;
  }

  public void setSynonym(Synonym synonym)
  {
    this.synonym = synonym;
  }

}
