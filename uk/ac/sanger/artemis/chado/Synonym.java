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
public class Synonym
{
  // synonym
  private Integer synonym_id;
  private String name;
  private Cvterm cvterm;
  private String synonym_sgml;
  
  public Integer getSynonym_id()
  {
    return synonym_id;
  }

  public void setSynonym_id(Integer synonym_id)
  {
    this.synonym_id = synonym_id;
  }

  public String getName()
  {
    return name;
  }
  
  public void setName(String name)
  {
    this.name = name;
  }
  
  public Cvterm getCvterm()
  {
    return cvterm;
  }

  public void setCvterm(Cvterm cvterm)
  {
    this.cvterm = cvterm;
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