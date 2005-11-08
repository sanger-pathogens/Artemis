/* Organism.java
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

package uk.ac.sanger.ibatis;


public class Organism
{

  private int id;
  private String abbreviation;
  private String genus;
  private String species;
  private String common_name;
  private String comment;

  public int getId()
  {
    return id;
  }

  public void setId(int id)
  {
    this.id = id;
  }

  public String getAbbreviation()
  {
    return abbreviation;
  }

  public void setAbbreviation(String abbreviation)
  {
    this.abbreviation = abbreviation;
  }

  public String getGenus()
  {
    return genus;
  }

  public void setGenus(String genus)
  {
    this.genus = genus;
  }

  public String getSpecies()
  {
    return species;
  }

  public void setSpecies(String species)
  {
    this.species = species;
  }

  public String getCommon_name()
  {
    return common_name;
  }

  public void setCommon_name(String common_name)
  {
    this.common_name = common_name;
  }

  public String getComment()
  {
    return comment;
  }

  public void setComment()
  {
    this.comment = comment;
  }

}
