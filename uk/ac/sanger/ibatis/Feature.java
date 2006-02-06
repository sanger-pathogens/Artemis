/* Feature.java
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

import java.util.Date;

public class Feature
{
 
  /** schema */
  private String schema;
  /** feature id */
  private int id;
  /** id of the parent feature */
  private int srcfeature_id;
  /** id of the parent feature */
  private String object_id;
  /** time last changed */
  private Date timelastmodified;
  /** +1 or -1 depending if on the forward or reverse */
  private int strand;
  /** start position */
  private int fmin;
  /** end position */
  private int fmax;
  /** sequence length */
  private int length;
  /** features unique name */
  private String uniquename = null;
  /** features name */
  private String name;
  /** features residues */
  private byte[] residues;
  /** feature type id */
  private long type_id;
  /** feature property type id */
  private long prop_type_id;
  /** feature description value */
  private String value;
  /** feature organism abbreviation */
  private String abbreviation;

 
  public int getId()
  {
    return id;
  } 

  public void setId(int id)
  {
    this.id = id;
  }

  public String getSchema()
  {
    return schema;
  }

  public void setSchema(String schema)
  {
    this.schema = schema;
  }

  public String getObject_id()
  {
    return object_id;
  }

  public void setObject_id(String object_id)
  {
    this.object_id = object_id;
  }

  public Date getTimelastmodified()
  {
    return timelastmodified;
  }

  public void setTimelastmodified(Date timelastmodified)
  {
    this.timelastmodified = timelastmodified;
  }

  public int getStrand()
  {
    return strand;
  }

  public void setStrand(int strand)
  {
    this.strand = strand;
  }

  public int getFmin()
  {
    return fmin;
  }

  public void setFmin(int fmin)
  {
    this.fmin = fmin;
  }

  public int getFmax()
  {
    return fmax;
  }

  public void setFmax(int fmax)
  {
    this.fmax = fmax;
  }

  public int getLength()
  {
    return length;
  }

  public void setLength(int length)
  {
    this.length = length;
  }

  public String getUniquename()
  {
    return uniquename;
  }

  public void setUniquename(String uniquename)
  {
    this.uniquename = uniquename;
  }

  public String getName()
  {
    return name;
  }

  public void setName(String name)
  {
    this.name = name;
  }

  public byte[] getResidues()
  {
    return residues;
  } 

  public void setResidues(byte[] residues)
  {
    this.residues = residues;
  }

  public long getType_id()
  {
    return type_id;
  }

  public void setType_id(long type_id)
  {
    this.type_id = type_id;
  }

  public long getProp_type_id()
  {
    return prop_type_id;
  }

  public void setProp_type_id(long prop_type_id)
  {
    this.prop_type_id = prop_type_id;
  }

  public String getValue()
  {
    return value;
  }

  public void setValue(String value)
  {
    this.value = value;
  }

  public String getAbbreviation()
  {
    return abbreviation;
  }

  public void setAbbreviation(String abbreviation)
  {
    this.abbreviation = abbreviation;
  }

  public int getSrcfeature_id()
  {
    return srcfeature_id;
  }

  public void setSrcfeature_id(int srcfeature_id)
  {
    this.srcfeature_id = srcfeature_id;
  }

}
