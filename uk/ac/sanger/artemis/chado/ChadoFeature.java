/* ChadoFeature.java
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


package uk.ac.sanger.artemis.chado;

import java.util.Date;
import java.util.Hashtable;
import java.util.Vector;

/**
*
* Chado feature representation of feature, featureloc, featureprop,
* feature_relationship tables
*
**/
public class ChadoFeature
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
  /** hashtable of qualifiers */
  private Hashtable qualifiers;
  /** phase */
  private int phase;

  /**
   *
   * Get the feature_id.
   *
   */
  public int getId()
  {
    return id;
  } 

  /**
   *
   * Set the feature_id.
   *
   */
  public void setId(int id)
  {
    this.id = id;
  }

  /**
   *
   * Get the postgres schema.
   *
   */
  public String getSchema()
  {
    return schema;
  }

  /**
   *
   * Set the postgres schema.
   *
   */
  public void setSchema(String schema)
  {
    this.schema = schema;
  }

  /**
   *
   * Get the parent feature_id.
   *
   */
  public String getObject_id()
  {
    return object_id;
  }

  /**
   *
   * Set the parent feature_id.
   *
   */
  public void setObject_id(String object_id)
  {
    this.object_id = object_id;
  }

  /**
   *
   * Get the last time feature was modified.
   *
   */
  public Date getTimelastmodified()
  {
    return timelastmodified;
  }

  /**
   *
   * Set the last time feature was modified.
   *
   */
  public void setTimelastmodified(Date timelastmodified)
  {
    this.timelastmodified = timelastmodified;
  }

  /**
   *
   * Get the strand.
   *
   */
  public int getStrand()
  {
    return strand;
  }

  /**
   *
   * Set the strand.
   *
   */
  public void setStrand(int strand)
  {
    this.strand = strand;
  }

  /**
   *
   * Get the fmin value.
   *
   */
  public int getFmin()
  {
    return fmin;
  }

  /**
   *
   * Set the fmin value.
   *
   */
  public void setFmin(int fmin)
  {
    this.fmin = fmin;
  }

  /**
   *
   * Get the fmax value.
   *
   */
  public int getFmax()
  {
    return fmax;
  }

  /**
   *
   * Set the fmax value.
   *
   */
  public void setFmax(int fmax)
  {
    this.fmax = fmax;
  }

  /**
   *
   * Get sequence length.
   *
   */
  public int getLength()
  {
    return length;
  }

  /**
   *
   * Set sequence length.
   *
   */
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

  public int getPhase()
  {
    return phase;  
  }

  public void setPhase(int phase)
  {
    this.phase = phase;
  }

  public void addQualifier(long prop_type_id, String value)
  {
    if(qualifiers == null)
      qualifiers = new Hashtable();
     
    final Long type_id = new Long(prop_type_id);
    if(qualifiers.containsKey(type_id))
    {
      Vector v = (Vector)qualifiers.get(type_id);
      v.add(value);
      qualifiers.put(type_id, v);
    }
    else
    {
      Vector v = new Vector();
      v.add(value);
      qualifiers.put(type_id, v);
    }
  }

  public Hashtable getQualifiers()
  {
    return qualifiers;
  }


}
