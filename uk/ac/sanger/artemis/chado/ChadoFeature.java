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
  /** feature organism identifier */
  private int organism_id;
  /** hashtable of qualifiers */
  private Hashtable qualifiers;
  /** phase */
  private int phase;

  /**
   *
   * Get the feature_id.
   * @return	the feature_id
   *
   */
  public int getId()
  {
    return id;
  } 

  /**
   *
   * Set the feature_id.
   * @param id	set the feature_id
   *
   */
  public void setId(int id)
  {
    this.id = id;
  }

  /**
   *
   * Get the postgres schema.
   * @return 	the postgres schema
   *
   */
  public String getSchema()
  {
    return schema;
  }

  /**
   *
   * Set the postgres schema.
   * @param schema	the postgres schema
   *
   */
  public void setSchema(String schema)
  {
    this.schema = schema;
  }

  /**
   *
   * Get the parent feature_id.
   * @return	the parent feature_id
   *
   */
  public String getObject_id()
  {
    return object_id;
  }

  /**
   *
   * Set the parent feature_id.
   * @param object_id	the parent feature_id
   *
   */
  public void setObject_id(String object_id)
  {
    this.object_id = object_id;
  }

  /**
   *
   * Get the last time feature was modified.
   * @return	the last time feature was modified
   *
   */
  public Date getTimelastmodified()
  {
    return timelastmodified;
  }

  /**
   *
   * Set the last time feature was modified.
   * @param	the last time the feature was modified
   *
   */
  public void setTimelastmodified(Date timelastmodified)
  {
    this.timelastmodified = timelastmodified;
  }

  /**
   *
   * Get the strand. The orientation/directionality of the location. 
   * Should be 0, -1 or +1.
   * @return	the strand
   *
   */
  public int getStrand()
  {
    return strand;
  }

  /**
   *
   * Set the strand. The orientation/directionality of the location. 
   * Should be 0, -1 or +1.
   * @param strand	the strand
   *
   */
  public void setStrand(int strand)
  {
    this.strand = strand;
  }

  /**
   *
   * Get the fmin value. The leftmost/minimal boundary in the linear range 
   * represented by the featureloc. To convert this to the eftmost position 
   * in a base-oriented system (e.g. GFF) add 1 to fmin.
   * @return	the fmin value
   *
   */
  public int getFmin()
  {
    return fmin;
  }

  /**
   *
   * Set the fmin value. The leftmost/minimal boundary in the linear range 
   * represented by the featureloc. To convert this to the eftmost position 
   * in a base-oriented system (e.g. GFF) add 1 to fmin.
   * @param fmin	the fmin value
   *
   */
  public void setFmin(int fmin)
  {
    this.fmin = fmin;
  }

  /**
   *
   * Get the fmax value. The rightmost/maximal boundary in the linear range 
   * represented by the featureloc. No conversion is required to go from fmax 
   * to the rightmost coordinate in a base-oriented system that counts from 1 
   * (eg GFF).
   * @return	the fmax value
   *
   */
  public int getFmax()
  {
    return fmax;
  }

  /**
   *
   * Set the fmax value. The rightmost/maximal boundary in the linear range 
   * represented by the featureloc. No conversion is required to go from fmax 
   * to the rightmost coordinate in a base-oriented system that counts from 1 
   * (eg GFF).
   * @param fmax	the fmax value
   *
   */
  public void setFmax(int fmax)
  {
    this.fmax = fmax;
  }

  /**
   *
   * Get sequence length. The length of the residue feature.
   * @return	the length of the residue feature
   *
   */
  public int getLength()
  {
    return length;
  }

  /**
   *
   * Set sequence length. The length of the residue feature.
   * @param length	the length of the residue feature
   *
   */
  public void setLength(int length)
  {
    this.length = length;
  }

  /**
   *
   * Get the unique name for a feature; may not be necessarily be 
   * particularly human-readable, although this is prefered. This name 
   * must be unique for this type of feature within this organism.
   * @return	the unique name for a feature.
   *
   */
  public String getUniquename()
  {
    return uniquename;
  }

  /**
   *
   * Set the unique name for a feature; may not be necessarily be 
   * particularly human-readable, although this is prefered. This name 
   * must be unique for this type of feature within this organism.
   * @param uniquename	the unique name for a feature.
   *
   */
  public void setUniquename(String uniquename)
  {
    this.uniquename = uniquename;
  }

  /**
   * 
   * The optional human-readable common name for a feature, for display
   * purposes.
   * @return 	the name for a feature
   *
   */
  public String getName()
  {
    return name;
  }

  /**
   * 
   * Set the optional human-readable common name for a feature, for display 
   * purposes.
   * @param name	the name for a feature
   *
   */
  public void setName(String name)
  {
    this.name = name;
  }

  /**
   *
   * A sequence of alphabetic characters representing biological residues
   * (nucleic acids, amino acids). 
   * @return     the feature residues
   *
   */
  public byte[] getResidues()
  {
    return residues;
  } 

  /**
   *
   * Set the sequence of alphabetic characters representing biological residues 
   * (nucleic acids, amino acids).
   * @param residues	the feature residues
   *
   */
  public void setResidues(byte[] residues)
  {
    this.residues = residues;
  }

  /**
   *
   * A required reference to a table:cvterm giving the feature type. 
   * This will typically be a Sequence Ontology identifier. 
   * @return	the feature SO identifier
   *
   */
  public long getType_id()
  {
    return type_id;
  }

  /**
   *
   * A required reference to a table:cvterm giving the feature type. 
   * This will typically be a Sequence Ontology identifier.
   * @param type_id	the feature SO identifier
   *
   */
  public void setType_id(long type_id)
  {
    this.type_id = type_id;
  }

  /**
   *
   * The name of the property/slot is a cvterm. The meaning of the property 
   * is defined in that cvterm. Certain properties will only apply to certain 
   * feature types; this will be handled by the Sequence Ontology
   * @return	the type identifier for the name of the property 
   *
   */
  public long getProp_type_id()
  {
    return prop_type_id;
  }

  /**
   *
   * Set the name of the property/slot is a cvterm. The meaning of the property 
   * is defined in that cvterm. Certain properties will only apply to certain 
   * feature types; this will be handled by the Sequence Ontology
   * @param prop_type_id	the type identifier for the name of the property
   *
   */
  public void setProp_type_id(long prop_type_id)
  {
    this.prop_type_id = prop_type_id;
  }

  /**
   *
   * Get the value of the property, represented as text.
   * @return	the value of the property
   *
   */
  public String getValue()
  {
    return value;
  }

  /**
   *
   * Set the value of the property, represented as text. 
   * @param value	the value of the property
   *
   */
  public void setValue(String value)
  {
    this.value = value;
  }

  /**
   *
   * Organism abbreviation.
   * @return	the organism abbreviation	
   *
   */
  public String getAbbreviation()
  {
    return abbreviation;
  }

  /**
   *
   * Organism abbreviation.
   * @return abbreviation	the organism abbreviation
   *
   */
  public void setAbbreviation(String abbreviation)
  {
    this.abbreviation = abbreviation;
  }

  /**
   *
   * Organism identifier.
   * @return     the organism id
   *
   */
  public int getOrganism_id()
  {
    return organism_id;
  }

  /**
   *
   * Organism identifier.
   * @param organism_id    the organism id
   *
   */
  public void setOrganism_id(final int organism_id)
  {
    this.organism_id = organism_id;
  }


  /**
   *
   * The source feature which the location in featureloc is relative to.
   * @return    the source feature
   *
   */
  public int getSrcfeature_id()
  {
    return srcfeature_id;
  }

  /**
   *
   * The source feature which the location in featureloc is relative to.
   * @param srcfeature_id	the source feature
   *
   */
  public void setSrcfeature_id(int srcfeature_id)
  {
    this.srcfeature_id = srcfeature_id;
  }

  /**
   *
   * The phase of translation wrt srcfeature_id. Values are 0,1,2. 
   * @return	the phase
   *
   */
  public int getPhase()
  {
    return phase;  
  }

  /**
   *
   * The phase of translation wrt srcfeature_id. Values are 0,1,2. 
   * @param phase	the phase
   *
   */
  public void setPhase(int phase)
  {
    this.phase = phase;
  }

  /**
   *
   * Used in merging the qualfiers to store them as a <code>Hashtable</code> of
   * the cvterm type_id of the property name and the values as a <code>Vector</code>.
   * @param	the cvterm type_id of the property name
   * @param	the property value	
   *
   */
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

  /**
   *
   * Get the qualfiers to store them as a <code>Hashtable</code> of cvterm
   * type_id (of the property name) and the values as a <code>Vector</code>.
   * @return	the qualifiers as a <code>Hashtable</code>
   *
   */
  public Hashtable getQualifiers()
  {
    return qualifiers;
  }


}
