/* ChadoFeature.java
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

import java.sql.Timestamp;
import java.util.Hashtable;
import java.util.Vector;
import java.util.List;

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
  /** time last changed */
  private Timestamp timelastmodified;
  /** sequence length */
  private int length;
  /** features unique name */
  private String uniquename = null;
  /** features name */
  private String name;
  /** features residues */
  private byte[] residues;
  /** feature cvterm */
  private Cvterm cvterm;
  /** feature property */
  private ChadoFeatureProp featureprop;
  /** feature location */
  private ChadoFeatureLoc featureloc;
  /** feature relationship */
  private ChadoFeatureRelationship feature_relationship;
  /** feature organism */
  private ChadoOrganism organism;
  /** merged featureprops */
  private Hashtable qualifiers;
  /** list of ChadoFeatureProp */
  private List featurepropList;

  /**
   * Get the feature_id.
   * @return	the feature_id
   */
  public int getId()
  {
    return id;
  } 

  /**
   * Set the feature_id.
   * @param id	set the feature_id
   */
  public void setId(int id)
  {
    this.id = id;
  }

  /**
   * Get the postgres schema.
   * @return 	the postgres schema
   */
  public String getSchema()
  {
    return schema;
  }

  /**
   * Set the postgres schema.
   * @param schema	the postgres schema*
   */
  public void setSchema(String schema)
  {
    this.schema = schema;
  }

  /**
   * Get the last time feature was modified.
   * @return	the last time feature was modified
   */
  public Timestamp getTimelastmodified()
  {
    return timelastmodified;
  }

  /**
   * Set the last time feature was modified.
   * @param	the last time the feature was modified
   */
  public void setTimelastmodified(Timestamp timelastmodified)
  {
    this.timelastmodified = timelastmodified;
  }

  /**
   * Get sequence length. The length of the residue feature.
   * @return	the length of the residue feature
   */
  public int getLength()
  {
    return length;
  }

  /**
   * Set sequence length. The length of the residue feature.
   * @param length	the length of the residue feature
   */
  public void setLength(int length)
  {
    this.length = length;
  }

  /**
   * Get the unique name for a feature; may not be necessarily be 
   * particularly human-readable, although this is prefered. This name 
   * must be unique for this type of feature within this organism.
   * @return	the unique name for a feature.
   */
  public String getUniquename()
  {
    return uniquename;
  }

  /**
   * Set the unique name for a feature; may not be necessarily be 
   * particularly human-readable, although this is prefered. This name 
   * must be unique for this type of feature within this organism.
   * @param uniquename	the unique name for a feature.
   */
  public void setUniquename(String uniquename)
  {
    this.uniquename = uniquename;
  }

  /**
   * The optional human-readable common name for a feature, for display
   * purposes.
   * @return 	the name for a feature
   */
  public String getName()
  {
    return name;
  }

  /** 
   * Set the optional human-readable common name for a feature, for display 
   * purposes.
   * @param name	the name for a feature
   */
  public void setName(String name)
  {
    this.name = name;
  }

  /**
   * A sequence of alphabetic characters representing biological residues
   * (nucleic acids, amino acids). 
   * @return     the feature residues
   */
  public byte[] getResidues()
  {
    return residues;
  } 

  /**
   * Set the sequence of alphabetic characters representing biological residues 
   * (nucleic acids, amino acids).
   * @param residues	the feature residues
   */
  public void setResidues(byte[] residues)
  {
    this.residues = residues;
  }

  /**
  * A reference to a table:cvterm giving the feature type. 
  * This will typically be a Sequence Ontology identifier. 
  * @return the feature SO cvterm
  */
  public Cvterm getCvterm()
  {
    return cvterm;
  }

  /**
  * A reference to a table:cvterm giving the feature type. 
  * This will typically be a Sequence Ontology identifier.
  * @param cvterm  the feature SO cvterm
  */
  public void setCvterm(Cvterm cvterm)
  {
    this.cvterm = cvterm;
  }
  
  /**
   * Get the the feature property.
   * @return the <code>ChadoFeatureProp</code>
   */
  public ChadoFeatureProp getFeatureprop()
  {
    return featureprop;
  }

  /**
   * Set the feature property. 
   * @param featureprop  the feature property
   */
  public void setFeatureprop(ChadoFeatureProp featureprop)
  {
    this.featureprop = featureprop;
  }

  /**
   * Reference to the featureloc table. 
   * @return featureloc the <code>ChadoFeatureLoc</code> object.
   */
  public ChadoFeatureLoc getFeatureloc()
  {
    return featureloc;
  }

  /**
   * Reference to the featureloc table.
   * @param featureloc  the feature location
   */
  public void setFeatureloc(ChadoFeatureLoc featureloc)
  {
    this.featureloc = featureloc;
  }
  
  /**
   * Reference to the feature_relationship table.
   * @return feature_relationship the <code>ChadoFeatureRelationship</code>
   *         object.
   */
  public ChadoFeatureRelationship getFeature_relationship()
  {
    return feature_relationship;
  }

  /**
   * Reference to the feature_relationship table.
   * @param feature_relationship the <code>ChadoFeatureRelationship</code>
   *        object.
   */
  public void setFeature_relationship(
      ChadoFeatureRelationship feature_relationship)
  {
    this.feature_relationship = feature_relationship;
  }

  /**
   * Reference to the organism table for this feature.
   * @return  the <code>ChadoOrganism</code> object
   */
  public ChadoOrganism getOrganism()
  {
    return organism;
  }

  /**
   * Reference to the organism table for this feature.
   * @param organism the <code>ChadoOrganism</code> object
   */
  public void setOrganism(ChadoOrganism organism)
  {
    this.organism = organism;
  }

  /**
   * A list of feature properties.
   * @return  featurepropList a <code>List</code> of featureprop's
   *          for a feature.
   */
  public List getFeaturepropList()
  {
    return featurepropList;
  }

  /**
   * A list of feature properties.
   * @param featurepropList a <code>List</code> of featureprop's
   *        for a feature.
   */
  public void setFeaturepropList(List featurepropList)
  {
    this.featurepropList = featurepropList;
    
    for(int i=0; i<featurepropList.size(); i++)
    {
      ChadoFeatureProp featureprop = (ChadoFeatureProp)featurepropList.get(i);
      addQualifier(featureprop.getCvterm().getId(), featureprop);
    }
  }
  
  /**
   * Used in merging the qualifiers to store them as a <code>Hashtable</code> of
   * the cvterm type_id (of the property name) and the property values as a 
   * <code>Vector</code>.
   * @param	the cvterm type_id of the property name
   * @param	the property value	
   */
  public void addQualifier(long prop_type_id, ChadoFeatureProp featprop)
  {
    if(qualifiers == null)
      qualifiers = new Hashtable();
     
    final Long type_id = new Long(prop_type_id);
    if(qualifiers.containsKey(type_id))
    {
      Vector v = (Vector)qualifiers.get(type_id);
      v.add(featprop);
      qualifiers.put(type_id, v);
    }
    else
    {
      Vector v = new Vector();
      v.add(featprop);
      qualifiers.put(type_id, v);
    }
  }

  /**
   * Get the qualifiers which are stored as a <code>Hashtable</code> of cvterm
   * type_id (of the property name) and the property values as a <code>Vector</code>.
   * @return	the qualifiers as a <code>Hashtable</code>
   */
  public Hashtable getQualifiers()
  {
    return qualifiers;
  }

}
