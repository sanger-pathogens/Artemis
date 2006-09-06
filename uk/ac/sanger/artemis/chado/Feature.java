/* Feature.java
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
import java.util.Collection;


/**
 * Chado feature representation of the feature table and with references
 * to the featureloc, featureprop, feature_relationship tables.
 */
public class Feature
{
 
  /** schema */
  private String schema;
  /** feature id */
  private int featureId;
  /** time last changed */
  private Timestamp timelastmodified;
  /** sequence length */
  private int seqLen;
  /** features unique name */
  private String uniqueName = null;
  /** features name */
  private String name;
  /** features residues */
  private byte[] residues;
  /** anotated or result of automated analysis */
  private boolean analysis;
  /** whether this feature has been obsoleted */
  private boolean obsolete;
  /** feature cvTerm */
  private CvTerm cvTerm;
  /** feature property */
  private FeatureProp featureprop;
  /** feature location (for a given srcfeature) */
  private FeatureLoc featureloc;
  /** feature relationship */
  private FeatureRelationship feature_relationship;
  /** feature organism */
  private Organism organism;
  /** optional primary public stable identifier */
  private FeatureDbXRef featureDbXRef;
  /** merged featureprops */
  private Hashtable qualifiers;
  /** list of FeatureProp */
  private Collection featureProps;
  /** list of FeatureRelationship children */
  private Collection featureRelationshipsForObjectId;
  /** list of FeatureRelationship parent */
  private Collection featureRelationshipsForSubjectId;
  /** list of feature dbxrefs (FeatureDbxref) */
  private Collection featureDbXRefs;
  /** list of feature locations for a feature_id */
  private Collection featurelocsForFeatureId;
  /** used by getResidueFeatures */
  private Collection cvTermIds;
  /** list of feature synonyms */
  private Collection featureSynonymsForFeatureId;



  /**
   * Get the feature_id.
   * @return	the feature_id
   */
  public int getFeatureId()
  {
    return featureId;
  } 

  /**
   * Set the feature_id.
   * @param id	set the feature_id
   */
  public void setFeatureId(int featureId)
  {
    this.featureId = featureId;
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
  public int getSeqLen()
  {
    return seqLen;
  }

  /**
   * Set sequence length. The length of the residue feature.
   * @param length	the length of the residue feature
   */
  public void setSeqLen(int seqLen)
  {
    this.seqLen = seqLen;
  }

  /**
   * Get the unique name for a feature; may not be necessarily be 
   * particularly human-readable, although this is prefered. This name 
   * must be unique for this type of feature within this organism.
   * @return	the unique name for a feature.
   */
  public String getUniqueName()
  {
    return uniqueName;
  }

  /**
   * Set the unique name for a feature; may not be necessarily be 
   * particularly human-readable, although this is prefered. This name 
   * must be unique for this type of feature within this organism.
   * @param uniquename	the unique name for a feature.
   */
  public void setUniqueName(String uniqueName)
  {
    this.uniqueName = uniqueName;
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
   * Detemine if this is anotated or result of automated analysis.
   * @return
   */
  public boolean isAnalysis()
  {
    return analysis;
  }

  /**
   * Detemine if this is anotated or result of automated analysis.
   * @param analysis
   */
  public void setAnalysis(boolean analysis)
  {
    this.analysis = analysis;
  }

  /**
   * Determine if this feature has been made obsolete.
   * @return
   */
  public boolean isObsolete()
  {
    return obsolete;
  }

  /**
   * Determine if this feature has been made obsolete.
   * @param obsolete
   */
  public void setObsolete(boolean obsolete)
  {
    this.obsolete = obsolete;
  }
  
  /**
  * A reference to a table:cvTerm giving the feature type. 
  * This will typically be a Sequence Ontology identifier. 
  * @return the feature SO cvTerm
  */
  public CvTerm getCvTerm()
  {
    return cvTerm;
  }

  /**
  * A reference to a table:cvTerm giving the feature type. 
  * This will typically be a Sequence Ontology identifier.
  * @param cvTerm  the feature SO cvTerm
  */
  public void setCvTerm(CvTerm cvTerm)
  {
    this.cvTerm = cvTerm;
  }
  
  /**
   * Get the the feature property.
   * @return the <code>FeatureProp</code>
   */
  public FeatureProp getFeatureprop()
  {
    return featureprop;
  }

  /**
   * Set the feature property. 
   * @param featureprop  the feature property
   */
  public void setFeatureprop(FeatureProp featureprop)
  {
    this.featureprop = featureprop;
  }

  /**
   * Reference to the featureloc table. 
   * @return featureloc the <code>FeatureLoc</code> object.
   */
  public FeatureLoc getFeatureloc()
  {
    return featureloc;
  }

  /**
   * Reference to the featureloc table.
   * @param featureloc  the feature location
   */
  public void setFeatureloc(FeatureLoc featureloc)
  {
    this.featureloc = featureloc;
  }
  
  /**
   * Reference to the feature_relationship table.
   * @return feature_relationship the <code>FeatureRelationship</code>
   *         object.
   */
  public FeatureRelationship getFeature_relationship()
  {
    return feature_relationship;
  }

  /**
   * Reference to the feature_relationship table.
   * @param feature_relationship the <code>FeatureRelationship</code>
   *        object.
   */
  public void setFeature_relationship(
      FeatureRelationship feature_relationship)
  {
    this.feature_relationship = feature_relationship;
  }

  /**
   * Reference to the organism table for this feature.
   * @return  the <code>Organism</code> object
   */
  public Organism getOrganism()
  {
    return organism;
  }

  /**
   * Reference to the organism table for this feature.
   * @param organism the <code>Organism</code> object
   */
  public void setOrganism(Organism organism)
  {
    this.organism = organism;
  }

  
  /**
   * An optional primary public stable identifier for this feature. 
   * Secondary identifiers and external dbxrefs go in table:feature_dbxref
   * @return
   */
  public FeatureDbXRef getFeatureDbXRef()
  {
    return featureDbXRef;
  }

  /**
   * An optional primary public stable identifier for this feature. 
   * Secondary identifiers and external dbxrefs go in table:feature_dbxref
   * @param dbxref
   */
  public void setFeatureDbXRef(FeatureDbXRef featureDbXRef)
  {
    this.featureDbXRef = featureDbXRef;
  }

  public Collection getCvTermIds()
  {
    return cvTermIds;
  }

  public void setCvTermIds(Collection cvTermIds)
  {
    this.cvTermIds = cvTermIds;
  }

  /**
   * A list of feature properties.
   * @return  featurepropList a <code>List</code> of featureprop's
   *          for a feature.
   */
  public Collection getFeatureProps()
  {
    return featureProps;
  }

  /**
   * A list of feature properties.
   * @param featurepropList a <code>List</code> of featureprop's
   *        for a feature.
   */
  public void setFeatureProps(Collection featureProps)
  {
    this.featureProps = featureProps;
    
    for(int i=0; i<featureProps.size(); i++)
    {
      FeatureProp featureprop = (FeatureProp)(((List)featureProps).get(i));
      addQualifier(featureprop.getCvTerm().getCvTermId(), featureprop);
    }
  }
  
  /**
   * Get a list of <code>FeatureRelationship</code> children
   * @return
   */
  public Collection getFeatureRelationshipsForObjectId()
  {
    return featureRelationshipsForObjectId;
  }

  /**
   * Set a list of <code>FeatureRelationship</code> children
   * @param featureRelationshipsForObjectId
   */
  public void setFeatureRelationshipsForObjectId(
      Collection featureRelationshipsForObjectId)
  {
    this.featureRelationshipsForObjectId = featureRelationshipsForObjectId;
  }

  /**
   * Get a list of <code>FeatureRelationship</code> parent
   * @return
   */
  public Collection getFeatureRelationshipsForSubjectId()
  {
    return featureRelationshipsForSubjectId;
  }

  /**
   * Set a list of <code>FeatureRelationship</code> parent
   * @param featureRelationshipsForSubjectId
   */
  public void setFeatureRelationshipsForSubjectId(
      Collection featureRelationshipsForSubjectId)
  {
    this.featureRelationshipsForSubjectId = featureRelationshipsForSubjectId;
  }
  
  /**
   * Get a list of feature dbxrefs (<code>FeatureDbxref</code>)
   * @return
   */
  public Collection getFeatureDbXRefs()
  {
    return featureDbXRefs;
  }

  /**
   * Set a list of feature dbxrefs (<code>FeatureDbxref</code>)
   * @param featureDbxrefs
   */
  public void setFeatureDbXRefs(Collection featureDbXRefs)
  {
    this.featureDbXRefs = featureDbXRefs;
  }

  /**
   * Get list of feature locations for a feature_id
   * @return
   */
  public Collection getFeaturelocsForFeatureId()
  {
    return featurelocsForFeatureId;
  }

  /**
   * Set a list of feature locations for a feature_id
   * @param featurelocsForFeatureId
   */
  public void setFeaturelocsForFeatureId(Collection featurelocsForFeatureId)
  {
    this.featurelocsForFeatureId = featurelocsForFeatureId;
  }
  
  public Collection getFeatureSynonymsForFeatureId()
  {
    return featureSynonymsForFeatureId;
  }

  public void setFeatureSynonymsForFeatureId(Collection featureSynonymsForFeatureId)
  {
    this.featureSynonymsForFeatureId = featureSynonymsForFeatureId;
  }
  
  
  
  /**
   * Used in merging the qualifiers to store them as a <code>Hashtable</code> of
   * the cvTerm type_id (of the property name) and the property values as a 
   * <code>Vector</code>.
   * @param	the cvTerm type_id of the property name
   * @param	the property value	
   */
  public void addQualifier(long prop_type_id, FeatureProp featprop)
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
   * Get the qualifiers which are stored as a <code>Hashtable</code> of cvTerm
   * type_id (of the property name) and the property values as a <code>Vector</code>.
   * @return	the qualifiers as a <code>Hashtable</code>
   */
  public Hashtable getQualifiers()
  {
    return qualifiers;
  }

  /**
   * Utility for finding a feature location from a List that corresponds
   * to a particular srcfeature.
   * @param locs            List of FeatureLoc
   * @param srcfeature_id   srcfeature id
   * @return
   */
  public static FeatureLoc getFeatureLoc(List locs, int srcfeature_id)
  {
    for(int i=0; i<locs.size(); i++)
    {
      FeatureLoc loc = (FeatureLoc)locs.get(i);
      if(loc.getSrcfeature_id() == srcfeature_id)
        return loc;
    }
    return null;
  }
 
}
