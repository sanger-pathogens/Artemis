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

import java.util.Collection;
import java.util.List;
import java.util.Vector;

import org.gmod.schema.sequence.FeatureLoc;
import org.gmod.schema.sequence.FeatureProp;
import org.gmod.schema.sequence.FeatureRelationship;

public class Feature extends org.gmod.schema.sequence.Feature
{
  private String schema;
  /** used by getResidueFeatures */
  private Collection cvTermIds;
  private FeatureLoc featureloc;
  private FeatureProp featureprop;
  private FeatureRelationship feature_relationship;
  
  /**
   * Get the postgres schema.
   * @return    the postgres schema
   */
  public String getSchema() 
  {
    return schema;
  }

  /**
   * Set the postgres schema.
   * @param schema      the postgres schema*
   */
  public void setSchema(String schema)
  {
    this.schema = schema;
  }

  /**
   * @return
   */
  public Collection getCvTermIds() 
  {
    return cvTermIds;
  }

  /**
   * @param cvTermIds
   */
  public void setCvTermIds(Collection cvTermIds) 
  {
    this.cvTermIds = cvTermIds;
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

  public FeatureProp getFeatureprop()
  {
    return featureprop;
  }

  public void setFeatureprop(FeatureProp featureprop)
  {
    this.featureprop = featureprop;
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

  public static FeatureLoc getFeatureLoc(List locs, int srcfeature_id)
  {
    for(int i=0; i<locs.size(); i++)
    {
      FeatureLoc loc = (FeatureLoc)locs.get(i);
      if(loc.getFeatureBySrcFeatureId().getFeatureId() == srcfeature_id)
        return loc;
    }
    return null;
  }
  
  /**
   * Used in merging the qualifiers to store them as a <code>List</code> of
   * <code>FeatureProp</code>.
   * @param the FeatureProp 
   */
  public void addFeatureProp(final FeatureProp featureprop)
  {
    if(getFeatureProps() == null || getFeatureProps().size() == 0)
      setFeatureProps(new Vector());

    getFeatureProps().add(featureprop);
  }

}