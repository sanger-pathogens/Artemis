/* 
 *
 * created: 2007
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2007  Genome Research Limited
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
import java.util.Iterator;
import java.util.List;
import java.util.Vector;

import org.gmod.schema.cv.CvTerm;
import org.gmod.schema.sequence.Feature;
import org.gmod.schema.sequence.FeatureRelationship;

import uk.ac.sanger.artemis.io.DocumentEntry;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.LazyQualifierValue;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.Document;
import uk.ac.sanger.artemis.util.StringVector;


public class ClusterLazyQualifierValue implements LazyQualifierValue
{
  /** force complete loading of the data */
  private boolean forceLoad = false;
  /** data loaded */
  private boolean lazyLoaded = false;
  private String value;
  private GFFStreamFeature feature;

  
  /**
   * Qualifier object to handle lazy loading of cluster/ortholog/paralog data
   * @param matchFeature
   * @param featureId
   */
  public ClusterLazyQualifierValue(final String value, final GFFStreamFeature feature)
  {
    this.value = value;
    this.feature = feature;
  }

  public String getString()
  {
    if(forceLoad && !lazyLoaded)
      return getHardString();
    else
      return value;
  }

  public boolean isLazyLoaded()
  {
    return lazyLoaded;
  }

  public void setForceLoad(boolean forceLoad)
  {
    this.forceLoad = forceLoad;
  }

  private String getHardString()
  {
    lazyLoaded = true;
    
    final String featureId = (String) feature.getQualifierByName("feature_id").getValues().get(0);
    final List featureIds = new Vector();
    
    StringVector strings = StringVector.getStrings(value, ";");
    String f_id[] = ArtemisUtils.getString(strings, "object_id=").split("=");
    featureIds.add( Integer.valueOf(f_id[1]) );
    
    String rank = ArtemisUtils.getString(strings, "rank");
    
    final Document document = ((DocumentEntry)feature.getEntry()).getDocument();
    List clusters = ((DatabaseDocument)document).getClustersByFeatureIds(featureIds);
    
    value = "";
    for(int i=0;i<clusters.size(); i++)
    {
      final Feature clusterFeature = (Feature)clusters.get(i);
      final Collection subjects = clusterFeature.getFeatureRelationshipsForSubjectId();
      final Iterator it = subjects.iterator();
      int cnt = 0;
      while(it.hasNext())
      {
        FeatureRelationship fr = (FeatureRelationship)it.next();
        Feature subjectFeature = fr.getFeatureBySubjectId();
        
        if(subjectFeature.getFeatureId() != Integer.parseInt(featureId))
        {
          if(!value.equals(""))
            value = value.concat(", ");
          
          value = value.concat(subjectFeature.getOrganism().getCommonName()+":");
          
          String geneName = subjectFeature.getUniqueName();
          if(!subjectFeature.getCvTerm().getName().equals("gene"))
          {
            Feature parent = getParentFeature(subjectFeature);
            if(parent.getCvTerm().getName().equals("gene"))
              geneName = parent.getUniqueName();
            else if(parent != null)
            {
              parent = getParentFeature(parent);
              if(parent.getCvTerm().getName().equals("gene"))
                geneName = parent.getUniqueName();
            }
          }
          
          value = value.concat(geneName+" link="+
              subjectFeature.getUniqueName());
          
          cnt++;
        }
      }
      if(cnt > 1)
        value = value.concat("; cluster_name="+clusterFeature.getUniqueName());
      else
        value = value.concat("; match_name="+clusterFeature.getUniqueName());
    }
    value = value.concat("; "+rank);

    return value;
  }
  
  /**
   * Given a chado Feature find the parent from its feature_relationship
   * @param childFeature
   * @return
   */
  private Feature getParentFeature(final Feature childFeature)
  {
    Collection featureRelationships = childFeature.getFeatureRelationshipsForSubjectId();
    Iterator it = featureRelationships.iterator();
    while(it.hasNext())
    {
      FeatureRelationship frSubject = (FeatureRelationship)it.next();
      CvTerm cvTerm = frSubject.getCvTerm();
      
      if(cvTerm.getName().equals("derives_from") ||
         cvTerm.getName().indexOf("part_of")>-1)
        return frSubject.getFeatureByObjectId();
     
    }
    return null;
  }
  

}