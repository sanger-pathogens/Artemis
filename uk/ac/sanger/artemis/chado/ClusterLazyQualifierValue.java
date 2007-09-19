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
import java.util.Hashtable;
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
  private List clusters;
  
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

  /**
   * This speeds up loading for a list of ortho/paralogs
   * @param values List of ClusterLazyQualifierValue
   * @param feature
   */
  public static void setClusterFromValueList(final List values, final GFFStreamFeature feature)
  {
    final List clusterFeatureIds = new Vector();
    final Hashtable hash = new Hashtable(values.size());
    for(int i=0; i<values.size(); i++)
    {
      ClusterLazyQualifierValue  lazyValue = (ClusterLazyQualifierValue)values.get(i);
      StringVector strings = StringVector.getStrings(lazyValue.getValue(), ";");
      String f_id[] = ArtemisUtils.getString(strings, "object_id=").split("=");
      Integer clusterFeatureId = Integer.valueOf(f_id[1]);
      clusterFeatureIds.add(clusterFeatureId);
      lazyValue.initCluster();
      hash.put(clusterFeatureId, lazyValue);
    }
    
    final Document document = ((DocumentEntry)feature.getEntry()).getDocument();
    List allClusters = ((DatabaseDocument)document).getClustersByFeatureIds(clusterFeatureIds);
    
   
    // 
    // get parent gene
    /*
    List subjectIds = new Vector();
    for(int i=0;i<allClusters.size(); i++)
    {
      final Feature clusterFeature = (Feature)allClusters.get(i);
      final Collection subjects = clusterFeature.getFeatureRelationshipsForSubjectId();
      final Iterator it = subjects.iterator();
      while(it.hasNext())
      {
        FeatureRelationship fr = (FeatureRelationship)it.next();
        Feature subjectFeature = fr.getFeatureBySubjectId();
        subjectIds.add(new Integer(subjectFeature.getFeatureId()));
      }
    }
    
    List geneFeatures = 
      ((DatabaseDocument)document).getParentFeaturesByChildFeatureIds(subjectIds);
    final Hashtable genes = new Hashtable(geneFeatures.size());
    subjectIds = new Vector();
    for(int i=0; i<geneFeatures.size(); i++)
    {
      FeatureRelationship gene = (FeatureRelationship)geneFeatures.get(i);

      if(gene.getFeatureByObjectId().getCvTerm().getName().equals("gene") ||
         gene.getFeatureByObjectId().getCvTerm().getName().equals("pseudogene"))
      {
        genes.put(new Integer(gene.getFeatureBySubjectId().getFeatureId()), 
                              gene.getFeatureByObjectId().getUniqueName());
      }
      else
      {
        Integer objectId = new Integer(gene.getFeatureByObjectId().getFeatureId());
        subjectIds.add(objectId);
        genes.put(new Integer(gene.getFeatureBySubjectId().getFeatureId()), objectId);
      }
    }
    
    geneFeatures = 
      ((DatabaseDocument)document).getParentFeaturesByChildFeatureIds(subjectIds);
    for(int i=0; i<geneFeatures.size(); i++)
    {
      FeatureRelationship gene = (FeatureRelationship)geneFeatures.get(i);

      if(gene.getFeatureByObjectId().getCvTerm().getName().equals("gene") ||
         gene.getFeatureByObjectId().getCvTerm().getName().equals("pseudogene"))
      {
        Integer subjectId = new Integer(gene.getFeatureBySubjectId().getFeatureId());
        
        if(genes.containsValue(subjectId))
        {
          Enumeration keys = genes.keys();
          while(keys.hasMoreElements())
          {
            Integer key = (Integer)keys.nextElement();
            Object val = genes.get(key);
            if(val instanceof Integer && subjectId.equals(val))
              genes.put(key, gene.getFeatureByObjectId().getUniqueName());
          }
        }
      }
    }
    
    Enumeration keys = genes.keys();
    while(keys.hasMoreElements())
    {
      Integer key = (Integer)keys.nextElement();
      Object val = genes.get(key);
      if(val instanceof String)
      {
        System.out.println(key.intValue()+"  "+val);
      }
    }
    */
    
    for(int i=0;i<allClusters.size(); i++)
    {
      final Feature clusterFeature = (Feature)allClusters.get(i);
      ClusterLazyQualifierValue lazyValue = 
        (ClusterLazyQualifierValue)hash.get(new Integer(clusterFeature.getFeatureId()));
      lazyValue.addToCluster(clusterFeature);
    }
  }
  
  private void addToCluster(final Feature clusterFeature)
  {
    if(clusters == null)
      clusters = new Vector();
    clusters.add(clusterFeature);
  }
  
  private void initCluster()
  {
    clusters = new Vector();
  }
  
  private String getHardString()
  {
    lazyLoaded = true;
    
    final String featureId = (String) feature.getQualifierByName("feature_id").getValues().get(0);
     
    StringVector strings = StringVector.getStrings(value, ";");   
    String rank = ArtemisUtils.getString(strings, "rank");
    
    //
    // should already have the clusters loaded by calling setClusterFromValueList()
    if(clusters == null)
    {
      final List featureIds = new Vector();
      String f_id[] = ArtemisUtils.getString(strings, "object_id=").split("=");
      featureIds.add( Integer.valueOf(f_id[1]) );
      final Document document = ((DocumentEntry)feature.getEntry()).getDocument();
      clusters = ((DatabaseDocument)document).getClustersByFeatureIds(featureIds);
      for(int i=0;i<clusters.size(); i++)
        System.out.println("*********NOT PRELOADED "+
            ((Feature)clusters.get(i)).getUniqueName() );
    }

    value = "";
    for(int i=0;i<clusters.size(); i++)
    {
      final Feature clusterFeature = (Feature)clusters.get(i);
      
      if(clusterFeature.getCvTerm().getName().indexOf("match") < 0)
        continue;

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
          if(!subjectFeature.getCvTerm().getName().equals("gene") ||
             !subjectFeature.getCvTerm().getName().equals("pseudogene"))
          {
            Feature parent = getParentFeature(subjectFeature);

            if(parent.getCvTerm().getName().equals("gene") ||
               parent.getCvTerm().getName().equals("pseudogene"))
              geneName = parent.getUniqueName();
            else if(parent != null)
            {
              parent = getParentFeature(parent);
              if(parent.getCvTerm().getName().equals("gene") ||
                 parent.getCvTerm().getName().equals("pseudogene"))
                geneName = parent.getUniqueName();
            }
          }
          
          value = value.concat(geneName+" link="+
              subjectFeature.getUniqueName());
          
          value = value.concat(" type="+fr.getCvTerm().getName());
          
          cnt++;
        }
      }
      if(cnt > 1)
        value = value.concat("; cluster_name="+clusterFeature.getUniqueName());
      else
        value = value.concat("; match_name="+clusterFeature.getUniqueName());
    }
    if(value.equals(""))
      return value;
    
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


  public String getValue()
  {
    return value;
  }
  

}