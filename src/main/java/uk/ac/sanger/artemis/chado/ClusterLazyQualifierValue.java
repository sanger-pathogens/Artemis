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

import org.gmod.schema.analysis.AnalysisFeature;
import org.gmod.schema.cv.CvTerm;
import org.gmod.schema.sequence.Feature;
import org.gmod.schema.sequence.FeatureCvTerm;
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
  /** include gene name */
  private boolean loadGeneName = true;
  
  private String value;
  private String name;
  private GFFStreamFeature feature;
  private List clusters;
  
  /**
   * Qualifier object to handle lazy loading of cluster/ortholog/paralog data
   * @param matchFeature
   * @param featureId
   */
  public ClusterLazyQualifierValue(final String value, 
                                   final String name,
                                   final GFFStreamFeature feature)
  {
    this.value = value;
    this.name = name;
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
   * This speeds up loading for a list of ortho/paralogs/clusters
   * @param values List of ClusterLazyQualifierValue
   * @param feature
   */
  public static void setClusterFromValueList(final List values, final DatabaseDocument document)
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
      
      final Vector v;
      if(hash.containsKey(clusterFeatureId))
        v = (Vector)hash.get(clusterFeatureId);
      else
        v = new Vector();
        
      v.add(lazyValue);
      hash.put(clusterFeatureId, v);
    }
    
    //final Document document = ((DocumentEntry)feature.getEntry()).getDocument();
    List allClusters = ((DatabaseDocument)document).getClustersByFeatureIds(clusterFeatureIds);
    List<Integer> featureIds = new Vector<Integer>();
    for(int i=0;i<allClusters.size(); i++)
    {
      final Feature clusterFeature = (Feature)allClusters.get(i);
      Collection subjects = clusterFeature.getFeatureRelationshipsForSubjectId();
      Iterator it = subjects.iterator();
      while(it.hasNext())
      {
        FeatureRelationship fr = (FeatureRelationship)it.next();
        featureIds.add(fr.getFeatureBySubjectId().getFeatureId());
      }
    }
    
    // bulk load parent features
    List<FeatureRelationship> parentFrs = 
      ((DatabaseDocument)document).getParentFeaturesByChildFeatureIds(featureIds);
    featureIds.clear();
    Hashtable<Integer, FeatureRelationship> hashFR = new Hashtable<Integer, FeatureRelationship>();
    for(int i=0; i<parentFrs.size(); i++)
    {
      FeatureRelationship frBulk = parentFrs.get(i);
      featureIds.add(frBulk.getFeatureByObjectId().getFeatureId());
      hashFR.put(frBulk.getFeatureBySubjectId().getFeatureId(), frBulk);
    }
    
    // bulk load grand parent features (to find gene)
    parentFrs = 
      ((DatabaseDocument)document).getParentFeaturesByChildFeatureIds(featureIds);
    for(int i=0; i<parentFrs.size(); i++)
    {
      FeatureRelationship frBulk = parentFrs.get(i);
      hashFR.put(frBulk.getFeatureBySubjectId().getFeatureId(), frBulk);
    }
    
    for(int i=0;i<allClusters.size(); i++)
    {
      final Feature clusterFeature = (Feature)allClusters.get(i);
      Collection<FeatureRelationship> sbjts = clusterFeature.getFeatureRelationshipsForSubjectId();
      List<FeatureRelationship> newSubject = new Vector<FeatureRelationship>();
      Iterator<FeatureRelationship> it = sbjts.iterator();
      while(it.hasNext())
      {
        FeatureRelationship fr = it.next();
        FeatureRelationship frBulk = hashFR.get(fr.getFeatureBySubjectId().getFeatureId());
        if(frBulk == null)
          continue;
        List<FeatureRelationship> frsForSbjtId = new Vector<FeatureRelationship>();

        Feature parent = frBulk.getFeatureByObjectId();
        if(!parent.getCvTerm().getName().equals("gene") ||
           !parent.getCvTerm().getName().equals("pseudogene") &&
           hashFR.contains(frBulk.getFeatureByObjectId().getFeatureId()))
        {
          frBulk = hashFR.get(frBulk.getFeatureByObjectId().getFeatureId());
        }
        frsForSbjtId.add(frBulk);
        
        fr.getFeatureBySubjectId().setFeatureRelationshipsForSubjectId(frsForSbjtId);
        newSubject.add(fr);
        
        //System.out.println(fr.getFeatureBySubjectId().getUniqueName()+" "+
        //                   frBulk.getFeatureByObjectId().getUniqueName());
      }
      clusterFeature.setFeatureRelationshipsForSubjectId(newSubject);
      
      final Vector v = (Vector)hash.get(new Integer(clusterFeature.getFeatureId()));
      for(int j=0; j<v.size(); j++)
      {
        ClusterLazyQualifierValue lazyValue = (ClusterLazyQualifierValue)v.get(j);
        lazyValue.addToCluster(clusterFeature);
      }
    }

    clusterFeatureIds.clear();
    featureIds.clear();
    parentFrs.clear();
    hashFR.clear();
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
  
  private synchronized String getHardString()
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
      
      if(f_id.length < 2)
        return value;
      
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
      {
        // this looks like an ortho/paralog stored without a match feature
        // and just a feature relationship between features
        if(clusterFeature.getCvTerm().getName().equals("polypeptide") ||
           clusterFeature.getCvTerm().getName().equals("gene"))
        {
          final DatabaseDocument document = 
            (DatabaseDocument)feature.getDocumentEntry().getDocument();
          Feature matchFeature = 
            document.getFeatureByUniquename(clusterFeature.getUniqueName());

          value = value.concat(matchFeature.getOrganism().getCommonName()+":");
          if(loadGeneName)
          {          
            String geneName = getGeneName(matchFeature);
            value = value.concat(geneName+" ");
          }
  
          value = value.concat("link="+matchFeature.getUniqueName());
          value = value.concat(" type="+name);
  
          String product = getProduct(matchFeature);
          if(product != null)
            value = value.concat(";product="+product);
          
          value = value.concat(";match_name="+matchFeature.getUniqueName());
        }
        
        continue;
      }

      final Collection subjects = clusterFeature.getFeatureRelationshipsForSubjectId();
      Iterator it = subjects.iterator();
      int cnt = 0;
      
      while(it.hasNext())
      {
        FeatureRelationship fr = (FeatureRelationship)it.next();
        Feature subjectFeature = fr.getFeatureBySubjectId();
        
        if(subjectFeature.getFeatureId() != Integer.parseInt(featureId))
          cnt++;
      }
      
      it = subjects.iterator();
      while(it.hasNext())
      {
        FeatureRelationship fr = (FeatureRelationship)it.next();
        Feature subjectFeature = fr.getFeatureBySubjectId();
        
        if(subjectFeature.getFeatureId() != Integer.parseInt(featureId))
        { 
          if(!value.equals(""))
            value = value.concat(", ");
          
          value = value.concat(subjectFeature.getOrganism().getCommonName()+":");
          
          if(loadGeneName)
          {
            try
            {
              String geneName = getGeneName(subjectFeature);
              value = value.concat(geneName+" ");
            }
            catch(NullPointerException npe)
            {
              System.err.println("Cannot get the gene name of "+
                    subjectFeature.getUniqueName());
            }
          }
          
          value = value.concat("link="+
              subjectFeature.getUniqueName());
          
          value = value.concat(" type="+fr.getCvTerm().getName());
          
          // if not a cluster
          if(cnt < 2 &&
             fr.getCvTerm().getName().equals("orthologous_to"))
          {
            String product = getProduct(subjectFeature);
            if(product != null)
              value = value.concat("; product="+product);
          }
        }
      }
      if(cnt > 1 || 
         (clusterFeature.getUniqueName().startsWith("CLUSTER_") && cnt > 0))
      {
        value = value.concat("; cluster_name="+clusterFeature.getUniqueName());
        
        Collection analysisFeatures = clusterFeature.getAnalysisFeatures();
        Iterator itAnalysis = analysisFeatures.iterator();
        if(itAnalysis.hasNext())
        {
          AnalysisFeature analysisFeature = (AnalysisFeature) itAnalysis.next();
          value = value.concat("; program="+analysisFeature.getAnalysis().getProgram());
        }
        //else
        //  value = value.concat("; program=unknown");
      }
      else if(cnt > 0)
        value = value.concat("; match_name="+clusterFeature.getUniqueName());
    }
    
    if(value.equals(""))
      return value;
    
    value = value.concat("; "+rank);

    //
    //
    loadGeneName = false;

    return value;
  }
  
  private String getProduct(final Feature subjectFeature)
  {
    Collection featureCvTerms = subjectFeature.getFeatureCvTerms();
    
    if(featureCvTerms != null)
    {
      Iterator itFCT = featureCvTerms.iterator();
      while(itFCT.hasNext())
      {
        FeatureCvTerm featureCvTerm = (FeatureCvTerm) itFCT.next();
        CvTerm cvTerm = featureCvTerm.getCvTerm();
        if(cvTerm.getCv().getName().equals(
            uk.ac.sanger.artemis.chado.ChadoTransactionManager.PRODUCT_CV))
          return featureCvTerm.getCvTerm().getName();
      }
    }
    return null;
  }
  
  private String getGeneName(final Feature subjectFeature)
  {
    String geneName = subjectFeature.getUniqueName();
    if(!subjectFeature.getCvTerm().getName().equals("gene") ||
       !subjectFeature.getCvTerm().getName().equals("pseudogene"))
    {
      Feature parent = getParentFeature(subjectFeature);  
      try
      {
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
      catch(NullPointerException npe)
      {
        System.err.println(geneName+" parent not found"); 
      }
    }
    return geneName;
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

  public void setLoadGeneName(boolean loadGeneName)
  {
    this.loadGeneName = loadGeneName;
  }
}