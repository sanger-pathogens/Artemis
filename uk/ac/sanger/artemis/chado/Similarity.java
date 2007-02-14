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
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Vector;

import org.gmod.schema.analysis.AnalysisFeature;
import org.gmod.schema.sequence.FeatureDbXRef;
import org.gmod.schema.sequence.FeatureProp;
import org.gmod.schema.sequence.Feature;

import uk.ac.sanger.artemis.io.LazyQualifierValue;

public class Similarity implements LazyQualifierValue
{
  /** match feature associated with the similarity */
  private Feature matchFeature;
  /** feature_id of the query feature */
  private int featureId;
  /** force complete loading of the data */
  private boolean forceLoad = false;
  
  /**
   * Qualifier object to handle lazy loading of similarity data
   * @param matchFeature
   * @param featureId
   */
  public Similarity(final Feature matchFeature, final int featureId)
  {
    this.matchFeature = matchFeature;
    this.featureId    = featureId;
  }

  /**
   * Handle the loading of the data into a String 
   */
  public String getString()
  {
    if(forceLoad)
      return getHardString();
    else
      return getSoftString();
  }
  
  public String getHardString()
  {
    StringBuffer buff = new StringBuffer();
    
    Collection featureLocs = matchFeature.getFeatureLocsForFeatureId();
    Iterator it2 = featureLocs.iterator();

    Collection analysisFeatures = matchFeature.getAnalysisFeatures();
    Iterator it3 = analysisFeatures.iterator();
    AnalysisFeature analysisFeature = (AnalysisFeature) it3.next();

    buff.append(analysisFeature.getAnalysis().getProgram()+";");

    org.gmod.schema.sequence.Feature subject = null;
    org.gmod.schema.sequence.FeatureLoc queryLoc   = null;
    org.gmod.schema.sequence.FeatureLoc subjectLoc   = null;
    
    while(it2.hasNext())
    {
      org.gmod.schema.sequence.FeatureLoc featureLoc = 
        (org.gmod.schema.sequence.FeatureLoc)it2.next();

      org.gmod.schema.sequence.Feature queryOrSubject = featureLoc
          .getFeatureBySrcFeatureId();

      if(queryOrSubject.getFeatureId() != featureId)
      {
        subject = queryOrSubject;
        subjectLoc = featureLoc;
      }
      else 
      {
        queryLoc = featureLoc;
      }
    }

    if(subject.getDbXRef() != null)
    {
      buff.append(subject.getDbXRef().getDb().getName() + ":");
      buff.append(subject.getDbXRef().getAccession());
    }

    Collection dbXRefs = subject.getFeatureDbXRefs();
    if(dbXRefs != null && dbXRefs.size() > 0)
    {
      buff.append(" (");
      Iterator it4 = dbXRefs.iterator();
      while(it4.hasNext())
      {
        FeatureDbXRef featureDbXRef = (FeatureDbXRef) it4.next();
        buff.append(featureDbXRef.getDbXRef().getDb().getName() + ":");
        buff.append(featureDbXRef.getDbXRef().getAccession());
        if(it4.hasNext())
          buff.append(",");
      }
      buff.append(")");
    }
    buff.append("; ");

    List featureProps = new Vector(subject.getFeatureProps());
    Collections.sort(featureProps, new FeaturePropComparator());
    
    for(int i=0; i<featureProps.size(); i++)
    {
      FeatureProp featureProp = (FeatureProp)featureProps.get(i);
      
      if(featureProp.getValue() != null)
        buff.append(featureProp.getValue().trim());
      buff.append("; ");
    }

    buff.append("length "+subject.getSeqLen());
    
    if(matchFeature.getCvTerm().getName().equals("protein_match"))
      buff.append(" aa; ");
    else
      buff.append("; ");
    
    if(analysisFeature.getIdentity() != null)
      buff.append("id="+analysisFeature.getIdentity()+"%; ");
    if(analysisFeature.getSignificance() != null)
      buff.append("E()="+analysisFeature.getSignificance()+"; ");
    if(analysisFeature.getRawScore() != null)
      buff.append("score="+analysisFeature.getRawScore()+"; ");
    
    if(queryLoc != null)
    {
      int fmin = queryLoc.getFmin().intValue()+1;
      buff.append("query "+fmin+"-"+queryLoc.getFmax());
      if(matchFeature.getCvTerm().getName().equals("protein_match"))
        buff.append(" aa; ");
      else
        buff.append("; ");
    }
    
    if(subjectLoc != null)
    {
      int fmin = subjectLoc.getFmin().intValue()+1;
      buff.append("subject "+fmin+"-"+subjectLoc.getFmax());
      if(matchFeature.getCvTerm().getName().equals("protein_match"))
        buff.append(" aa; ");
      else
        buff.append("; ");
    }
    
    if(matchFeature.getFeatureProps() != null)
    {
      featureProps = new Vector(matchFeature.getFeatureProps());
      Collections.sort(featureProps, new FeaturePropComparator());
      
      for(int i=0; i<featureProps.size(); i++)
      {
        FeatureProp featureProp = (FeatureProp)featureProps.get(i);
        buff.append(featureProp.getCvTerm().getName()+"="+featureProp.getValue());
        if(i < featureProps.size()-1)
          buff.append("; ");
      }
    }
    
    return new String(buff);
  }
  
  
  public String getSoftString()
  {
    return new String("LAZY LOADING...;");
  }

  public boolean isForceLoad()
  {
    return forceLoad;
  }

  public void setForceLoad(boolean forceLoad)
  {
    this.forceLoad = forceLoad;
  }


  class FeaturePropComparator
        implements Comparator
  {

    public int compare(Object o1, Object o2)
    {
      int rank1 = ((FeatureProp)o1).getRank();
      int rank2 = ((FeatureProp)o2).getRank();
      
      return rank1-rank2;
    }

  }


  public Feature getMatchFeature()
  {
    return matchFeature;
  }

}
