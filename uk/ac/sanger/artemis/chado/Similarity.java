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
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;

import org.gmod.schema.analysis.AnalysisFeature;
import org.gmod.schema.sequence.FeatureDbXRef;
import org.gmod.schema.sequence.FeatureProp;

import uk.ac.sanger.artemis.io.GFFStreamFeature;

public class Similarity
{
  public static String getSimilarityString(final GFFStreamFeature gffFeature)
  {
    HashSet similarities = gffFeature.getSimilarityFeatures();
    StringBuffer buff = new StringBuffer();
    
    if(similarities != null)
    {
      int f_id = 0;
      try
      {
        f_id = Integer.parseInt((String)gffFeature.getQualifierByName(
            "feature_id").getValues().get(0));
      }
      catch(Exception e)
      {
        e.printStackTrace();
      }

      Iterator it = similarities.iterator();
      while(it.hasNext())
      {
        org.gmod.schema.sequence.Feature matchFeature = 
          (org.gmod.schema.sequence.Feature)it.next();

        
        Collection featureLocs = matchFeature.getFeatureLocsForFeatureId();
        Iterator it2 = featureLocs.iterator();

        Collection analysisFeatures = matchFeature.getAnalysisFeatures();
        Iterator it3 = analysisFeatures.iterator();
        AnalysisFeature analysisFeature = (AnalysisFeature) it3.next();

        buff.append("/similarity=" + analysisFeature.getAnalysis().getProgram()
            + "; ");

        org.gmod.schema.sequence.Feature subject = null;
        while(it2.hasNext())
        {
          org.gmod.schema.sequence.FeatureLoc featureLoc = 
            (org.gmod.schema.sequence.FeatureLoc)it2.next();

          org.gmod.schema.sequence.Feature queryOrSubject = featureLoc
              .getFeatureBySrcFeatureId();

          if(queryOrSubject.getFeatureId() != f_id)
            subject = queryOrSubject;
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

        List featureProps = (List)subject.getFeatureProps();
        
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
        
        if(matchFeature.getFeatureProps() != null)
        {
          featureProps = (List)matchFeature.getFeatureProps();
          for(int i=0; i<featureProps.size(); i++)
          {
            FeatureProp featureProp = (FeatureProp)featureProps.get(i);
            buff.append(featureProp.getCvTerm().getName()+"="+featureProp.getValue());
            if(i<featureProps.size()-1)
              buff.append("; ");
          }
        }
        
        if(it.hasNext())
          buff.append("\n");
      }
    }
    
    //uk.ac.sanger.artemis.components.Splash.logger4j.debug(buff);
    return new String(buff);
  }
  
}