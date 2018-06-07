/* FeatureLocLazyQualifierValue
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
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Vector;

import org.gmod.schema.analysis.AnalysisFeature;
import org.gmod.schema.general.DbXRef;
import org.gmod.schema.sequence.FeatureDbXRef;
import org.gmod.schema.sequence.FeatureLoc;
import org.gmod.schema.sequence.FeatureProp;
import org.gmod.schema.sequence.Feature;

import uk.ac.sanger.artemis.io.LazyQualifierValue;
import uk.ac.sanger.artemis.util.DatabaseDocument;

public class FeatureLocLazyQualifierValue implements LazyQualifierValue
{
  /** match feature associated with the similarity */
  private Feature matchFeature;
  /** feature_id of the query feature */
  private int featureId;
  /** force complete loading of the data */
  private boolean forceLoad = false;
  /** data loaded */
  private boolean lazyLoaded = false;
  
  /**
   * Qualifier object to handle lazy loading of properties that
   * are featureloc'ed to the feature being read in. e.g. similarity,
   * polypeptide_domain, protein predictions (TMHMM, signal_peptide).
   * @param matchFeature
   * @param featureId
   */
  public FeatureLocLazyQualifierValue(final Feature matchFeature, final int featureId)
  {
    this.matchFeature = matchFeature;
    this.featureId    = featureId;
  }

  /**
   * Bulk retrieval of lazy properties (used to speed up writing to files)
   * @param similarity  a <code>List</code> of Similarity qualifier values
   * @param doc         the Document to which these features belong
   */
  public static void bulkRetrieve(final List similarity,
                                  final DatabaseDocument doc)
  {
    final Iterator it = similarity.iterator();
    final Hashtable featureLocHash = new Hashtable(similarity.size()*2);
    final Hashtable matchFeatures = new Hashtable(similarity.size());
    
    while(it.hasNext())
    {
      FeatureLocLazyQualifierValue thisSimilarity = (FeatureLocLazyQualifierValue)it.next();
      Feature thisMatchFeature = thisSimilarity.getMatchFeature();

      Collection featureLocs = thisMatchFeature.getFeatureLocsForFeatureId();
      Iterator it2 = featureLocs.iterator();
      
      while(it2.hasNext())
      {
        org.gmod.schema.sequence.FeatureLoc featureLoc = 
          (org.gmod.schema.sequence.FeatureLoc)it2.next();

        if(featureLoc.getSrcFeatureId() <= 0)
          continue;
        
        final Integer srcFeatureId = new Integer(featureLoc.getSrcFeatureId());
        List locs;
        if(featureLocHash.containsKey(srcFeatureId))
          locs = (Vector)featureLocHash.get(srcFeatureId);
        else
          locs = new Vector();
        
        locs.add(featureLoc);
        featureLocHash.put(srcFeatureId, locs);
      }
     
      matchFeatures.put(new Integer(thisMatchFeature.getFeatureId()), thisMatchFeature);
    }
    
    final List queryAndSubjectFeatureIds = new Vector(featureLocHash.keySet());
    //
    // bulk load the subject and query features
    //
    final List sims = doc.getFeaturesByListOfIds(queryAndSubjectFeatureIds);

    for(int i=0; i<sims.size();i++)
    {
      Feature srcFeature = (Feature)sims.get(i);
      Integer srcFeatureId = new Integer(srcFeature.getFeatureId());
      if(featureLocHash.containsKey(srcFeatureId))
      {
        Vector locs = (Vector)featureLocHash.get(srcFeatureId);
        for(int j=0;j<locs.size();j++)
        {
          FeatureLoc featureLoc = (FeatureLoc)locs.get(j);
          featureLoc.setFeatureBySrcFeatureId(srcFeature);
        }
      }
    }
    sims.clear();
    
    //
    // bulk load subject / query feature_dbxref's
    //
    final List featureDbXRefs = doc.getFeatureDbXRefsByFeatureId(queryAndSubjectFeatureIds);
    
    for(int i=0;i<featureDbXRefs.size();i++)
    {
      Feature srcFeature = (Feature)featureDbXRefs.get(i);
      Integer srcFeatureId = new Integer(srcFeature.getFeatureId());
      if(featureLocHash.containsKey(srcFeatureId))
      {
        Vector locs = (Vector)featureLocHash.get(srcFeatureId);
        for(int j=0;j<locs.size();j++)
        {
          FeatureLoc featureLoc = (FeatureLoc)locs.get(j);
          featureLoc.getFeatureBySrcFeatureId().setFeatureDbXRefs(
              srcFeature.getFeatureDbXRefs());
        }
      }
    }
    featureDbXRefs.clear();
    
    //
    // bulk load the match feature properties
    //
    final List matchFeaturesWithProps = doc.getFeaturePropByFeatureIds(
                                        new Vector(matchFeatures.keySet()) );
    
    for(int i=0; i<matchFeaturesWithProps.size(); i++)
    {
      Feature thisMatch = (Feature)matchFeaturesWithProps.get(i);
      Integer featureId = new Integer(thisMatch.getFeatureId());
      if(matchFeatures.containsKey(featureId))
      {
        Feature storedMatch = ((Feature)matchFeatures.get(featureId));
        storedMatch.setFeatureProps(thisMatch.getFeatureProps());
        //storedMatch.setDbXRef(thisMatch.getDbXRef());
      }
    }
    matchFeaturesWithProps.clear();
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
  
  /**
   * This returns the completed value, loading any lazy properties
   * @return
   */
  private String getHardString()
  {
    final StringBuffer buff = new StringBuffer();
    
    Collection featureLocs = matchFeature.getFeatureLocsForFeatureId();
    Iterator it2 = featureLocs.iterator();

    Collection analysisFeatures = matchFeature.getAnalysisFeatures();
    Iterator it3 = analysisFeatures.iterator();
    AnalysisFeature analysisFeature = null;
    
    if(it3.hasNext())  // attached analysisfeature 
    {
      analysisFeature = (AnalysisFeature) it3.next();
      buff.append(analysisFeature.getAnalysis().getProgram()+";");
    }
    
    if(analysisFeature == null ||
       matchFeature.getCvTerm().getName().equals("polypeptide_domain"))
    {
      // predictions and polypeptide_domains have dbxrefs
      buff.append(getMatchFeatureDbXRefs());
    }
    

    org.gmod.schema.sequence.Feature subject = null;
    org.gmod.schema.sequence.FeatureLoc queryLoc   = null;
    org.gmod.schema.sequence.FeatureLoc subjectLoc = null;
    
    while(it2.hasNext())
    {
      org.gmod.schema.sequence.FeatureLoc featureLoc = 
        (org.gmod.schema.sequence.FeatureLoc)it2.next();

      if(featureLoc.getSrcFeatureId() <= 0)
        continue;
      
      org.gmod.schema.sequence.Feature queryOrSubject = 
        featureLoc.getFeatureBySrcFeatureId();

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

    if(subject != null)
    {
      if(subject.getDbXRef() != null)
      {
        buff.append(subject.getDbXRef().getDb().getName() + ":");
        buff.append(subject.getDbXRef().getAccession());
      }

      Collection dbXRefs = subject.getFeatureDbXRefs();
      
      if(dbXRefs != null && dbXRefs.size() > 0)
      {
        final StringBuffer buffDbXRefs = new StringBuffer();
        Iterator it4 = dbXRefs.iterator();
        while(it4.hasNext())
        {
          FeatureDbXRef featureDbXRef = (FeatureDbXRef) it4.next();
          featureDbXRef.getDbXRef();
          try
          {
            buffDbXRefs.append(featureDbXRef.getDbXRef().getDb().getName() + ":");
            buffDbXRefs.append(featureDbXRef.getDbXRef().getAccession());
            if(it4.hasNext())
              buffDbXRefs.append(",");
          }
          catch(NullPointerException npe){}
        }
        
        if(buffDbXRefs.length() > 0)
        {
          buff.append(" (");
          buff.append(buffDbXRefs);
          buff.append(")");
        }
      }
      buff.append(";");

      List featureProps = new Vector(subject.getFeatureProps());
      Collections.sort(featureProps, new FeaturePropComparator());

      for(int i = 0; i < featureProps.size(); i++)
      {
        FeatureProp featureProp = (FeatureProp) featureProps.get(i);

        if(featureProp.getValue() != null)
          buff.append(featureProp.getValue().trim());
        buff.append(";");
      }

      buff.append("length " + subject.getSeqLen());
    }
    
    if(matchFeature.getCvTerm().getName().equals("protein_match"))
      buff.append(" aa; ");
    else
      buff.append(";");
    
    if(analysisFeature != null && analysisFeature.getIdentity() != null)
      buff.append("id="+analysisFeature.getIdentity()+"%;");
    if(analysisFeature != null && analysisFeature.getSignificance() != null)
      buff.append("E()="+analysisFeature.getSignificance()+";");
    if(analysisFeature != null && analysisFeature.getRawScore() != null)
      buff.append("score="+analysisFeature.getRawScore()+";");
    
    if(queryLoc != null && queryLoc.getFmin().intValue() > -1)
    {
      final int fmin;
      if(queryLoc.getFmin().compareTo(queryLoc.getFmax()) == 0)
        fmin = queryLoc.getFmin().intValue();
      else
        fmin = queryLoc.getFmin().intValue()+1;
      buff.append("query "+fmin+"-"+queryLoc.getFmax());
      if(matchFeature.getCvTerm().getName().equals("protein_match"))
        buff.append(" aa;");
      else
        buff.append(";");
    }
    
    if(subjectLoc != null && subjectLoc.getFmin().intValue() > -1)
    {
      int fmin = subjectLoc.getFmin().intValue()+1;
      buff.append("subject "+fmin+"-"+subjectLoc.getFmax());
      if(matchFeature.getCvTerm().getName().equals("protein_match"))
        buff.append(" aa;");
      else
        buff.append(";");
    }
    
    if(matchFeature.getFeatureProps() != null)
    {
      List featureProps = new Vector(matchFeature.getFeatureProps());
      Collections.sort(featureProps, new FeaturePropComparator());
      
      for(int i=0; i<featureProps.size(); i++)
      {
        FeatureProp featureProp = (FeatureProp)featureProps.get(i);
        
        final String cvTermName;
        if(featureProp.getCvTerm().getName() == null ||
           featureProp.getCvTerm().getName().equals("null"))
          cvTermName = DatabaseDocument.getCvTermByCvTermId(
              featureProp.getCvTerm().getCvTermId(), null).getName();
        else
          cvTermName = featureProp.getCvTerm().getName();

        buff.append(cvTermName+"="+featureProp.getValue());
        if(i < featureProps.size()-1)
          buff.append(";");
      }
    }
    
    lazyLoaded = true;
    return new String(buff);
  }
  
  /**
   * Get dbxrefs associated with the match feature
   * @return
   */
  private String getMatchFeatureDbXRefs()
  {
    final StringBuffer dbXRefs = new StringBuffer();
    final Collection featureDbXRefs = matchFeature.getFeatureDbXRefs();
    final Iterator it3 = featureDbXRefs.iterator();
    if(it3 == null)
      return "";
    
    while(it3.hasNext())
    {
      DbXRef dbXRef = ((FeatureDbXRef) it3.next()).getDbXRef();
      dbXRefs.append(dbXRef.getDb().getName()+":"+
                               dbXRef.getAccession());
      
      if( dbXRef.getDescription() != null && 
         !dbXRef.getDescription().equals("") )
        dbXRefs.append(" :\t"+dbXRef.getDescription());
      dbXRefs.append(";");
    }
    
    final DbXRef dbXRef = matchFeature.getDbXRef();
    
    if(dbXRef != null)
    {
      dbXRefs.append(dbXRef.getDb().getName()+":"+
                   dbXRef.getAccession());
    
      /*if( dbXRef.getDescription() != null && 
         !dbXRef.getDescription().equals("") )
       dbXRefs.append(" :\t"+dbXRef.getDescription()+";");*/
      dbXRefs.append(";");
    }
    
    return dbXRefs.toString();
  }
  
  private String getSoftString()
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

  public Feature getMatchFeature()
  {
    return matchFeature;
  }

  public boolean isLazyLoaded()
  {
    return lazyLoaded;
  }
  
  class FeaturePropComparator implements Comparator
  {
    public int compare(Object o1, Object o2)
    {
      int rank1 = ((FeatureProp)o1).getRank();
      int rank2 = ((FeatureProp)o2).getRank();
      return rank1-rank2;
    }
  }

}
