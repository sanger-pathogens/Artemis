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

import org.gmod.schema.analysis.Analysis;
import org.gmod.schema.analysis.AnalysisFeature;
import org.gmod.schema.cv.CvTerm;
import org.gmod.schema.general.Db;
import org.gmod.schema.general.DbXRef;
import org.gmod.schema.sequence.FeatureDbXRef;
import org.gmod.schema.sequence.FeatureLoc;
import org.gmod.schema.sequence.FeatureProp;
import org.gmod.schema.sequence.Feature;

import uk.ac.sanger.artemis.components.Splash;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.LazyQualifierValue;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.StringVector;

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

  public static AnalysisFeature getAnalysisFeature(final String uniquename,
      final String qualifier_string, final GFFStreamFeature feature)
  {
    int queryFeatureId = Integer.parseInt((String) feature.getQualifierByName(
        "feature_id").getValues().get(0));

    AnalysisFeature analysisFeature = new AnalysisFeature();
    Analysis analysis = new Analysis();

    org.gmod.schema.sequence.Feature queryFeature = new org.gmod.schema.sequence.Feature();
    org.gmod.schema.sequence.Feature subjectFeature = new org.gmod.schema.sequence.Feature();
    org.gmod.schema.sequence.Feature matchFeature = new org.gmod.schema.sequence.Feature();
    FeatureDbXRef featureDbXRef = new FeatureDbXRef();

    subjectFeature.setCvTerm(getCvTerm("similarity_region"));

    queryFeature.setUniqueName(uniquename);
    queryFeature.setFeatureId(queryFeatureId);
    matchFeature.setUniqueName("MATCH_" + uniquename);

    analysisFeature.setAnalysis(analysis);
    analysisFeature.setFeature(matchFeature);

    List analysisFeatures = new Vector();
    analysisFeatures.add(analysisFeature);
    matchFeature.setAnalysisFeatures(analysisFeatures);

    // algorithm
    // StringTokenizer tok = new StringTokenizer(qualifier_string, ";");
    StringVector qualifier_strings = StringVector.getStrings(qualifier_string,
        ";");

    // String method = tok.nextToken();
    analysis.setProgram((String) qualifier_strings.get(0));

    // primary dbxref
    DbXRef dbXRef_1 = new DbXRef();
    Db db_1 = new Db();
    dbXRef_1.setDb(db_1);
    String value = (String) qualifier_strings.get(1);
    String values[] = value.split(" ");
    int ind = values[0].indexOf(':');
    final String primary_name = values[0].substring(0, ind);
    db_1.setName(primary_name);
    dbXRef_1.setAccession(values[0].substring(ind + 1));
    Splash.logger4j.debug("Primary dbXRef  " + db_1.getName() + ":"
        + dbXRef_1.getAccession());
    subjectFeature.setDbXRef(dbXRef_1);
    subjectFeature
        .setUniqueName(db_1.getName() + ":" + dbXRef_1.getAccession());

    if(primary_name.equalsIgnoreCase("UniProt"))
      matchFeature.setCvTerm(getCvTerm("protein_match"));
    else
      matchFeature.setCvTerm(getCvTerm("nucleotide_match"));

    // secondary dbxref
    if(values.length > 1)
    {
      DbXRef dbXRef_2 = new DbXRef();
      Db db_2 = new Db();
      dbXRef_2.setDb(db_2);

      values[1] = values[1].replaceAll("^\\W", "");
      values[1] = values[1].replaceAll("\\W$", "");

      ind = values[1].indexOf(':');
      db_2.setName(values[1].substring(0, ind));
      dbXRef_2.setAccession(values[1].substring(ind + 1));
      Splash.logger4j.debug("Secondary dbXRef  " + db_2.getName() + " "
          + dbXRef_2.getAccession());
      featureDbXRef.setDbXRef(dbXRef_2);
      featureDbXRef.setFeature(subjectFeature);
      List featureDbXRefs = new Vector();
      featureDbXRefs.add(featureDbXRef);
      subjectFeature.setFeatureDbXRefs(featureDbXRefs);
    }

    // organism
    final String organismStr = (String) qualifier_strings.get(2);
    if(!organismStr.equals(""))
    {
      FeatureProp featureProp = new FeatureProp();
      featureProp.setCvTerm(getCvTerm("organism"));
      featureProp.setValue(organismStr);
      featureProp.setRank(0);
      subjectFeature.addFeatureProp(featureProp);
    }

    // product
    final String product = (String) qualifier_strings.get(3);
    if(!product.equals(""))
    {
      FeatureProp featureProp = new FeatureProp();
      featureProp.setCvTerm(getCvTerm("product"));
      featureProp.setValue(product);
      featureProp.setRank(1);
      subjectFeature.addFeatureProp(featureProp);
    }

    // gene
    final String gene = (String) qualifier_strings.get(4);
    if(!gene.equals(""))
    {
      FeatureProp featureProp = new FeatureProp();
      featureProp.setCvTerm(getCvTerm("gene"));
      featureProp.setValue(gene);
      featureProp.setRank(2);
      subjectFeature.addFeatureProp(featureProp);
    }

    // length
    String length = getString(qualifier_strings, "length");
    if(!length.equals(""))
    {
      if(length.startsWith("length=") || length.startsWith("length "))
        length = length.substring(7);
      if(length.endsWith("aa"))
        length = length.substring(0, length.length() - 2).trim();
      subjectFeature.setSeqLen(new Integer(length));
    }

    // percentage identity
    String id = getString(qualifier_strings, "id");
    if(!id.equals(""))
    {
      if(id.startsWith("id="))
        id = id.substring(3);
      if(id.endsWith("%"))
        id = id.substring(0, id.length() - 1);
      analysisFeature.setIdentity(new Double(id));
    }

    // ungapped id
    String ungappedId = getString(qualifier_strings, "ungapped id=");
    if(!ungappedId.equals(""))
    {
      if(ungappedId.startsWith("ungapped id="))
        ungappedId = ungappedId.substring(12);
      if(ungappedId.endsWith("%"))
        ungappedId = ungappedId.substring(0, ungappedId.length() - 1);
      FeatureProp featureProp = new FeatureProp();
      featureProp.setCvTerm(getCvTerm("ungapped id"));
      featureProp.setValue(ungappedId);
      matchFeature.addFeatureProp(featureProp);
    }

    // e-value
    String evalue = getString(qualifier_strings, "E()=");
    if(!evalue.equals(""))
    {
      if(evalue.startsWith("E()="))
        evalue = evalue.substring(4);
      analysisFeature.setSignificance(new Double(evalue));
    }

    // score
    String score = getString(qualifier_strings, "score=");
    if(!score.equals(""))
    {
      if(score.startsWith("score="))
        score = score.substring(6);
      analysisFeature.setRawScore(new Double(score));
    }

    // overlap
    String overlap = getString(qualifier_strings, "overlap");
    if(!overlap.equals(""))
    {
      if(overlap.startsWith("overlap="))
        overlap = overlap.substring(8);

      FeatureProp featureProp = new FeatureProp();
      featureProp.setCvTerm(getCvTerm("overlap"));
      featureProp.setValue(overlap);
      matchFeature.addFeatureProp(featureProp);
    }

    Short strand;
    if(feature.getLocation().isComplement())
      strand = new Short("-1");
    else
      strand = new Short("1");

    // query location
    String queryLoc = getString(qualifier_strings, "query");
    if(!queryLoc.equals(""))
    {
      FeatureLoc featureLoc = new FeatureLoc();
      String locs[] = queryLoc.split(" ");
      locs = locs[1].split("-");
      int fmin = Integer.parseInt(locs[0]) - 1;
      featureLoc.setFmin(new Integer(fmin));
      int fmax = Integer.parseInt(locs[1]);
      featureLoc.setFmax(new Integer(fmax));
      featureLoc.setRank(1);
      featureLoc.setStrand(strand);
      featureLoc.setFeatureBySrcFeatureId(queryFeature);
      matchFeature.addFeatureLocsForFeatureId(featureLoc);
    }

    // subject location
    String subjectLoc = getString(qualifier_strings, "subject");
    if(!subjectLoc.equals(""))
    {
      FeatureLoc featureLoc = new FeatureLoc();
      String locs[] = subjectLoc.split(" ");
      locs = locs[1].split("-");
      int fmin = Integer.parseInt(locs[0]) - 1;
      featureLoc.setFmin(new Integer(fmin));
      int fmax = Integer.parseInt(locs[1]);
      featureLoc.setFmax(new Integer(fmax));
      featureLoc.setRank(0);
      featureLoc.setStrand(strand);
      featureLoc.setFeatureBySrcFeatureId(subjectFeature);
      matchFeature.addFeatureLocsForFeatureId(featureLoc);
    }

    //Similarity sim = new Similarity(matchFeature, queryFeatureId);
    //System.out.println(sim.getHardString());
    return analysisFeature;
  }

  private static String getString(final StringVector sv, final String name)
  {
    for(int i = 0; i < sv.size(); i++)
    {
      String value = (String) sv.get(i);
      if(value.trim().startsWith(name))
        return value.trim();
    }
    return "";
  }

  /**
   * Get CvTerm that have been cached
   * @param cvTermName
   * @return
   */
  private static CvTerm getCvTerm(String cvTermName)
  {
    if(cvTermName.startsWith("\""))
      cvTermName = cvTermName.substring(1, cvTermName.length() - 1);

    CvTerm cvTerm = DatabaseDocument.getCvTermByCvTermName(cvTermName);

    if(cvTerm != null)
    {
      Splash.logger4j.debug("USE CvTerm from cache, CvTermId=" + cvTermName
          + "  -> " + cvTerm.getCvTermId() + " " + cvTerm.getName() + ":"
          + cvTerm.getCv().getName());
    }
    else
    {
      Splash.logger4j.warn("CvTerm not found in cache = " + cvTermName);
      cvTerm = new CvTerm();
      cvTerm.setName(cvTermName);
    }
    return cvTerm;
  }

  public Feature getMatchFeature()
  {
    return matchFeature;
  }

}
