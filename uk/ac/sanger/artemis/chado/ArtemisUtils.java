/* ArtemisUtil
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
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Vector;

import org.gmod.schema.analysis.Analysis;
import org.gmod.schema.analysis.AnalysisFeature;
import org.gmod.schema.cv.CvTerm;
import org.gmod.schema.general.Db;
import org.gmod.schema.general.DbXRef;
import org.gmod.schema.organism.Organism;
import org.gmod.schema.sequence.Feature;
import org.gmod.schema.sequence.FeatureCvTerm;
import org.gmod.schema.sequence.FeatureCvTermDbXRef;
import org.gmod.schema.sequence.FeatureCvTermProp;
import org.gmod.schema.sequence.FeatureDbXRef;
import org.gmod.schema.sequence.FeatureLoc;
import org.gmod.schema.sequence.FeatureProp;
import org.gmod.schema.sequence.FeatureRelationship;

import uk.ac.sanger.artemis.io.DatabaseDocumentEntry;
import uk.ac.sanger.artemis.io.GFF3Encoder;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.StringVector;


public class ArtemisUtils
{
  private static org.apache.log4j.Logger logger4j = 
    org.apache.log4j.Logger.getLogger(ArtemisUtils.class);
  
  protected static String getCurrentSchema()
  {
    String schema = System.getProperty("chado");
    int index  = schema.indexOf("?");
    int index2 = schema.indexOf("user=");
    if(index2 < 0)
      schema = schema.substring(index+1);
    else
      schema = schema.substring(index2+5);
    return schema;
  }
  
  /**
   * Insert featureCvTerm with an appropriate rank
   * @param dao
   * @param featureCvTerm
   */
  public static void inserFeatureCvTerm(GmodDAO dao, FeatureCvTerm featureCvTerm)
  {
    List<FeatureCvTerm> featureCvTerms = dao.getFeatureCvTermsByFeature(featureCvTerm.getFeature());
    logger4j.debug("In inserFeatureCvTerm() inserting");
   
    boolean isNotProduct = true;
    try
    {
      // product feature_cvterm.rank is used to define alternative products
      // where there are multiple products, so the rank may already be set
      if(featureCvTerm.getCvTerm().getCv().getName().equals(
          DatabaseDocument.PRODUCTS_TAG_CVNAME))
        isNotProduct = false;
    } catch(Exception e){}
    
    if(isNotProduct || featureCvTerm.getRank() == 0)
    {
      int rank = 0;
      for(FeatureCvTerm this_feature_cvterm: featureCvTerms)
      {
        if(this_feature_cvterm.getCvTerm().getName().equals( 
           featureCvTerm.getCvTerm().getName() )  &&
           this_feature_cvterm.getCvTerm().getCv().getName().equals( 
               featureCvTerm.getCvTerm().getCv().getName() ))
        {
          rank++;
        }
      }
      featureCvTerm.setRank(rank);
    }
    dao.persist(featureCvTerm);
  }
  
  /**
   * Delete featureCvTerm and update associated feature_cvterm.rank's 
   * if appropriate
   * @param dao
   * @param featureCvTerm
   */
  public static void deleteFeatureCvTerm(GmodDAO dao, FeatureCvTerm featureCvTerm)
  {
    List featureCvTerms = dao.getFeatureCvTermsByFeature(featureCvTerm.getFeature());
   
    logger4j.debug("In deleteFeatureCvTerm() deleting one of "+featureCvTerms.size());
    
    List featureCvTermDbXRefs = new Vector();
    
    if(featureCvTerm.getFeatureCvTermDbXRefs() != null &&
       featureCvTerm.getFeatureCvTermDbXRefs().size() > 0)
      featureCvTermDbXRefs = (List)featureCvTerm.getFeatureCvTermDbXRefs();
    
    List featureCvTermProps = new Vector();
    
    if(featureCvTerm.getFeatureCvTermProps()!= null &&
        featureCvTerm.getFeatureCvTermProps().size() > 0)
      featureCvTermProps = (List)featureCvTerm.getFeatureCvTermProps();
    
    // delete feature_cvterm and update ranks if appropriate
    FeatureCvTerm deleteme = null;
    Vector rankable = null;
    
    logger4j.debug("In deleteFeatureCvTerm() looking to delete ");
    for(int i=0; i<featureCvTerms.size(); i++)
    {
      FeatureCvTerm this_feature_cvterm = (FeatureCvTerm)featureCvTerms.get(i);
      
      if(this_feature_cvterm.getCvTerm().getName().equals( 
         featureCvTerm.getCvTerm().getName() )  &&
         this_feature_cvterm.getCvTerm().getCv().getName().equals( 
             featureCvTerm.getCvTerm().getCv().getName() ))
      {
         logger4j.debug("Found CvTerm.name "+featureCvTerm.getCvTerm().getName());
         Collection this_featureCvTermDbXRefs = this_feature_cvterm.getFeatureCvTermDbXRefs();
         Collection this_featureCvTermProps   = this_feature_cvterm.getFeatureCvTermProps();
         
         if(this_featureCvTermDbXRefs == null)
           this_featureCvTermDbXRefs = new Vector();
         
         
         if(this_featureCvTermDbXRefs.size() != featureCvTermDbXRefs.size() ||
            featureCvTermProps.size() != this_featureCvTermProps.size())
         {
           if(rankable == null)
             rankable = new Vector();
           
           rankable.add(this_feature_cvterm);
           
           logger4j.debug("FeatureCvTermDbXRefs not the same - ignore "+
               this_featureCvTermDbXRefs.size()+" != "+featureCvTermDbXRefs.size() + " || "+
               featureCvTermProps.size()+" != "+this_featureCvTermProps.size());
           continue;
         }
         
         boolean found = true;
         Iterator it = this_featureCvTermDbXRefs.iterator();
         while(it.hasNext())
         {
           FeatureCvTermDbXRef fcd = (FeatureCvTermDbXRef)it.next();

           if(!containsFeatureCvTermDbXRef(fcd, featureCvTermDbXRefs))
           {
             found = false;
             break;
           }
         }
         
         it = this_featureCvTermProps.iterator();
         while(it.hasNext())
         {
           FeatureCvTermProp fcp   = (FeatureCvTermProp)it.next();
           if(!containsFeatureCvTermProp(fcp, featureCvTermProps))
           {
             logger4j.debug(fcp.getCvTerm().getName()+" "+fcp.getValue());
             
             found = false;
             break;
           }
         }
         
         if(!found)
         {
           if(rankable == null)
             rankable = new Vector();
           
           rankable.add(this_feature_cvterm);
           continue;
         }
         
         deleteme = this_feature_cvterm;
      }
    }
    if(deleteme != null)
      dao.delete(deleteme);
    else
      logger4j.debug("FeatureCvTerm not found");
    
    if(rankable != null)
    {
      // feature_cvterm.rank may need updating for those stored here
      for(int i=0; i<rankable.size(); i++)
      {
        FeatureCvTerm fc = (FeatureCvTerm)rankable.get(i);
        
        if(fc.getRank() == i)
          continue;
        
        logger4j.debug("UPDATE rank for "+ fc.getCvTerm().getCv().getName() + "   rank = " +
                              fc.getRank()+" -> "+i);
        fc.setRank(i);
        dao.merge(fc);
      }
    }
  }
  
  /**
   * Return true if the list contains a given feature_cvterm_dbxref
   * @param fcd
   * @param featureCvTermDbXRefs
   * @return
   */
  private static boolean containsFeatureCvTermDbXRef(FeatureCvTermDbXRef fcd, 
                                              List featureCvTermDbXRefs)
  {
    for(int i=0; i<featureCvTermDbXRefs.size(); i++)
    {
      FeatureCvTermDbXRef this_fcd = (FeatureCvTermDbXRef)featureCvTermDbXRefs.get(i);
      if( this_fcd.getDbXRef().getAccession().equals( fcd.getDbXRef().getAccession() ) &&
          this_fcd.getDbXRef().getDb().getName().equals( fcd.getDbXRef().getDb().getName() ))
        return true;
    }
    return false;
  }
  
  /**
   * Return true if the list contains a given feature_cvterm_prop
   * @param fcp
   * @param featureCvTermProps
   * @return
   */
  private static boolean containsFeatureCvTermProp(FeatureCvTermProp fcp, 
                                            List featureCvTermProps)
  {
    for(int i = 0; i < featureCvTermProps.size(); i++)
    {
      FeatureCvTermProp this_fcp = (FeatureCvTermProp)featureCvTermProps.get(i);
      if(this_fcp.getValue().equals(fcp.getValue()))
        return true;

      if(this_fcp.getValue().equals(GFF3Encoder.decode(fcp.getValue())))
        return true;
    }
    return false;
  }
  
  //
  //
  //
  
  
  /**
   * Build an AnalysisFeature object for a /similarity qualifier
   * @param uniqueName
   * @param qualifier_string
   * @param feature
   * @return
   */
  protected static AnalysisFeature getAnalysisFeature(String uniqueName,
               String qualifier_string, final GFFStreamFeature feature)
  {
    // - allow for case with no featureId!
    int queryFeatureId;
    if(feature.getQualifierByName("feature_id") == null)
      queryFeatureId = -1;
    else
      queryFeatureId = Integer.parseInt((String) feature.getQualifierByName(
                                          "feature_id").getValues().get(0));

    org.gmod.schema.sequence.Feature queryFeature = new org.gmod.schema.sequence.Feature();
    org.gmod.schema.sequence.Feature subjectFeature = new org.gmod.schema.sequence.Feature();
    org.gmod.schema.sequence.Feature matchFeature = new org.gmod.schema.sequence.Feature();
    FeatureDbXRef featureDbXRef = new FeatureDbXRef();
    
    final Hashtable rangeHash = feature.getSegmentRangeStore();
    if(rangeHash != null)
    {
      final Enumeration id_keys= rangeHash.keys();
      uniqueName = (String)id_keys.nextElement();
      
      DatabaseDocument doc;
      
      try
      {
        doc = (DatabaseDocument)feature.getDocumentEntry().getDocument();
      }
      catch(Exception e)
      {
        doc = (DatabaseDocument) (
            (DatabaseDocumentEntry)((uk.ac.sanger.artemis.Feature)feature.getUserData()).
            getEntry().getEMBLEntry()).getDocument();
      }
      FeatureLoc loc = new FeatureLoc();
      Feature srcFeature = new Feature();
      srcFeature.setFeatureId(Integer.parseInt(doc.getSrcFeatureId()));
      loc.setFeatureBySrcFeatureId(srcFeature);
      queryFeature.setFeatureLoc(loc);
      logger4j.debug(uniqueName);
    }
    
    AnalysisFeature analysisFeature = new AnalysisFeature();
    Analysis analysis = new Analysis();

    subjectFeature.setCvTerm(getCvTerm("region"));  // similarity_region

    queryFeature.setUniqueName(uniqueName);
    queryFeature.setFeatureId(queryFeatureId);
    matchFeature.setUniqueName("MATCH_" + uniqueName);

    analysisFeature.setAnalysis(analysis);
    analysisFeature.setFeature(matchFeature);

    List analysisFeatures = new Vector();
    analysisFeatures.add(analysisFeature);
    matchFeature.setAnalysisFeatures(analysisFeatures);

    // algorithm
    // StringTokenizer tok = new StringTokenizer(qualifier_string, ";");
    while(qualifier_string.startsWith("\""))
      qualifier_string = qualifier_string.substring(1);
    while(qualifier_string.endsWith("\""))
      qualifier_string = qualifier_string.substring(0,qualifier_string.length()-1);
 
    final StringVector qualifier_strings = StringVector.getStrings(qualifier_string,
        ";");

    analysis.setProgram((String) qualifier_strings.get(0));
    analysis.setProgramVersion("1.0");

    // primary dbxref
    DbXRef dbXRef_1 = new DbXRef();
    Db db_1 = new Db();
    dbXRef_1.setDb(db_1);
    String value = ((String) qualifier_strings.get(1)).trim();
    
    if(value.startsWith("with="))
      value = value.substring(5);
    String values[] = value.split(" ");

    int ind = values[0].indexOf(':');
    final String primary_name = values[0].substring(0, ind);
    db_1.setName(primary_name);
    dbXRef_1.setAccession(values[0].substring(ind + 1));
    logger4j.debug("Primary dbXRef  " + db_1.getName() + ":"
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
      logger4j.debug("Secondary dbXRef  " + db_2.getName() + " "
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
    if(qualifier_strings.size() > 4)
    {
      final String gene = (String) qualifier_strings.get(4);
      if(!gene.equals(""))
      {
        FeatureProp featureProp = new FeatureProp();
        featureProp.setCvTerm(getCvTerm("gene"));
        featureProp.setValue(gene);
        featureProp.setRank(2);
        subjectFeature.addFeatureProp(featureProp);
      }
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
      
      int index = id.indexOf(" ");
      if(index > -1)
        id = id.substring(index);
      
      analysisFeature.setIdentity(new Double(id));
    }

    // ungapped id
    String ungappedId = getString(qualifier_strings, "ungapped id");
    if(!ungappedId.equals(""))
    {
      if(ungappedId.startsWith("ungapped id"))
        ungappedId = ungappedId.substring(11);
      if(ungappedId.startsWith("="))
        ungappedId = ungappedId.substring(1);

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
    if(overlap.equals(""))
      overlap = getStringEndsWith(qualifier_strings, "overlap");
    
    if(!overlap.equals(""))
    {
      if(overlap.startsWith("overlap="))
        overlap = overlap.substring(8);

      FeatureProp featureProp = new FeatureProp();
      featureProp.setCvTerm(getCvTerm("overlap"));
      featureProp.setValue(overlap);
      matchFeature.addFeatureProp(featureProp);
    }

    final Short strand;
    if(feature.getLocation().isComplement())
      strand = new Short("-1");
    else
      strand = new Short("1");

    // query location
    String queryLoc = getString(qualifier_strings, "query");
    final FeatureLoc featureLoc = getFeatureLoc(queryLoc, queryFeature, strand, 0);
    matchFeature.addFeatureLocsForFeatureId(featureLoc);

    // subject location
    String subjectLoc = getString(qualifier_strings, "subject");
    final FeatureLoc subjectFeatureLoc = getFeatureLoc(subjectLoc, subjectFeature, strand, 1);
    matchFeature.addFeatureLocsForFeatureId(subjectFeatureLoc);

    return analysisFeature;
  }
  
  
  
  
  
  
  /**
   * Build an AnalysisFeature object for a orthologous_to/paralogous_to qualifier
   * @param uniqueName
   * @param qualifier_string
   * @param feature
   * @return
   */
  protected static Feature getMatchFeatureWithFeatureRelations(final String uniqueName,
      final String qualifierName, String qualifierString, final GFFStreamFeature feature,
      final boolean isDelete)
  {
    final int queryFeatureId;
    
    if(feature.getQualifierByName("feature_id") != null)
      queryFeatureId = Integer.parseInt((String) feature.getQualifierByName(
                                        "feature_id").getValues().get(0));
    else
      queryFeatureId = 0;

    Feature thisFeature     = new Feature();
    Feature orthologFeature = new Feature();
    Feature matchFeature    = new Feature();

    thisFeature.setUniqueName(uniqueName);
    thisFeature.setFeatureId(queryFeatureId);  

    while(qualifierString.startsWith("\""))
      qualifierString = qualifierString.substring(1);
    while(qualifierString.endsWith("\""))
      qualifierString = qualifierString.substring(0,qualifierString.length()-1);
 
    final StringVector qualifierStrings = StringVector.getStrings(qualifierString,
        ";");

    String ortholog = (String)qualifierStrings.get(0);
    String orthologStr[] = ortholog.split(":", 2);
    
    final Organism organism = new Organism();
    organism.setCommonName(orthologStr[0]);
    
    int index1 = orthologStr[1].indexOf("link=")+5;
    int index2 = orthologStr[1].indexOf("type=");
    String link = orthologStr[1].substring(index1, index2).trim();

    orthologFeature.setUniqueName(link);
    orthologFeature.setOrganism(organism);
    
    String matchName   = getString(qualifierStrings, "match_name");
    String clusterName = getString(qualifierStrings, "cluster_name");
    if(!matchName.equals(""))
    {
      // there is no match feature just a direct feature_relationship link
      if(link.equals(matchName.substring(11)))
      {
        if(isDelete)
          return null;
        else
          matchFeature.setUniqueName("ORTHO_PARA_" + uniqueName + "_" + link);
      }
      else
        matchFeature.setUniqueName(matchName.substring(11));
    }
    else if(!clusterName.equals(""))
      matchFeature.setUniqueName(clusterName.substring(13));
    else
      matchFeature.setUniqueName("ORTHO_PARA_" + uniqueName + "_" + link);
    matchFeature.setCvTerm(getCvTerm("protein_match"));
    matchFeature.setOrganism(organism);
    
    // rank
    String rank = getString(qualifierStrings, "rank");
    if(!rank.equals(""))
    {
      if(rank.startsWith("rank=") || rank.startsWith("rank "))
        rank = rank.substring(5);
    }
    
    FeatureRelationship thisFeatureRelationship = new FeatureRelationship();
    thisFeatureRelationship.setFeatureByObjectId(matchFeature);
    thisFeatureRelationship.setFeatureBySubjectId(thisFeature);
    thisFeatureRelationship.setCvTerm(getCvTerm(qualifierName));
    if(!rank.equals(""))
      thisFeatureRelationship.setRank(Integer.parseInt(rank));
    
    List featureRelationshipsForSubjectId = new Vector();
    featureRelationshipsForSubjectId.add(thisFeatureRelationship);
    
    FeatureRelationship thisOrthologRelationship = new FeatureRelationship();
    thisOrthologRelationship.setFeatureByObjectId(matchFeature);
    thisOrthologRelationship.setFeatureBySubjectId(orthologFeature);
    thisOrthologRelationship.setCvTerm(getCvTerm(qualifierName));
    featureRelationshipsForSubjectId.add(thisOrthologRelationship);
    matchFeature.setFeatureRelationshipsForSubjectId(featureRelationshipsForSubjectId);
    
    return matchFeature;
  }
  
  
  
  
  /**
   * Return a FeatureLoc for a Feature
   * @param location
   * @param feature
   * @param strand
   * @param rank
   * @return
   */
  private static FeatureLoc getFeatureLoc(
      final String location,
      final org.gmod.schema.sequence.Feature feature,
      final Short strand,
      final int rank)
  {
    final FeatureLoc featureLoc = new FeatureLoc();
    
    if(!location.equals(""))
    {  
      String locs[] = location.split(" ");
      locs = locs[1].split("-");
      int fmin = Integer.parseInt(locs[0]) - 1;
      featureLoc.setFmin(new Integer(fmin));
      int fmax = Integer.parseInt(locs[1]);
      featureLoc.setFmax(new Integer(fmax));
    }
    else
    {
      featureLoc.setFmin(new Integer(-1));
      featureLoc.setFmax(new Integer(-1));
    }
    featureLoc.setRank(rank);
    featureLoc.setStrand(strand);
    featureLoc.setFeatureBySrcFeatureId(feature);
    return featureLoc;
  }

  public static String getString(final StringVector sv, final String name)
  {
    for(int i = 0; i < sv.size(); i++)
    {
      String value = (String) sv.get(i);
      if(value.trim().startsWith(name))
        return value.trim();
    }
    return "";
  }
  
  private static String getStringEndsWith(final StringVector sv, final String name)
  {
    for(int i = 0; i < sv.size(); i++)
    {
      String value = (String) sv.get(i);
      if(value.trim().endsWith(name))
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
      logger4j.debug("USE CvTerm from cache, CvTermId=" + cvTermName
          + "  -> " + cvTerm.getCvTermId() + " " +  
          cvTerm.getCv().getName() + ":" + cvTerm.getName());
    }
    else
    {
      logger4j.warn("CvTerm not found in cache = " + cvTermName);
      cvTerm = new CvTerm();
      cvTerm.setName(cvTermName);
    }
    return cvTerm;
  }
}