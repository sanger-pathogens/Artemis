/* IBatisDAO                                                                                                 /* IBatisDAO
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

import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.Hashtable;
import java.util.Vector;
import java.sql.*;

import org.gmod.schema.sequence.Feature;
import org.gmod.schema.sequence.FeatureCvTerm;
import org.gmod.schema.sequence.FeatureCvTermDbXRef;
import org.gmod.schema.sequence.FeatureCvTermProp;
import org.gmod.schema.sequence.FeatureCvTermPub;
import org.gmod.schema.sequence.FeatureDbXRef;
import org.gmod.schema.sequence.FeatureProp;
import org.gmod.schema.sequence.FeaturePub;
import org.gmod.schema.sequence.Synonym;
import org.gmod.schema.sequence.FeatureLoc;
import org.gmod.schema.sequence.FeatureRelationship;
import org.gmod.schema.sequence.FeatureSynonym;
import org.gmod.schema.general.Db;
import org.gmod.schema.general.DbXRef;
import org.gmod.schema.organism.Organism;
import org.gmod.schema.pub.Pub;
import org.gmod.schema.pub.PubDbXRef;
import org.gmod.schema.analysis.AnalysisFeature;
import org.gmod.schema.cv.Cv;
import org.gmod.schema.cv.CvTerm;

import javax.sql.DataSource;
import javax.swing.JPasswordField;

/**
 *
 * iBATIS implemetation of the <code>DAO</code> data
 * access interface.
 *
 */
public class IBatisDAO extends GmodDAO
{
  private static org.apache.log4j.Logger logger4j = 
    org.apache.log4j.Logger.getLogger(IBatisDAO.class);
  private SqlMapClientWrapper sqlMap;
  private boolean featureCvTermRank        = true;
  private boolean featureCvTermRankChecked = false;
  
  /**
   * Define a iBatis data access object. This uses <code>DbSqlConfig</code>
   * to read the configuration in. The system property <quote>chado</quote>
   * can be used to define the database location <i>e.g.</i>
   * -Dchado=host:port/database?user
   */
  public IBatisDAO(final JPasswordField pfield)
  {
    SqlMapClientWrapper sqlMap = new SqlMapClientWrapper(pfield);
    
    this.sqlMap = sqlMap;
    
/*    DbSqlConfig sql_config = new DbSqlConfig();
    sql_config.init(pfield);
    this.sqlMap = sql_config.getSqlMapInstance();*/
  }
  
  public void close() throws SQLException
  {
    sqlMap.close();
  }
  
  public DataSource getDataSource() throws SQLException
  {
    return sqlMap.getSqlMap().getDataSource();
  }

  /**
   * Test to see if this is an old chado database version
   * without the feature_cvterm.rank column.
   * @return true if the feature_cvterm.rank column exists.
   */
  boolean isFeatureCvTermRank()
  {
  	if(featureCvTermRankChecked)
  		return featureCvTermRank;
  	
    String schema = ArtemisUtils.getCurrentSchema();   
    // check column names
    List list = sqlMap.queryForList("getFeatureCvTermColumnsForASchema", schema);   
    featureCvTermRank = false;
    for(int i=0; i<list.size(); i++)
    {
      if( ((String)list.get(i)).equals("rank") )
      {  
      	featureCvTermRank = true;
        break;
      }
    }
    
    featureCvTermRankChecked = true;
    return featureCvTermRank;
  }
  
  //////
  ////// GeneralDaoI
  //////
  //////
  public List getDbs()
  {
    return sqlMap.queryForList("getDbs", null);
  }
  
  //////
  ////// SequenceDaoI
  //////
  //////
  
  public Feature getLazyFeatureNoResiduesById(final Integer featureId)
  {
    return (Feature)sqlMap.queryForObject("getLazyFeatureNoResiduesById",featureId);
  }
  
  public List getFeatureDbXRefsByFeatureId(final List featureIds)
  {
    return sqlMap.queryForList("getFeatureDbXRefsByFeatureId",featureIds);
  }
  
  public List getResiduesByUniqueName(final String uniqueName)
  {
    return sqlMap.queryForList("getResiduesByUniqueName",uniqueName);
  }
  
  public List getResidueFeatures(final Integer organismId)
  {

    CvTerm cvTerm = getCvTermByNameAndCvName("top_level_seq", "genedb_misc");
    if(cvTerm  != null)
      return sqlMap.queryForList("getTopLevelFeatures",organismId);
    else
      return sqlMap.queryForList("getResidueFeatures",organismId);
  }
  
  public List getResidueFeaturesByOrganismCommonName(final String commonName)
  {
    return sqlMap.queryForList("getResidueFeaturesByOrganismCommonName",commonName);
  }
  
  public List getFeatureLocsByListOfIds(final Collection featureIds)
  {
    return sqlMap.queryForList("getFeatureLocsByListOfIds",featureIds);
  }
  
  public List getParentFeaturesByChildFeatureIds(final List featureIds)
  {
    return sqlMap.queryForList("getParentFeaturesByChildFeatureIds",featureIds);
  }
  
  public List getFeatureCvTermsBySrcFeature(Feature srcFeature)
  {
  	if(isFeatureCvTermRank())
    {
      logger4j.debug("USE getFeatureCvTermsBySrcFeature()");
      return sqlMap.queryForList("getFeatureCvTermsBySrcFeature",
          srcFeature);
    }
    else
    {
      logger4j.debug("USE getFeatureCvTermsNoRankBySrcFeature()");
      return sqlMap.queryForList("getFeatureCvTermsNoRankBySrcFeature",
          srcFeature);
    }
  }
  
  public List getFeaturePubsBySrcFeature(final Feature srcFeature)
  {
    return sqlMap.queryForList("getFeaturePubsBySrcFeature",
                               srcFeature);
  }
  
  
  public List getFeaturePubsByFeature(final Feature feature)
  {
    return sqlMap.queryForList("getFeaturePubsByFeature",
                               feature);
  }

  public List getClustersByFeatureIds(final List featureIds)
  {
    return sqlMap.queryForList("getLazyClustersByFeatureIds", featureIds);
  }
  
  /**
   * Returns matches for all features on a given srcfeature
   */
  public List getSimilarityMatches(final Integer srcFeatureId)
  {
    return sqlMap.queryForList("getLazySimilarityMatches", srcFeatureId);
  }
  
  /**
   * Returns matches for a list of feature_id's
   */
  public List getSimilarityMatchesByFeatureIds(final List featureIds)
  {
    return sqlMap.queryForList("getLazySimilarityMatchesByFeatureIds", featureIds);
  }
  
  public List getFeaturePropByFeatureIds(final List featureIds)
  {
    return sqlMap.queryForList("getFeaturePropByFeatureIds", featureIds);
  }
  
  /**
   * Return the feature corresponding to this feature_id 
   * 
   * @param id the systematic id
   * @return the Feature, or null
   */
  public Feature getFeatureById(int id) 
  {
    org.gmod.schema.sequence.Feature feature = 
      new org.gmod.schema.sequence.Feature();
    feature.setFeatureId(id);
    
    if(isFeatureCvTermRank())
      return (Feature)sqlMap.queryForObject("getLazyFeature", feature);
    else
    	return (Feature)sqlMap.queryForObject("getLazyFeatureNoFeatureCvTermRank", feature);
  }
  
  public List getFeaturesByListOfIds(final List featureIds)
  {
    return sqlMap.queryForList("getFeaturesByListOfIds", featureIds);
  }
  
  public List getFeaturesByUniqueName(String uniquename) 
  {
    org.gmod.schema.sequence.Feature feature = 
      new org.gmod.schema.sequence.Feature();
    feature.setUniqueName(uniquename);
    
    if(isFeatureCvTermRank())
      return sqlMap.queryForList("getLazyFeature", feature);
    else
    	return sqlMap.queryForList("getLazyFeatureNoFeatureCvTermRank", feature);
  }
  
  /**
   * Return a <code>Feature</code> based on its uniquename and type_id
   */
  public Feature getFeatureByUniqueName(String uniquename, String featureType) 
  {
    org.gmod.schema.sequence.Feature feature = 
      new org.gmod.schema.sequence.Feature();
    feature.setUniqueName(uniquename);
    
    CvTerm cvTerm = new CvTerm();
    cvTerm.setName(featureType);
    feature.setCvTerm(cvTerm);
    
    if(isFeatureCvTermRank())
      return (Feature)sqlMap.queryForObject("getLazyFeatureExact", feature);
    else
    	return (Feature)sqlMap.queryForObject("getLazyFeatureExactNoFeatureCvTermRank", feature);
  }
   
  
  /**
   * Return a list of features with any current (ie non-obsolete) name or synonym  
   * @param name the lookup name
   * @return a (possibly empty) List<Feature> of children with this current name
   */
  public List getFeaturesByAnyCurrentName(String name) 
  {
    final Synonym alias = new Synonym();
    alias.setName(name);
    
    List feature_synonym_list = 
      sqlMap.queryForList("getFeatureSynonymsByName", alias);
    
    Feature feature = new Feature();
    feature.setUniqueName(name);
    feature.setFeatureSynonyms(feature_synonym_list);

    if(isFeatureCvTermRank())
      return sqlMap.queryForList("getLazyFeature", feature);
    else
    	return sqlMap.queryForList("getLazyFeatureNoFeatureCvTermRank", feature);
  }
  
  /**
   * Return a list of features with this name or synonym (including obsolete names)
   *  
   * @param name the lookup name
   * @return a (possibly empty) List<Feature> of children with this name
   */
  public List getFeaturesByAnyName(String name, String featureType)
  {
    return null;
  }
  
  /**
   * Return a list of features located on a source Feature, within a given range
   *  
   * @param min the minimum (interbase) coordinate
   * @param max the maximum (interbase) coordinate
   * @param strand 
   * @param parent the source feature
   * @param type 
   * @return a List<Feature> which ??? this range
   */
  public List getFeaturesByRange(int min, int max, int strand,
      Feature parent, String type)
  {
    Feature feature = new Feature();
    FeatureLoc featureLoc = new FeatureLoc();
    featureLoc.setFmin(new Integer(min));
    featureLoc.setFmax(new Integer(max));
    feature.setFeatureLoc(featureLoc);
    featureLoc.setFeatureBySrcFeatureId(parent);
    return sqlMap.queryForList("getFeature", feature); 
  }
  
  /**
   * This can be used to get individual features or children.
   * If Feature.featureloc.srcfeature_id is set this is used
   * to return the children of that srcfeature_id.
   * @param feature  the feature to query
   * @return    the <code>List</code> of child <code>Feature</code> objects
   */
  public List getFeaturesByLocatedOnFeature(final Feature feature)
  { 
  	if(isFeatureCvTermRank())
      return sqlMap.queryForList("getFeature", feature);
  	else
  		return sqlMap.queryForList("getFeatureNoFeatureCvTermRank", feature);
  }
  
  /**
   * Return the FeatureCvTerm that links a given Feature and CvTerm, 
   * with a given value of 'not'
   * 
   * @param feature the Feature to test the link for
   * @param cvTerm the CvTerm to test the link for
   * @param not test for the not flag in the FeatureCvTerm 
   * @return the Feature, or null
   */
  public FeatureCvTerm getFeatureCvTermByFeatureAndCvTerm(
          Feature feature,
          CvTerm cvTerm, boolean not)
  {
    return null;
  }
  
  /**
   * Return a list of FeatureCvterm's for a Feature, or a list
   * of all FeatureCvTerm's if Feature is null.
   * @param feature the Feature to retrieve associated FeatureCvTerm's
   * @return the FeatureCvTerm's
   */
  public List getFeatureCvTermsByFeature(Feature feature)
  {  
    if(isFeatureCvTermRank())
      return
        sqlMap.queryForList("getFeatureCvTermsByFeature", feature);
    else
      return
        sqlMap.queryForList("getFeatureCvTermsNoRankByFeature", feature);
  }
  
  /**
   * Return a synonym of the given name and type if it exists
   * @param name the name to lookup
   * @param type the type of the Synonym
   * @return a Synonym, or null  
   */
  public Synonym getSynonymByNameAndCvTerm(
      String name, CvTerm type) 
  {
    Synonym synonym = new Synonym();
    synonym.setName(name);
    synonym.setCvTerm(type);

    return (Synonym)sqlMap.queryForObject("getSynonymByNameAndType", 
           synonym);
  }
  
  /**
   * Return a list of FeatureSynonyms which link a given Feature
   * and Synonym 
   * @param feature the test Feature
   * @param synonym the test Synonym
   * @return a (possibly empty) List<FeatureSynonym>
   */
  public List getFeatureSynonymsByFeatureAndSynonym(
      Feature feature, Synonym synonym)
  {
    return
      sqlMap.queryForList("getFeatureSynonymsByName", synonym);
  }
  
  public List getFeatureSynonymsBySrcFeature(final Feature srcFeature)
  {
    return
      sqlMap.queryForList("getFeatureSynonymsBySrcFeature", srcFeature);
  }
  
  public List getFeatureSynonymsByFeatureIds(final List featuresIds)
  {
    return
      sqlMap.queryForList("getFeatureSynonymsByFeatureIds", featuresIds);
  }
  
  /**
   * Return all the FeatureDbXRefs for a given feature, <b>specified by name</b>, or all if 
   * <code>null</code> is passed
   * 
   * @param uniqueName the uniquename of a Feature, or null for all FeatureDbXRefs
   * @return a (possibly empty) List<FeatureDbXRefI> 
   */
  public List getFeatureDbXRefsByFeatureUniquename(final String uniqueName)
  {
    Feature feature = new Feature();
    feature.setUniqueName(uniqueName);
      
    return sqlMap.queryForList("getFeatureDbXRef", feature);  
  }
  
  public List getFeatureDbXRefsBySrcFeature(final Feature srcFeature)
  {
    return sqlMap.queryForList("getFeatureDbXRefsBySrcFeature", srcFeature);
  }
  
  
  /**
   * Return the list of FeatureSynonyms for a given Feature, <b>specified by name</b>, or all if 
   * <code>null</code> is passed
   * 
   * @param uniqueName the uniquename of a Feature, or null for all
   * @return a (possibly empty) List<FeatureSynonymI> of matching synonyms
   */
  public List getFeatureSynonymsByFeatureUniquename(final String uniqueName)
  {
    Feature feature = new Feature();
    feature.setUniqueName(uniqueName);
    
    return sqlMap.queryForList("getFeatureSynonymsByUniquename", feature); 
  }

  public List getAllFeatureSynonymsAsFeature()
  {
    return sqlMap.queryForList("getAllFeatureSynonymsAsFeature", null);
  }
  
  /**
   * Return the list of Features for a given GO number 
   * 
   * 
   * @param go the GO number
   * @return a (possibly empty) List<Feature> of matching genes
   */
  public List getFeatureByGO(final String go)
  {
    return null;
  }
  
  /**
   * Return a list of features contained in this organisms with this name or synonym (including obsolete names). The 
   * name can contain an SQL wildcard (%) 
   *  
   * @param name the lookup name
   * @param featureType the type of feature to return eg "gene"
   * @param organisms the list of organisms
   * @return a (possibly empty) List<Feature> of children with this name
   */
  public List getFeaturesByAnyNameAndOrganism(String nl,List ids,String featureType)
  {
    return null; 
  }
  
  /**
   * Return a list of features that have this particular cvterm 
   * 
   *  
   * @param cvTermName the CvTerm name
   * @return a (possibly empty) List<Feature> of children
   */
  public List getFeaturesByCvTermName(String cvTermName)
  {
    return null;
  }
  
  /**
   * Return a list of top-level features 
   * @return a (possibly empty) List<Feature> of children
   */
  public List getTopLevelFeatures()
  {
    return null;
  }
  
  
  /**
   * Get a list of all FeatureCvTermDbXRef's for a Feature, or a list
   * of all FeatureCvTermDbXRef's if Feature is null.
   * @param feature the Feature to retrieve associated FeatureCvTermDbXRef's
   * @return the FeatureCvTermDbXRef's
   */
  public List getFeatureCvTermDbXRefByFeature(Feature feature)
  {
    return sqlMap.queryForList("getFeatureCvTermDbXRefByFeature", feature);
  }
  
  /**
   * Get a list of all FeatureCvTermDbXRef's for a Feature, or a list
   * of all FeatureCvTermDbXRef's if Feature is null.
   * @param feature the Feature to retrieve associated FeatureCvTermDbXRef's
   * @return the FeatureCvTermDbXRef's
   */
  public List getFeatureCvTermDbXRefBySrcFeature(final Feature srcFeature)
  {
    return sqlMap.queryForList("getFeatureCvTermDbXRefBySrcFeature",
        srcFeature);
  }
  
  
  
  /**
   * Get a list of all FeatureCvTermPub's for a Feature, or a list
   * of all FeatureCvTermPub's if Feature is null.
   * @param feature the Feature to retrieve associated FeatureCvTermPub's
   * @return the FeatureCvTermPub's
   */
  public List getFeatureCvTermPubByFeature(Feature feature)
  {
    return sqlMap.queryForList("getFeatureCvTermPubByFeature", feature);
  }
  
  
  /**
   * Get a list of all FeatureCvTermPub's for a Feature, or a list
   * of all FeatureCvTermPub's if Feature is null.
   * @param feature the Feature to retrieve associated FeatureCvTermPub's
   * @return the FeatureCvTermPub's
   */
  public List getFeatureCvTermPubBySrcFeature(final Feature srcFeature)
  {
    return sqlMap.queryForList("getFeatureCvTermPubBySrcFeature", 
                               srcFeature);
  }
  
  
  public List getProducts()
  {
    return null;  
  }

  
  //////
  ////// SchemaDaoI
  //////
  //////

  /**
   * Return a list of chado features with residues
   * with residues.
   * @return    the <code>List</code> of <code>Feature</code> objects
   */
  public List getResidueFeatures()
  { 
    CvTerm cvTerm = getCvTermByNameAndCvName("top_level_seq", "genedb_misc");
    if(cvTerm  != null)
      return sqlMap.queryForList("getTopLevelFeatures",null);
    else
      return sqlMap.queryForList("getResidueFeatures",null);
  }

  public List getOrganismsContainingSrcFeatures()
  {
    CvTerm cvTerm = getCvTermByNameAndCvName("populated", "genedb_misc");
    if(cvTerm  != null)
      return sqlMap.queryForList("getTopLevelOrganisms", null);
    else  
      return sqlMap.queryForList("getOrganismsContainingSrcFeatures", null);
  }
  
  /**
   *
   * For a schema return the type_id's with residues.
   * @param schema      schema/organism name or null
   * @return    the <code>List</code> of type_id's as <code>String</code>
   *            objects
   */
  public List getResidueType(final String schema)
  { 
    return sqlMap.queryForList("getResidueType", schema);
  }

  /**
   *
   * Get available schemas (as a <code>List</code> of <code>String</code>       
   * objects).
   * @return    the available schemas
   */
  public List getSchema()
  {
    return sqlMap.queryForList("getSchema", null);
  }

  //////
  ////// CvDaoI
  //////
  //////
  public List getAllCvs()
  {
    return sqlMap.queryForList("getAllCvs", null);
  }
  
  /**
   * Retrieve a named CvTerm from a given Cv
   * 
   * @param cvTermName the name of the cvterm
   * @param cv the controlled vocabulary this cvterm is part of
   * @return a (possibly empty) list of matching cvterms
   */
  public List getCvTermByNameInCv(String cvTermName, Cv cv)
  {
    CvTerm cvTerm = new CvTerm();
    cvTerm.setName(cvTermName);
    cvTerm.setCv(cv);
    
    return sqlMap.queryForList("getCvterm", cvTerm);
  }
  
  /**
   * Get the full list of cvterm_id and name as a <code>List</code> of 
   * <code>CvTerm</code> objects.
   * @return    the full list of cvterm_id and name
   */
  public List getCvTerms()
  {
    return sqlMap.queryForList("getCvterm", null);
  }
  
  /**
   * Retrieve a named CvTerm from a given Cv
   * 
   * @param cvTermName the name of the cvterm
   * @param name the controlled vocabulary name this cvterm could be part of
   * @return a (possibly empty) cvterm
   */
  public CvTerm getCvTermByNameAndCvName(String cvTermName, String name)
  {
    Cv cv = new Cv();
    cv.setName(name);
    CvTerm cvTerm = new CvTerm();
    cvTerm.setName(cvTermName);
    cvTerm.setCv(cv);
    
    return (CvTerm)sqlMap.queryForObject("getCvterm", cvTerm);
  }
  
  
  public CvTerm getCvTermById(final int cvTermId)
  {
    return (CvTerm)sqlMap.queryForObject("getCvtermByCvTermId", new Integer(cvTermId));
  }

  //////
  ////// OrganismDaoI
  //////
  //////
  
  public List getOrganisms()
  {
    return sqlMap.queryForList("getOrganism", null);
  }
  
  public Organism getOrganismByCommonName(String commonName)
  {
    Organism o = new Organism();
    o.setCommonName(commonName);
    return (Organism)sqlMap.queryForObject("getOrganism", o);
  }
  
  //////
  ////// PubDaoI
  //////
  //////
  
  /**
   * Get a list of all PubDbXRef's
   * @return list of PubDbXRef's
   */
  public List getPubDbXRef()
  {
    return sqlMap.queryForList("getPubDbXRef", null);
  }
  
  
  //
  // Graph data
  public Graph getGraph(final Integer graphId)
  {
    //return (Graph) sqlMap.queryForObject("getGraph", graphId);
    return null;
  }
  
  
  public List getGraphs(final Integer featureId)
  {
    return sqlMap.queryForList("getGraphs", featureId);
  }
  
  public List getTableColumns(String tableName)
  {
    return sqlMap.queryForList("getTableColumns", tableName);
  }
  
//
// WRITE BACK
//
  
  /**
   * Merge (update) an already persistent object back to the database (at the end of 
   * the current transaction, or depending upon flush mode). This method is defined in 
   * all the DAOs. It's recommended to call it through an appropriate one eg SequenceDaoI
   *  for FeatureI 
   * 
   * @param o The object to merge
   */
  public void merge(Object o) 
  {
    if(o instanceof FeatureLoc)
      updateFeatureLoc((FeatureLoc)o);
    else if(o instanceof Feature)
    {
      /*Feature f = (Feature)o;
      if(f.getSeqLen() > 0 &&
          f.getCvTerm().getName().equals("region"))
       {
         // insert sequence region
         sqlMap.insert("updateRegionSequence", o);
       } */
      
      if(o instanceof FeatureForUpdatingResidues)
      {
        if(! ((FeatureForUpdatingResidues)o).isResidueUpdate() )
          sqlMap.update("updateFeatureLocByChangingSequence", o);
        sqlMap.update("updateFeatureResidues", o);
      }
      else
        sqlMap.update("updateFeature", o);
    }
    else if(o instanceof FeatureProp)
      sqlMap.update("updateFeatureProp", o);
    else if(o instanceof FeatureRelationship)
      sqlMap.update("updateFeatureRelationshipsForSubjectId", o);
    else if(o instanceof FeatureCvTerm)
      sqlMap.update("updateFeatureCvTerm", o);
    
  }
  
  
  /**
   * Save the object to the database (at the end of the current transaction, 
   * or depending upon flush mode). This method is defined in all the DAOs. 
   * It's recommended to call it through an appropriate one eg SequenceDaoI 
   * for FeatureI 
   * @param o The object to store
   */
  public void persist(Object o)
  {
    if(o instanceof FeatureProp)
      sqlMap.insert("insertFeatureProp", o);
    else if(o instanceof Feature)
      insertFeature((Feature)o);
    else if(o instanceof FeatureDbXRef)
      insertFeatureDbXRef((FeatureDbXRef)o);
    else if(o instanceof FeatureSynonym)
      insertFeatureAlias((FeatureSynonym)o);
    else if(o instanceof FeatureCvTerm)
      insertAllFeatureCvTerm((FeatureCvTerm)o);
    else if(o instanceof FeaturePub)
      insertFeaturePub((FeaturePub)o);
    else if(o instanceof FeatureRelationship)
      insertFeatureRelationship((FeatureRelationship)o);
    else if(o instanceof AnalysisFeature)
      insertAnalysisFeature((AnalysisFeature)o);
    else if(o instanceof CvTerm)
      insertCvTerm((CvTerm)o);
  }
  
  
  /**
   * Remove the object from the database (at the end of the current transaction, 
   * or depending upon flush mode). This method is defined in all the DAOs. 
   * It's recommended to call it through an appropriate one eg SequenceDaoI for 
   * FeatureI 
   * @param o The object to delete
   */
  public void delete(Object o)
  {
    if(o instanceof Feature)
      sqlMap.delete("deleteFeature", o);
    else if(o instanceof FeatureProp)
      sqlMap.delete("deleteFeatureProp", o);
    else if(o instanceof FeatureDbXRef)
      sqlMap.delete("deleteFeatureDbXRef", o);
    else if(o instanceof FeatureSynonym)
      deleteFeatureSynonym((FeatureSynonym)o);
    else if(o instanceof FeatureCvTerm)
      sqlMap.delete("deleteFeatureCvTerm", o);
    else if(o instanceof FeaturePub)
      sqlMap.delete("deleteFeaturePub", o);
    else if(o instanceof FeatureRelationship)
      sqlMap.delete("deleteFeatureRelationship", o);
    else if(o instanceof AnalysisFeature)
      deleteAnalysisFeature((AnalysisFeature)o);
  }

  /**
   * Update a featureLoc. This checks if there are multiple featureLocs
   * associated with this feature and will update all of them.
   * @param featureLoc
   */
  private void updateFeatureLoc(FeatureLoc featureLoc)
  {
    List featureLocs = sqlMap.queryForList("getFeatureLoc", 
        featureLoc.getFeatureByFeatureId().getFeatureId());
    
    if(featureLocs.size() > 1)
    {
      // the feature appears to have multiple featureLocs
      int srcFeatureId = featureLoc.getFeatureBySrcFeatureId().getFeatureId();
      int newFmin = featureLoc.getFmin();
      int newFmax = featureLoc.getFmax();
      int diffFmin = 0;
      int diffFmax = 0;

      for(int i=0; i<featureLocs.size(); i++)
      {
        FeatureLoc thisFeatureLoc = (FeatureLoc) featureLocs.get(i);
        int thisSrcFeatureId = 
          thisFeatureLoc.getFeatureBySrcFeatureId().getFeatureId();
        if(srcFeatureId == thisSrcFeatureId)
        {
          diffFmin = newFmin-thisFeatureLoc.getFmin();
          diffFmax = newFmax-thisFeatureLoc.getFmax();
          break;
        }
      }
      
      for(int i=0; i<featureLocs.size(); i++)
      {
        FeatureLoc thisFeatureLoc = (FeatureLoc) featureLocs.get(i);
        int thisSrcFeatureId = 
          thisFeatureLoc.getFeatureBySrcFeatureId().getFeatureId();
        
        if(srcFeatureId != thisSrcFeatureId)
        {
          // check to see if the featureloc fmin/fmax are marked as partial
          // if they are then ignore as they are likely to be the end of
          // a contig
          boolean changed = false;
          if(!thisFeatureLoc.isFminPartial())
          {
            thisFeatureLoc.setFmin(thisFeatureLoc.getFmin()+diffFmin);
            changed = true;
          }
          if(!thisFeatureLoc.isFmaxPartial())
          {
            thisFeatureLoc.setFmax(thisFeatureLoc.getFmax()+diffFmax);
            changed = true;
          }
          if(changed)
            sqlMap.update("updateFeatureLoc", thisFeatureLoc);
        }
      }
    }
    sqlMap.update("updateFeatureLoc", featureLoc);
  }

  public List getFeatureLocsByFeatureId(int featureId)
  {
    return sqlMap.queryForList("getFeatureLoc",featureId);
  }
  
  public List getFeatureLocsBySrcFeatureId(int srcFeatureId)
  {
    return sqlMap.queryForList("getFeatureLocBySrcFeatureId",
                                srcFeatureId);
  }
  
  /**
   * Insert a feature into the database defined by the <code>Feature</code>.
   * @param feature   the feature to insert
   */
  private void insertFeature
                    (final Feature feature)
  {
    final Organism organism;
    if(feature.getFeatureId() > 0 || feature.getFeatureLoc() != null)
    {
      Integer organism_id = (Integer)sqlMap.queryForObject(
              "getOrganismIdBySrcFeatureIdOrFeatureId", feature);
      organism = new Organism();
      organism.setOrganismId(organism_id.intValue()); 
    }
    else
      organism = (Organism)sqlMap.queryForObject("getOrganism", feature.getOrganism());
    
    feature.setOrganism(organism);
    
    //
    // insert feature into feature table
 
    sqlMap.insert("insertFeature", feature);

    //
    // get the current feature_id sequence value
    int feature_id = ((Integer)sqlMap.queryForObject("currval", 
                              "feature_feature_id_seq")).intValue();

    //
    // insert feature location into featureloc
    if(feature.getFeatureLoc() != null)
    {
      feature.setFeatureId(feature_id);
      FeatureLoc featureLoc = feature.getFeatureLoc();
      featureLoc.setFeatureByFeatureId(feature);
      sqlMap.insert("insertFeatureLoc", featureLoc);
    }
    
    // insert feature relationships
    if(feature.getFeatureRelationshipsForSubjectId() != null)
    {
      List parents = new Vector(
          feature.getFeatureRelationshipsForSubjectId());
      
      for(int i=0; i<parents.size(); i++)
      {
        FeatureRelationship feature_relationship =
               (FeatureRelationship)parents.get(i);
        sqlMap.insert("insertFeatureRelationship", feature_relationship);
      }
    }
    
    if(feature.getFeatureProps() != null)
    {
      List featureProps = new Vector(feature.getFeatureProps());
      
      for(int i=0; i<featureProps.size(); i++)
        sqlMap.insert("insertFeatureProp", (FeatureProp)featureProps.get(i));
    }
  }

  private void insertFeatureRelationship(FeatureRelationship feature_relationship)
  {
    sqlMap.insert("insertFeatureRelationship", feature_relationship);
  }
  
  /**
   * Insert a feature_cvterm and associated feature_cvtermprop's,
   * feature_cvterm_dbxref's and feature_cvterm_pub.
   * @param feature_cvterm
   */
  private void insertFeaturePub(final FeaturePub featurePub)
  {
    // get the pub_id and create a new Pub if necessary
    if(featurePub.getPub() != null)
      featurePub.setPub( loadPub(featurePub.getPub()) );
    
    sqlMap.insert("insertFeaturePub", featurePub);
  }
  
  /**
   * Insert a feature_dbxref for a feature.
   * @param feature_dbxref    the <code>FeatureDbXRef</code>
   */
  private void insertFeatureDbXRef(final FeatureDbXRef feature_dbxref)
  {
    DbXRef dbXRef = loadDbXRef(feature_dbxref.getDbXRef());
    feature_dbxref.setDbXRef(dbXRef);
    
    //  get the feature id's  
    List features = getFeaturesByUniqueName(
        feature_dbxref.getFeature().getUniqueName());
    
    feature_dbxref.getFeature().setFeatureId( 
        ((Feature)features.get(0)).getFeatureId() );

    sqlMap.insert("insertFeatureDbXRef", feature_dbxref);
  }
  
  /**
   * Insert analysis_feature, match feature and associated query and subject features
   * @param analysisFeature
   */
  private void insertAnalysisFeature(AnalysisFeature analysisFeature)
  {
    Feature matchFeature   = analysisFeature.getFeature();
    Feature queryFeature   = null;
    Feature subjectFeature = null;
    FeatureLoc queryLoc    = null;
    FeatureLoc subjectLoc  = null;
    
    Collection featureLocs = matchFeature.getFeatureLocsForFeatureId();
    Iterator it = featureLocs.iterator();
    while(it.hasNext())
    {
      FeatureLoc featureLoc = (FeatureLoc)it.next();

      if(featureLoc.getFeatureBySrcFeatureId().getFeatureId() > 0 ||
         featureLoc.getFeatureBySrcFeatureId().getFeatureId() == -1)
      {
        queryFeature = featureLoc.getFeatureBySrcFeatureId();
        queryLoc = featureLoc;
      }
      else
      {
        subjectFeature = featureLoc.getFeatureBySrcFeatureId();
        subjectLoc = featureLoc;
      }
    }
    
    if(queryFeature.getFeatureId() == -1)
      queryFeature = (Feature)getFeaturesByUniqueName(queryFeature.getUniqueName()).get(0);
    
    Integer organism_id = (Integer)sqlMap.queryForObject(
        "getOrganismIdBySrcFeatureIdOrFeatureId", queryFeature);
    Organism organism = new Organism();
    organism.setOrganismId(organism_id.intValue());
    matchFeature.setOrganism(organism);
    subjectFeature.setOrganism(organism);
    
    // insert match feature
    int value = 1;
    List matches = getFeaturesByUniqueName(matchFeature.getUniqueName()+"_%");
    
    for(int i=0; i<matches.size(); i++)
    {
      String name = ((Feature)matches.get(i)).getUniqueName();
      int index = name.lastIndexOf('_');
      name = name.substring(index+1);
      if( Integer.parseInt(name) >= value )
        value = Integer.parseInt(name)+1;
    }
    
    matchFeature.setUniqueName( matchFeature.getUniqueName()+"_"+value);
    sqlMap.insert("insertFeature", matchFeature);
    int matchFeatureId = ((Integer)sqlMap.queryForObject("currval", 
                               "feature_feature_id_seq")).intValue();
    matchFeature.setFeatureId(matchFeatureId);
    
    // insert primary dbXRef
    DbXRef primaryDbXRef = subjectFeature.getDbXRef();
    primaryDbXRef = loadDbXRef(primaryDbXRef);
    subjectFeature.setDbXRef(primaryDbXRef);
    
    // insert subject feature
    Feature f = getFeatureByUniqueName(subjectFeature.getUniqueName(), 
                                subjectFeature.getCvTerm().getName());
    if(f != null)
      subjectFeature.setFeatureId(f.getFeatureId());
    else
    {
      sqlMap.insert("insertFeature", subjectFeature);
      int subjectFeatureId = ((Integer)sqlMap.queryForObject("currval", 
                                      "feature_feature_id_seq")).intValue();
      subjectFeature.setFeatureId(subjectFeatureId);
    }
    
    // insert match featureLocs
    queryLoc.setFeatureBySrcFeatureId(queryFeature);
    queryLoc.setFeatureByFeatureId(matchFeature);
    sqlMap.insert("insertFeatureLoc", queryLoc);
    
    subjectLoc.setFeatureBySrcFeatureId(subjectFeature);
    subjectLoc.setFeatureByFeatureId(matchFeature);
    sqlMap.insert("insertFeatureLoc", subjectLoc);
    
    
    // insert analysis
    Integer analysisId =  (Integer) sqlMap.insert("insertAnalysis", analysisFeature.getAnalysis());
    analysisFeature.getAnalysis().setAnalysisId(analysisId.intValue());
    
    // insert analysis feature
    matchFeature.setFeatureId(matchFeatureId);
    analysisFeature.setFeature(matchFeature);
    sqlMap.insert("insertAnalysisFeature", analysisFeature);
    
    // if subject feature not inserted add dbXRefs and featureProps
    if(f == null)
    {
      // insert dbXRefs
      Collection featureCvTermDbXRefs = subjectFeature.getFeatureDbXRefs();

      if(featureCvTermDbXRefs != null)
      {
        it = featureCvTermDbXRefs.iterator();

        while(it.hasNext())
          insertFeatureDbXRef((FeatureDbXRef) it.next());
      }

      // insert subject featureprops
      Collection featureProps = subjectFeature.getFeatureProps();
      if(featureProps != null)
      {
        it = featureProps.iterator();

        while(it.hasNext())
        {
          FeatureProp featureProp = (FeatureProp) it.next();
          featureProp.setFeature(subjectFeature);
          sqlMap.insert("insertFeatureProp", featureProp);
        }
      }
    }
    
    // insert match featureprops
    Collection featureProps = matchFeature.getFeatureProps();
    if(featureProps != null)
    {
      it = featureProps.iterator();

      while(it.hasNext())
      {
        FeatureProp featureProp = (FeatureProp) it.next();
        featureProp.setFeature(matchFeature);
        sqlMap.insert("insertFeatureProp", featureProp);
      }
    }
  }
  
  /**
   * Insert a cvterm
   */
  private void insertCvTerm(final CvTerm cvTerm)
  {
    Cv cv = (Cv)sqlMap.queryForObject("getCvByName", cvTerm.getCv().getName());
    cvTerm.setCv(cv);
    
    DbXRef dbXRef = cvTerm.getDbXRef();
    dbXRef = loadDbXRef(dbXRef);
    cvTerm.setDbXRef(dbXRef);
    
    sqlMap.insert("insertCvTerm", cvTerm);
  }
  
  /**
   * Insert a feature_synonym for a feature.
   * @param feature_synonym    the <code>FeatureSynonym</code>
   */
  private void insertFeatureAlias(final FeatureSynonym feature_synonym)
  {
    Object synonym  = 
      sqlMap.queryForObject("getSynonymByNameAndType", 
          feature_synonym.getSynonym());
    
    if(synonym == null)
    {
      logger4j.debug("synonym "+feature_synonym.getSynonym()+" needs inserting");
      // create a new synonym name     
      sqlMap.insert("insertAlias", feature_synonym);
      
      synonym  =
        sqlMap.queryForObject("getSynonymByNameAndType", 
            feature_synonym.getSynonym());
    }
    
    feature_synonym.setSynonym((Synonym)synonym);
    sqlMap.insert("insertFeatureAlias", feature_synonym);
  }
  
  
  protected void insertFeatureCvTerm(final FeatureCvTerm feature_cvterm)
  {
    sqlMap.insert("insertFeatureCvTerm", feature_cvterm);
  }
  
  protected int getCurrval(String seq_id)
  {
    return ((Integer)sqlMap.queryForObject("currval", 
                                           seq_id)).intValue();
  }
  
  protected void insertFeatureCvTermProp(FeatureCvTermProp featureCvTermProp)
  {
    sqlMap.insert("insertFeatureCvTermProp", featureCvTermProp);
  }
  
  protected void insertFeatureCvTermPub(FeatureCvTermPub featureCvTermPub)
  {
    sqlMap.insert("insertFeatureCvTermPub", featureCvTermPub);
  }
  
  protected void insertFeatureCvTermDbXRef(FeatureCvTermDbXRef featureCvTermDbXRef)
  {
    sqlMap.insert("insertFeatureCvTermDbXRef", featureCvTermDbXRef);
  }
  
  private void deleteAnalysisFeature(AnalysisFeature analysisFeature)
  {
    //
    Feature matchFeature = analysisFeature.getFeature();
    List featureLocs = new Vector(matchFeature.getFeatureLocsForFeatureId());
    
    Feature subjectFeature = null;
    Feature queryFeature   = null;
    int nsubject = 0;
    int nquery   = 0;
    
    for(int i=0; i<featureLocs.size(); i++)
    {
      FeatureLoc featureLoc = (FeatureLoc)featureLocs.get(i);

      if(featureLoc.getFeatureBySrcFeatureId().getFeatureId() > 0 ||
         featureLoc.getFeatureBySrcFeatureId().getFeatureId() == -1)
      {
        queryFeature = featureLoc.getFeatureBySrcFeatureId();
        nquery = i;
      }
      else
      {
        subjectFeature = featureLoc.getFeatureBySrcFeatureId();
        nsubject = i;
      }
    }
    
    if(queryFeature.getFeatureId() == -1)
    {
      Feature feat = 
        (Feature)getFeaturesByUniqueName(queryFeature.getUniqueName()).get(0);
      queryFeature.setFeatureId(feat.getFeatureId());
      ((FeatureLoc)featureLocs.get(nquery)).setFeatureBySrcFeatureId(queryFeature);
    }
    
    /*
    if(((FeatureLoc)featureLocs.get(0)).getFeatureBySrcFeatureId().getDbXRef() != null)
    {
      subjectFeature = ((FeatureLoc)featureLocs.get(0)).getFeatureBySrcFeatureId();
      nsubject = 0;
    }
    else
    {
      subjectFeature = ((FeatureLoc)featureLocs.get(1)).getFeatureBySrcFeatureId();
      nsubject = 1;
    }
    */
    
    // look for SWALL: and UNIPROT:
    subjectFeature.setUniqueName( "%"+subjectFeature.getDbXRef().getAccession() );
    List subjectFeatures = sqlMap.queryForList("getLazyFeature", subjectFeature);
    Integer matchFeatureId = null;
    
    for(int i=0; i<subjectFeatures.size(); i++)
    {
      ((FeatureLoc)featureLocs.get(nsubject)).setFeatureBySrcFeatureId(
                                              (Feature)subjectFeatures.get(i));
      matchFeature.setFeatureLocsForFeatureId(featureLocs);
      List result = sqlMap.queryForList("getFeatureIdBySrcFeatureId", matchFeature);
      if(result != null)
      {
        matchFeatureId = (Integer)result.get(0);
        break;
      }
    }
    
    matchFeature.setFeatureId(matchFeatureId.intValue());
    sqlMap.delete("deleteFeatureById", matchFeature);
  }
  
  /**
   * Delete a feature_synonym for a feature.
   * @param feature_synonym     the <code>FeatureSynonym</code>
   */
  private int deleteFeatureSynonym(FeatureSynonym feature_synonym)
  {
    List feature_synonym_list = 
      sqlMap.queryForList("getFeatureSynonymsByName", feature_synonym.getSynonym());
    
    feature_synonym.setSynonym( 
        ((FeatureSynonym)feature_synonym_list.get(0)).getSynonym() );
     
    // check this name is not used some where else, 
    // i.e. in more than one row
    if(feature_synonym_list.size() > 1)
      return sqlMap.delete("deleteFeatureAlias", feature_synonym);
    else
      return sqlMap.delete("deleteAlias", feature_synonym);
  }
  
  protected Integer getDbId(Db db)
  {
    Integer dbId = (Integer)sqlMap.queryForObject("getDbId", db);
    
    if(dbId == null)
    {
      List dbIds = sqlMap.queryForList("getDbIdIgnoreCase", db);
      if(dbIds.size() > 0)
        dbId = (Integer)dbIds.get(0);
    }
    
    return dbId;
  }
  
  protected Integer getDbXRefId(DbXRef dbXRef)
  {
    return (Integer)sqlMap.queryForObject("getDbXRefId", dbXRef);
  }
  
  protected void insertDbXRef(DbXRef dbXRef)
  {
    sqlMap.insert("insertDbXRef", dbXRef);
  }
  
  protected Pub getPubByUniqueName(Pub pub)
  {
    return (Pub)sqlMap.queryForObject("getPubByUniqueName", pub);
  }
  
  protected void insertPub(Pub pub)
  {
    sqlMap.insert("insertPub", pub);
  }
  
  protected void insertPubDbXRef(PubDbXRef pubDbXRef)
  {
    sqlMap.insert("insertPubDbXRef", pubDbXRef);
  }

  public void startTransaction() throws SQLException
  { 
    sqlMap.startTransaction();
  }
  
  public void endTransaction() throws SQLException
  { 
    sqlMap.endTransaction();
  }
  
  public void commitTransaction() throws SQLException
  {
    sqlMap.commitTransaction();
  }


  /**
   * Takes a list and creates a <code>Hashtable</code> with the keys
   * being the feature_id and the value a <code>Vector</code> of the dbxrefs.
   * @param list  a <code>List</code> of <code>DbXRef</code> objects.
   * @return a <code>Hashtable</code> of dbxrefs.
   */
  public static Hashtable<Integer, List<String>> mergeDbXRef(final List<FeatureDbXRef> list)
  {
    Hashtable<Integer, List<String>> dbxrefHash = new Hashtable<Integer, List<String>>();
    for(FeatureDbXRef dbxref : list)
    {
      Integer feature_id = new Integer(dbxref.getFeature().getFeatureId());
      String value = dbxref.getDbXRef().getDb().getName() + ":" + 
                     dbxref.getDbXRef().getAccession();
      if(dbxrefHash.containsKey(feature_id))
      {
        List<String> v = dbxrefHash.get(feature_id);
        v.add(value);
        dbxrefHash.put(feature_id, v);
      }  
      else
      {
        List<String> v = new Vector<String>();
        v.add(value);
        dbxrefHash.put(feature_id, v);
      }
    }
    return dbxrefHash;
  }
}
