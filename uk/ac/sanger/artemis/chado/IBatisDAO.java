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

import java.util.List;
import java.util.Hashtable;
import java.util.Vector;
import java.sql.*;

import org.gmod.schema.sequence.Feature;
import org.gmod.schema.sequence.FeatureCvTerm;
import org.gmod.schema.sequence.FeatureDbXRef;
import org.gmod.schema.sequence.FeatureProp;
import org.gmod.schema.sequence.Synonym;
import org.gmod.schema.sequence.FeatureLoc;
import org.gmod.schema.sequence.FeatureRelationship;
import org.gmod.schema.sequence.FeatureSynonym;
import org.gmod.schema.general.DbXRef;
import org.gmod.schema.organism.Organism;
import org.gmod.schema.cv.CvTerm;

import javax.swing.JPasswordField;

/**
 *
 * iBATIS implemetation of the <code>DAO</code> data
 * access interface.
 *
 */
public class IBatisDAO extends GmodDAO
{
  private SqlMapClientWrapper sqlMap;
  
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


  //////
  ////// SequenceDaoI
  //////
  //////

  
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
    return getLazyFeature(feature);
  }
  
  public Feature getFeatureByUniqueName(String uniquename) 
  {
    org.gmod.schema.sequence.Feature feature = 
      new org.gmod.schema.sequence.Feature();
    feature.setUniqueName(uniquename);
    return getLazyFeature(feature);
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

    List features = sqlMap.queryForList("getLazyFeature", feature);
    
    return features;
    
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
    return null; 
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
    List feature_list = sqlMap.queryForList("getFeature", feature);

    // merge same features in the list
    //return mergeList(feature_list);
    return feature_list;
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
  
  public List getFeatureCvTermsByFeature(Feature feature)
  {
    return
      sqlMap.queryForList("getFeatureCvTermsByFeature", feature);
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
   * 
   *  
   * @return a (possibly empty) List<Feature> of children
   */
  public List getTopLevelFeatures()
  {
    return null;
  }
  
  
  public List getFeatureCvTermDbXRef(Feature feature)
  {
    return sqlMap.queryForList("getFeatureCvTermDbXRef", feature);
  }
  
  public List getProducts()
  {
    return null;  
  }
  
  //////
  //////
  
  
  /**
   * Get the properties of a feature.
   * @param uniquename  the unique name of the feature
   * @return  the <code>List</code> of <code>Feature</code>
   */
  private Feature getLazyFeature(
      final org.gmod.schema.sequence.Feature feature)
  { 
    return (Feature)sqlMap.queryForObject("getLazyFeature", feature);
  }
  
  //////
  ////// SchemaDaoI
  //////
  //////

  /**
   * Given a list of distict cvterm_id/type_id's of feature types
   * that have residues (from getResidueType()) in the given schema 
   * and the schema name return a list of chado features in the schema
   * with residues.
   * @param cvTermIds   list of cvterm_id/type_id's
   * @param schema      schema/organism name or null
   * @return    the <code>List</code> of <code>Feature</code> objects
   */
  public List getResidueFeatures(List cvTermIds, 
                                 final String schema)
  { 
    java.util.Map map = new java.util.HashMap();
    map.put("cvTermIds", cvTermIds);
    return sqlMap.queryForList("getResidueFeatures",
                               map);
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
    return null;
  }

  //////
  ////// OrganismDaoI
  //////
  //////
  
  public List getOrganisms()
  {
    return sqlMap.queryForList("getOrganism", null);
  }
  
  //////
  ////// PubDaoI
  //////
  //////
  
  public List getPubDbXRef()
  {
    return sqlMap.queryForList("getPubDbXRef", null);
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
      sqlMap.update("updateFeatureLoc", o);
    else if(o instanceof Feature)
      sqlMap.update("updateFeature", o);
    else if(o instanceof FeatureProp)
      sqlMap.update("updateFeatureProp", o);
    else if(o instanceof FeatureRelationship)
      sqlMap.update("updateFeatureRelationshipsForSubjectId", o);
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
  }


  /**
   * Insert a feature into the database defined by the <code>Feature</code>.
   * @param feature   the feature to insert
   */
  private void insertFeature
                    (final Feature feature)
  {
    Integer organism_id = (Integer)sqlMap.queryForObject("getOrganismID", feature);

    //
    // insert feature into feature table
    Organism organism = new Organism();
    organism.setOrganismId(organism_id.intValue());
    feature.setOrganism(organism);  
    sqlMap.insert("insertFeature", feature);

    //
    // get the current feature_id sequence value
    int feature_id = ((Integer)sqlMap.queryForObject("currval", 
                              "feature_feature_id_seq")).intValue();

    //
    // insert feature location into featureloc
    feature.setFeatureId(feature_id);
    sqlMap.insert("insertFeatureLoc", feature);
    
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
  }

  /**
   * Insert a feature_dbxref for a feature.
   * @param feature_dbxref    the <code>FeatureDbXRef</code>
   */
  private void insertFeatureDbXRef(final FeatureDbXRef feature_dbxref)
  {
    Integer db_id = (Integer)sqlMap.queryForObject("getDbId", 
                         feature_dbxref.getDbXRef().getDb());
    
    System.out.print(db_id.intValue());

    
    if(db_id == null)
      throw new RuntimeException("No database called "+
          feature_dbxref.getDbXRef().getDb().getName()+" found (for "+
          feature_dbxref.getFeature().getUniqueName()+
          ") check the spelling!");
    
    feature_dbxref.getDbXRef().setDbXRefId(db_id.intValue());
    
    Integer dbxref_id = 
      (Integer)sqlMap.queryForObject("getDbXRefId", feature_dbxref.getDbXRef());
    if(dbxref_id == null)
    {
      // create a new accession entry in dbxref
      sqlMap.insert("insertDbXRef", feature_dbxref.getDbXRef());
      // now get the new dbxref_id
      dbxref_id = (Integer)sqlMap.queryForObject("getDbXRefId", 
          feature_dbxref.getDbXRef());
    }
    
    DbXRef dbXRef = new DbXRef();
    dbXRef.setDbXRefId(dbxref_id.intValue());
    feature_dbxref.setDbXRef(dbXRef);
    
    //  get the feature id's  
    Feature feature = (Feature)getFeatureByUniqueName(
        feature_dbxref.getFeature().getUniqueName());
    feature_dbxref.getFeature().setFeatureId( feature.getFeatureId() );

    sqlMap.insert("insertFeatureDbXRef", feature_dbxref);
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
      System.out.println("HERE2");
      // create a new synonym name     
      sqlMap.insert("insertAlias", feature_synonym);
      
      synonym  =
        sqlMap.queryForObject("getSynonymByNameAndType", 
            feature_synonym.getSynonym());
    }
    
    System.out.println("HERE "+((Synonym)synonym).getSynonymId()+" "+((Synonym)synonym).getName());
    feature_synonym.setSynonym((Synonym)synonym);
    System.out.println("HERE "+((Synonym)synonym).getSynonymId());
    System.out.println("\n----------------------------------> "+feature_synonym.getSynonym().getSynonymId()+"\n");
    //feature_synonym.setFeatureSynonymId(synonym.getSynonymId());
    sqlMap.insert("insertFeatureAlias", feature_synonym);
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
  public static Hashtable mergeDbXRef(final List list)
  {
    Hashtable dbxrefHash = new Hashtable();
    for(int i = 0; i < list.size(); i++)
    {
      FeatureDbXRef dbxref = (FeatureDbXRef)list.get(i);
      Integer feature_id = new Integer(dbxref.getFeature().getFeatureId());
      String value = dbxref.getDbXRef().getDb().getName() + ":" + 
                     dbxref.getDbXRef().getAccession();
      if(dbxrefHash.containsKey(feature_id))
      {
        Vector v = (Vector)dbxrefHash.get(feature_id);
        v.add(value);
        dbxrefHash.put(feature_id, v);
      }  
      else
      {
        Vector v = new Vector();
        v.add(value);
        dbxrefHash.put(feature_id, v);
      }
    }
    return dbxrefHash;
  }

}
