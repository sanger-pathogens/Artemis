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

import com.ibatis.sqlmap.client.SqlMapClient;

import java.util.List;
import java.util.Hashtable;
import java.util.Vector;
import java.sql.*;

import javax.swing.JPasswordField;

/**
 *
 * iBATIS implemetation of the <code>DAO</code> data
 * access interface.
 *
 */
public class IBatisDAO implements ChadoDAO
{
  private SqlMapClient sqlMap;
  
  /**
   * Define a iBatis data access object. This uses <code>DbSqlConfig</code>
   * to read the configuration in. The system property <quote>chado</quote>
   * can be used to define the database location <i>e.g.</i>
   * -Dchado=host:port/database?user
   */
  public IBatisDAO(final JPasswordField pfield)
  {
    DbSqlConfig sql_config = new DbSqlConfig();
    sql_config.init(pfield);
    this.sqlMap = sql_config.getSqlMapInstance();
  }

  
  /**
   * Return the feature corresponding to this feature_id 
   * 
   * @param id the systematic id
   * @return the Feature, or null
   * @throws SQLException 
   */
  public Feature getFeatureById(int id) 
                      throws SQLException
  {
    Feature feature = new Feature();
    feature.setId(id);
    return getLazyFeature(feature);
  }
  
  public Feature getFeatureByUniqueName(String uniquename) 
                      throws SQLException
  {
    Feature feature = new Feature();
    feature.setUniquename(uniquename);
    return getLazyFeature(feature);
  }
   

  /**
   * This can be used to get individual features or children.
   * If Feature.featureloc.srcfeature_id is set this is used
   * to return the children of that srcfeature_id.
   * @param feature  the feature to query
   * @return    the <code>List</code> of child <code>Feature</code> objects
   * @throws SQLException
   */
  public List getFeaturesByLocatedOnFeature(final Feature feature)
               throws SQLException
  { 
    List feature_list = sqlMap.queryForList("getFeature", feature);

    // merge same features in the list
    return mergeList(feature_list);
  }

  /**
   * Return a list of features with any current (ie non-obsolete) name or synonym  
   * @param name the lookup name
   * @return a (possibly empty) List<Feature> of children with this current name
   * @throws SQLException 
   */
  public List getFeaturesByAnyCurrentName(String name) 
              throws SQLException
  {
    final Synonym alias = new Synonym();
    alias.setName(name);
    
    List feature_synonym_list = 
      sqlMap.queryForList("getFeatureSynonymsByName", alias);
    
    Feature feature = new Feature();
    feature.setUniquename(name);
    feature.setFeatureSynonymsForFeatureId(feature_synonym_list);

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
   * Get the properties of a feature.
   * @param uniquename  the unique name of the feature
   * @return  the <code>List</code> of <code>Feature</code>
   * @throws SQLException
   */
  private Feature getLazyFeature(final Feature feature)
                             throws SQLException
  { 
    return (Feature)sqlMap.queryForObject("getLazyFeature", feature);
  }

  /**
   * Given a list of distict cvterm_id/type_id's of feature types
   * that have residues (from getResidueType()) in the given schema 
   * and the schema name return a list of chado features in the schema
   * with residues.
   * @param cvterm_ids list of cvterm_id/type_id's
   * @param schema      schema/organism name or null
   * @return    the <code>List</code> of <code>Feature</code> objects
   * @throws SQLException
   */
  public List getResidueFeatures(List cvterm_ids, 
                                 final String schema)
                     throws SQLException
  { 
    Feature feature = new Feature();
    feature.setSchema(schema);
    feature.setFeatureCvterms(cvterm_ids);

    return sqlMap.queryForList("getResidueFeatures",
                                feature);
  }

  /**
   *
   * For a schema return the type_id's with residues.
   * @param schema      schema/organism name or null
   * @return    the <code>List</code> of type_id's as <code>String</code>
   *            objects
   * @throws SQLException
   */
  public List getResidueType(final String schema)
                     throws SQLException
  { 
    return sqlMap.queryForList("getResidueType", schema);
  }

  /**
   *
   * Get available schemas (as a <code>List</code> of <code>String</code>       
   * objects).
   * @return    the available schemas
   * @throws SQLException
   */
  public List getSchema()
                throws SQLException
  {
    return sqlMap.queryForList("getSchema", null);
  }

  /**
   * Get the full list of cvterm_id and name as a <code>List</code> of 
   * <code>Cvterm</code> objects.
   * @return    the full list of cvterm_id and name
   * @throws SQLException
   */
  public List getCvterm()
              throws SQLException
  {
    return sqlMap.queryForList("getCvterm", null);
  }

  
  /**
   * Get dbxref for a feature.
   * @param uniquename  the unique name for the feature. If set to NULL
   *                    all <code>FeatureDbxref</code> are returned.
   * @return a <code>List</code> of feature_dbxrefs.
   * @throws SQLException
   */
  public List getFeatureDbxrefByUniquename(final String uniquename)
              throws SQLException
  {
    Feature feature = new Feature();
    feature.setUniquename(uniquename);
    
    return sqlMap.queryForList("getFeatureDbxref", feature);  
  }
   
  /**
   * Return a list of FeatureSynonyms for a uniquename
   * @param uniquename  the unique name for the feature. If set to NULL
   *                    all <code>FeatureSynonym</code> are returned.
   * @return
   * @throws SQLException
   */
  public List getFeatureSynonymsByUniquename(final String uniquename)
         throws SQLException
  {
    Feature feature = new Feature();
    feature.setUniquename(uniquename);
    
    return sqlMap.queryForList("getFeatureSynonymsByUniquename", feature);  
  }
  
  /**
   * Return a synonym of the given name and type if it exists
   * @param name the name to lookup
   * @param type the type of the Synonym
   * @return a Synonym, or null  
   */
  public Synonym getSynonymByNameAndCvTerm(
      String name, Cvterm type) 
      throws SQLException
  {
    Synonym synonym = new Synonym();
    synonym.setName(name);
    synonym.setCvterm(type);

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
      throws SQLException
  {
    return
      sqlMap.queryForList("getFeatureSynonymsByName", synonym);
  }


  
  /**
   *
   * @param name cvterm name
   * @param cv_name ontology name (e.g. gene, sequence)
   * @throws SQLException
   */
  public Cvterm getCvtermID(String name, String cv_name)
                throws SQLException
  { 
    Cvterm cvterm   = new Cvterm();
    Cv cv = new Cv();
    cv.setName(cv_name);
    cvterm.setCv(cv);
    cvterm.setName(name);

    return (Cvterm)sqlMap.queryForObject("getCvterm", cvterm);
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
   * @throws SQLException 
   */
  public void merge(Object o) throws SQLException
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
   * @throws SQLException 
   */
  public void persist(Object o) throws SQLException
  {
    if(o instanceof FeatureProp)
      sqlMap.insert("insertFeatureProp", o);
    else if(o instanceof Feature)
      insertFeature((Feature)o);
    else if(o instanceof FeatureDbxref)
      insertFeatureDbxref((FeatureDbxref)o);
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
  public void delete(Object o) throws SQLException
  {
    if(o instanceof Feature)
      sqlMap.delete("deleteFeature", o);
    else if(o instanceof FeatureProp)
      sqlMap.delete("deleteFeatureProp", o);
    else if(o instanceof FeatureDbxref)
      sqlMap.delete("deleteFeatureDbxref", o);
    else if(o instanceof FeatureSynonym)
      deleteFeatureSynonym((FeatureSynonym)o);
  }


  /**
   * Insert a feature into the database defined by the <code>Feature</code>.
   * @param feature   the feature to insert
   * @throws SQLException
   */
  private void insertFeature
                    (final Feature feature)
                     throws SQLException
  {
    Integer organism_id = (Integer)sqlMap.queryForObject("getOrganismID", feature);

    //
    // insert feature into feature table
    Organism organism = new Organism();
    organism.setId(organism_id.intValue());
    feature.setOrganism(organism);  
    sqlMap.insert("insertFeature", feature);

    //
    // get the current feature_id sequence value
    int feature_id = ((Integer)sqlMap.queryForObject("currval", 
                              "feature_feature_id_seq")).intValue();

    //
    // insert feature location into featureloc
    feature.setId(feature_id);
    sqlMap.insert("insertFeatureLoc", feature);
    
    // insert feature relationships
    if(feature.getFeatureRelationshipsForSubjectId() != null)
    {
      List parents = feature.getFeatureRelationshipsForSubjectId();
      
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
   * @param feature_dbxref    the <code>FeatureDbxref</code>
   * @throws SQLException
   */
  private void insertFeatureDbxref(final FeatureDbxref feature_dbxref)
                     throws SQLException
  {
    Integer db_id = (Integer)sqlMap.queryForObject("getDbId", 
                         feature_dbxref.getDbxref().getDb());
    
    System.out.print(db_id.intValue());

    
    if(db_id == null)
      throw new SQLException("No database called "+
          feature_dbxref.getDbxref().getDb().getName()+" found (for "+
          feature_dbxref.getFeature().getUniquename()+
          ") check the spelling!");
    
    feature_dbxref.getDbxref().setDb_id(db_id.intValue());
    
    Integer dbxref_id = 
      (Integer)sqlMap.queryForObject("getDbxrefId", feature_dbxref.getDbxref());
    if(dbxref_id == null)
    {
      // create a new accession entry in dbxref
      sqlMap.insert("insertDbxref", feature_dbxref.getDbxref());
      // now get the new dbxref_id
      dbxref_id = (Integer)sqlMap.queryForObject("getDbxrefId", 
          feature_dbxref.getDbxref());
    }
    
    feature_dbxref.setDbxref_id(dbxref_id.intValue());
    
    //  get the feature id's  
    Feature feature = getFeatureByUniqueName(
        feature_dbxref.getFeature().getUniquename());
    feature_dbxref.getFeature().setId( feature.getId() );

    sqlMap.insert("insertFeatureDbxref", feature_dbxref);
  }
  
  
  /**
   * Insert a feature_synonym for a feature.
   * @param feature_synonym    the <code>FeatureSynonym</code>
   * @throws SQLException
   */
  private void insertFeatureAlias(final FeatureSynonym feature_synonym)
                     throws SQLException
  {
    Synonym synonym  = 
      (Synonym)sqlMap.queryForObject("getSynonymByNameAndType", 
          feature_synonym.getSynonym());
    
    if(synonym == null)
    {
      // create a new synonym name     
      sqlMap.insert("insertAlias", feature_synonym);
      
      synonym  =
        (Synonym)sqlMap.queryForObject("getSynonymByNameAndType", 
            feature_synonym.getSynonym());
    }
    
    feature_synonym.setSynonym_id(synonym.getSynonym_id());
    sqlMap.insert("insertFeatureAlias", feature_synonym);
  }
  
  /**
   * Delete a feature_synonym for a feature.
   * @param feature_synonym     the <code>FeatureSynonym</code>
   * @throws SQLException
   */
  private int deleteFeatureSynonym(final FeatureSynonym feature_synonym)
                     throws SQLException
  {
    List feature_synonym_list = 
      sqlMap.queryForList("getFeatureSynonymsByName", feature_synonym.getSynonym());
    
    final FeatureSynonym synonym = 
      (FeatureSynonym)feature_synonym_list.get(0); 
    feature_synonym.setSynonym_id(synonym.getSynonym().getSynonym_id());
    
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
   * Takes a list and creates a new one merging all feature objects
   * within it with the same feature and stores the qualifiers/attributes
   *  as a hash
   * @param list of feature objects
   * @return list of flattened/merged feature objects
   */
  protected static List mergeList(final List list)
  {
    // merge same features in the list
    int feature_size  = list.size();
    final List flatten_list = new Vector();
    Feature featNext  = null;

    for(int i = 0; i < feature_size; i++)
    {
      Feature feat = (Feature)list.get(i);
      String name  = feat.getUniquename();

      feat.addQualifier(feat.getFeatureprop().getCvterm().getCvtermId(),
                        feat.getFeatureprop());

      if(i < feature_size - 1)
        featNext = (Feature)list.get(i + 1);
      else
        featNext = null;
      
      // merge next line if part of the same feature
      while(featNext != null && featNext.getUniquename().equals(name))
      {
        feat.addQualifier(featNext.getFeatureprop().getCvterm().getCvtermId(),
                          featNext.getFeatureprop());
        i++;

        if(i < feature_size - 1)
          featNext = (Feature)list.get(i + 1);
        else
          break;
      }

      flatten_list.add(feat);
    }

    return flatten_list;
  }

  /**
   * Takes a list and creates a <code>Hashtable</code> with the keys
   * being the feature_id and the value a <code>Vector</code> of the dbxrefs.
   * @param list  a <code>List</code> of <code>Dbxref</code> objects.
   * @return a <code>Hashtable</code> of dbxrefs.
   */
  public static Hashtable mergeDbxref(final List list)
  {
    Hashtable dbxrefHash = new Hashtable();
    for(int i = 0; i < list.size(); i++)
    {
      FeatureDbxref dbxref = (FeatureDbxref)list.get(i);
      Integer feature_id = new Integer(dbxref.getFeature().getId());
      String value = dbxref.getDbxref().getDb().getName() + ":" + 
                     dbxref.getDbxref().getAccession();
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
