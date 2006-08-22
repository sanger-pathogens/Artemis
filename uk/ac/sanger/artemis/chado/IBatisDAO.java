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
 * iBATIS implemetation of the <code>ChadoDAO</code> data
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
  public ChadoFeature getFeatureById(int id) 
                      throws SQLException
  {
    ChadoFeature feature = new ChadoFeature();
    feature.setId(id);
    return getLazyFeature(feature);
  }
  
  public ChadoFeature getFeatureByUniqueName(String uniquename) 
                      throws SQLException
  {
    ChadoFeature feature = new ChadoFeature();
    feature.setUniquename(uniquename);
    return getLazyFeature(feature);
  }
   

  /**
   * This can be used to get individual features or children.
   * If ChadoFeature.featureloc.srcfeature_id is set this is used
   * to return the children of that srcfeature_id.
   * @param feature  the feature to query
   * @return    the <code>List</code> of child <code>ChadoFeature</code> objects
   * @throws SQLException
   */
  public List getFeaturesByLocatedOnFeature(final ChadoFeature feature)
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
    final ChadoSynonym alias = new ChadoSynonym();
    alias.setName(name);
    
    List feature_synonym_list = 
      sqlMap.queryForList("getFeatureSynonymsByName", alias);
    
    ChadoFeature feature = new ChadoFeature();
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
   * @return  the <code>List</code> of <code>ChadoFeature</code>
   * @throws SQLException
   */
  private ChadoFeature getLazyFeature(final ChadoFeature feature)
                             throws SQLException
  { 
    return (ChadoFeature)sqlMap.queryForObject("getLazyFeature", feature);
  }

  /**
   * Given a list of distict cvterm_id/type_id's of feature types
   * that have residues (from getResidueType()) in the given schema 
   * and the schema name return a list of chado features in the schema
   * with residues.
   * @param cvterm_ids list of cvterm_id/type_id's
   * @param schema      schema/organism name or null
   * @return    the <code>List</code> of <code>ChadoFeature</code> objects
   * @throws SQLException
   */
  public List getResidueFeatures(List cvterm_ids, 
                                 final String schema)
                     throws SQLException
  { 
    ChadoFeature feature = new ChadoFeature();
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
   *                    all <code>ChadoFeatureDbxref</code> are returned.
   * @return a <code>List</code> of feature_dbxrefs.
   * @throws SQLException
   */
  public List getFeatureDbxrefByUniquename(final String uniquename)
              throws SQLException
  {
    ChadoFeature feature = new ChadoFeature();
    feature.setUniquename(uniquename);
    
    return sqlMap.queryForList("getFeatureDbxref", feature);  
  }
   
  /**
   * Return a list of ChadoFeatureSynonyms for a uniquename
   * @param uniquename  the unique name for the feature. If set to NULL
   *                    all <code>ChadoFeatureSynonym</code> are returned.
   * @return
   * @throws SQLException
   */
  public List getFeatureSynonymsByUniquename(final String uniquename)
         throws SQLException
  {
    ChadoFeature feature = new ChadoFeature();
    feature.setUniquename(uniquename);
    
    return sqlMap.queryForList("getFeatureSynonymsByUniquename", feature);  
  }
  
  /**
   * Return a synonym of the given name and type if it exists
   * @param name the name to lookup
   * @param type the type of the Synonym
   * @return a Synonym, or null  
   */
  public ChadoSynonym getSynonymByNameAndCvTerm(
      String name, ChadoCvterm type) 
      throws SQLException
  {
    ChadoSynonym synonym = new ChadoSynonym();
    synonym.setName(name);
    synonym.setCvterm(type);

    return (ChadoSynonym)sqlMap.queryForObject("getSynonymByNameAndType", 
           synonym);
  }
  
  /**
   * Return a list of ChadoFeatureSynonyms which link a given Feature
   * and Synonym 
   * @param feature the test Feature
   * @param synonym the test Synonym
   * @return a (possibly empty) List<FeatureSynonym>
   */
  public List getFeatureSynonymsByFeatureAndSynonym(
      ChadoFeature feature, ChadoSynonym synonym)
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
  public ChadoCvterm getCvtermID(String name, String cv_name)
                throws SQLException
  { 
    ChadoCvterm cvterm   = new ChadoCvterm();
    ChadoCv cv = new ChadoCv();
    cv.setName(cv_name);
    cvterm.setCv(cv);
    cvterm.setName(name);

    return (ChadoCvterm)sqlMap.queryForObject("getCvterm", cvterm);
  }

//
// WRITE BACK
//
  /**
   *
   * Update attributes defined by the <code>ChadoTransaction</code>.
   * @param tsn         the <code>ChadoTransaction</code>
   * @return	number of rows changed
   * @throws SQLException
   */
  public int updateAttributes
                    (final ChadoTransaction tsn)
                     throws SQLException 
  { 
    return sqlMap.update("updateAttributes", tsn);
  }

  /**
   * Insert attributes defined by the <code>ChadoTransaction</code>.
   * @param tsn         the <code>ChadoTransaction</code>
   * @throws SQLException
   */
  public void insertAttributes
                    (final ChadoTransaction tsn)
                     throws SQLException
  {
    // get the feature id's
    List feature_ids = sqlMap.queryForList("getFeatureID", tsn);

    for(int i=0; i<feature_ids.size(); i++)
    {
      tsn.setFeature_id( ((ChadoFeature)feature_ids.get(i)).getId() );
      sqlMap.insert("insertAttributes", tsn);
    }
  }

  /**
   * Delete attributes defined by the <code>ChadoTransaction</code>.
   * @param tsn         the <code>ChadoTransaction</code>
   * @throws SQLException
   */
  public void deleteAttributes
                    (final ChadoTransaction tsn)
                     throws SQLException
  {
    // get the feature id's
    List feature_ids = sqlMap.queryForList("getFeatureID", tsn);

    for(int i=0; i<feature_ids.size(); i++)
    {
      tsn.setFeature_id( ((ChadoFeature)feature_ids.get(i)).getId() );
      sqlMap.delete("deleteAttributes", tsn);
    }
  }

  /**
   * Insert a feature into the database defined by the <code>ChadoTransaction</code>.
   * @param tsn                 the <code>ChadoTransaction</code>
   * @parma srcfeature_id       the parent feature identifier
   * @throws SQLException
   */
  public void insertFeature
                    (final ChadoTransaction tsn,
                     final String srcfeature_id)
                     throws SQLException
  {
    // get the organism id from the srcfeature_id 
    ChadoFeature feature = new ChadoFeature();
    ChadoFeatureLoc featureloc = new ChadoFeatureLoc();
    feature.setFeatureloc(featureloc);
    
    featureloc.setSrcfeature_id(Integer.parseInt(srcfeature_id));
    
    Integer organism_id = (Integer)sqlMap.queryForObject("getOrganismID", feature);

    //
    // insert feature into feature table
    ChadoFeature chadoFeature = tsn.getChadoFeature();
    ChadoOrganism organism = new ChadoOrganism();
    organism.setId(organism_id.intValue());
    chadoFeature.setOrganism(organism);  
    sqlMap.insert("insertFeature", chadoFeature);

    //
    // get the current feature_id sequence value
    int feature_id = ((Integer)sqlMap.queryForObject("currval", 
                              "feature_feature_id_seq")).intValue();

    //
    // insert feature location into featureloc
    chadoFeature.getFeatureloc().setSrcfeature_id(Integer.parseInt(srcfeature_id));
    chadoFeature.setId(feature_id);
    sqlMap.insert("insertFeatureLoc", chadoFeature);
    
    // insert feature relationship
    if(tsn.getParents() != null || 
       tsn.getDerives_from() != null)
    {
      List feature_ids = sqlMap.queryForList("getFeatureID", tsn);
      
      for(int i=0; i<feature_ids.size(); i++)
      {
        ChadoFeatureRelationship feature_relationship =
              new ChadoFeatureRelationship();
        feature_relationship.setObject_id( ((ChadoFeature)feature_ids.get(i)).getId() );
        feature_relationship.setSubject_id(feature_id);
        
        ChadoCvterm cvterm = new ChadoCvterm();
        if(tsn.getParents().contains( 
            ((ChadoFeature)feature_ids.get(i)).getUniquename() ))
          cvterm.setName("part_of");
        else
          cvterm.setName("derives_from");
          
        feature_relationship.setCvterm(cvterm);
        chadoFeature.setFeature_relationship(feature_relationship);
        sqlMap.insert("insertFeatureRelationship", chadoFeature);
      }
   
    }
  }


  /**
   * Delete a feature from the database defined by the <code>ChadoTransaction</code>.
   * @param tsn         the <code>ChadoTransaction</code>
   * @return    number of rows deleted
   * @throws SQLException
   */
  public int deleteFeature
                    (final ChadoTransaction tsn)
                     throws SQLException
  {
    ChadoFeature chadoFeature = new ChadoFeature();
    chadoFeature.setUniquename(tsn.getUniqueName());

    return sqlMap.delete("deleteFeature", chadoFeature);
  }

  /**
   * Insert a dbxref for a feature.
   * @param tsn           the <code>ChadoTransaction</code>
   * @return    number of rows changed
   * @throws SQLException
   */
  public int insertFeatureDbxref(final ChadoTransaction tsn)
                     throws SQLException
  {
    ChadoFeatureDbxref dbxref = tsn.getFeatureDbxref();
    
    Integer db_id = (Integer)sqlMap.queryForObject("getDbId", 
                                 dbxref.getDbxref().getDb());
    if(db_id == null)
      throw new SQLException("No database called "+
                             dbxref.getDbxref().getDb().getName()+" found (for "+
                             tsn.getUniqueName()+
                             ") check the spelling!");
    
    dbxref.getDbxref().setDb_id(db_id.intValue());
    
    Integer dbxref_id = (Integer)sqlMap.queryForObject("getDbxrefId", dbxref.getDbxref());
    if(dbxref_id == null)
    {
      // create a new accession entry in dbxref
      sqlMap.insert("insertDbxref", dbxref.getDbxref());
      // now get the new dbxref_id
      dbxref_id = (Integer)sqlMap.queryForObject("getDbxrefId", dbxref.getDbxref());
    }
    
    dbxref.setDbxref_id(dbxref_id.intValue());
    
    //  get the feature id's
    List feature_ids = sqlMap.queryForList("getFeatureID", tsn);
    dbxref.setFeature_id( ((ChadoFeature)feature_ids.get(0)).getId() );
    
    sqlMap.insert("insertFeatureDbxref", dbxref);

    return 1;
  }
  
  /**
   * Delete a dbxref for a feature.
   * @param tsn           the <code>ChadoTransaction</code>
   * @return    number of rows changed
   * @throws SQLException
   */
  public int deleteFeatureDbxref(final ChadoTransaction tsn)
                     throws SQLException
  {
    ChadoFeatureDbxref dbxref = tsn.getFeatureDbxref();
    
    // get the feature id's
    List feature_ids = sqlMap.queryForList("getFeatureID", tsn);
    
    dbxref.setFeature_id( ((ChadoFeature)feature_ids.get(0)).getId() );
    return sqlMap.delete("deleteFeatureDbxref", dbxref);
  }
  
  /**
   * Insert a synonym for a feature.
   * @param tsn           the <code>ChadoTransaction</code>
   * @return    number of rows changed
   * @throws SQLException
   */
  public int insertFeatureAlias(final ChadoTransaction tsn)
                     throws SQLException
  {
    final ChadoFeatureSynonym alias = tsn.getAlias();
    
    ChadoSynonym synonym  = 
      (ChadoSynonym)sqlMap.queryForObject("getSynonymByNameAndType", 
          alias.getSynonym());
    
    if(synonym == null)
    {
      // create a new synonym name     
      sqlMap.insert("insertAlias", alias);
      
      synonym  =
        (ChadoSynonym)sqlMap.queryForObject("getSynonymByNameAndType", 
            alias.getSynonym());
    }
    
    alias.setSynonym_id(synonym.getSynonym_id());
    sqlMap.insert("insertFeatureAlias", alias);
    return 1;
  }
  
  /**
   * Delete a synonym for a feature.
   * @param tsn           the <code>ChadoTransaction</code>
   * @return    number of rows changed
   * @throws SQLException
   */
  public int deleteFeatureAlias(final ChadoTransaction tsn)
                     throws SQLException
  {
    final ChadoFeatureSynonym alias = tsn.getAlias();
     
    List feature_synonym_list = 
      sqlMap.queryForList("getFeatureSynonymsByName", alias.getSynonym());
    
    final ChadoFeatureSynonym synonym = 
      (ChadoFeatureSynonym)feature_synonym_list.get(0); 
    alias.setSynonym_id(synonym.getSynonym().getSynonym_id());
    
    // check this name is not used some where else, 
    // i.e. in more than one row
    if(feature_synonym_list.size() > 1)
      return sqlMap.delete("deleteFeatureAlias", alias);
    else
      return sqlMap.delete("deleteAlias", alias);
  }
  
  /**
   * Update feature_relationship for a feature.
   * @param tsn           the <code>ChadoTransaction</code>
   * @return    number of rows changed
   * @throws SQLException
   */
  public void updateFeatureRelationshipsForSubjectId(
      final ChadoTransaction tsn)
                     throws SQLException
  { 
    sqlMap.update("updateFeatureRelationshipsForSubjectId", tsn);
  }
  
  /**
   * Write the time a feature was last modified
   * @param uniquename  the unique name of the feature
   * @param timestamp   the time stamp to use, 
   *                    if NULL use CURRENT_TIMESTAMP
   * @return  number of rows changed
   * @throws SQLException
   */
  public int writeTimeLastModified
                    (final String uniquename,
                     final Timestamp timestamp)
                     throws SQLException
  {
    ChadoTransaction tsn = new ChadoTransaction(ChadoTransaction.UPDATE,
                                                uniquename, "feature", 
                                                null, null);
    if(timestamp == null)
      tsn.addProperty("timelastmodified", "CURRENT_TIMESTAMP");
    else
      tsn.addProperty("timelastmodified", "'"+ timestamp.toString() + "'");
    
    return updateAttributes(tsn);
  }

  /**
   * Write the time a feature was last accessed
   * @param uniquename  the unique name of the feature
   * @param timestamp   the time stamp to use, 
   *                    if NULL use CURRENT_TIMESTAMP
   * @return  number of rows changed
   * @throws SQLException
   */
  public int writeTimeAccessioned
                    (final String uniquename,
                     final Timestamp timestamp)
                     throws SQLException
  {
    ChadoTransaction tsn = new ChadoTransaction(ChadoTransaction.UPDATE,
                                                uniquename, "feature", 
                                                null, null);
    if(timestamp == null)
      tsn.addProperty("timelastmodified", "CURRENT_TIMESTAMP");
    else
      tsn.addProperty("timelastmodified", "'"+ timestamp.toString() + "'");
    
    return updateAttributes(tsn);
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
    ChadoFeature featNext  = null;

    for(int i = 0; i < feature_size; i++)
    {
      ChadoFeature feat = (ChadoFeature)list.get(i);
      String name  = feat.getUniquename();

      feat.addQualifier(feat.getFeatureprop().getCvterm().getCvtermId(),
                        feat.getFeatureprop());

      if(i < feature_size - 1)
        featNext = (ChadoFeature)list.get(i + 1);

      // merge next line if part of the same feature
      while(featNext != null && featNext.getUniquename().equals(name))
      {
        feat.addQualifier(featNext.getFeatureprop().getCvterm().getCvtermId(),
                          featNext.getFeatureprop());
        i++;

        if(i < feature_size - 1)
          featNext = (ChadoFeature)list.get(i + 1);
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
      ChadoFeatureDbxref dbxref = (ChadoFeatureDbxref)list.get(i);
      Integer feature_id = new Integer(dbxref.getFeature_id());
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

