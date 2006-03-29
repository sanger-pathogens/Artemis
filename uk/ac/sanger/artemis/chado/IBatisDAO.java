/* IBatisDAO
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
  /**
   *
   * Define a iBatis data access object. This uses <code>DbSqlConfig</code>
   * to read the configuration in. The system property <quote>chado</quote>
   * can be used to define the database location <i>e.g.</i>
   * -Dchado=host:port/database?user
   *
   */
  public IBatisDAO(final JPasswordField pfield)
  {
    DbSqlConfig.init(pfield);
  }

  /**
   *
   * Get the feature name given a feature_id and schema.
   * @param feature_id  id of feature to query
   * @param schema      schema/organism name or null
   * @return    the feature name
   */
  public String getFeatureName(final ChadoFeature feature)
                throws SQLException
  {
    SqlMapClient sqlMap = DbSqlConfig.getSqlMapInstance();
    return (String)sqlMap.queryForObject("getFeatureName", feature);
  }

  /**
   *
   * Get the residues of a feature.
   * @param feature_id  id of feature to query
   * @param schema      schema/organism name or null
   * @return    the <code>ChadoFeature</code> with the residues
   * @throws SQLException
   */
  public String getFeatureName(final int feature_id,
                               final String schema)
                       throws SQLException
  {
    ChadoFeature feature = new ChadoFeature();
    feature.setId(feature_id);
    if(schema != null)
      feature.setSchema(schema);
    return getFeatureName(feature);
  }

  /**
   *
   * Get child feature properties for a given parent
   * feature to be able to construct a GFF like feature.
   * @param parentFeatureID  the id of parent feature to query
   * @param schema           the schema/organism name or null
   * @return    the <code>List</code> of child <code>ChadoFeature</code> objects
   * @throws SQLException
   */
  public List getGff(final int feature_id,
                     final String schema)
                     throws SQLException
  {
    SqlMapClient sqlMap = DbSqlConfig.getSqlMapInstance();
    ChadoFeature feature = new ChadoFeature();
    feature.setId(feature_id);
    if(schema != null)
      feature.setSchema(schema);

    List feature_list = sqlMap.queryForList("getGffLine", feature);

    // merge same features in the list
    return JdbcDAO.mergeList(feature_list);
  }

  /**
   *
   * Get the residues of a feature.
   * @param feature_id  id of feature to query
   * @param schema      schema/organism name or null
   * @throws SQLException
   */
  public ChadoFeature getSequence(final int feature_id,
                             final String schema)
                        throws SQLException
  {
    SqlMapClient sqlMap = DbSqlConfig.getSqlMapInstance();
    ChadoFeature feature = new ChadoFeature();
    feature.setId(feature_id);
    if(schema != null)
      feature.setSchema(schema);
    return (ChadoFeature)sqlMap.queryForObject("getSequence",
                                           feature);
  }

  /**
   *
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
    SqlMapClient sqlMap = DbSqlConfig.getSqlMapInstance();
    SchemaCVList schema_CVlist = new SchemaCVList();
    schema_CVlist.setSchema(schema);
    schema_CVlist.setCvlist(cvterm_ids);

    return sqlMap.queryForList("getSchemaResidueFeatures",
                                schema_CVlist);
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
    SqlMapClient sqlMap = DbSqlConfig.getSqlMapInstance();
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
    SqlMapClient sqlMap = DbSqlConfig.getSqlMapInstance();
    return sqlMap.queryForList("getSchema", null);
  }

  /**
   *
   * Get the full list of cvterm_id and name as a <code>List</code> of 
   * <code>Cvterm</code> objects.
   * @return    the full list of cvterm_id and name
   * @throws SQLException
   */
  public List getCvterm()
              throws SQLException
  {
    SqlMapClient sqlMap = DbSqlConfig.getSqlMapInstance();
    return sqlMap.queryForList("getCvterm", null);
  }

  /**
   * 
   * Get dbxref for a feature.
   * @param schema      the postgres schema name
   * @param uniquename  the unique name for the feature. If set to NULL
   *                    all <code>Dbxref</code> are returned.
   * @return a <code>List</code> of <code>Dbxref</code>
   * @throws SQLException
   */
  public List getDbxref(final String schema, final String uniquename)
              throws SQLException
  {
    ChadoFeature feature = new ChadoFeature();
    feature.setSchema(schema);
    
    SqlMapClient sqlMap = DbSqlConfig.getSqlMapInstance();
    return sqlMap.queryForList("getDbxref", feature);  
  }
  
  /**
   *
   * @param name cvterm name
   * @param cv_name ontology name (e.g. gene, sequence)
   * @throws SQLException
   */
  public static Cvterm getCvtermID(String name, String cv_name)
                throws SQLException
  {
    SqlMapClient sqlMap = DbSqlConfig.getSqlMapInstance();
    Cvterm cvterm   = new Cvterm();
    cvterm.setName(name);
    cvterm.setCv_name(cv_name);
    return (Cvterm)sqlMap.queryForObject("getCvterm", cvterm);
  }

//
// WRITE BACK
//
  /**
   *
   * Update attributes defined by the <code>ChadoTransaction</code>.
   * @param schema      schema/organism name or null
   * @param tsn         the <code>ChadoTransaction</code>
   * @return	number of rows changed
   * @throws SQLException
   */
  public int updateAttributes
                    (final String schema, final ChadoTransaction tsn)
                     throws SQLException 
  {
    SqlMapClient sqlMap = DbSqlConfig.getSqlMapInstance();
    tsn.setSchema(schema);
    return sqlMap.update("updateAttributes", tsn);
  }

  /**
   *
   * Insert attributes defined by the <code>ChadoTransaction</code>.
   * @param schema      schema/organism name or null
   * @param tsn         the <code>ChadoTransaction</code>
   * @throws SQLException
   */
  public void insertAttributes
                    (final String schema, final ChadoTransaction tsn)
                     throws SQLException
  {
    SqlMapClient sqlMap = DbSqlConfig.getSqlMapInstance();
    tsn.setSchema(schema);

    // get the feature id's
    List feature_ids = sqlMap.queryForList("getFeatureID", tsn);

    for(int i=0; i<feature_ids.size(); i++)
    {
      tsn.setFeature_id( ((Integer)feature_ids.get(i)).intValue() );
      sqlMap.insert("insertAttributes", tsn);
    }
  }

  /**
   *
   * Delete attributes defined by the <code>ChadoTransaction</code>.
   * @param schema      schema/organism name or null
   * @param tsn         the <code>ChadoTransaction</code>
   * @throws SQLException
   */
  public void deleteAttributes
                    (final String schema, final ChadoTransaction tsn)
                     throws SQLException
  {
    SqlMapClient sqlMap = DbSqlConfig.getSqlMapInstance();
    tsn.setSchema(schema);
  
    // get the feature id's
    List feature_ids = sqlMap.queryForList("getFeatureID", tsn);

    for(int i=0; i<feature_ids.size(); i++)
    {
      tsn.setFeature_id( ((Integer)feature_ids.get(i)).intValue() );
      sqlMap.delete("deleteAttributes", tsn);
    }
  }

  /**
   *
   * Insert a feature into the database defined by the <code>ChadoTransaction</code>.
   * @param schema              schema/organism name or null
   * @param tsn                 the <code>ChadoTransaction</code>
   * @parma srcfeature_id       the parent feature identifier
   * @throws SQLException
   */
  public void insertFeature
                    (final String schema, final ChadoTransaction tsn,
                     final String srcfeature_id)
                     throws SQLException
  {
    // get the organism id from the srcfeature_id 
    ChadoFeature feature = new ChadoFeature();
    feature.setSchema(schema);
    feature.setSrcfeature_id(Integer.parseInt(srcfeature_id));
    SqlMapClient sqlMap = DbSqlConfig.getSqlMapInstance();
    Integer organism_id = (Integer)sqlMap.queryForObject("getOrganismID", feature);

    //
    // insert feature into feature table
    ChadoFeature chadoFeature = tsn.getChadoFeature();
    chadoFeature.setSchema(schema);
    chadoFeature.setOrganism_id(organism_id.intValue());  
    sqlMap.insert("insertFeature", chadoFeature);

    //
    // get the current feature_id sequence value
    int feature_id = ((Integer)sqlMap.queryForObject("currval", 
                              schema+".feature_feature_id_seq")).intValue();

    //
    // insert feature location into featureloc
    chadoFeature.setSrcfeature_id(Integer.parseInt(srcfeature_id));
    chadoFeature.setId(feature_id);
    sqlMap.insert("insertFeatureLoc", chadoFeature);
  }


  /**
   *
   * Delete a feature from the database defined by the <code>ChadoTransaction</code>.
   * @param schema      schema/organism name or null
   * @param tsn         the <code>ChadoTransaction</code>
   * @return    number of rows deleted
   * @throws SQLException
   */
  public int deleteFeature
                    (final String schema, final ChadoTransaction tsn)
                     throws SQLException
  {
    SqlMapClient sqlMap = DbSqlConfig.getSqlMapInstance();
    ChadoFeature chadoFeature = new ChadoFeature();
    chadoFeature.setSchema(schema);
    chadoFeature.setUniquename(tsn.getUniqueName());

    return sqlMap.delete("deleteFeature", chadoFeature);
  }

  /**
   *
   * Write the time a feature was last modified
   * @param schema      schema/organism name or null
   * @param uniquename  the unique name of the feature
   * @return  number of rows changed
   * @throws SQLException
   */
  public int writeTimeLastModified
                    (final String schema, final String uniquename)
                     throws SQLException
  {
    ChadoTransaction tsn = new ChadoTransaction(ChadoTransaction.UPDATE,
                                                uniquename, "feature");
    tsn.addProperty("timelastmodified", "CURRENT_TIMESTAMP");
    return updateAttributes(schema, tsn);
  }

  /**
   *
   * Write the time a feature was last accessed
   * @param schema      schema/organism name or null
   * @param uniquename  the unique name of the feature
   * @return  number of rows changed
   * @throws SQLException
   */
  public int writeTimeAccessioned
                    (final String schema, final String uniquename)
                     throws SQLException
  {
    ChadoTransaction tsn = new ChadoTransaction(ChadoTransaction.UPDATE,
                                                uniquename, "feature");
    tsn.addProperty("timeaccessioned", "CURRENT_TIMESTAMP");
    return updateAttributes(schema, tsn);
  }
}

