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

public class IBatisDAO implements ChadoDAO
{
  public IBatisDAO(final JPasswordField pfield)
  {
    DbSqlConfig.init(pfield);
  }

  /**
   *
   * Get feature name given the feature_id and schema
   *
   */
  public String getFeatureName(final Feature feature)
                throws SQLException
  {
    SqlMapClient sqlMap = DbSqlConfig.getSqlMapInstance();
    return (String)sqlMap.queryForObject("getFeatureName", feature);
  }

  /**
   *
   * Get feature name given the feature_id and schema.
   * @param feature_id  id of feature to query
   * @param schema      schema/organism name or null
   *
   */
  public String getFeatureName(final int feature_id,
                               final String schema)
                       throws SQLException
  {
    Feature feature = new Feature();
    feature.setId(feature_id);
    if(schema != null)
      feature.setSchema(schema);
    return getFeatureName(feature);
  }

  /**
   *
   * Get child feature properties for a given parent
   * feature to be able to construct a GFF like feature.
   *
   * @param feature_id  id of parent feature to query
   * @param schema      schema/organism name or null
   *
   */
  public List getGff(final int feature_id,
                     final String schema)
                     throws SQLException
  {
    SqlMapClient sqlMap = DbSqlConfig.getSqlMapInstance();
    Feature feature = new Feature();
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
   *
   */
  public Feature getSequence(final int feature_id,
                             final String schema)
                        throws SQLException
  {
    SqlMapClient sqlMap = DbSqlConfig.getSqlMapInstance();
    Feature feature = new Feature();
    feature.setId(feature_id);
    if(schema != null)
      feature.setSchema(schema);
    return (Feature)sqlMap.queryForObject("getSequence",
                                           feature);
  }

  /**
   *
   * Given a list of distict cvterm_id/type_id's of feature types
   * that have residues in the given schema and the schema name
   * return a list of features in the schema with residues.
   * @param cvterm_ids list of cvterm_id/type_id's
   * @param schema      schema/organism name or null
   *
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
   * For a schema return the type_id's with residues
   * @param schema      schema/organism name or null
   * @return list of type_id's
   *
   */
  public List getResidueType(final String schema)
                     throws SQLException
  {
    SqlMapClient sqlMap = DbSqlConfig.getSqlMapInstance();
    return sqlMap.queryForList("getResidueType", schema);
  }

  /**
   *
   * Get available schemas (as a List of Feature objects)
   *
   */
  public List getSchema()
                throws SQLException
  {
    SqlMapClient sqlMap = DbSqlConfig.getSqlMapInstance();
    return sqlMap.queryForList("getSchema", null);
  }

  public List getCvterm()
              throws SQLException
  {
    SqlMapClient sqlMap = DbSqlConfig.getSqlMapInstance();
    return sqlMap.queryForList("getCvterm", null);
  }

  /**
   *
   * @param name cvterm name
   * @param cv_name ontology name (e.g. gene, sequence)
   *
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
   * @param schema schema to update.
   *
   */
  public void updateAttributes
                    (final String schema, final ChadoTransaction tsn)
                     throws SQLException 
  {
    SqlMapClient sqlMap = DbSqlConfig.getSqlMapInstance();
    tsn.setSchema(schema);
    sqlMap.update("updateAttributes", tsn);
  }

  /**
   *
   * @param schema schema to update.
   *
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
   * @param schema schema to update.
   *
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

}

