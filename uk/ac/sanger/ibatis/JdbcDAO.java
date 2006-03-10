/* 
 *
 * created: 2006
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2006  Genome Research Limited
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

package uk.ac.sanger.ibatis;

import javax.swing.JPasswordField;
import java.sql.*;
import java.io.*;
import java.util.List;
import java.util.Vector;

public class JdbcDAO 
             implements ChadoDAO
{

  private String sqlLog = System.getProperty("user.home") +
                          System.getProperty("file.separator") +
                          "art_sql_debug.log";
  private Connection conn;

  public JdbcDAO(final String location, final JPasswordField pfield)
         throws java.sql.SQLException, java.net.ConnectException
  {
    if(pfield == null || pfield.getPassword().length == 0)
      conn = DriverManager.getConnection(location);

    // assume we have a password
    final int index = location.indexOf("?user=");
    conn = DriverManager.getConnection(location.substring(0, index),
                                       location.substring(index + 6),
                                       new String(pfield.getPassword()));
  }

  
  public Connection getConnection()
  {
    return conn;
  }

  public Feature getSequence(final int feature_id,
                             final String schema)
                        throws SQLException
  {
    Statement st = conn.createStatement();
    String sql = "SELECT name, residues from " + schema +
                 ".feature where feature_id = '" + feature_id + "'";

    appendToLogFile(sql, sqlLog);

    ResultSet rs = st.executeQuery(sql);
    rs.next();

    Feature feature = new Feature();
    feature.setName(rs.getString("name"));
    feature.setResidues(rs.getBytes("residues"));
    return feature;
  }

  
  public String getFeatureName(final int feature_id,
                               final String schema)
                       throws SQLException
  {
    Statement st = conn.createStatement();

    String sql = "SELECT name FROM " + schema + ".feature WHERE feature_id= " +
                  feature_id;
    appendToLogFile(sql, sqlLog);
    ResultSet rs = st.executeQuery(sql);
    rs.next();
    return rs.getString("name");
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
  public List getGff(final int parentFeatureID,
                     final String schema)
                     throws SQLException
  {
    Statement st = conn.createStatement(ResultSet.TYPE_SCROLL_INSENSITIVE,
                                        ResultSet.CONCUR_UPDATABLE);

    String sql = "SELECT timelastmodified, f.feature_id, object_id, "
        + "strand, fmin, fmax, uniquename, f.type_id, "
        + schema + ".featureprop.type_id AS prop_type_id, featureprop.value"
        + " FROM  "
        + schema + ".featureloc fl, "
        + schema + ".feature f"
        + " LEFT JOIN "
        + schema + ".feature_relationship fr ON "
        + "fr.subject_id="
        + "f.feature_id"
        + " LEFT JOIN "
        + schema + ".featureprop ON "
        + schema + ".featureprop.feature_id="
        + "f.feature_id"
        + " WHERE srcfeature_id = "
        + parentFeatureID + " AND "
        + "fl.feature_id="
        + "f.feature_id"
        + " AND ("
        + "fl.rank="
        + "fr.rank OR "
        + "fr.rank IS NULL)"
        + " ORDER BY "
        + "f.type_id, uniquename";

    appendToLogFile(sql, sqlLog);
    ResultSet rs = st.executeQuery(sql);

    List list = new Vector();
    while(rs.next())
    {
      Feature feature = new Feature();
      feature.setFmin( rs.getInt("fmin") );
      feature.setFmax( rs.getInt("fmax") );
      feature.setType_id( rs.getLong("type_id") );
      feature.setProp_type_id( rs.getLong("prop_type_id") );
      feature.setStrand( rs.getInt("strand") );
      feature.setUniquename( rs.getString("uniquename") );
      feature.setTimelastmodified( rs.getDate("timelastmodified") );
      feature.setId( rs.getInt("feature_id") );
      feature.setObject_id( rs.getString("object_id") );
      feature.setValue( rs.getString("value"));

      list.add(feature);
    }

    // merge same features in the list
    return mergeList(list);
  }

  /**
   *
   * Takes a list and creates a new one merging all feature objects
   * within it with the same feature and stores the qualifiers/attributes
   *  as a hash
   * @param list of feature objects
   * @return list of flattened/merged feature objects
   * 
   */
  protected static List mergeList(List list)
  {
    // merge same features in the list
    int feature_size  = list.size();
    List flatten_list = new Vector();
    Feature featNext  = null;

    for(int i = 0; i < feature_size; i++)
    {
      Feature feat = (Feature)list.get(i);
      String name  = feat.getUniquename();

      feat.addQualifier(feat.getProp_type_id(),
                        feat.getValue());

      if(i < feature_size - 1)
        featNext = (Feature)list.get(i + 1);

      // merge next line if part of the same feature
      while(featNext != null && featNext.getUniquename().equals(name))
      {
        feat.addQualifier(featNext.getProp_type_id(),
                          featNext.getValue());
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
    String sql = new String(
            "SELECT abbreviation, name, feature_id, type_id FROM organism, "+
            schema + ".feature WHERE (");

    for(int j = 0; j < cvterm_ids.size(); j++)
    {
      sql = sql + " type_id = " + (String)cvterm_ids.get(j);
      if(j < cvterm_ids.size() - 1)
        sql = sql + " OR ";
    }

    sql = sql + ") and organism.organism_id=" + schema
            + ".feature.organism_id " + "and residues notnull "
            + "ORDER BY abbreviation";

    appendToLogFile(sql, sqlLog);

    Statement st = conn.createStatement();
    ResultSet rs = st.executeQuery(sql);
    List list = new Vector();
    while(rs.next())
    {
      Feature feature = new Feature();
      feature.setId( rs.getInt("feature_id") );
      feature.setAbbreviation( rs.getString("abbreviation") );
      feature.setName( rs.getString("name") );
      feature.setType_id( rs.getLong("type_id") );
      
      list.add(feature);
    }
    return list;
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
    String sql = "SELECT DISTINCT type_id FROM " + schema +
                     ".feature WHERE residues notnull";
    appendToLogFile(sql, sqlLog);

    List cvterm_ids = new Vector();
    Statement st = conn.createStatement();
    ResultSet rs = st.executeQuery(sql);

    while(rs.next())
      cvterm_ids.add(rs.getString("type_id"));

    return cvterm_ids;
  }

  
  /**
   * 
   * Get available schemas (as a List of Feature objects)
   *
   */
  public List getSchema()
                throws SQLException
  {
    Statement st = conn.createStatement();

    String query = "SELECT schema_name FROM information_schema.schemata "+
                   "WHERE schema_name=schema_owner ORDER BY schema_name";
    appendToLogFile(query, sqlLog);

    ResultSet rs = st.executeQuery(query);
    List schemas = new Vector();

    while(rs.next())
      schemas.add(rs.getString("schema_name"));

    return schemas;
  } 


  public List getCvterm()
              throws SQLException
  {
    String sql = "SELECT cvterm.cvterm_id, cvterm.name " +
                 "FROM cvterm, cv WHERE cv.cv_id = cvterm.cv_id";

    appendToLogFile(sql, sqlLog);
    Statement st = conn.createStatement();
    ResultSet rs = st.executeQuery(sql);
    List cvterms = new Vector();

    while(rs.next())
    {
      Cvterm cvterm = new Cvterm();
      cvterm.setId( rs.getLong("cvterm_id") );
      cvterm.setName( rs.getString("name") );
      cvterms.add(cvterm);
    }

    return cvterms;
  }


  /**
   *
   * Appends a log entry to the log file
   *
   * @param logEntry
   *          entry to add to log file
   * @param logFileName
   *          log file name
   *
   */
  private void appendToLogFile(String logEntry, String logFileName)
  {
    if(System.getProperty("debug") == null)
      return;

    BufferedWriter bw = null;
    try
    {
      String dat = new java.util.Date().toString();
      bw = new BufferedWriter(new FileWriter(logFileName, true));
      bw.write(dat + ":: " + logEntry);
      bw.newLine();
      bw.flush();
    }
    catch(Exception ioe)
    {
      System.out.println("Error writing to log file " + logFileName);
      ioe.printStackTrace();
    }
    finally
    // always close the file
    {
      if(bw != null)
        try
        {
          bw.close();
        }
        catch(IOException ioe2)
        {
        }
    }
  }


}
