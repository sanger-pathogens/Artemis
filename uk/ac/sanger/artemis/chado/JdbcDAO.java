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

package uk.ac.sanger.artemis.chado;

import javax.swing.JPasswordField;
import java.sql.*;
import java.io.*;
import java.util.List;
import java.util.Vector;

/**
 *
 * Java Database Connectivity (JDBC) implemetation of the
 * <code>ChadoDAO</code> data access interface.
 *
 */
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

  /**
   *
   * Get the residues of a feature.
   * @param feature_id  id of feature to query
   * @param schema      schema/organism name or null
   * @return    the <code>ChadoFeature</code> with the residues
   *
   */
  public ChadoFeature getSequence(final int feature_id,
                             final String schema)
                        throws SQLException
  {
    Statement st = conn.createStatement();
    String sql = "SELECT name, residues from " + schema +
                 ".feature where feature_id = '" + feature_id + "'";

    appendToLogFile(sql, sqlLog);

    ResultSet rs = st.executeQuery(sql);
    rs.next();

    ChadoFeature feature = new ChadoFeature();
    feature.setName(rs.getString("name"));
    feature.setResidues(rs.getBytes("residues"));
    return feature;
  }

  /**
   *
   * Get the feature name given a feature_id and schema.
   * @param feature_id  id of feature to query
   * @param schema      schema/organism name or null
   * @return    the feature name
   */
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
   * @param parentFeatureID  the id of parent feature to query
   * @param schema           the schema/organism name or null
   * @return    the <code>List</code> of <code>ChadoFeature</code> objects
   *
   */
  public List getGff(final int parentFeatureID,
                     final String schema)
                     throws SQLException
  {
    Statement st = conn.createStatement(ResultSet.TYPE_SCROLL_INSENSITIVE,
                                        ResultSet.CONCUR_UPDATABLE);

    String sql = "SELECT timelastmodified, f.feature_id, object_id, "
        + "fl.strand, fmin, fmax, uniquename, f.type_id, "
        + "fp.type_id AS prop_type_id, fp.value, fl.phase"
        + " FROM  "
        + schema + ".feature f"
        + " LEFT JOIN " + schema + ".feature_relationship fr ON "
                        + "fr.subject_id=" + "f.feature_id"
        + " LEFT JOIN " + schema + ".featureprop fp ON "
                        + "fp.feature_id=" + "f.feature_id"
        + " LEFT JOIN " + schema + ".featureloc fl ON "
                        + "f.feature_id=" + "fl.feature_id"
        + " WHERE srcfeature_id = "
        + parentFeatureID 
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
      ChadoFeature feature = new ChadoFeature();
      feature.setFmin( rs.getInt("fmin") );
      feature.setFmax( rs.getInt("fmax") );
      feature.setType_id( rs.getLong("type_id") );
      feature.setProp_type_id( rs.getLong("prop_type_id") );
      feature.setStrand( rs.getInt("strand") );
      
      int phase = rs.getInt("phase");
      if(rs.wasNull())
        feature.setPhase(10);
      else 
        feature.setPhase( rs.getInt("phase") );

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
    ChadoFeature featNext  = null;

    for(int i = 0; i < feature_size; i++)
    {
      ChadoFeature feat = (ChadoFeature)list.get(i);
      String name  = feat.getUniquename();

      feat.addQualifier(feat.getProp_type_id(),
                        feat.getValue());

      if(i < feature_size - 1)
        featNext = (ChadoFeature)list.get(i + 1);

      // merge next line if part of the same feature
      while(featNext != null && featNext.getUniquename().equals(name))
      {
        feat.addQualifier(featNext.getProp_type_id(),
                          featNext.getValue());
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
   *
   * Given a list of distict cvterm_id/type_id's of feature types
   * that have residues (from getResidueType()) in the given schema 
   * and the schema name return a list of chado features in the schema
   * with residues.
   * @param cvterm_ids list of cvterm_id/type_id's
   * @param schema      schema/organism name or null
   * @return    the <code>List</code> of <code>ChadoFeature</code> objects
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
      ChadoFeature feature = new ChadoFeature();
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
   * For a schema return the type_id's with residues.
   * @param schema      schema/organism name or null
   * @return	the <code>List</code> of type_id's as <code>String</code>
   *            objects
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
   * Get available schemas (as a <code>List</code> of <code>String</code>       
   * objects).
   * @return	the available schemas
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

  /**
   *
   * Get the full list of cvterm_id and name as a <code>List</code> of 
   * <code>Cvterm</code> objects.
   * @return    the full list of cvterm_id and name
   *
   */
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

//
// WRITE 
//
  /**
   *
   * Update attributes defined by the <code>ChadoTransaction</code>.
   * @param schema      schema/organism name or null
   * @param tsn         the <code>ChadoTransaction</code>
   *
   */
  public void updateAttributes
                    (final String schema, final ChadoTransaction tsn)
                     throws SQLException
  {
    StringBuffer sqlBuff = new StringBuffer();
    String chadoTable    = tsn.getChadoTable();
    sqlBuff.append("UPDATE "+schema+"."+chadoTable);
    sqlBuff.append(" SET ");

    List properties = tsn.getProperties();
    for(int i=0; i<properties.size(); i++)
    {
      sqlBuff.append((String)properties.get(i));
      if(i < properties.size()-1)
        sqlBuff.append(" , ");
    }

    sqlBuff.append(" FROM "+schema+".feature");
    sqlBuff.append(" WHERE "+schema+".feature.feature_id="+
                             schema+"."+chadoTable+".feature_id AND (");

    List uniquename = tsn.getUniquename();
    for(int i=0; i<uniquename.size(); i++)
    {
      sqlBuff.append(" "+schema+"."+"feature.uniquename='" + 
                     (String)uniquename.get(i) +"' ");
      if(i < uniquename.size()-1)
        sqlBuff.append("OR");
    }

    sqlBuff.append(")");

    List constraints = tsn.getConstraint();
    if(constraints != null)
    {
      for(int i=0; i<constraints.size(); i++)
      {
        sqlBuff.append(" AND ");
        // looks like specifying table, so include schema
        String constraint = (String)constraints.get(i);
        int index;
        if( (index = constraint.indexOf(".")) > -1 &&
            constraint.indexOf("=") > index)
         sqlBuff.append(schema+".");
        sqlBuff.append(constraint);
      }
    }

    Statement st = conn.createStatement();
    int rowCount = st.executeUpdate(new String(sqlBuff));

    if(rowCount > 0)
    {
      for(int i=0; i<uniquename.size(); i++)
        writeTimeLastModified(schema, (String)uniquename.get(i));
    }
  }

  /**
   *
   * Insert attributes defined by the <code>ChadoTransaction</code>.
   * @param schema      schema/organism name or null
   * @param tsn         the <code>ChadoTransaction</code>
   *
   */
  public void insertAttributes
                    (final String schema, final ChadoTransaction tsn)
                     throws SQLException
  {
    StringBuffer sqlBuff;

    List uniquename = tsn.getUniquename();
    String chadoTable   = tsn.getChadoTable();
    for(int i=0; i<uniquename.size(); i++)
    {
      sqlBuff = new StringBuffer();
      sqlBuff.append("INSERT INTO "+schema+"."+chadoTable);
      StringBuffer sqlKeys   = new StringBuffer();
      StringBuffer sqlValues = new StringBuffer();

      sqlKeys.append("feature_id , ");
      sqlValues.append("(SELECT feature_id FROM "+schema+".feature WHERE uniquename='"+
                         (String)uniquename.get(i)+"') , ");

      String name;
      String property;
      List propertiesName  = tsn.getPropertiesName();
      List propertiesValue = tsn.getPropertiesValue();
      for(int j=0; j<propertiesName.size(); j++)
      {
        name = (String)propertiesName.get(j);
        sqlKeys.append(name);
        sqlValues.append((String)propertiesValue.get(j));
        if(j < propertiesName.size()-1)
        {
          sqlKeys.append(" , ");
          sqlValues.append(" , ");
        }
      }
  
      sqlBuff.append(" ( "+sqlKeys.toString()+" ) ");
      sqlBuff.append(" values ");
      sqlBuff.append(" ( "+sqlValues.toString()+" )");

      appendToLogFile(new String(sqlBuff), sqlLog);

      Statement st = conn.createStatement();
      int rowCount = st.executeUpdate(new String(sqlBuff));
      System.out.println(rowCount+" row(s) inserted");
    }
  }

  /**
   *
   * Delete attributes defined by the <code>ChadoTransaction</code>.
   * @param schema      schema/organism name or null
   * @param tsn         the <code>ChadoTransaction</code>
   *
   */
  public void deleteAttributes
                    (final String schema, final ChadoTransaction tsn)
                     throws SQLException
  {
    StringBuffer sqlBuff;

    List uniquename = tsn.getUniquename();
    String chadoTable   = tsn.getChadoTable();
    for(int i=0; i<uniquename.size(); i++)
    {
      sqlBuff = new StringBuffer();

      sqlBuff.append("DELETE FROM "+schema+"."+chadoTable+" WHERE ");

      String name;
      String value;

      List constraints = tsn.getConstraint();
      for(int j=0; j<constraints.size(); j++)
        sqlBuff.append((String)constraints.get(j)+" AND ");
      
      sqlBuff.append("feature_id=(SELECT feature_id FROM "+
                     schema+".feature WHERE uniquename='"+
                     (String)uniquename.get(i)+"')");

      appendToLogFile(new String(sqlBuff), sqlLog);

      Statement st = conn.createStatement();
      int rowCount = st.executeUpdate(new String(sqlBuff));
      System.out.println(rowCount+" row(s) deleted");
    }
  }


  /**
   *
   * Insert a feature into the database defined by the <code>ChadoTransaction</code>.
   * @param schema              schema/organism name or null
   * @param tsn                 the <code>ChadoTransaction</code>
   * @parma srcfeature_id       the parent feature identifier
   *
   */
  public void insertFeature
                    (final String schema, final ChadoTransaction tsn, 
                     final String srcfeature_id)
                     throws SQLException
  {
    //
    // get the organism_id
    Statement st = conn.createStatement();
    String str_sql = "SELECT organism_id from " + schema +
                 ".feature where feature_id = '" + srcfeature_id + "'";

    System.out.println(str_sql);
    appendToLogFile(str_sql, sqlLog);

    ResultSet rs = st.executeQuery(str_sql);
    rs.next();

    ChadoFeature feature = new ChadoFeature();
    final int organism_id = rs.getInt("organism_id");

    //
    // insert feature into feature table
    ChadoFeature chadoFeature = tsn.getChadoFeature();
    // insert new feature into feature table
    StringBuffer sql_buff = new StringBuffer();
    sql_buff.append("INSERT INTO "+schema+".feature (");
    sql_buff.append(" feature_id ,");
    sql_buff.append(" organism_id ,");
    sql_buff.append(" name ,");
    sql_buff.append(" uniquename ,");
    sql_buff.append(" type_id");
    sql_buff.append(" ) VALUES ( ");
    sql_buff.append("nextval('"+schema+".feature_feature_id_seq') , ");
    sql_buff.append(organism_id+" , ");
    sql_buff.append("'"+chadoFeature.getUniquename()+"'"+" , ");
    sql_buff.append("'"+chadoFeature.getUniquename()+"'"+" , ");
    sql_buff.append(Long.toString(chadoFeature.getType_id()));
    sql_buff.append(" )");

    System.out.println(new String(sql_buff));
    st = conn.createStatement();
    int rowCount = st.executeUpdate(new String(sql_buff));

    //
    // get the current feature_id sequence value
    String sql = "SELECT currval('"+schema+".feature_feature_id_seq')";
    System.out.println(sql);
    
    rs = st.executeQuery(sql);
    rs.next();
    final int feature_id = rs.getInt("currval");

//  System.out.println("SELECT currval('"+schema+".featureprop_featureprop_id_seq')");

    //
    // insert feature location into featureloc
    sql_buff = new StringBuffer();
    sql_buff.append("INSERT INTO "+schema+".featureloc (");
    sql_buff.append(" featureloc_id ,");
    sql_buff.append(" feature_id ,");
    sql_buff.append(" srcfeature_id ,");
    sql_buff.append(" fmin ,");
    sql_buff.append(" fmax ");
    sql_buff.append(" ) VALUES ( ");
    sql_buff.append("nextval('"+schema+".featureloc_featureloc_id_seq') , ");
    sql_buff.append(feature_id+" , ");
    sql_buff.append(srcfeature_id+" , ");
    sql_buff.append(chadoFeature.getFmin()+" , ");
    sql_buff.append(chadoFeature.getFmax());
    sql_buff.append(" )");

    System.out.println(new String(sql_buff));
    st = conn.createStatement();
    rowCount = st.executeUpdate(new String(sql_buff));
  }

  /**
   *
   * Delete a feature from the database defined by the <code>ChadoTransaction</code>.
   * @param schema      schema/organism name or null
   * @param tsn         the <code>ChadoTransaction</code>
   *
   */
  public void deleteFeature
                    (final String schema, final ChadoTransaction tsn)
                     throws SQLException
  {
    String sql = "DELETE FROM "+schema+".feature WHERE uniquename='"+
                 tsn.getUniqueName()+"'";

    Statement st = conn.createStatement();
    int rowCount = st.executeUpdate(sql);
    System.out.println(sql);
  }

  
  /**
   *
   * Write the time a feature was last modified
   * @param schema      schema/organism name or null
   * @param uniquename  the unique name of the feature
   *
   */
  public void writeTimeLastModified
                    (final String schema, final String uniquename)
                     throws SQLException
  {
    String sql = "UPDATE "+schema+
                 ".feature SET timelastmodified=CURRENT_TIMESTAMP WHERE uniquename='"+
                 uniquename+"'";
    Statement st = conn.createStatement();
    int rowCount = st.executeUpdate(sql);
  }

  /**
   *
   * Write the time a feature was last accessed
   * @param schema      schema/organism name or null
   * @param uniquename  the unique name of the feature
   *
   */
  public void writeTimeAccessioned
                    (final String schema, final String uniquename)
                     throws SQLException
  {
    String sql = "UPDATE "+schema+
                 ".feature SET timeaccessioned=CURRENT_TIMESTAMP WHERE uniquename='"+
                 uniquename+"'";
    Statement st = conn.createStatement();
    int rowCount = st.executeUpdate(sql);
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
