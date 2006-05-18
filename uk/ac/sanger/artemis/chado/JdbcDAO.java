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
import java.util.Hashtable;

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

  /**
   * Define a JDBC data access object and establish a <code>Connection</code>.
   * @param location	the database location <i>e.g.</i> 
   *                    jdbc:postgresql://localhost:2997/chado?user=tjc
   * @param pfield	the password for this connection
   * @throws SQLException
   */
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


  /**
   * Get the residues of a feature.
   * @param feature_id  id of feature to query
   * @param schema      schema/organism name or null
   * @return    the <code>ChadoFeature</code> with the residues
   * @throws SQLException
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
   * Get the feature name given a feature_id and schema.
   * @param feature_id  id of feature to query
   * @param schema      schema/organism name or null
   * @return    the feature name
   * @throws SQLException
   */
  public String getFeatureName(final int feature_id,
                               final String schema)
                       throws SQLException
  {
    Statement st = conn.createStatement();

    String sql = "SELECT name FROM " + schema + 
                 ".feature WHERE feature_id= " +
                 feature_id;
    appendToLogFile(sql, sqlLog);
    ResultSet rs = st.executeQuery(sql);
    rs.next();
    return rs.getString("name");
  }

  /**
   * Get child feature properties for a given parent
   * feature to be able to construct a GFF like feature.
   * @param parentFeatureID  the id of parent feature to query
   * @param schema           the schema/organism name or null
   * @return    the <code>List</code> of child <code>ChadoFeature</code> objects
   * @throws SQLException
   */
  public List getGff(final int parentFeatureID,
                     final String schema)
                     throws SQLException
  {
    return getFeatureQuery(null, parentFeatureID, schema);
  }

  /**
   * Get the properties of a feature.
   * @param uniquename  the unique name of the feature
   * @param schema_list the <code>List</code> of schemas to search
   * @return  the <code>List</code> of <code>ChadoFeature</code>
   * @throws SQLException
   */
  public List getFeature(final String uniquename,
                         List schema_list)
                         throws SQLException
  {
    List list = new Vector();
    for(int i=0; i<schema_list.size(); i++)
    {
      String schema = (String)schema_list.get(i);
      List feat_list = getFeatureQuery(uniquename, -1, schema);
      list.addAll(feat_list);
    }
     
    return list;
  }
  
  /**
   * Get the properties of a feature.
   * @param uniquename  the unique name of the feature
   * @param parentFeatureID  the id of parent feature to query
   * @param schema_list the <code>List</code> of schemas to search
   * @return  the <code>List</code> of <code>ChadoFeature</code>
   * @throws SQLException
   */
  private List getFeatureQuery(final String uniquename,
                               final int parentFeatureID,
                               final String schema)
                               throws SQLException
  {
    Statement st = conn.createStatement(ResultSet.TYPE_SCROLL_INSENSITIVE,
                                        ResultSet.CONCUR_UPDATABLE);

    String sql = "SELECT timelastmodified, f.feature_id, object_id,"
               + " fl.strand, fmin, fmax, uniquename, f.type_id,"
               + " fp.type_id AS prop_type_id, fp.value, fl.phase,"
               + " f.organism_id, abbreviation, genus, species, common_name, comment"
               + " FROM  "     + schema + ".feature f"
               + " LEFT JOIN " + schema + ".feature_relationship fr ON "
                                        + "fr.subject_id=" + "f.feature_id"
               + " LEFT JOIN " + schema + ".featureprop fp ON "
                                        + "fp.feature_id=" + "f.feature_id"
               + " LEFT JOIN " + schema + ".featureloc fl ON "
                                        + "f.feature_id=" + "fl.feature_id"
               + " LEFT JOIN organism ON organism.organism_id=f.organism_id"                                
               + " WHERE ";
    
    
    if(uniquename != null)
      sql = sql + "uniquename LIKE '" + uniquename +"'";
    
    if(parentFeatureID > -1)
      sql = sql + "srcfeature_id = " + parentFeatureID;
    
    sql = sql  + " AND (fl.rank=fr.rank OR fr.rank IS NULL)"
               + " ORDER BY f.type_id, uniquename";
    
    appendToLogFile(sql, sqlLog);
    ResultSet rs = st.executeQuery(sql);
    
    final List list = new Vector();
    while(rs.next())
    {
      ChadoFeature feature = new ChadoFeature();
      
      ChadoFeatureLoc featureloc = new ChadoFeatureLoc();
      featureloc.setFmin( rs.getInt("fmin") );
      featureloc.setFmax( rs.getInt("fmax") );
      featureloc.setStrand( rs.getInt("strand") );
      
      int phase = rs.getInt("phase");
      if(rs.wasNull())
        featureloc.setPhase(10);
      else 
        featureloc.setPhase(phase);

      feature.setFeatureloc(featureloc);
      feature.setCvterm(new ChadoCvterm());
      feature.getCvterm().setId( rs.getLong("type_id") );
      
      // feature properties
      ChadoFeatureProp featureprop = new ChadoFeatureProp();
      ChadoCvterm cvterm = new ChadoCvterm();
      cvterm.setId(rs.getLong("prop_type_id"));
      featureprop.setCvterm(cvterm);
      featureprop.setValue(rs.getString("value"));
      feature.setFeatureprop(featureprop);

      feature.setSchema(schema);

      feature.setUniquename( rs.getString("uniquename") );
      feature.setTimelastmodified( rs.getTimestamp("timelastmodified") );
      feature.setId( rs.getInt("feature_id") );
      
      // feature relationship
      ChadoFeatureRelationship feature_relationship = new ChadoFeatureRelationship();
      feature_relationship.setObject_id( rs.getInt("object_id") );
      feature.setFeature_relationship(feature_relationship);
  
      // feature organism
      ChadoOrganism organism = new ChadoOrganism();
      organism.setAbbreviation(rs.getString("abbreviation"));
      organism.setComment(rs.getString("comment"));
      organism.setCommon_name(rs.getString("common_name"));
      organism.setGenus(rs.getString("genus"));
      organism.setId(rs.getInt("organism_id"));
      organism.setSpecies(rs.getString("species"));
      feature.setOrganism(organism);
      
      list.add(feature);
    }
    // merge same features in the list
    return mergeList(list);
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

      feat.addQualifier(feat.getFeatureprop().getCvterm().getId(),
                        feat.getFeatureprop());

      if(i < feature_size - 1)
        featNext = (ChadoFeature)list.get(i + 1);

      // merge next line if part of the same feature
      while(featNext != null && featNext.getUniquename().equals(name))
      {
        feat.addQualifier(featNext.getFeatureprop().getCvterm().getId(),
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
  protected static Hashtable mergeDbxref(final List list)
  {
    Hashtable dbxrefHash = new Hashtable();
    for(int i = 0; i < list.size(); i++)
    {
      ChadoDbxref dbxref = (ChadoDbxref)list.get(i);
      Integer feature_id = new Integer(dbxref.getFeature_id());
      String value = dbxref.getName() + ":" + dbxref.getAccession();
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
      ChadoOrganism organism = new ChadoOrganism();
      organism.setAbbreviation( rs.getString("abbreviation") );
      
      feature.setOrganism(organism);
      feature.setId( rs.getInt("feature_id") );
      feature.setName( rs.getString("name") );
      feature.setCvterm(new ChadoCvterm());
      feature.getCvterm().setId( rs.getLong("type_id") );
      
      list.add(feature);
    }
    return list;
  }
 

  /**
   * For a schema return the type_id's with residues.
   * @param schema      schema/organism name or null
   * @return	the <code>List</code> of type_id's as <code>String</code>
   *            objects
   * @throws SQLException
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
   * Get available schemas (as a <code>List</code> of <code>String</code>       
   * objects).
   * @return	the available schemas
   * @throws SQLException
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
   * Get the full list of cvterm_id and name as a <code>List</code> of 
   * <code>Cvterm</code> objects.
   * @return    the full list of cvterm_id and name
   * @throws SQLException
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
      ChadoCvterm cvterm = new ChadoCvterm();
      cvterm.setId( rs.getLong("cvterm_id") );
      cvterm.setName( rs.getString("name") );
      cvterms.add(cvterm);
    }

    return cvterms;
  }

  /**
   * Get dbxref for a feature.
   * @param schema      the postgres schema name
   * @param uniquename  the unique name for the feature. If set to NULL
   *                    all <code>Dbxref</code> are returned.
   * @return a <code>Hashtable</code> of dbxrefs.
   * @throws SQLException
   */
  public Hashtable getDbxref(final String schema, final String uniquename)
              throws SQLException
  {
    String sql = "SELECT db.name, dbx.accession, f.feature_id FROM "+
                 schema+".feature_dbxref dbx_f "+
                 "LEFT JOIN dbxref dbx ON dbx.dbxref_id=dbx_f.dbxref_id "+
                 "LEFT JOIN db ON db.db_id=dbx.db_id "+
                 "LEFT JOIN "+schema+".feature f ON dbx_f.feature_id=f.feature_id ";
    
    if(uniquename != null)
      sql = sql + "WHERE f.uniquename='"+uniquename+"'";
    
    appendToLogFile(sql, sqlLog);
    Statement st = conn.createStatement();
    ResultSet rs = st.executeQuery(sql);
    List dbxrefs = new Vector();

    while(rs.next())
    {
      ChadoDbxref dbxref = new ChadoDbxref();
      dbxref.setName( rs.getString("name") );
      dbxref.setAccession( rs.getString("accession") );
      dbxref.setFeature_id( rs.getInt("feature_id") );
      dbxrefs.add(dbxref);
    }

    return mergeDbxref(dbxrefs);
  }
  
  /**
   * Get dbxref for a feature.
   * @param schema      the postgres schema name
   * @param uniquename  the unique name for the feature. If set to NULL
   *                    all synonyms are returned.
   * @return a <code>Hashtable</code> of synonym values with the 
   *         feature_id as the key.
   * @throws SQLException
   */
  public Hashtable getAlias(final String schema, final String uniquename)
              throws SQLException
  {
    String sql = "SELECT s.name, f.feature_id, cvterm.name AS cvterm_name FROM "+
                                                   schema+".feature_synonym fs "+
                 "LEFT JOIN "+schema+".feature f ON f.feature_id=fs.feature_id "+
                 "LEFT JOIN "+schema+".synonym s ON fs.synonym_id=s.synonym_id "+
                 "LEFT JOIN cvterm ON s.type_id=cvterm_id";
    
    if(uniquename != null)
      sql = sql + " WHERE uniquename='"+uniquename+"'";
    
    appendToLogFile(sql, sqlLog);
    Statement st = conn.createStatement();
    ResultSet rs = st.executeQuery(sql);
    Hashtable synonym = new Hashtable();
    Integer feature_id;
    Vector value;
    ChadoFeatureSynonym alias;
    while(rs.next())
    {
      feature_id = new Integer(rs.getInt("feature_id"));
      if(synonym.containsKey(feature_id))
        value = (Vector)synonym.get(feature_id);
      else
        value = new Vector();
      
      alias = new ChadoFeatureSynonym();
      ChadoCvterm cvterm = new ChadoCvterm();
      cvterm.setName(rs.getString("cvterm_name"));
      
      ChadoSynonym syn = new ChadoSynonym();
      syn.setName( rs.getString("name") );
      syn.setCvterm(cvterm);
      
      alias.setSynonym(syn);
      value.add(alias);
      synonym.put(feature_id, value);
    }
    
    return synonym;
  }
  
  /**
   * Get the time a feature was last modified.
   * @param schema      schema/organism name or null
   * @param uniquename  the unique name of the feature
   * @return  number of rows changed
   * @throws SQLException
   */
  public Timestamp getTimeLastModified
                   (final String schema, final String uniquename)
                   throws SQLException
  {
    StringBuffer sqlBuff = new StringBuffer("SELECT timelastmodified FROM ");
    
    sqlBuff.append(schema);
    sqlBuff.append(".feature WHERE ");
    sqlBuff.append(schema);
    sqlBuff.append(".feature.uniquename='");
    sqlBuff.append(uniquename);
    sqlBuff.append("'");
    
    String sql = sqlBuff.toString();
    appendToLogFile(sql, sqlLog);
    Statement st = conn.createStatement();
    ResultSet rs = st.executeQuery(sql);
    
    boolean result = rs.next();
    
    if(result)
      return rs.getTimestamp("timelastmodified");
    return null;
  }
  
  
//
// WRITE 
//
  /**
   * Update attributes defined by the <code>ChadoTransaction</code>.
   * @param schema      schema/organism name or null
   * @param tsn         the <code>ChadoTransaction</code>
   * @return    number of rows changed
   * @throws SQLException
   */
  public int updateAttributes
                    (final String schema, final ChadoTransaction tsn)
                     throws SQLException
  {
    List uniquename = tsn.getUniquename();
    
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


    for(int i=0; i<uniquename.size(); i++)
    {
      sqlBuff.append(" "+schema+".feature.uniquename='" + 
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

    String sql = sqlBuff.toString();

    appendToLogFile(sql, sqlLog);
    Statement st = conn.createStatement();
    return st.executeUpdate(sql);
  }

  /**
   * Insert attributes defined by the <code>ChadoTransaction</code>.
   * @param schema      schema/organism name or null
   * @param tsn         the <code>ChadoTransaction</code>
   * @throws SQLException
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
   * Delete attributes defined by the <code>ChadoTransaction</code>.
   * @param schema      schema/organism name or null
   * @param tsn         the <code>ChadoTransaction</code>
   * @throws SQLException
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
    //
    // get the organism_id
    Statement st = conn.createStatement();
    String sql = "SELECT organism_id from " + schema +
                 ".feature where feature_id = '" + srcfeature_id + "'";

    appendToLogFile(sql, sqlLog);
    ResultSet rs = st.executeQuery(sql);
    rs.next();

    final int organism_id = rs.getInt("organism_id");

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
    sql_buff.append("'"+chadoFeature.getName()+"'"+" , ");
    sql_buff.append("'"+chadoFeature.getUniquename()+"'"+" , ");
    sql_buff.append(Long.toString(chadoFeature.getCvterm().getId()));
    sql_buff.append(" )");

    sql = new String(sql_buff);
    appendToLogFile(sql, sqlLog);
    st = conn.createStatement();
    int rowCount = st.executeUpdate(sql);

    //
    // get the current feature_id sequence value
    sql = "SELECT currval('"+schema+".feature_feature_id_seq')";
    appendToLogFile(sql, sqlLog);
    
    rs = st.executeQuery(sql);
    rs.next();
    final int feature_id = rs.getInt("currval");
    
    //
    // insert feature location into featureloc
    sql_buff = new StringBuffer();
    sql_buff.append("INSERT INTO "+schema+".featureloc (");
    sql_buff.append(" featureloc_id ,");
    sql_buff.append(" feature_id ,");
    sql_buff.append(" srcfeature_id ,");
    sql_buff.append(" fmin ,");
    sql_buff.append(" fmax ,");
    sql_buff.append(" strand ,");
    sql_buff.append(" phase ");
    sql_buff.append(" ) VALUES ( ");
    sql_buff.append("nextval('"+schema+".featureloc_featureloc_id_seq') , ");
    sql_buff.append(feature_id+" , ");
    sql_buff.append(srcfeature_id+" , ");
    sql_buff.append(chadoFeature.getFeatureloc().getFmin()+" , ");
    sql_buff.append(chadoFeature.getFeatureloc().getFmax()+" , ");
    sql_buff.append(chadoFeature.getFeatureloc().getStrand()+" , ");
    sql_buff.append(chadoFeature.getFeatureloc().getPhase());
    sql_buff.append(" )");

    sql = new String(sql_buff);
    appendToLogFile(sql, sqlLog);
    st = conn.createStatement();
    rowCount = st.executeUpdate(sql);
  }

  /**
   * Delete a feature from the database defined by the 
   * <code>ChadoTransaction</code>.
   * @param schema      schema/organism name or null
   * @param tsn         the <code>ChadoTransaction</code>
   * @return	number of rows deleted
   * @throws SQLException
   */
  public int deleteFeature
                    (final String schema, final ChadoTransaction tsn)
                     throws SQLException
  {
    String sql = "DELETE FROM "+schema+".feature WHERE uniquename='"+
                 tsn.getUniqueName()+"'";
    appendToLogFile(sql, sqlLog);
    
    Statement st = conn.createStatement();
    return st.executeUpdate(sql);
  }

  /**
   * Insert a dbxref for a feature.
   * @param schema        schema/organism name or null
   * @param tsn           the <code>ChadoTransaction</code>
   * @return    number of rows changed
   * @throws SQLException
   */
  public int insertFeatureDbxref(final String schema, final ChadoTransaction tsn)
                     throws SQLException
  {
    final ChadoDbxref dbxref = tsn.getFeatureDbxref();
    final String uniquename  = tsn.getUniqueName();
    
    // find database id
    String sql = "SELECT db_id FROM db WHERE name='"+dbxref.getName()+"'";
    
    Statement st   = conn.createStatement();
    ResultSet rs   = st.executeQuery(sql);
    boolean exists = rs.next();
    
    if(!exists)
      throw new SQLException("No database called "+
                             dbxref.getName()+
                             " found (for "+uniquename+
                             ") check the spelling!");

    final int db_id = rs.getInt("db_id");
    // find if accession exists already
    String sqlDbxrefId = "SELECT dbxref_id FROM dbxref WHERE accession='"+
                          dbxref.getAccession()+"' AND db_id="+db_id;
    appendToLogFile(sqlDbxrefId, sqlLog);
    rs     = st.executeQuery(sqlDbxrefId);
    exists = rs.next();
    
    if(!exists)
    {
      // create a new accession entry in dbxref
      sql = "INSERT INTO dbxref ( db_id, accession ) "+
            "VALUES ("+db_id+", "+dbxref.getAccession()+" )";
      appendToLogFile(sql, sqlLog);
      int rowCount = st.executeUpdate(new String(sql));
      
      // now get the new dbxref_id
      appendToLogFile(sqlDbxrefId, sqlLog);
      rs = st.executeQuery(sqlDbxrefId);
      rs.next();
    }
    
    final int dbxref_id = rs.getInt("dbxref_id"); 
    sql = "INSERT INTO "+schema+".feature_dbxref "+
          "(feature_id, dbxref_id, is_current)"+
          " VALUES "+
          "( (SELECT feature_id FROM "+schema+
             ".feature WHERE  uniquename='"+uniquename+"'), "+
          dbxref_id+", "+ Boolean.toString(dbxref.isCurrent())+")";
    appendToLogFile(sql, sqlLog);
    return st.executeUpdate(new String(sql));
  }
  
  /**
   * Delete a dbxref for a feature.
   * @param schema        schema/organism name or null
   * @param tsn           the <code>ChadoTransaction</code>
   * @return    number of rows changed
   * @throws SQLException
   */
  public int deleteFeatureDbxref(final String schema, final ChadoTransaction tsn)
                     throws SQLException
  {
    final ChadoDbxref dbxref = tsn.getFeatureDbxref();
    final String uniquename = tsn.getUniqueName();
    
    final String sql = 
      "DELETE FROM "+schema+".feature_dbxref "+
      "WHERE dbxref_id="+
      "(SELECT dbxref_id FROM dbxref WHERE accession='"+dbxref.getAccession()+"' "+
             "AND db_id=(SELECT db_id FROM db WHERE name='"+dbxref.getName()+"'))"+
      "AND feature_id=(SELECT feature_id FROM "+schema+
             ".feature WHERE  uniquename='"+uniquename+"')";
    
    Statement st = conn.createStatement();
    return st.executeUpdate(sql);
  }
  
  /**
   * Insert a synonym for a feature.
   * @param schema        schema/organism name or null
   * @param tsn           the <code>ChadoTransaction</code>
   * @return    number of rows changed
   * @throws SQLException
   */
  public int insertFeatureAlias(final String schema, final ChadoTransaction tsn)
                     throws SQLException
  {
    final ChadoFeatureSynonym alias  = tsn.getAlias();
    final String uniquename   = alias.getUniquename();
    final String synonym_name = alias.getSynonym().getName();
      
    String sql;
     
    String sqlAliasId = "SELECT synonym_id FROM "+
                        schema+".synonym WHERE synonym.name='"+synonym_name+"'";

    appendToLogFile(sqlAliasId, sqlLog);
    Statement st   = conn.createStatement(ResultSet.TYPE_SCROLL_SENSITIVE,
                                          ResultSet.CONCUR_UPDATABLE);
    ResultSet rs   = st.executeQuery(sqlAliasId);
    boolean exists = rs.next();
    
    if(!exists)
    {
      // create a new synonym name     
      String type_id = Long.toString(alias.getSynonym().getCvterm().getId());
      
      sql = "INSERT INTO "+schema+
            ".synonym (name, type_id, synonym_sgml) values ( '"+
            synonym_name+"',"+type_id+",'"+synonym_name+"')" ;

      st.executeUpdate(sql);
      appendToLogFile(sql, sqlLog);
      
      rs = st.executeQuery(sqlAliasId);
      rs.next();
      appendToLogFile(sqlAliasId, sqlLog);
    }
    
    final int synonym_id = rs.getInt("synonym_id"); 
    sql = "INSERT INTO "+schema+
           ".feature_synonym ( synonym_id, feature_id, pub_id )"+
           " values ( "+synonym_id+" ,"+
               "(SELECT feature_id FROM "+schema+
               ".feature WHERE  uniquename='"+uniquename+"'), "+
               " 1)";
 
    appendToLogFile(sql, sqlLog);
    return st.executeUpdate(sql);
  }
  
  /**
   * Delete a synonym for a feature.
   * @param schema        schema/organism name or null
   * @param tsn           the <code>ChadoTransaction</code>
   * @return    number of rows changed
   * @throws SQLException
   */
  public int deleteFeatureAlias(final String schema, final ChadoTransaction tsn)
                     throws SQLException
  {
    final ChadoFeatureSynonym alias = tsn.getAlias();
    final String uniquename   = alias.getUniquename();
    final String synonym_name = alias.getSynonym().getName();
    String sql = "SELECT synonym_id FROM "+schema+".feature_synonym WHERE "+ 
                 "synonym_id=(SELECT synonym_id FROM "+schema+".synonym WHERE "+
                 "synonym.name='"+synonym_name+"')";
    
    appendToLogFile(sql, sqlLog);
    Statement st   = conn.createStatement(ResultSet.TYPE_SCROLL_INSENSITIVE,
                                          ResultSet.CONCUR_READ_ONLY);
    ResultSet rs   = st.executeQuery(sql);
    rs.last();
    int nrows = rs.getRow();
    final int synonym_id = rs.getInt("synonym_id"); 
    
    // check this name is not used some where else, 
    // i.e. in more than one row
    if(nrows>1)
    {
      sql = "DELETE FROM "+schema+".feature_synonym WHERE "+
            "synonym_id="+synonym_id+" AND "+
            "feature_id=(SELECT feature_id FROM "+schema+
                     ".feature WHERE  uniquename='"+uniquename+"')";
    }
    else
      sql = "DELETE FROM "+schema+".synonym WHERE synonym_id="+synonym_id;
    
    st   = conn.createStatement();
    return st.executeUpdate(sql);
  }
  
  /**
   * Write the time a feature was last modified
   * @param schema      schema/organism name or null
   * @param uniquename  the unique name of the feature
   * @param timestamp   the time stamp to use, 
   *                    if NULL use CURRENT_TIMESTAMP
   * @return  number of rows changed
   * @throws SQLException
   */
  public int writeTimeLastModified
                    (final String schema, final String uniquename,
                     final Timestamp timestamp)
                     throws SQLException
  {
    String sql = "UPDATE "+schema+
                 ".feature SET timelastmodified=";
    
    if(timestamp == null)
      sql = sql +"CURRENT_TIMESTAMP";
    else
      sql = sql +"?";
      
    sql = sql + " WHERE uniquename= ?";
    
    PreparedStatement pstmt = conn.prepareStatement(sql);

    int param = 1;
    if(timestamp != null)
    {
      pstmt.setTimestamp(param, timestamp);
      param++;
    }
    pstmt.setString(param, uniquename);
    
    int rowCount = pstmt.executeUpdate();
    return rowCount;
  }

  /**
   * Write the time a feature was last accessed
   * @param schema      schema/organism name or null
   * @param uniquename  the unique name of the feature
   * @param timestamp   the time stamp to use, 
   *                    if NULL use CURRENT_TIMESTAMP
   * @throws SQLException
   */
  public int writeTimeAccessioned
                    (final String schema, final String uniquename,
                     final Timestamp timestamp)
                     throws SQLException
  {   
    String sql = "UPDATE "+schema+
                 ".feature SET timeaccessioned=";

    if(timestamp == null)
      sql = sql +"CURRENT_TIMESTAMP";
    else
      sql = sql +"?";

    sql = sql + " WHERE uniquename= ?";

    PreparedStatement pstmt = conn.prepareStatement(sql);

    int param = 1;
    if(timestamp != null)
    {
      pstmt.setTimestamp(param, timestamp);
      param++;
    }
    pstmt.setString(param, uniquename);

    int rowCount = pstmt.executeUpdate();
    return rowCount;
  }

  /**
   * Appends a log entry to the log file
   * @param logEntry    entry to add to log file
   * @param logFileName log file name
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
