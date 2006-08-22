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

  /**
   * Return a features with this systematic id
   *  
   * @param name the systematic id
   * @return the Feature, or null
   * @throws SQLException 
   */
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
    return getFeatureQuery(null, 
                 feature.getFeatureloc().getSrcfeature_id(), -1);
  }

  /**
   * Return a list of features with any current (ie non-obsolete) name or synonym
   *  
   * @param name the lookup name
   * @return a (possibly empty) List<Feature> of children with this current name
   * @throws SQLException 
   */
  public List getFeaturesByAnyCurrentName(String name) 
              throws SQLException
  {
    ChadoFeature feature = new ChadoFeature();
    feature.setUniquename(name);
    return getFeatureQuery(name, -1, -1);
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
    List list = getFeatureQuery(feature.getUniquename(), 
                                -1, feature.getId());
    return (ChadoFeature)list.get(0);
  }
  
  /**
   * Get the properties of a feature.
   * @param uniquename  the unique name of the feature
   * @param parentFeatureID  the id of parent feature to query
   * @return  the <code>List</code> of <code>ChadoFeature</code>
   * @throws SQLException
   */
  private List getFeatureQuery(final String uniquename,
                               final int parentFeatureID, 
                               final int feature_id)
                               throws SQLException
  {
    Statement st = conn.createStatement(ResultSet.TYPE_SCROLL_INSENSITIVE,
                                        ResultSet.CONCUR_UPDATABLE);

    String sql = "SELECT timelastmodified, f.feature_id, residues,"
               + " fl.strand, fmin, fmax, uniquename, f.type_id,"
               + " fp.type_id AS prop_type_id, fp.value, fl.phase,"
               + " f.organism_id, abbreviation, genus, species, common_name, comment,"
               + " fr.object_id, fr.type_id AS relation_type_id, fr.rank"
               + " FROM  feature f"
               + " LEFT JOIN feature_relationship fr ON "
                                        + "fr.subject_id=" + "f.feature_id"
               + " LEFT JOIN featureprop fp ON "
                                        + "fp.feature_id=" + "f.feature_id"
               + " LEFT JOIN featureloc fl ON "
                                        + "f.feature_id=" + "fl.feature_id"
               + " LEFT JOIN organism ON organism.organism_id=f.organism_id"                                
               + " WHERE ";
    
    
    if(uniquename != null)
      sql = sql + "uniquename LIKE '" + uniquename +"'";
    
    if(parentFeatureID > -1)
      sql = sql + "srcfeature_id = " + parentFeatureID;
    
    if(feature_id > -1)
      sql = sql + "f.feature_id = " + feature_id;
    
    sql = sql  //+ " AND (fl.rank=fr.rank OR fr.rank IS NULL)"
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
      
      feature.setResidues( rs.getBytes("residues") );

      feature.setFeatureloc(featureloc);
      feature.setCvterm(new ChadoCvterm());
      feature.getCvterm().setCvtermId( rs.getLong("type_id") );
      
      // feature properties
      ChadoFeatureProp featureprop = new ChadoFeatureProp();
      ChadoCvterm cvterm = new ChadoCvterm();
      cvterm.setCvtermId(rs.getLong("prop_type_id"));
      featureprop.setCvterm(cvterm);
      featureprop.setValue(rs.getString("value"));
      feature.setFeatureprop(featureprop);

      feature.setUniquename( rs.getString("uniquename") );
      feature.setTimelastmodified( rs.getTimestamp("timelastmodified") );
      feature.setId( rs.getInt("feature_id") );
      
      // feature relationship
      ChadoFeatureRelationship feature_relationship = new ChadoFeatureRelationship();
      cvterm = new ChadoCvterm();
      cvterm.setCvtermId(rs.getLong("relation_type_id"));
      feature_relationship.setCvterm(cvterm);
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
    return IBatisDAO.mergeList(list);
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
            "SELECT abbreviation, uniquename, name, feature_id, type_id FROM organism, "+
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
      feature.setUniquename( rs.getString("uniquename") );
      feature.setCvterm(new ChadoCvterm());
      feature.getCvterm().setCvtermId( rs.getLong("type_id") );
      
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
    String sql = "SELECT DISTINCT type_id FROM " +schema+
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
      cvterm.setCvtermId( rs.getLong("cvterm_id") );
      cvterm.setName( rs.getString("name") );
      cvterms.add(cvterm);
    }

    return cvterms;
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
    String sql = "SELECT db.name, dbx.accession, f.feature_id FROM "+
                 "feature_dbxref dbx_f "+
                 "LEFT JOIN dbxref dbx ON dbx.dbxref_id=dbx_f.dbxref_id "+
                 "LEFT JOIN db ON db.db_id=dbx.db_id "+
                 "LEFT JOIN feature f ON dbx_f.feature_id=f.feature_id ";
    
    if(uniquename != null)
      sql = sql + "WHERE f.uniquename='"+uniquename+"'";
    
    appendToLogFile(sql, sqlLog);
    Statement st = conn.createStatement();
    ResultSet rs = st.executeQuery(sql);
    List dbxrefs = new Vector();

    while(rs.next())
    {
      ChadoFeatureDbxref feature_dbxref = new ChadoFeatureDbxref();
      ChadoDbxref dbxref = new ChadoDbxref();
      ChadoDb db = new ChadoDb();
      db.setName( rs.getString("name") );
      dbxref.setAccession( rs.getString("accession") );
      dbxref.setDb(db);
      feature_dbxref.setDbxref(dbxref);
      feature_dbxref.setFeature_id( rs.getInt("feature_id") );
      dbxrefs.add(feature_dbxref);
    }

    return dbxrefs;
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
    String sql = "SELECT fs.*, s.name, s.type_id FROM "+
    "feature_synonym fs "+
    "LEFT JOIN feature f ON f.feature_id=fs.feature_id "+
    "LEFT JOIN synonym s ON fs.synonym_id=s.synonym_id ";
  
    if(uniquename != null)
      sql = sql + " WHERE uniquename='" + uniquename + "'";

    appendToLogFile(sql, sqlLog);
    Statement st = conn.createStatement();
    ResultSet rs = st.executeQuery(sql);
    List synonym = new Vector();
    ChadoFeatureSynonym alias;
    while(rs.next())
    {
      alias = new ChadoFeatureSynonym();
      ChadoCvterm cvterm = new ChadoCvterm();
      cvterm.setCvtermId(rs.getLong("type_id"));
      ChadoSynonym syn = new ChadoSynonym();
      syn.setName(rs.getString("name"));
      syn.setCvterm(cvterm);

      alias.setSynonym(syn);
      alias.setFeature_id(new Integer(rs.getInt("feature_id")));
      alias.setPub_id(new Integer(rs.getInt("pub_id")));
      alias.setInternal(rs.getBoolean("is_internal"));
      alias.setCurrent(rs.getBoolean("is_current"));
      synonym.add(alias);
    }

    return synonym;
  }
  
  public List getFeatureSynonymsByFeatureAndSynonym(
      ChadoFeature feature, ChadoSynonym synonym)
  {
    return null;
  }


  public ChadoSynonym getSynonymByNameAndCvTerm(
      String name, ChadoCvterm type)
  {
    return null;
  }
  
//
// WRITE 
//
  /**
   * Update attributes defined by the <code>ChadoTransaction</code>.
   * @param tsn         the <code>ChadoTransaction</code>
   * @return    number of rows changed
   * @throws SQLException
   */
  public int updateAttributes
                    (final ChadoTransaction tsn)
                     throws SQLException
  {
    List uniquename = tsn.getUniquename();
    
    StringBuffer sqlBuff = new StringBuffer();
    String chadoTable    = tsn.getChadoTable();
    sqlBuff.append("UPDATE "+chadoTable);
    sqlBuff.append(" SET ");

    List properties = tsn.getProperties();
    for(int i=0; i<properties.size(); i++)
    {
      sqlBuff.append((String)properties.get(i));
      if(i < properties.size()-1)
        sqlBuff.append(" , ");
    }

    sqlBuff.append(" FROM feature");
    sqlBuff.append(" WHERE feature.feature_id="+
                           chadoTable+".feature_id AND (");


    for(int i=0; i<uniquename.size(); i++)
    {
      sqlBuff.append(" feature.uniquename='" + 
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
        /*
        int index;
        if( (index = constraint.indexOf(".")) > -1 &&
            constraint.indexOf("=") > index)
         sqlBuff.append(schema+".");
        */
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
   * @param tsn         the <code>ChadoTransaction</code>
   * @throws SQLException
   */
  public void insertAttributes
                    (final ChadoTransaction tsn)
                     throws SQLException
  {
    StringBuffer sqlBuff;

    List uniquename = tsn.getUniquename();
    String chadoTable   = tsn.getChadoTable();
    for(int i=0; i<uniquename.size(); i++)
    {
      sqlBuff = new StringBuffer();
      sqlBuff.append("INSERT INTO "+chadoTable);
      StringBuffer sqlKeys   = new StringBuffer();
      StringBuffer sqlValues = new StringBuffer();

      sqlKeys.append("feature_id , ");
      sqlValues.append("(SELECT feature_id FROM feature WHERE uniquename='"+
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
   * @param tsn         the <code>ChadoTransaction</code>
   * @throws SQLException
   */
  public void deleteAttributes
                    (final ChadoTransaction tsn)
                     throws SQLException
  {
    StringBuffer sqlBuff;

    List uniquename = tsn.getUniquename();
    String chadoTable   = tsn.getChadoTable();
    for(int i=0; i<uniquename.size(); i++)
    {
      sqlBuff = new StringBuffer();

      sqlBuff.append("DELETE FROM "+chadoTable+" WHERE ");

      List constraints = tsn.getConstraint();
      for(int j=0; j<constraints.size(); j++)
        sqlBuff.append((String)constraints.get(j)+" AND ");
      
      sqlBuff.append("feature_id=(SELECT feature_id FROM "+
                     "feature WHERE uniquename='"+
                     (String)uniquename.get(i)+"')");

      System.out.println(sqlBuff.toString());
      appendToLogFile(new String(sqlBuff), sqlLog);

      Statement st = conn.createStatement();
      int rowCount = st.executeUpdate(new String(sqlBuff));
      System.out.println(rowCount+" row(s) deleted");
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
    //
    // get the organism_id
    Statement st = conn.createStatement();
    String sql = "SELECT organism_id from " +
                 "feature where feature_id = '" + srcfeature_id + "'";

    appendToLogFile(sql, sqlLog);
    ResultSet rs = st.executeQuery(sql);
    rs.next();

    final int organism_id = rs.getInt("organism_id");

    ChadoFeature chadoFeature = tsn.getChadoFeature();
    // insert new feature into feature table
    StringBuffer sql_buff = new StringBuffer();
    sql_buff.append("INSERT INTO feature (");
    sql_buff.append(" feature_id ,");
    sql_buff.append(" organism_id ,");
    sql_buff.append(" name ,");
    sql_buff.append(" uniquename ,");
    sql_buff.append(" type_id");
    sql_buff.append(" ) VALUES ( ");
    sql_buff.append("nextval('feature_feature_id_seq') , ");
    sql_buff.append(organism_id+" , ");
    sql_buff.append("'"+chadoFeature.getName()+"'"+" , ");
    sql_buff.append("'"+chadoFeature.getUniquename()+"'"+" , ");
    sql_buff.append(Long.toString(chadoFeature.getCvterm().getCvtermId()));
    sql_buff.append(" )");

    sql = new String(sql_buff);
    appendToLogFile(sql, sqlLog);
    st = conn.createStatement();
    int rowCount = st.executeUpdate(sql);

    //
    // get the current feature_id sequence value
    sql = "SELECT currval('feature_feature_id_seq')";
    appendToLogFile(sql, sqlLog);
    
    rs = st.executeQuery(sql);
    rs.next();
    final int feature_id = rs.getInt("currval");
    
    //
    // insert feature location into featureloc
    sql_buff = new StringBuffer();
    sql_buff.append("INSERT INTO featureloc (");
    sql_buff.append(" featureloc_id ,");
    sql_buff.append(" feature_id ,");
    sql_buff.append(" srcfeature_id ,");
    sql_buff.append(" fmin ,");
    sql_buff.append(" fmax ,");
    sql_buff.append(" strand ,");
    sql_buff.append(" phase ");
    sql_buff.append(" ) VALUES ( ");
    sql_buff.append("nextval('featureloc_featureloc_id_seq') , ");
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
   * @param tsn         the <code>ChadoTransaction</code>
   * @return	number of rows deleted
   * @throws SQLException
   */
  public int deleteFeature
                    (final ChadoTransaction tsn)
                     throws SQLException
  {
    String sql = "DELETE FROM feature WHERE uniquename='"+
                 tsn.getUniqueName()+"'";
    appendToLogFile(sql, sqlLog);
    
    Statement st = conn.createStatement();
    return st.executeUpdate(sql);
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
    final ChadoFeatureDbxref dbxref = tsn.getFeatureDbxref();
    final String uniquename  = tsn.getUniqueName();
    
    // find database id
    String sql = "SELECT db_id FROM db WHERE name='"+
                 dbxref.getDbxref().getDb().getName()+"'";
    
    Statement st   = conn.createStatement();
    ResultSet rs   = st.executeQuery(sql);
    boolean exists = rs.next();
    
    if(!exists)
      throw new SQLException("No database called "+
                             dbxref.getDbxref().getDb().getName()+
                             " found (for "+uniquename+
                             ") check the spelling!");

    final int db_id = rs.getInt("db_id");
    // find if accession exists already
    String sqlDbxrefId = "SELECT dbxref_id FROM dbxref WHERE accession='"+
                          dbxref.getDbxref().getAccession()+"' AND db_id="+db_id;
    
    appendToLogFile(sqlDbxrefId, sqlLog);
    rs     = st.executeQuery(sqlDbxrefId);
    exists = rs.next();
    
    if(!exists)
    {
      // create a new accession entry in dbxref
      sql = "INSERT INTO dbxref ( db_id, accession ) "+
            "VALUES ("+db_id+", '"+dbxref.getDbxref().getAccession()+"' )";
      
      appendToLogFile(sql, sqlLog);
      int rowCount = st.executeUpdate(new String(sql));
      
      // now get the new dbxref_id
      appendToLogFile(sqlDbxrefId, sqlLog);
      rs = st.executeQuery(sqlDbxrefId);
      rs.next();
    }
    
    final int dbxref_id = rs.getInt("dbxref_id"); 
    sql = "INSERT INTO feature_dbxref "+
          "(feature_id, dbxref_id, is_current)"+
          " VALUES "+
          "( (SELECT feature_id FROM "+
             "feature WHERE  uniquename='"+uniquename+"'), "+
          dbxref_id+", "+ Boolean.toString(dbxref.isCurrent())+")";
    System.out.println(sql);
    appendToLogFile(sql, sqlLog);
    return st.executeUpdate(new String(sql));
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
    final ChadoFeatureDbxref dbxref = tsn.getFeatureDbxref();
    final String uniquename = tsn.getUniqueName();
    
    final String sql = 
      "DELETE FROM feature_dbxref "+
      "WHERE dbxref_id="+
      "(SELECT dbxref_id FROM dbxref WHERE accession='"+dbxref.getDbxref().getAccession()+"' "+
             "AND db_id=(SELECT db_id FROM db WHERE name='"+dbxref.getDbxref().getDb().getName()+"'))"+
      "AND feature_id=(SELECT feature_id FROM "+
             "feature WHERE  uniquename='"+uniquename+"')";
    
    Statement st = conn.createStatement();
    return st.executeUpdate(sql);
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
    final ChadoFeatureSynonym alias  = tsn.getAlias();
    final String uniquename   = alias.getUniquename();
    final String synonym_name = alias.getSynonym().getName();
      
    String sql;
     
    String sqlAliasId = "SELECT synonym_id FROM "+
                        "synonym WHERE synonym.name='"+synonym_name+"'";

    appendToLogFile(sqlAliasId, sqlLog);
    Statement st   = conn.createStatement(ResultSet.TYPE_SCROLL_SENSITIVE,
                                          ResultSet.CONCUR_UPDATABLE);
    ResultSet rs   = st.executeQuery(sqlAliasId);
    boolean exists = rs.next();
    
    if(!exists)
    {
      // create a new synonym name     
      String type_id = Long.toString(alias.getSynonym().getCvterm().getCvtermId());
      
      sql = "INSERT INTO "+
            "synonym (name, type_id, synonym_sgml) values ( '"+
            synonym_name+"',"+type_id+",'"+synonym_name+"')" ;

      st.executeUpdate(sql);
      appendToLogFile(sql, sqlLog);
      
      rs = st.executeQuery(sqlAliasId);
      rs.next();
      appendToLogFile(sqlAliasId, sqlLog);
    }
    
    final int synonym_id = rs.getInt("synonym_id"); 
    sql = "INSERT INTO "+
           "feature_synonym ( synonym_id, feature_id, pub_id )"+
           " values ( "+synonym_id+" ,"+
               "(SELECT feature_id FROM "+
               "feature WHERE  uniquename='"+uniquename+"'), "+
               " 1)";
 
    appendToLogFile(sql, sqlLog);
    return st.executeUpdate(sql);
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
    final String uniquename   = alias.getUniquename();
    final String synonym_name = alias.getSynonym().getName();
    String sql = "SELECT synonym_id FROM synonym WHERE "+
                 "synonym.name='"+synonym_name+"'";
    
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
      sql = "DELETE FROM feature_synonym WHERE "+
            "synonym_id="+synonym_id+" AND "+
            "feature_id=(SELECT feature_id FROM "+
                     "feature WHERE  uniquename='"+uniquename+"')";
    }
    else
      sql = "DELETE FROM synonym WHERE synonym_id="+synonym_id;
    
    st   = conn.createStatement();
    return st.executeUpdate(sql);
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
    final ChadoFeature chado_feature = tsn.getChadoFeature();
    final String parent = tsn.getParent_uniquename();
    
      
    StringBuffer sqlBuff = new StringBuffer();
    sqlBuff.append("UPDATE feature_relationship ");
    sqlBuff.append(" SET ");

    List properties = tsn.getProperties();
    for(int i=0; i<properties.size(); i++)
    {
      sqlBuff.append((String)properties.get(i));
      if(i < properties.size()-1)
        sqlBuff.append(" , ");
    }

    sqlBuff.append(" WHERE feature_relationship.subject_id="+
        "(SELECT feature_id FROM feature WHERE uniquename='"+
        chado_feature.getUniquename()+"') AND "+
        "feature_relationship.object_id="+
        "(SELECT feature_id FROM feature WHERE uniquename='"+
        parent+"')"); 
 
    String sql = sqlBuff.toString();

    System.out.println(sql);
    appendToLogFile(sql, sqlLog);
    Statement st = conn.createStatement();
    st.executeUpdate(sql);
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
    String sql = "UPDATE "+
                 "feature SET timelastmodified=";
    
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
   * @param uniquename  the unique name of the feature
   * @param timestamp   the time stamp to use, 
   *                    if NULL use CURRENT_TIMESTAMP
   * @throws SQLException
   */
  public int writeTimeAccessioned
                    (final String uniquename,
                     final Timestamp timestamp)
                     throws SQLException
  {   
    String sql = "UPDATE "+
                 "feature SET timeaccessioned=";

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
