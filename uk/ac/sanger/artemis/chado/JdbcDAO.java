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
   */
  public Feature getFeatureById(int id)
  {
    Feature feature = new Feature();
    feature.setId(id);
    return getLazyFeature(feature);
  }

  /**
   * Return a features with this systematic id
   *  
   * @param name the systematic id
   * @return the Feature, or null
   */
  public Feature getFeatureByUniqueName(String uniquename)
  {
    Feature feature = new Feature();
    feature.setUniquename(uniquename);
    feature.setId(-1);
    return getLazyFeature(feature);
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
    return getFeatureQuery(null, 
                 feature.getFeatureloc().getSrcfeature_id(), -1);
  }

  /**
   * Return a list of features with any current (ie non-obsolete) name or synonym
   *  
   * @param name the lookup name
   * @return a (possibly empty) List<Feature> of children with this current name 
   */
  public List getFeaturesByAnyCurrentName(String name) 
  {
    Feature feature = new Feature();
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
   * @return  the <code>List</code> of <code>Feature</code>
   */
  private Feature getLazyFeature(final Feature feature)
  {
    List list = getFeatureQuery(feature.getUniquename(), 
                                -1, feature.getId());
    return (Feature)list.get(0);
  }
  
  /**
   * Get the properties of a feature.
   * @param uniquename  the unique name of the feature
   * @param parentFeatureID  the id of parent feature to query
   * @return  the <code>List</code> of <code>Feature</code>
   */
  private List getFeatureQuery(final String uniquename,
                               final int parentFeatureID, 
                               final int feature_id)
  {
    final List list = new Vector();
    try
    {
      Statement st = conn.createStatement(ResultSet.TYPE_SCROLL_INSENSITIVE,
          ResultSet.CONCUR_UPDATABLE);

      String sql = "SELECT timelastmodified, f.feature_id, residues,"
          + " fl.strand, fmin, fmax, uniquename, f.type_id,"
          + " fp.type_id AS prop_type_id, fp.value, fl.phase,"
          + " f.organism_id, abbreviation, genus, species, common_name, comment,"
          + " fr.object_id, fr.type_id AS relation_type_id, fr.rank"
          + " FROM  feature f" + " LEFT JOIN feature_relationship fr ON "
          + "fr.subject_id=" + "f.feature_id" + " LEFT JOIN featureprop fp ON "
          + "fp.feature_id=" + "f.feature_id" + " LEFT JOIN featureloc fl ON "
          + "f.feature_id=" + "fl.feature_id"
          + " LEFT JOIN organism ON organism.organism_id=f.organism_id"
          + " WHERE ";

      if(uniquename != null)
        sql = sql + "uniquename LIKE '" + uniquename + "'";

      if(parentFeatureID > -1)
        sql = sql + "srcfeature_id = " + parentFeatureID;

      if(feature_id > -1)
        sql = sql + "f.feature_id = " + feature_id;

      sql = sql // + " AND (fl.rank=fr.rank OR fr.rank IS NULL)"
          + " ORDER BY f.type_id, uniquename";

      appendToLogFile(sql, sqlLog);
      ResultSet rs = st.executeQuery(sql);

      while(rs.next())
      {
        Feature feature = new Feature();

        FeatureLoc featureloc = new FeatureLoc();
        featureloc.setFmin(rs.getInt("fmin"));
        featureloc.setFmax(rs.getInt("fmax"));
        featureloc.setStrand(rs.getInt("strand"));

        int phase = rs.getInt("phase");
        if(rs.wasNull())
          featureloc.setPhase(10);
        else
          featureloc.setPhase(phase);

        feature.setResidues(rs.getBytes("residues"));

        feature.setFeatureloc(featureloc);
        feature.setCvterm(new Cvterm());
        feature.getCvterm().setCvtermId(rs.getLong("type_id"));

        // feature properties
        FeatureProp featureprop = new FeatureProp();
        Cvterm cvterm = new Cvterm();
        cvterm.setCvtermId(rs.getLong("prop_type_id"));
        featureprop.setCvterm(cvterm);
        featureprop.setValue(rs.getString("value"));
        feature.setFeatureprop(featureprop);

        feature.setUniquename(rs.getString("uniquename"));
        feature.setTimelastmodified(rs.getTimestamp("timelastmodified"));
        feature.setId(rs.getInt("feature_id"));

        // feature relationship
        FeatureRelationship feature_relationship = new FeatureRelationship();
        cvterm = new Cvterm();
        cvterm.setCvtermId(rs.getLong("relation_type_id"));
        feature_relationship.setCvterm(cvterm);

        Feature object = new Feature();
        object.setId(rs.getInt("object_id"));
        feature_relationship.setObject(object);
        feature.setFeature_relationship(feature_relationship);

        // feature organism
        Organism organism = new Organism();
        organism.setAbbreviation(rs.getString("abbreviation"));
        organism.setComment(rs.getString("comment"));
        organism.setCommon_name(rs.getString("common_name"));
        organism.setGenus(rs.getString("genus"));
        organism.setId(rs.getInt("organism_id"));
        organism.setSpecies(rs.getString("species"));
        feature.setOrganism(organism);

        list.add(feature);
      }
    }
    catch(SQLException sqle)
    {
      throw new RuntimeException(sqle);
    }
    // merge same features in the list
    return IBatisDAO.mergeList(list);
  }
  

  /**
   * Given a list of distict cvterm_id/type_id's of feature types that have
   * residues (from getResidueType()) in the given schema and the schema name
   * return a list of chado features in the schema with residues.
   * 
   * @param cvterm_ids
   *          list of cvterm_id/type_id's
   * @param schema
   *          schema/organism name or null
   * @return the <code>List</code> of <code>Feature</code> objects
   */
  public List getResidueFeatures(List cvterm_ids, 
                                 final String schema)
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

    List list = new Vector();
    try
    {
      Statement st = conn.createStatement();
      ResultSet rs = st.executeQuery(sql);
      while(rs.next())
      {
        Feature feature = new Feature();
        Organism organism = new Organism();
        organism.setAbbreviation(rs.getString("abbreviation"));

        feature.setOrganism(organism);
        feature.setId(rs.getInt("feature_id"));
        feature.setName(rs.getString("name"));
        feature.setUniquename(rs.getString("uniquename"));
        feature.setCvterm(new Cvterm());
        feature.getCvterm().setCvtermId(rs.getLong("type_id"));

        list.add(feature);
      }
    }
    catch(SQLException sqle)
    {
      throw new RuntimeException(sqle);
    }
    return list;
  }
 

  /**
   * For a schema return the type_id's with residues.
   * @param schema      schema/organism name or null
   * @return	the <code>List</code> of type_id's as <code>String</code>
   *            objects
   */
  public List getResidueType(final String schema)
  {
    String sql = "SELECT DISTINCT type_id FROM " +schema+
                 ".feature WHERE residues notnull";
    appendToLogFile(sql, sqlLog);

    List cvterm_ids = new Vector();
    try
    {
      Statement st = conn.createStatement();
      ResultSet rs = st.executeQuery(sql);

      while(rs.next())
        cvterm_ids.add(rs.getString("type_id"));
    }
    catch(SQLException sqle)
    {
      throw new RuntimeException(sqle);
    }
    
    return cvterm_ids;
  }

  /**
   * Get available schemas (as a <code>List</code> of <code>String</code>       
   * objects).
   * @return	the available schemas
   */
  public List getSchema()
  {
    List schemas = new Vector();
    try
    {
      Statement st = conn.createStatement();

      String query = "SELECT schema_name FROM information_schema.schemata "
          + "WHERE schema_name=schema_owner ORDER BY schema_name";
      appendToLogFile(query, sqlLog);

      ResultSet rs = st.executeQuery(query);

      while(rs.next())
        schemas.add(rs.getString("schema_name"));
    }
    catch(SQLException sqle)
    {
      throw new RuntimeException(sqle);
    }
    
    return schemas;
  } 

  /**
   * Get the full list of cvterm_id and name as a <code>List</code> of 
   * <code>Cvterm</code> objects.
   * @return    the full list of cvterm_id and name
   */
  public List getCvterm()
  {
    String sql = "SELECT cvterm.cvterm_id, cvterm.name " +
                 "FROM cvterm, cv WHERE cv.cv_id = cvterm.cv_id";

    appendToLogFile(sql, sqlLog);
    
    try
    {
      Statement st = conn.createStatement();
      ResultSet rs = st.executeQuery(sql);
      List cvterms = new Vector();

      while(rs.next())
      {
        Cvterm cvterm = new Cvterm();
        cvterm.setCvtermId(rs.getLong("cvterm_id"));
        cvterm.setName(rs.getString("name"));
        cvterms.add(cvterm);
      }
      return cvterms;
    }
    catch(SQLException sqle)
    {
      throw new RuntimeException(sqle);
    }

  }

  /**
   * Get dbxref for a feature.
   * @param uniquename  the unique name for the feature. If set to NULL
   *                    all <code>FeatureDbxref</code> are returned.
   * @return a <code>List</code> of feature_dbxrefs.
   */
  public List getFeatureDbxrefByUniquename(final String uniquename)
  {
    String sql = "SELECT db.name, dbx.accession, f.feature_id FROM "+
                 "feature_dbxref dbx_f "+
                 "LEFT JOIN dbxref dbx ON dbx.dbxref_id=dbx_f.dbxref_id "+
                 "LEFT JOIN db ON db.db_id=dbx.db_id "+
                 "LEFT JOIN feature f ON dbx_f.feature_id=f.feature_id ";
    
    if(uniquename != null)
      sql = sql + "WHERE f.uniquename='"+uniquename+"'";
    
    appendToLogFile(sql, sqlLog);
    
    try
    {
      Statement st = conn.createStatement();
      ResultSet rs = st.executeQuery(sql);
      List dbxrefs = new Vector();

      while(rs.next())
      {
        FeatureDbxref feature_dbxref = new FeatureDbxref();
        Dbxref dbxref = new Dbxref();
        Db db = new Db();
        db.setName(rs.getString("name"));
        dbxref.setAccession(rs.getString("accession"));
        dbxref.setDb(db);
        Feature feat = new Feature();
        feat.setId(rs.getInt("feature_id"));
        feature_dbxref.setDbxref(dbxref);
        feature_dbxref.setFeature(feat);
        dbxrefs.add(feature_dbxref);
      }

      return dbxrefs;
    }
    catch(SQLException sqle)
    {
      throw new RuntimeException(sqle);
    }
  }
  
  
  /**
   * Return a list of FeatureSynonyms for a uniquename
   * @param uniquename  the unique name for the feature. If set to NULL
   *                    all <code>FeatureSynonym</code> are returned.
   * @return
   */
  public List getFeatureSynonymsByUniquename(final String uniquename)
  {
    String sql = "SELECT fs.*, s.name, s.type_id FROM "+
    "feature_synonym fs "+
    "LEFT JOIN feature f ON f.feature_id=fs.feature_id "+
    "LEFT JOIN synonym s ON fs.synonym_id=s.synonym_id ";
  
    if(uniquename != null)
      sql = sql + " WHERE uniquename='" + uniquename + "'";

    appendToLogFile(sql, sqlLog);
    
    try
    {
      Statement st = conn.createStatement();
      ResultSet rs = st.executeQuery(sql);
      List synonym = new Vector();
      FeatureSynonym alias;
      while(rs.next())
      {
        alias = new FeatureSynonym();
        Cvterm cvterm = new Cvterm();
        cvterm.setCvtermId(rs.getLong("type_id"));
        Synonym syn = new Synonym();
        syn.setName(rs.getString("name"));
        syn.setCvterm(cvterm);
        Feature feat = new Feature();
        feat.setId(rs.getInt("feature_id"));

        alias.setSynonym(syn);
        alias.setFeature(feat);
        alias.setPub_id(new Integer(rs.getInt("pub_id")));
        alias.setInternal(rs.getBoolean("is_internal"));
        alias.setCurrent(rs.getBoolean("is_current"));
        synonym.add(alias);
      }

      return synonym;
    }
    catch(SQLException sqle)
    {
      throw new RuntimeException(sqle);
    }
  }
  
  public List getFeatureSynonymsByFeatureAndSynonym(
      Feature feature, Synonym synonym)
  {
    return null;
  }


  public Synonym getSynonymByNameAndCvTerm(
      String name, Cvterm type)
  {
    return null;
  }
  
//
// WRITE 
//
  
  
  /**
   * Merge (update) an already persistent object back to the database 
   * (at the end of the current transaction, or depending upon flush mode). 
   * This method is defined in all the DAOs. It's recommended to call it 
   * through an appropriate one eg SequenceDaoI for FeatureI 
   * @param o The object to merge
   */
  public void merge(Object o) 
  {
    if(o instanceof FeatureLoc)
      updateFeatureLoc((FeatureLoc)o);
    else if(o instanceof Feature)
      updateFeature((Feature)o);
    else if(o instanceof FeatureProp)
      updateFeatureProp((FeatureProp)o);
    else if(o instanceof FeatureRelationship)
      updateFeatureRelationship((FeatureRelationship)o);
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
      insertFeatureProp((FeatureProp)o);
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
  public void delete(Object o)
  {
    if(o instanceof Feature)
      deleteFeature((Feature)o);
    else if(o instanceof FeatureProp)
      deleteFeatureProp((FeatureProp)o);
    else if(o instanceof FeatureDbxref)
      deleteFeatureDbxref((FeatureDbxref)o);
    else if(o instanceof FeatureSynonym)
      deleteFeatureSynonym((FeatureSynonym)o);
  }
  
  
  /**
   * Update a feature location with the give <code>FeatureLoc</code>
   * object.
   * @param featureloc  the new <code>FeatureLoc</code> object.
   */
  private void updateFeatureLoc(FeatureLoc featureloc) 
  {
    final String sql = "UPDATE featureloc SET fmin=?, fmax=?, rank=?, strand=?, phase=? "+
                       "WHERE feature_id=(SELECT feature_id FROM feature WHERE uniquename=?)";
    try
    {
    PreparedStatement pstmt = conn.prepareStatement(sql);
      pstmt.setInt(1, featureloc.getFmin());
      pstmt.setInt(2, featureloc.getFmax());
      pstmt.setInt(3, featureloc.getRank());
      pstmt.setInt(4, featureloc.getStrand());
      pstmt.setInt(5, featureloc.getPhase());
      pstmt.setString(6, featureloc.getFeature().getUniquename());
      appendToLogFile(sql, sqlLog);

      pstmt.executeUpdate();
    }
    catch(SQLException sqle)
    {
      throw new RuntimeException(sqle);
    }
  }
  
  /**
   * Update a feature with a given <code>Feature</code> object.
   * @param feature the new <code>Feature</code> object.
   */
  private void updateFeature(Feature feature)
  {
    String sql = "UPDATE feature SET uniquename=?";

    if(feature.getCvterm() != null)
      sql = sql+", type_id=?";
    
    if(feature.getTimelastmodified() != null)
      sql = sql+", timelastmodified=?";

    sql = sql+"WHERE feature_id=?";
    
    try
    {
      PreparedStatement pstmt = conn.prepareStatement(sql);
      pstmt.setString(1, feature.getUniquename());
      int param = 2;
      if(feature.getCvterm() != null)
      {
        pstmt.setLong(param, feature.getCvterm().getCvtermId());
        param++;
      }

      if(feature.getTimelastmodified() != null)
      {
        pstmt.setTimestamp(param, feature.getTimelastmodified());
        param++;
      }

      pstmt.setInt(param, feature.getId());
      appendToLogFile(sql, sqlLog);
      pstmt.executeUpdate();
    }
    catch(SQLException sqle)
    {
      throw new RuntimeException(sqle);
    }
  }
  
  /**
   * Update a feature property with a given <code>FeatureProp</code>
   * object.
   * @param featureprop the new <code>FeatureProp</code> object.
   */
  private void updateFeatureProp(FeatureProp featureprop)
  { 
    String sql = "UPDATE featureprop SET value=? "+
                 "WHERE rank=? AND type_id=? AND "+
                 "feature_id=(SELECT feature_id FROM feature WHERE uniquename=?)";
    try
    {
      PreparedStatement pstmt = conn.prepareStatement(sql);
      pstmt.setString(1, featureprop.getValue());
      pstmt.setInt(2, featureprop.getRank());
      pstmt.setLong(3, featureprop.getCvterm().getCvtermId());
      pstmt.setString(4, featureprop.getFeature().getUniquename());

      appendToLogFile(sql, sqlLog);
      pstmt.executeUpdate();
    }
    catch(SQLException sqle)
    {
      throw new RuntimeException(sqle);
    }
  }
  

  /**
   * Insert attributes defined by the <code>FeatureProp</code>.
   * @param featureprop     the new <code>FeatureProp</code>
   */
  private void insertFeatureProp
                    (final FeatureProp featureprop)
  {
    StringBuffer sqlBuff = new StringBuffer();

    sqlBuff.append("INSERT INTO featureprop");
    sqlBuff.append(" ( feature_id, type_id, value, rank ) ");
    sqlBuff.append("VALUES ");
    
    sqlBuff.append("((SELECT feature_id FROM feature WHERE uniquename=");
    sqlBuff.append("'"+ featureprop.getFeature().getUniquename()+"')," );
    sqlBuff.append(featureprop.getCvterm().getCvtermId()+", ");
    sqlBuff.append(featureprop.getValue()+",");
    sqlBuff.append(featureprop.getRank());

    appendToLogFile(new String(sqlBuff), sqlLog);

    try
    {
      Statement st = conn.createStatement();
      int rowCount = st.executeUpdate(new String(sqlBuff));
      System.out.println(rowCount + " row(s) inserted");
    }
    catch(SQLException sqle)
    {
      throw new RuntimeException(sqle);
    }
  }

  /**
   * Delete attributes defined by the <code>FeatureProp</code>.
   * @param featureprop      the new <code>FeatureProp</code>
   */
  private void deleteFeatureProp
                    (final FeatureProp featureprop)
  {
    StringBuffer sqlBuff = new StringBuffer();
    String uniquename = featureprop.getFeature().getUniquename();

    sqlBuff.append("DELETE FROM featureprop WHERE ");

    if(uniquename != null)
      sqlBuff.append("feature_id="+
          "(SELECT feature_id FROM feature WHERE uniquename='"+
          uniquename+"') AND ");
    
    if(featureprop.getRank() > -1)
      sqlBuff.append("rank="+featureprop.getRank()+" AND ");
    
    if(featureprop.getValue() != null)
      sqlBuff.append("value="+featureprop.getValue()+" AND ");
    
    sqlBuff.append("type_id="+featureprop.getCvterm().getCvtermId());
      
    appendToLogFile(new String(sqlBuff), sqlLog);

    try
    {
      Statement st = conn.createStatement();
      int rowCount = st.executeUpdate(new String(sqlBuff));
      System.out.println(rowCount+" row(s) deleted");
    }
    catch(SQLException sqle)
    {
      throw new RuntimeException(sqle);
    }
  }


  /**
   * Insert a feature into the database defined by the <code>Feature</code>.
   * @param feature    the new <code>Feature</code>
   */
  private void insertFeature
                    (final Feature feature)
  {
    try
    {
      //
      // get the organism_id
      Statement st = conn.createStatement();
      String sql = "SELECT organism_id from " + "feature where feature_id = '"
          + feature.getFeatureloc().getSrcfeature_id() + "'";

      appendToLogFile(sql, sqlLog);
      ResultSet rs = st.executeQuery(sql);
      rs.next();

      final int organism_id = rs.getInt("organism_id");

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
      sql_buff.append(organism_id + " , ");
      sql_buff.append("'" + feature.getName() + "'" + " , ");
      sql_buff.append("'" + feature.getUniquename() + "'" + " , ");
      sql_buff.append(Long.toString(feature.getCvterm().getCvtermId()));
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
      sql_buff.append(feature_id + " , ");
      sql_buff.append(feature.getFeatureloc().getSrcfeature_id() + " , ");
      sql_buff.append(feature.getFeatureloc().getFmin() + " , ");
      sql_buff.append(feature.getFeatureloc().getFmax() + " , ");
      sql_buff.append(feature.getFeatureloc().getStrand() + " , ");
      sql_buff.append(feature.getFeatureloc().getPhase());
      sql_buff.append(" )");

      sql = new String(sql_buff);
      appendToLogFile(sql, sqlLog);
      st = conn.createStatement();
      rowCount = st.executeUpdate(sql);

      // insert feature relationships
      if(feature.getFeatureRelationshipsForSubjectId() != null)
      {
        List parents = feature.getFeatureRelationshipsForSubjectId();
        for(int i = 0; i < parents.size(); i++)
        {
          // insert feature_relationship
          FeatureRelationship feature_relationship = (FeatureRelationship) parents
              .get(i);

          sql_buff = new StringBuffer();
          sql_buff.append("INSERT INTO feature_relationship ");
          sql_buff.append("( subject_id, object_id, type_id ) ");
          sql_buff.append("VALUES ");
          sql_buff
              .append("( (SELECT feature_id FROM feature WHERE uniquename='");
          sql_buff.append(feature_relationship.getSubject().getUniquename()
              + "'), ");
          sql_buff.append("(SELECT feature_id FROM feature WHERE uniquename='");
          sql_buff.append(feature_relationship.getObject().getUniquename()
              + "'), ");
          sql_buff.append(feature_relationship.getCvterm().getCvtermId() + ")");

          sql = new String(sql_buff);
          appendToLogFile(sql, sqlLog);
          st = conn.createStatement();
          rowCount = st.executeUpdate(sql);

        }
      }
    }
    catch(SQLException sqle)
    {
      throw new RuntimeException(sqle);
    }

  }

  /**
   * Delete a feature from the database defined by the 
   * <code>Feature</code>.
   * @param feature  the new <code>Feature</code>
   */
  private int deleteFeature
                    (final Feature feature)
  {
    try
    {
      String sql = "DELETE FROM feature WHERE uniquename='"
          + feature.getUniquename() + "'";
      appendToLogFile(sql, sqlLog);

      Statement st = conn.createStatement();
      return st.executeUpdate(sql);
    }
    catch(SQLException sqle)
    {
      throw new RuntimeException(sqle);
    }
  }

  /**
   * Insert a feature_dbxref for a feature.
   * @param feature_dbxref    the new <code>FeatureDbxref</code>
   */
  private void insertFeatureDbxref(final FeatureDbxref feature_dbxref)
  {   
    // find database id
    String sql = "SELECT db_id FROM db WHERE name='"+
                 feature_dbxref.getDbxref().getDb().getName()+"'";
    
    try
    {
      Statement st = conn.createStatement();
      ResultSet rs = st.executeQuery(sql);
      boolean exists = rs.next();

      if(!exists)
        throw new SQLException("No database called "
            + feature_dbxref.getDbxref().getDb().getName() + " found (for "
            + feature_dbxref.getFeature().getUniquename()
            + ") check the spelling!");

      final int db_id = rs.getInt("db_id");
      // find if accession exists already
      String sqlDbxrefId = "SELECT dbxref_id FROM dbxref WHERE accession='"
          + feature_dbxref.getDbxref().getAccession() + "' AND db_id=" + db_id;

      appendToLogFile(sqlDbxrefId, sqlLog);
      rs = st.executeQuery(sqlDbxrefId);
      exists = rs.next();

      if(!exists)
      {
        // create a new accession entry in dbxref
        sql = "INSERT INTO dbxref ( db_id, accession ) " + "VALUES (" + db_id
            + ", '" + feature_dbxref.getDbxref().getAccession() + "' )";

        appendToLogFile(sql, sqlLog);
        int rowCount = st.executeUpdate(new String(sql));

        // now get the new dbxref_id
        appendToLogFile(sqlDbxrefId, sqlLog);
        rs = st.executeQuery(sqlDbxrefId);
        rs.next();
      }

      final int dbxref_id = rs.getInt("dbxref_id");
      sql = "INSERT INTO feature_dbxref "
          + "(feature_id, dbxref_id, is_current)" + " VALUES "
          + "( (SELECT feature_id FROM " + "feature WHERE  uniquename='"
          + feature_dbxref.getFeature().getUniquename() + "'), " + dbxref_id
          + ", " + Boolean.toString(feature_dbxref.isCurrent()) + ")";
      System.out.println(sql);
      appendToLogFile(sql, sqlLog);
      st.executeUpdate(new String(sql));
    }
    catch(SQLException sqle)
    {
      throw new RuntimeException(sqle);
    }
  }
  
  /**
   * Delete a feature_dbxref for a feature.
   * @param feature_dbxref  the  new <code>FeatureDbxref</code>
   */
  private void deleteFeatureDbxref(final FeatureDbxref feature_dbxref)
  {
    final String uniquename = feature_dbxref.getFeature().getUniquename();
    
    final String sql = 
      "DELETE FROM feature_dbxref "+
      "WHERE dbxref_id="+
      "(SELECT dbxref_id FROM dbxref WHERE accession='"+
         feature_dbxref.getDbxref().getAccession()+"' "+
      "AND db_id=(SELECT db_id FROM db WHERE name='"+
         feature_dbxref.getDbxref().getDb().getName()+"'))"+
      "AND feature_id=(SELECT feature_id FROM "+
             "feature WHERE  uniquename='"+uniquename+"')";
    
    try
    {
      Statement st = conn.createStatement();
      st.executeUpdate(sql);
    }
    catch(SQLException sqle)
    {
      throw new RuntimeException(sqle);
    }
  }
  
  /**
   * Insert a feature_synonym for a feature.
   * @param feature_synonym  the new <code>FeatureSynonym</code>
   */
  private void insertFeatureAlias(final FeatureSynonym feature_synonym)
  {
    final String uniquename   = feature_synonym.getFeature().getUniquename();
    final String synonym_name = feature_synonym.getSynonym().getName();
      
    String sql;
     
    String sqlAliasId = "SELECT synonym_id FROM "+
                        "synonym WHERE synonym.name='"+synonym_name+"'";

    appendToLogFile(sqlAliasId, sqlLog);
    
    try
    {
      Statement st = conn.createStatement(ResultSet.TYPE_SCROLL_SENSITIVE,
          ResultSet.CONCUR_UPDATABLE);
      ResultSet rs = st.executeQuery(sqlAliasId);
      boolean exists = rs.next();

      if(!exists)
      {
        // create a new synonym name
        String type_id = Long.toString(feature_synonym.getSynonym().getCvterm()
            .getCvtermId());

        sql = "INSERT INTO "
            + "synonym (name, type_id, synonym_sgml) values ( '" + synonym_name
            + "'," + type_id + ",'" + synonym_name + "')";

        st.executeUpdate(sql);
        appendToLogFile(sql, sqlLog);

        rs = st.executeQuery(sqlAliasId);
        rs.next();
        appendToLogFile(sqlAliasId, sqlLog);
      }

      final int synonym_id = rs.getInt("synonym_id");
      sql = "INSERT INTO "
          + "feature_synonym ( synonym_id, feature_id, pub_id )" + " values ( "
          + synonym_id + " ," + "(SELECT feature_id FROM "
          + "feature WHERE  uniquename='" + uniquename + "'), " + " 1)";

      appendToLogFile(sql, sqlLog);
      st.executeUpdate(sql);
    }
    catch(SQLException sqle)
    {
      throw new RuntimeException(sqle);
    }
  }
  
  /**
   * Delete a feature_synonym for a feature.
   * @param feature_synonym  the new <code>FeatureSynonym</code>
   * @return    number of rows changed
   */
  private void deleteFeatureSynonym(final FeatureSynonym feature_synonym)
  {
    final String uniquename   = feature_synonym.getFeature().getUniquename();
    final String synonym_name = feature_synonym.getSynonym().getName();
    String sql = "SELECT synonym_id FROM synonym WHERE "+
                 "synonym.name='"+synonym_name+"'";
    
    appendToLogFile(sql, sqlLog);
    
    try
    {
      Statement st = conn.createStatement(ResultSet.TYPE_SCROLL_INSENSITIVE,
          ResultSet.CONCUR_READ_ONLY);
      ResultSet rs = st.executeQuery(sql);
      rs.last();
      int nrows = rs.getRow();
      final int synonym_id = rs.getInt("synonym_id");

      // check this name is not used some where else,
      // i.e. in more than one row
      if(nrows > 1)
      {
        sql = "DELETE FROM feature_synonym WHERE " + "synonym_id=" + synonym_id
            + " AND " + "feature_id=(SELECT feature_id FROM "
            + "feature WHERE  uniquename='" + uniquename + "')";
      }
      else
        sql = "DELETE FROM synonym WHERE synonym_id=" + synonym_id;

      st = conn.createStatement();
      st.executeUpdate(sql);
    }
    catch(SQLException sqle)
    {
      throw new RuntimeException(sqle);
    }
  }
  
  /**
   * Update feature_relationship for a feature.
   * @param feature_relationship  the <code>FeatureRelationship</code>
   */
  private void updateFeatureRelationship(
      final FeatureRelationship feature_relationship)
  {  
    StringBuffer sqlBuff = new StringBuffer();
    sqlBuff.append("UPDATE feature_relationship ");
    sqlBuff.append(" SET ");
    sqlBuff.append(" rank="+feature_relationship.getRank()+", ");
    sqlBuff.append(" type_id="+feature_relationship.getCvterm().getCvtermId());
    sqlBuff.append(" WHERE ");
    sqlBuff.append("subject_id=");
    sqlBuff.append("( SELECT feature_id FROM feature WHERE uniquename='");
    sqlBuff.append(feature_relationship.getSubject().getUniquename()+"' ) ");
    sqlBuff.append("AND ");
    sqlBuff.append("object_id=");
    sqlBuff.append("( SELECT feature_id FROM feature WHERE uniquename='");
    sqlBuff.append(feature_relationship.getObject().getUniquename()+"' ) ");
 
    String sql = sqlBuff.toString();

    System.out.println(sql);
    appendToLogFile(sql, sqlLog);
    try
    {
      Statement st = conn.createStatement();
      st.executeUpdate(sql);
    }
    catch(SQLException sqle)
    {
      throw new RuntimeException(sqle);
    }
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
        catch(IOException ioe2){}
    }
  }

}
