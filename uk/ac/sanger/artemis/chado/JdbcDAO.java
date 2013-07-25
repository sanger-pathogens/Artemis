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
import java.util.Collection;
import java.util.List;
import java.util.Vector;

import uk.ac.sanger.artemis.util.DatabaseLocationParser;

import org.gmod.schema.sequence.Feature;
import org.gmod.schema.sequence.FeatureCvTerm;
import org.gmod.schema.sequence.FeatureCvTermProp;
import org.gmod.schema.sequence.FeatureDbXRef;
import org.gmod.schema.sequence.FeatureProp;
import org.gmod.schema.sequence.Synonym;
import org.gmod.schema.sequence.FeatureLoc;
import org.gmod.schema.sequence.FeatureRelationship;
import org.gmod.schema.sequence.FeatureSynonym;
import org.gmod.schema.sequence.FeatureCvTermDbXRef;
import org.gmod.schema.sequence.FeatureCvTermPub;
import org.gmod.schema.general.Db;
import org.gmod.schema.general.DbXRef;
import org.gmod.schema.organism.Organism;
import org.gmod.schema.cv.CvTerm;
import org.gmod.schema.cv.Cv;
import org.gmod.schema.pub.Pub;
import org.gmod.schema.pub.PubDbXRef;



/**
 *
 * Java Database Connectivity (JDBC) implemetation of the
 * <code>ChadoDAO</code> data access interface.
 *
 */
public class JdbcDAO extends GmodDAO
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
    DatabaseLocationParser dlp = new DatabaseLocationParser(location);
    conn = DriverManager.getConnection(dlp.getConnectionString(),
                                       dlp.getUsername(),
                                       new String(pfield.getPassword()));
  }

  //////
  ////// GeneralDaoI
  //////
  //////
  public List getDbs()
  {
    return null;
  }
  
  //////
  ////// SequenceDaoI
  //////
  //////
  
  public List getOrganismsContainingSrcFeatures()
  {
    return null;
  }
  
  public Feature getLazyFeatureNoResiduesById(final Integer featureId)
  {
    return null;
  }
  
  public List getClustersByFeatureIds(List featureIds)
  {
    // TODO Auto-generated method stub
    return null;
  }
  
  public List getFeatureDbXRefsByFeatureId(List featureIds)
  {
    // TODO Auto-generated method stub
    return null;
  }
  
  public List getFeaturePropByFeatureIds(List featureIds)
  {
    // TODO Auto-generated method stub
    return null;
  }

  public List getFeaturesByListOfIds(List featureIds)
  {
    // TODO Auto-generated method stub
    return null;
  }
  
  public List getFeatureDbXRefsBySrcFeature(Feature srcFeature)
  {
    // TODO Auto-generated method stub
    return null;
  }

  public List getFeatureSynonymsBySrcFeature(Feature srcFeature)
  {
    // TODO Auto-generated method stub
    return null;
  }

  public List getFeatureSynonymsByFeatureIds(final List featuresIds)
  {
    return null;
  }
  
  public List getFeatureLocsByFeatureId(int featureId)
  {
    return null;
  }
  
  public List getFeatureLocsBySrcFeatureId(int srcFeatureId)
  {
    return null; 
  }
  
  public List getFeatureLocsByListOfIds(final Collection featureIds)
  {
    return null;
  }
  
  public List getResiduesByUniqueName(String uniqueName)
  {
    // TODO Auto-generated method stub
    return null;
  }
  
  public List getFeaturePubsBySrcFeature(Feature srcFeature)
  {
    // TODO Auto-generated method stub
    return null;
  }
  
  public List getFeaturePubsByFeature(final Feature feature)
  {
    return null;
  }
  
  public List getSimilarityMatches(final Integer srcFeatureId)
  {
    return null;
  }
  
  public List getSimilarityMatchesByFeatureIds(List featureIds)
  {
    return null;
  }
  
  public List getFeatureCvTermDbXRefBySrcFeature(Feature srcFeature)
  {
    // TODO Auto-generated method stub
    return null;
  }

  public List getFeatureCvTermPubBySrcFeature(Feature srcFeature)
  {
    // TODO Auto-generated method stub
    return null;
  }

  public List getFeatureCvTermsBySrcFeature(Feature srcFeature)
  {
    // TODO Auto-generated method stub
    return null;
  }

  public List getResidueFeatures(Integer organismId)
  {
    // TODO Auto-generated method stub
    return null;
  }
  
  public List getResidueFeaturesByOrganismCommonName(final String commonName)
  {
    return null;
  }
  
  public List getParentFeaturesByChildFeatureIds(final List featureIds)
  {
    return null;
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
    feature.setFeatureId(id);
    return getLazyFeature(feature);
  }

  /**
   * Return a features with this systematic id
   *  
   * @param name the systematic id
   * @return the Feature, or null
   */
  public Feature getFeatureByUniqueName(final String uniquename, final String featureType)
  {
    Feature feature = new Feature();
    feature.setUniqueName(uniquename);
    feature.setFeatureId(-1);
    
    CvTerm cvTerm = new CvTerm();
    cvTerm.setName(featureType);
    feature.setCvTerm(cvTerm);
    
    return getLazyFeature(feature);
  }
  
  /**
   * Return a features with this systematic id
   *  
   * @param name the systematic id
   * @return the Feature, or null
   */
  public List getFeaturesByUniqueName(String uniquename)
  {
    return getFeatureQuery(uniquename, -1, -1, null);
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
    feature.setUniqueName(name);
    
    // getFeatureSynonymsByName() needs implementing
    //List feature_synonym_list = getFeatureSynonymsByName();
    
    return getFeatureQuery(name, -1, -1, null);
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
      org.gmod.schema.sequence.Feature parent, String type)
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
    return getFeatureQuery(null, 
                 feature.getFeatureLoc().getFeatureBySrcFeatureId().getFeatureId(), -1, null);
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
    String sqlTest = "SELECT pg_attribute.attname "+
                     "FROM pg_attribute, pg_class, pg_namespace "+
                     "WHERE pg_namespace.oid=pg_class.relnamespace AND "+
                     "attrelid=pg_class.oid AND "+
                     "relname='feature_cvterm' AND "+
                     "attnum > 0 AND "+
                     "nspname='"+ArtemisUtils.getCurrentSchema()+"'";
    appendToLogFile(sqlTest, sqlLog);
    
    boolean fcRank = false;
    try
    {
      Statement st = conn.createStatement();
      ResultSet rs = st.executeQuery(sqlTest);
      while(rs.next())
      {
        if(rs.getString("attname").equals("rank"))
          fcRank = true;
      }
      
    }
    catch(SQLException sqle)
    {
      throw new RuntimeException(sqle);
    }
    
    String sql = "SELECT fc.*, "+
     "fcp.type_id, fcp.value, fcp.rank AS fcp_rank, "+
     "cvterm.name AS cvterm_name, cv.name AS cv_name, "+
     "pub.pub_id, pub.uniquename, "+
     "db.name, dbxref.accession "+
     "FROM feature_cvterm fc "+
     "LEFT JOIN feature_cvtermprop fcp ON fc.feature_cvterm_id=fcp.feature_cvterm_id "+
     "LEFT JOIN cvterm ON cvterm.cvterm_id=fc.cvterm_id "+
     "LEFT JOIN cv ON cvterm.cv_id=cv.cv_id "+
     "LEFT JOIN pub ON fc.pub_id=pub.pub_id "+
     "LEFT JOIN dbxref ON cvterm.dbxref_id=dbxref.dbxref_id "+
     "LEFT JOIN db ON dbxref.db_id=db.db_id ";

    
    if(feature != null && feature.getUniqueName() != null)
      sql = sql + " WHERE "+
        "feature_id=(SELECT feature_id FROM feature WHERE uniquename='"+
        feature.getUniqueName()+"')";
    
    if(fcRank)
      sql = sql + " ORDER BY fc.feature_cvterm_id, fc.rank, type_id, fcp_rank";
    else
      sql = sql + " ORDER BY fc.feature_cvterm_id, type_id, fcp_rank";
    
    appendToLogFile(sql, sqlLog);

    try
    {
      Statement st = conn.createStatement(ResultSet.TYPE_SCROLL_INSENSITIVE,
          ResultSet.CONCUR_UPDATABLE);
      ResultSet rs = st.executeQuery(sql);
      List featureCvTerms = new Vector();

      while(rs.next())
      {
        int feature_id = rs.getInt("feature_id");
        Feature this_feature = new Feature();
        this_feature.setFeatureId(feature_id);
  
        CvTerm cvterm = new CvTerm();
        cvterm.setCvTermId(rs.getInt("cvterm_id"));
        cvterm.setName(rs.getString("cvterm_name"));
        Cv cv = new Cv();
        cv.setName(rs.getString("cv_name"));
        cvterm.setCv(cv);
        
        DbXRef dbxref = new DbXRef();
        dbxref.setAccession(rs.getString("accession"));   
        Db db = new Db();
        db.setName(rs.getString("name"));
        dbxref.setDb(db);
        
        cvterm.setDbXRef(dbxref);      
        
        Pub pub = new Pub();
        pub.setPubId(rs.getInt("pub_id"));
        pub.setUniqueName(rs.getString("uniquename"));
        
        
        int fc_rank = 0;
        
        if(fcRank)
          fc_rank = rs.getInt("rank");
        
        FeatureCvTerm feature_cvterm =
           new FeatureCvTerm(cvterm, this_feature, pub, rs.getBoolean("is_not"), fc_rank);
        
        // feature_cvtermprop's group by feature_cvterm_id
        List featureCvTermProps = new Vector();
        int next_fc_rank = -1;      
        int next_feature_cvterm_id = -1;  
        int feature_cvterm_id = rs.getInt("feature_cvterm_id");
        feature_cvterm.setFeatureCvTermId(feature_cvterm_id);
        
        do
        {
          FeatureCvTermProp featureProp = new FeatureCvTermProp();
          CvTerm featurePropCvTerm = new CvTerm();
          featurePropCvTerm.setCvTermId(rs.getInt("type_id"));
          featureProp.setCvTerm(featurePropCvTerm);
          featureProp.setValue(rs.getString("value"));
          featureProp.setRank(rs.getInt("fcp_rank"));
          
          featureCvTermProps.add(featureProp);
          
          if(rs.next())
          {
            next_feature_cvterm_id = rs.getInt("feature_cvterm_id");
            next_fc_rank = 0; 
            
            if(fcRank)
              next_fc_rank = rs.getInt("rank");
            if(feature_cvterm_id != next_feature_cvterm_id ||
               fc_rank != next_fc_rank)
              rs.previous();
          }
          else
            next_feature_cvterm_id = -1;
        } while(feature_cvterm_id == next_feature_cvterm_id &&
                fc_rank == next_fc_rank);

        feature_cvterm.setFeatureCvTermProps(featureCvTermProps);

        featureCvTerms.add(feature_cvterm);
      }

      return featureCvTerms;
    }
    catch(SQLException sqle)
    {
      throw new RuntimeException(sqle);
    }
  }
  
  public Synonym getSynonymByNameAndCvTerm(
      String name, CvTerm cvTerm)
  {
    String sql = "SELECT * FROM synonym WHERE ";
    
    if(name != null)
     sql = sql + "name="+name+" AND "; 
    if(cvTerm != null)
      sql = sql + "type_id="+cvTerm.getCvTermId()+" AND ";
    
    sql = sql + "synonym_id > 0";
    
    appendToLogFile(sql, sqlLog);

    try
    {
      Statement st = conn.createStatement();
      ResultSet rs = st.executeQuery(sql);
      
      CvTerm cvterm = new CvTerm();
      cvterm.setCvTermId(rs.getInt("type_id"));
      
      Synonym synonym = new Synonym(cvterm, rs.getString("name"), 
                                    rs.getString("synonym_sgml"));
      synonym.setSynonymId(rs.getInt("synonym_id"));
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
  
  /**
   * Return a list of FeatureSynonyms which link a given Feature
   * and Synonym 
   * @param feature the test Feature
   * @param synonym the test Synonym
   * @return a (possibly empty) List<FeatureSynonym>
   */
  public List getFeatureSynonymsByFeatureUniquename(
      Feature feature, Synonym synonym)
  {
    return null;
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
    String sql = "SELECT db.name, dbx.accession, dbx.version, dbx.description, "
        + "dbx_f.feature_id, dbx_f.is_current FROM "
        + "feature_dbxref dbx_f "
        + "LEFT JOIN dbxref dbx ON dbx.dbxref_id=dbx_f.dbxref_id "
        + "LEFT JOIN db ON db.db_id=dbx.db_id "
        + "LEFT JOIN feature f ON dbx_f.feature_id=f.feature_id ";

    if(uniqueName != null)
      sql = sql + "WHERE f.uniquename='" + uniqueName + "'";

    sql = sql + " ORDER BY f.type_id,  uniquename";
    
    appendToLogFile(sql, sqlLog);

    try
    {
      Statement st = conn.createStatement();
      ResultSet rs = st.executeQuery(sql);
      List dbxrefs = new Vector();

      while(rs.next())
      {
        FeatureDbXRef feature_dbxref = new FeatureDbXRef();
        DbXRef dbxref = new DbXRef();
        Db db = new Db();
        db.setName(rs.getString("name"));
        dbxref.setAccession(rs.getString("accession"));
        dbxref.setVersion(rs.getString("version"));
        dbxref.setDescription(rs.getString("description"));
        dbxref.setDb(db);
        Feature feat = new Feature();
        feat.setFeatureId(rs.getInt("feature_id"));
        feature_dbxref.setDbXRef(dbxref);
        feature_dbxref.setFeature(feat);
        feature_dbxref.setCurrent(rs.getBoolean("is_current"));
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
   * Return the list of FeatureSynonyms for a given Feature, <b>specified by name</b>, or all if 
   * <code>null</code> is passed
   * 
   * @param uniqueName the uniquename of a Feature, or null for all
   * @return a (possibly empty) List<FeatureSynonymI> of matching synonyms
   */
  public List getFeatureSynonymsByFeatureUniquename(final String uniqueName)
  {
    String sql = "SELECT fs.*, s.name, s.type_id FROM "+
    "feature_synonym fs "+
    "LEFT JOIN feature f ON f.feature_id=fs.feature_id "+
    "LEFT JOIN synonym s ON fs.synonym_id=s.synonym_id ";
  
    if(uniqueName != null)
      sql = sql + " WHERE uniquename='" + uniqueName + "'";

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
        CvTerm cvterm = new CvTerm();
        cvterm.setCvTermId(rs.getInt("type_id"));
        Synonym syn = new Synonym();
        syn.setName(rs.getString("name"));
        syn.setCvTerm(cvterm);
        Feature feat = new Feature();
        feat.setFeatureId(rs.getInt("feature_id"));

        alias.setSynonym(syn);
        alias.setFeature(feat);
        
        Pub pub = new Pub();
        pub.setPubId(rs.getInt("pub_id"));
        alias.setPub(pub);
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
  
  
  public List getAllFeatureSynonymsAsFeature()
  {
    String sql = "SELECT fs.*, s.name, s.type_id , s.synonym_id FROM "+
    "feature_synonym fs "+
    "LEFT JOIN synonym s ON fs.synonym_id=s.synonym_id ORDER BY feature_id";

    appendToLogFile(sql, sqlLog);
    
    try
    {
      Statement st = conn.createStatement(ResultSet.TYPE_SCROLL_INSENSITIVE,
          ResultSet.CONCUR_UPDATABLE);
      ResultSet rs = st.executeQuery(sql);
      List features = new Vector();
      FeatureSynonym alias;
      while(rs.next())
      {
        int feature_id = rs.getInt("feature_id");
        
        Feature feature = new Feature();
        feature.setFeatureId(feature_id);
        java.util.Collection synonyms = feature.getFeatureSynonyms();
        if(synonyms == null || synonyms.size() == 0)
          synonyms = new Vector();
        int next_feature_id = -1;
        
        do
        {
          alias = new FeatureSynonym();
          CvTerm cvterm = new CvTerm();
          cvterm.setCvTermId(rs.getInt("type_id"));
          Synonym syn = new Synonym();
          syn.setName(rs.getString("name"));
          syn.setCvTerm(cvterm);
          Feature feat = new Feature();
          feat.setFeatureId(rs.getInt("feature_id"));

          alias.setSynonym(syn);
          alias.setFeature(feat);
        
          Pub pub = new Pub();
          pub.setPubId(rs.getInt("pub_id"));
          alias.setPub(pub);
          alias.setInternal(rs.getBoolean("is_internal"));
          alias.setCurrent(rs.getBoolean("is_current"));
          synonyms.add(alias);
          
          if(rs.next())  
          {
            next_feature_id = rs.getInt("feature_id");
            if(feature_id != next_feature_id)
              rs.previous();
          }
          else
            next_feature_id = -1;
          feature.setFeatureSynonyms(synonyms);
        }
        while(feature_id == next_feature_id);

        
        features.add(feature);
      }

      return features;
    }
    catch(SQLException sqle)
    {
      throw new RuntimeException(sqle);
    }
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
  
  public List getFeatureCvTermDbXRefByFeature(Feature feature)
  { 
    String sql = "SELECT fcd.feature_cvterm_id, dbx.*, db.name "+
      "FROM feature_cvterm_dbxref fcd "+
      "LEFT JOIN dbxref dbx ON dbx.dbxref_id=fcd.dbxref_id "+
      "LEFT JOIN db ON db.db_id=dbx.db_id";

    if(feature != null && feature.getUniqueName() != null)
      sql = sql+ " " +
          "LEFT JOIN feature_cvterm fc ON fcd.feature_cvterm_id=fc.feature_cvterm_id "+
          "WHERE feature_id=(SELECT feature_id FROM feature WHERE uniquename='"+
          feature.getUniqueName()+"')";
    
    appendToLogFile(sql, sqlLog);
    
    try
    {
      Statement st = conn.createStatement();
      ResultSet rs = st.executeQuery(sql);
      List featureCvTermDbXRefs = new Vector();
      FeatureCvTermDbXRef featureCvTermDbXRef;
      while(rs.next())
      {
        featureCvTermDbXRef = new FeatureCvTermDbXRef();
        FeatureCvTerm featureCvTerm = new FeatureCvTerm();
        featureCvTerm.setFeatureCvTermId( rs.getInt("feature_cvterm_id") );
        featureCvTermDbXRef.setFeatureCvTerm(featureCvTerm);
        
        DbXRef dbXRef = new DbXRef();
        dbXRef.setAccession( rs.getString("accession") );
        dbXRef.setDescription( rs.getString("description") );
        dbXRef.setVersion( rs.getString("version") );
        
        Db db = new Db();
        db.setName( rs.getString("name") );
        
        dbXRef.setDb(db);
        
        featureCvTermDbXRef.setDbXRef(dbXRef);
        featureCvTermDbXRefs.add(featureCvTermDbXRef);
      }
      
      return featureCvTermDbXRefs;
    }
    catch(SQLException sqle)
    {
      throw new RuntimeException(sqle);
    }
  }
  
  public List getFeatureCvTermPubByFeature(Feature feature)
  {
    String sql = "SELECT fcp.feature_cvterm_id, pub.* " +
                 "FROM feature_cvterm_pub fcp " +
                 "LEFT JOIN pub ON fcp.pub_id=pub.pub_id ";
    
    if(feature != null && feature.getUniqueName() != null)
      sql = sql+ " " +
          "LEFT JOIN feature_cvterm fc ON fcp.feature_cvterm_id=fc.feature_cvterm_id "+
          "WHERE feature_id=(SELECT feature_id FROM feature WHERE uniquename='"+
          feature.getUniqueName()+"')";
    
    appendToLogFile(sql, sqlLog);
    
    try
    {
      Statement st = conn.createStatement();
      ResultSet rs = st.executeQuery(sql);
      List featureCvTermPubs = new Vector();
      FeatureCvTermPub featureCvTermPub;
      while(rs.next())
      {
        featureCvTermPub = new FeatureCvTermPub();
        FeatureCvTerm featureCvTerm = new FeatureCvTerm();
        featureCvTerm.setFeatureCvTermId( rs.getInt("feature_cvterm_id") );
        featureCvTermPub.setFeatureCvTerm(featureCvTerm);
        
        Pub pub = new Pub();
        pub.setUniqueName(rs.getString("uniquename"));
        
        featureCvTermPub.setPub(pub);
        featureCvTermPubs.add(featureCvTermPub);
      }
      
      return featureCvTermPubs;
    }
    catch(SQLException sqle)
    {
      throw new RuntimeException(sqle);
    }
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
  private Feature getLazyFeature(final Feature feature)
  {
    List list = getFeatureQuery(feature.getUniqueName(), 
                                -1, feature.getFeatureId(),
                                feature.getCvTerm());
    if(list == null || list.size() < 1)
      return null;
    
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
                               final int feature_id,
                               final CvTerm cvTerm)
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
          + " LEFT JOIN organism ON organism.organism_id=f.organism_id ";
      
      if(cvTerm != null && cvTerm.getName() != null)
        sql = sql + "LEFT JOIN cvterm ON f.type_id=cvterm.cvterm_id ";
      
      sql = sql + " WHERE ";

      if(uniquename != null)
        sql = sql + "uniquename LIKE '" + uniquename + "'";

      if(parentFeatureID > -1)
        sql = sql + "srcfeature_id = " + parentFeatureID;

      if(feature_id > -1)
        sql = sql + "f.feature_id = " + feature_id;

      if(cvTerm != null && cvTerm.getName() != null)
        sql = sql + " AND cvterm.name="+cvTerm.getName();
      
      sql = sql
          + " ORDER BY f.type_id, uniquename";

      appendToLogFile(sql, sqlLog);
      ResultSet rs = st.executeQuery(sql);
 
      while(rs.next())
      {
        int feat_id = rs.getInt("feature_id");
        
        Feature feature = new Feature();
           
        FeatureLoc featureloc = new FeatureLoc();
        featureloc.setFmin(new Integer(rs.getInt("fmin")));
        featureloc.setFmax(new Integer(rs.getInt("fmax")));
        featureloc.setStrand(new Short(rs.getShort("strand")));

        int phase = rs.getInt("phase");
        if(rs.wasNull())
          featureloc.setPhase(null);
        else
          featureloc.setPhase(new Integer(phase));

        feature.setResidues(rs.getBytes("residues"));

        feature.setFeatureLoc(featureloc);
        feature.setCvTerm(new CvTerm());
        feature.getCvTerm().setCvTermId(rs.getInt("type_id"));
        feature.setUniqueName(rs.getString("uniquename"));
        feature.setTimeLastModified(rs.getTimestamp("timelastmodified"));
        feature.setFeatureId(rs.getInt("feature_id"));

        // feature organism
        Organism organism = new Organism();
        organism.setAbbreviation(rs.getString("abbreviation"));
        organism.setComment(rs.getString("comment"));
        organism.setCommonName(rs.getString("common_name"));
        organism.setGenus(rs.getString("genus"));
        organism.setOrganismId(rs.getInt("organism_id"));
        organism.setSpecies(rs.getString("species"));
        feature.setOrganism(organism);

        boolean next = false;
        do
        {
          // feature properties
          int prop_type_id = rs.getInt("prop_type_id");
          
          if(prop_type_id != 0)
          {
            FeatureProp featureprop = new FeatureProp();
            CvTerm cvterm = new CvTerm();
            cvterm.setCvTermId(prop_type_id);
            featureprop.setCvTerm(cvterm);
            featureprop.setValue(rs.getString("value"));

            if(feature.getFeatureProps() == null
                || feature.getFeatureProps().size() == 0)
              feature.setFeatureProps(new Vector());
            feature.addFeatureProp(featureprop);
          }
          else 
            feature.setFeatureProps(new Vector(0));
          
          // feature relationship
          FeatureRelationship feature_relationship = new FeatureRelationship();
          CvTerm cvterm = new CvTerm();
          cvterm.setCvTermId(rs.getInt("relation_type_id"));
          feature_relationship.setCvTerm(cvterm);

          int obj_id = rs.getInt("object_id");
          
          if(obj_id != 0)
          {
            Feature object = new Feature();
            object.setFeatureId(obj_id);
            feature_relationship.setFeatureByObjectId(object);

            if(feature.getFeatureRelationshipsForSubjectId() == null
                || feature.getFeatureRelationshipsForSubjectId().size() == 0)
              feature.setFeatureRelationshipsForSubjectId(new Vector());

            feature.addFeatureRelationshipsForSubjectId(feature_relationship);
          }
          else
            feature.setFeatureRelationshipsForSubjectId(new Vector(0));
          
          if(!rs.isLast())
          {
            rs.next();
            if(rs.getInt("feature_id") == feat_id)
              next = true;
            else
            {
              rs.previous();
              next = false;
            }
          }
          else
            next = false;
        }while(next);
        
        
        list.add(feature);
      }
    }
    catch(SQLException sqle)
    {
      throw new RuntimeException(sqle);
    }
    // merge same features in the list
    //return mergeList(list);
    return list;
  }
  
  //////
  ////// SchemaDaoI
  //////
  //////

  /**
   * Return a list of chado features with residues.
   * @return the <code>List</code> of <code>Feature</code> objects
   */
  public List getResidueFeatures()
  {
    String sql = new String(
            "SELECT uniquename, name, feature_id, type_id FROM ");
    
    sql = sql + "feature WHERE residues notnull ";

    appendToLogFile(sql, sqlLog);

    List list = new Vector();
    try
    {
      Statement st = conn.createStatement();
      ResultSet rs = st.executeQuery(sql);
      while(rs.next())
      {
        Feature feature = new Feature();
        
        feature.setFeatureId(rs.getInt("feature_id"));
        feature.setName(rs.getString("name"));
        feature.setUniqueName(rs.getString("uniquename"));
        feature.setCvTerm(new CvTerm());
        feature.getCvTerm().setCvTermId(rs.getInt("type_id"));

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
    String sql = "SELECT DISTINCT type_id FROM ";
    
    if(schema != null || !schema.equals(""))
      sql = sql + schema +"." ;
    sql = sql + "feature WHERE residues notnull";
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

  //////
  ////// CvDaoI
  //////
  //////
  
  public List getCvTermByNameInCv(String cvTermName, Cv cv)
  {
    return null;
  }
  
  public List getAllCvs()
  {
    return null;
  }
  
  /**
   * Get the full list of cvterm_id and name as a <code>List</code> of 
   * <code>CvTerm</code> objects.
   * @return    the full list of cvterm_id and name
   */
  public List getCvTerms()
  {
    String sql = "SELECT cvterm.cvterm_id, cvterm.name as cvterm_name, cv.NAME as cv_name, accession " +
                 "FROM cvterm " +
                 "LEFT JOIN dbxref ON dbxref.dbxref_id=cvterm.dbxref_id LEFT JOIN cv ON cv.cv_id = cvterm.cv_id";

    appendToLogFile(sql, sqlLog);
    
    try
    {
      Statement st = conn.createStatement();
      ResultSet rs = st.executeQuery(sql);
      List cvterms = new Vector();

      while(rs.next())
      {
        CvTerm cvterm = new CvTerm();
        cvterm.setCvTermId(rs.getInt("cvterm_id"));
        cvterm.setName(rs.getString("cvterm_name"));
        Cv cv = new Cv();
        cv.setName(rs.getString("cv_name"));
        cvterm.setCv(cv);
        DbXRef dbXRef = new DbXRef();
        dbXRef.setAccession(rs.getString("accession"));
        cvterm.setDbXRef(dbXRef);
        
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
  
  public CvTerm getCvTermById(final int cvTermId)
  {
    // TODO Auto-generated method stub
    return null;
  }
  
  //////
  ////// OrganismDaoI
  //////
  //////
  
  public List getOrganisms()
  {
    String sql = "SELECT organism_id AS organismId, abbreviation, "+
      "genus, species, common_name AS commonName, comment "+ 
      "FROM organism ORDER BY commonName";
    
    appendToLogFile(sql, sqlLog);
    List organisms = new Vector();
    
    try
    {
      Statement st = conn.createStatement();
      ResultSet rs = st.executeQuery(sql);

      while(rs.next())
      {
        Organism organism = new Organism();
        organism.setOrganismId(rs.getInt("organismId"));
        organism.setAbbreviation(rs.getString("abbreviation"));
        organism.setGenus(rs.getString("genus"));
        organism.setSpecies(rs.getString("species"));
        organism.setCommonName(rs.getString("commonName"));
        organism.setComment(rs.getString("comment"));
        organisms.add(organism);
      }
    }
    catch(SQLException sqle)
    {
      throw new RuntimeException(sqle);
    }
    return organisms;
  }
  
  //////
  ////// PubDaoI
  //////
  //////
  
  public List getPubDbXRef()
  {
    String sql = "SELECT pub_id, pub_dbxref.dbxref_id, "+
       "accession, version, dbx.description AS dbx_description, "+ 
       "db.name, db.description, db.urlprefix, db.url FROM pub_dbxref "+
       "LEFT JOIN dbxref dbx ON pub_dbxref.dbxref_id=dbx.dbxref_id "+
       "LEFT JOIN db ON db.db_id=dbx.db_id";
    
    appendToLogFile(sql, sqlLog);
    List pubDbXRefs = new Vector();
    
    try
    {
      Statement st = conn.createStatement();
      ResultSet rs = st.executeQuery(sql);
      
      PubDbXRef pubDbXRef;
      while(rs.next())
      {
        pubDbXRef = new PubDbXRef();
        Pub pub = new Pub();
        pub.setPubId(rs.getInt("pub_id"));
        pubDbXRef.setPub(pub);
        
        DbXRef dbXRef = new DbXRef();
        dbXRef.setAccession(rs.getString("accession"));
        dbXRef.setDescription(rs.getString("description"));
        dbXRef.setVersion(rs.getString("version"));
        
        Db db = new Db();
        db.setName(rs.getString("name"));
        db.setDescription(rs.getString("description"));
        db.setUrl(rs.getString("url"));
        db.setUrlPrefix(rs.getString("urlPrefix"));
        
        dbXRef.setDb(db);
        pubDbXRef.setDbXRef(dbXRef);
        pubDbXRefs.add(pubDbXRef);
      }
    }
    catch(SQLException sqle)
    {
      throw new RuntimeException(sqle);
    }
    
    return pubDbXRefs;
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
    {
      if(o instanceof FeatureForUpdatingResidues)
        updateFeatureResidues((FeatureForUpdatingResidues)o);
      else
        updateFeature((Feature)o);
    }   
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
    else if(o instanceof FeatureDbXRef)
      insertFeatureDbXRef((FeatureDbXRef)o);
    else if(o instanceof FeatureSynonym)
      insertFeatureAlias((FeatureSynonym)o);
    else if(o instanceof FeatureCvTerm)
      insertAllFeatureCvTerm((FeatureCvTerm)o);
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
    else if(o instanceof FeatureDbXRef)
      deleteFeatureDbXRef((FeatureDbXRef)o);
    else if(o instanceof FeatureSynonym)
      deleteFeatureSynonym((FeatureSynonym)o);
    else if(o instanceof FeatureCvTerm)
      deleteFeatureCvTerm((FeatureCvTerm)o);
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
      pstmt.setInt(1, featureloc.getFmin().intValue());
      pstmt.setInt(2, featureloc.getFmax().intValue());
      pstmt.setInt(3, featureloc.getRank());
      pstmt.setShort(4, featureloc.getStrand().shortValue());
      
      if(featureloc.getPhase() != null)
        pstmt.setInt(5, featureloc.getPhase().intValue());
      else
        pstmt.setNull(5, Types.INTEGER);
      
      pstmt.setString(6, featureloc.getFeatureByFeatureId().getUniqueName());
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

    if(feature.getCvTerm() != null)
      sql = sql+", type_id=?";
    
    if(feature.getTimeLastModified() != null)
      sql = sql+", timelastmodified=?";

    sql = sql+"WHERE feature_id=?";
    
    try
    {
      PreparedStatement pstmt = conn.prepareStatement(sql);
      pstmt.setString(1, feature.getUniqueName());
      int param = 2;
      if(feature.getCvTerm() != null)
      {
        pstmt.setLong(param, feature.getCvTerm().getCvTermId());
        param++;
      }

      if(feature.getTimeLastModified() != null)
      {
        pstmt.setTimestamp(param, feature.getTimeLastModified());
        param++;
      }

      pstmt.setInt(param, feature.getFeatureId());
      appendToLogFile(sql, sqlLog);
      pstmt.executeUpdate();
    }
    catch(SQLException sqle)
    {
      throw new RuntimeException(sqle);
    }
  }
  
  
  /**
   * Update feature residues (inserting or deleting a dna sequence).
   * @param feature the new <code>FeatureForUpdatingResidues</code> object.
   */
  private void updateFeatureResidues(FeatureForUpdatingResidues feature)
  {
    String sql1 =
      "UPDATE featureloc SET ";
    
    if(feature.getNewSubSequence() != null)
      sql1 = sql1 + "fmin=fmin+" + feature.getLength() + 
                 " , fmax=fmax+" + feature.getLength();
    else
      sql1 = sql1 + "fmin=fmin-" + feature.getLength() + 
                 " , fmax=fmax-" + feature.getLength();

    sql1 = sql1 + " WHERE fmin >= " + feature.getStartBase() +
                  " AND srcfeature_id="+feature.getFeatureId();
    appendToLogFile(sql1, sqlLog);  
      
    String sql2 = " UPDATE feature SET "+
            "residues=substring(residues from 1 for "+ feature.getStartBase() + ") || ";
    
    if(feature.getNewSubSequence() != null)
      sql2 = sql2 + "'" + feature.getNewSubSequence() + "' || ";
    
    sql2 = sql2 + "substring(residues from "+ feature.getEndBase() + 
                                   " for "+ feature.getSeqLen() + "), "+
                                   "seqlen=" + feature.getSeqLen() +
                                   " WHERE feature_id="+feature.getFeatureId();

    appendToLogFile(sql2, sqlLog);

    try
    {
      Statement st = conn.createStatement();
      st.executeUpdate(sql1);
      st.executeUpdate(sql2);
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
      pstmt.setLong(3, featureprop.getCvTerm().getCvTermId());
      pstmt.setString(4, featureprop.getFeature().getUniqueName());

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
    
    sqlBuff.append("( (SELECT feature_id FROM feature WHERE uniquename=");
    sqlBuff.append("'"+ featureprop.getFeature().getUniqueName()+"')," );
    sqlBuff.append(featureprop.getCvTerm().getCvTermId()+", '");
    sqlBuff.append(featureprop.getValue()+"',");
    sqlBuff.append(featureprop.getRank()+" )");

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
    String uniquename = featureprop.getFeature().getUniqueName();

    sqlBuff.append("DELETE FROM featureprop WHERE ");

    if(uniquename != null)
      sqlBuff.append("feature_id="+
          "(SELECT feature_id FROM feature WHERE uniquename='"+
          uniquename+"') AND ");
    
    if(featureprop.getRank() > -1)
      sqlBuff.append("rank="+featureprop.getRank()+" AND ");
    
    if(featureprop.getValue() != null)
      sqlBuff.append("value="+featureprop.getValue()+" AND ");
    
    sqlBuff.append("type_id="+featureprop.getCvTerm().getCvTermId());
      
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
          + feature.getFeatureLoc().getFeatureBySrcFeatureId().getFeatureId() + "'";

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
      sql_buff.append("'" + feature.getUniqueName() + "'" + " , ");
      sql_buff.append(Long.toString(feature.getCvTerm().getCvTermId()));
      sql_buff.append(" )");

      sql = new String(sql_buff);
      appendToLogFile(sql, sqlLog);
      st = conn.createStatement();
      int rowCount = st.executeUpdate(sql);

      //
      // get the current feature_id sequence value
      final int feature_id = getCurrval("feature_feature_id_seq");

      //
      // insert feature location into featureloc
      sql_buff = new StringBuffer();
      sql_buff.append("INSERT INTO featureloc (");
      sql_buff.append(" featureloc_id ,");
      sql_buff.append(" feature_id ,");
      sql_buff.append(" srcfeature_id ,");
      sql_buff.append(" fmin ,");
      sql_buff.append(" fmax ,");
      sql_buff.append(" strand ");
      
      if(feature.getFeatureLoc().getPhase() != null)
        sql_buff.append(" , phase ");
      
      sql_buff.append(" ) VALUES ( ");
      sql_buff.append("nextval('featureloc_featureloc_id_seq') , ");
      sql_buff.append(feature_id + " , ");
      sql_buff.append(feature.getFeatureLoc().getFeatureBySrcFeatureId().getFeatureId() + " , ");
      sql_buff.append(feature.getFeatureLoc().getFmin() + " , ");
      sql_buff.append(feature.getFeatureLoc().getFmax() + " , ");
      sql_buff.append(feature.getFeatureLoc().getStrand());
      
      if(feature.getFeatureLoc().getPhase() != null)
        sql_buff.append(" , "+feature.getFeatureLoc().getPhase());
      
      sql_buff.append(" )");

      sql = new String(sql_buff);
      appendToLogFile(sql, sqlLog);
      st = conn.createStatement();
      rowCount = st.executeUpdate(sql);

      // insert feature relationships
      if(feature.getFeatureRelationshipsForSubjectId() != null)
      {
        List parents = (List)feature.getFeatureRelationshipsForSubjectId();
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
          sql_buff.append(feature_relationship.getFeatureBySubjectId().getUniqueName()
              + "'), ");
          sql_buff.append("(SELECT feature_id FROM feature WHERE uniquename='");
          sql_buff.append(feature_relationship.getFeatureByObjectId().getUniqueName()
              + "'), ");
          sql_buff.append(feature_relationship.getCvTerm().getCvTermId() + ")");

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
          + feature.getUniqueName() + "'";
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
   * @param feature_dbxref    the new <code>FeatureDbXRef</code>
   */
  private void insertFeatureDbXRef(final FeatureDbXRef feature_dbxref)
  {   
    // find database id
    String sql = "SELECT db_id FROM db WHERE name='"+
                 feature_dbxref.getDbXRef().getDb().getName()+"'";
    
    try
    {
      Statement st = conn.createStatement();
      ResultSet rs = st.executeQuery(sql);
      boolean exists = rs.next();

      if(!exists)
        throw new SQLException("No database called "
            + feature_dbxref.getDbXRef().getDb().getName() + " found (for "
            + feature_dbxref.getFeature().getUniqueName()
            + ") check the spelling!");

      final int db_id = rs.getInt("db_id");
      // find if accession exists already
      String sqlDbXRefId = "SELECT dbxref_id FROM dbxref WHERE accession='"
          + feature_dbxref.getDbXRef().getAccession() + "' AND db_id=" + db_id;

      appendToLogFile(sqlDbXRefId, sqlLog);
      rs = st.executeQuery(sqlDbXRefId);
      exists = rs.next();

      if(!exists)
      {
        // create a new accession entry in dbxref
        sql = "INSERT INTO dbxref ( db_id, accession ) " + "VALUES (" + db_id
            + ", '" + feature_dbxref.getDbXRef().getAccession() + "' )";

        appendToLogFile(sql, sqlLog);
        int rowCount = st.executeUpdate(new String(sql));

        // now get the new dbxref_id
        appendToLogFile(sqlDbXRefId, sqlLog);
        rs = st.executeQuery(sqlDbXRefId);
        rs.next();
      }

      final int dbxref_id = rs.getInt("dbxref_id");
      sql = "INSERT INTO feature_dbxref "
          + "(feature_id, dbxref_id, is_current)" + " VALUES "
          + "( (SELECT feature_id FROM " + "feature WHERE  uniquename='"
          + feature_dbxref.getFeature().getUniqueName() + "'), " + dbxref_id
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
   * @param feature_dbxref  the  new <code>FeatureDbXRef</code>
   */
  private void deleteFeatureDbXRef(final FeatureDbXRef feature_dbxref)
  {
    final String uniquename = feature_dbxref.getFeature().getUniqueName();
    
    final String sql = 
      "DELETE FROM feature_dbxref "+
      "WHERE dbxref_id="+
      "(SELECT dbxref_id FROM dbxref WHERE accession='"+
         feature_dbxref.getDbXRef().getAccession()+"' "+
      "AND db_id=(SELECT db_id FROM db WHERE name='"+
         feature_dbxref.getDbXRef().getDb().getName()+"'))"+
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
    final String uniquename   = feature_synonym.getFeature().getUniqueName();
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
        String type_id = Long.toString(feature_synonym.getSynonym().getCvTerm()
            .getCvTermId());

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
    final String uniquename   = feature_synonym.getFeature().getUniqueName();
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
   * Delete featureCvTerm and update associated feature_cvterm.rank's 
   * if appropriate
   * @param featureCvTerm
   */
  private void deleteFeatureCvTerm(FeatureCvTerm feature_cvterm)
  {
    final String sql = "DELETE FROM feature_cvterm WHERE feature_cvterm_id="+
          feature_cvterm.getFeatureCvTermId();
      
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
    sqlBuff.append(" type_id="+feature_relationship.getCvTerm().getCvTermId());
    sqlBuff.append(" WHERE ");
    sqlBuff.append("subject_id=");
    sqlBuff.append("( SELECT feature_id FROM feature WHERE uniquename='");
    sqlBuff.append(feature_relationship.getFeatureBySubjectId().getUniqueName()+"' ) ");
    sqlBuff.append("AND ");
    sqlBuff.append("object_id=");
    sqlBuff.append("( SELECT feature_id FROM feature WHERE uniquename='");
    sqlBuff.append(feature_relationship.getFeatureByObjectId().getUniqueName()+"' ) ");
 
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

  protected int getCurrval(String seq_id)
  {
    int currval;
    try
    {
      //
      // get the current feature_id sequence value
      String sql = "SELECT currval('"+seq_id+"')";
      appendToLogFile(sql, sqlLog);
      Statement st = conn.createStatement();
      ResultSet rs = st.executeQuery(sql);
      rs.next();
      currval = rs.getInt("currval");
    }
    catch(SQLException sqle)
    {
      throw new RuntimeException(sqle);
    }
    return currval;
  }

  protected Integer getDbId(Db db)
  {
    Integer db_id;
    try
    {
      String sql = "SELECT db_id FROM db WHERE name='"+db.getName()+"'";
      appendToLogFile(sql, sqlLog);
      Statement st = conn.createStatement();
      ResultSet rs = st.executeQuery(sql);
      rs.next();
      db_id = new Integer(rs.getInt("db_id"));
    }
    catch(SQLException sqle)
    {
      throw new RuntimeException(sqle);
    }

    return db_id;
  }

  protected Integer getDbXRefId(DbXRef dbXRef)
  {
    Integer dbxref_id;
    try
    {
      String sql = "SELECT dbxref_id FROM dbxref WHERE accession='"+
                    dbXRef.getAccession()+"' AND db_id="+dbXRef.getDb().getDbId();
      appendToLogFile(sql, sqlLog);
      Statement st = conn.createStatement();
      ResultSet rs = st.executeQuery(sql);
      rs.next();
      dbxref_id = new Integer(rs.getInt("dbxref_id"));
    }
    catch(SQLException sqle)
    {
      throw new RuntimeException(sqle);
    }

    return dbxref_id;
  }

  protected Pub getPubByUniqueName(Pub pub)
  {
    try
    {
      String sql = "SELECT * FROM pub WHERE uniquename='"+pub.getUniqueName()+"'";
      appendToLogFile(sql, sqlLog);
      Statement st = conn.createStatement();
      ResultSet rs = st.executeQuery(sql);
      rs.next();
      
      pub.setPubId( rs.getInt("pub_id") );
      //
    }
    catch(SQLException sqle)
    {
      throw new RuntimeException(sqle);
    }

    return pub;
  }

  protected void insertDbXRef(DbXRef dbXRef)
  {
    try
    {
      String sql = "INSERT INTO dbxref ( db_id, accession, version ) VALUES ("+
                   dbXRef.getDb().getDbId() +", '"+ 
                   dbXRef.getAccession()    +"', "+ 
                   dbXRef.getVersion()+")";
      Statement st = conn.createStatement();
      appendToLogFile(sql, sqlLog);
      st.executeUpdate(new String(sql));
    }
    catch(SQLException sqle)
    {
      throw new RuntimeException(sqle);
    }   
  }

  protected void insertFeatureCvTerm(FeatureCvTerm feature_cvterm)
  {
    try
    {
      String uniqueName = feature_cvterm.getFeature().getUniqueName();
      final String pubIdStr;
      
      if(feature_cvterm.getPub() != null)
      {
        if(feature_cvterm.getPub().getPubId() == 0)
          pubIdStr = "(SELECT pub_id FROM pub WHERE uniquename="+
                      feature_cvterm.getPub().getUniqueName()+")";
        else
          pubIdStr = Integer.toString(feature_cvterm.getPub().getPubId());
      }
      else
        pubIdStr = "0";
      
      String sql = "INSERT INTO feature_cvterm "+
        "( feature_cvterm_id, feature_id, cvterm_id, pub_id, rank, is_not ) "+
        "VALUES "+
        "( nextval('feature_cvterm_feature_cvterm_id_seq'), "+
          "(SELECT feature_id FROM feature WHERE uniquename='"+ uniqueName + "'), "+
          feature_cvterm.getCvTerm().getCvTermId() + " , "+
          pubIdStr + " , "+
          feature_cvterm.getRank() + " , "+
          feature_cvterm.isNot() + " ) ";
      appendToLogFile(sql, sqlLog);
      Statement st = conn.createStatement();
      st.executeUpdate(new String(sql));
    }
    catch(SQLException sqle)
    {
      throw new RuntimeException(sqle);
    }   
  }

  protected void insertFeatureCvTermDbXRef(FeatureCvTermDbXRef featureCvTermDbXRef)
  {
    try
    {
      String sql = "INSERT INTO feature_cvterm_dbxref "+
        "( feature_cvterm_id, dbxref_id ) "+
        "VALUES ( "+
        featureCvTermDbXRef.getFeatureCvTerm().getFeatureCvTermId() + ", " +
        featureCvTermDbXRef.getDbXRef().getDbXRefId()+" )";
      appendToLogFile(sql, sqlLog);
      Statement st = conn.createStatement();
      st.executeUpdate(new String(sql));
    }
    catch(SQLException sqle)
    {
      throw new RuntimeException(sqle);
    }   
  }

  protected void insertFeatureCvTermProp(FeatureCvTermProp featureCvTermProp)
  {
    try
    {
      String sql = "INSERT INTO feature_cvtermprop "+
        "( feature_cvtermprop_id, feature_cvterm_id, type_id, value, rank ) "+
        "VALUES ( "+
        "nextval('feature_cvtermprop_feature_cvtermprop_id_seq'), "+
        featureCvTermProp.getFeatureCvTerm().getFeatureCvTermId() + ", "+
        featureCvTermProp.getCvTerm().getCvTermId() + ", '"+
        featureCvTermProp.getValue() +"', "+
        featureCvTermProp.getRank() +")";
      appendToLogFile(sql, sqlLog);
      Statement st = conn.createStatement();
      st.executeUpdate(new String(sql));
    }
    catch(SQLException sqle)
    {
      throw new RuntimeException(sqle);
    }
  }

  protected void insertFeatureCvTermPub(FeatureCvTermPub featureCvTermPub)
  {
    try
    {
      String sql = "INSERT INTO feature_cvterm_pub "+
        "( feature_cvterm_id, pub_id ) "+
        "VALUES ( "+
        featureCvTermPub.getFeatureCvTerm().getFeatureCvTermId() + ", "+
        featureCvTermPub.getPub().getPubId() + ")";
      appendToLogFile(sql, sqlLog);
      Statement st = conn.createStatement();
      st.executeUpdate(new String(sql));
    }
    catch(SQLException sqle)
    {
      throw new RuntimeException(sqle);
    }
  }

  protected void insertPub(Pub pub)
  {
    try
    {
      String sql = "INSERT INTO pub ( ";
      
      if(pub.getTitle() != null)
        sql = sql + "title, ";
      
      /*if(pub.getVolumeTitle() != null)
        sql = sql + "volumetitle, ";
      
      if(pub.getVolume() != null)
        sql = sql + "volume, ";
      
      if(pub.getSeriesName() != null)
        sql = sql + "series_name, ";
      
      if(pub.getIssue() != null)
        sql = sql + "issue, ";
      
      if(pub.getPyear() != null)
        sql = sql + "pyear, ";
      
      if(pub.getPages() != null)
        sql = sql + "pages, ";
      
      if(pub.getMiniRef() != null)
        sql = sql + "miniref, ";

       sql = sql + "uniquename, type_id ";
       
       if(pub.isObsolete() != null)
         sql = sql + " , is_obsolete ";
       
       if(pub.getPublisher() != null)
         sql = sql + " , publisher ";
       
       if(pub.getPubPlace() != null)
         sql = sql + " , pubplace ";*/
       
       sql = sql + ") VALUES (";

       if(pub.getTitle() != null)
         sql = sql + pub.getTitle() + ", ";
       
       // FIX THIS
       //sql = sql + "'" + pub.getUniqueName() + "'," + 
       //                  pub.getCvTerm().getCvTermId()+ ")";
       
      appendToLogFile(sql, sqlLog);
      Statement st = conn.createStatement();
      st.executeUpdate(new String(sql));
    }
    catch(SQLException sqle)
    {
      throw new RuntimeException(sqle);
    }
  }
  
  protected void insertPubDbXRef(PubDbXRef pubDbXRef)
  {
    
  }
  
  public Organism getOrganismByCommonName(String commonName)
  {
    return null; 
  }

  public Graph getGraph(Integer graphId)
  {
    // TODO Auto-generated method stub
    return null;
  }

  public List getGraphs(Integer featureId)
  {
    // TODO Auto-generated method stub
    return null;
  }

  public List getFeaturesByUniqueNames(List arg0)
  {
    // TODO Auto-generated method stub
    return null;
  }

  public List getTableColumns(String tableName)
  {
    // TODO Auto-generated method stub
    return null;
  }
}
