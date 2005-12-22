/* DatabaseDocument.java
 *
 * created: 2005
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

package uk.ac.sanger.artemis.util;

import com.ibatis.sqlmap.client.SqlMapClient;
import uk.ac.sanger.ibatis.*;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.chado.ChadoTransaction;

import java.sql.*;
import java.io.*;
import java.util.Hashtable;
import java.util.Vector;
import java.util.Enumeration;
import java.util.Date;
import java.util.List;
import java.util.Iterator;
import javax.swing.JOptionPane;
import javax.swing.JPasswordField;

/**
 * Objects of this class are Documents created from a relational database.
 * 
 */

public class DatabaseDocument extends Document
{
  private String name = null;

  private String feature_id = "1";

  /** database schema */
  private String schema = "public";

  private static Hashtable cvterm;

  private InputStreamProgressListener progress_listener;

  private Hashtable db;

  private Vector organism;

  private Hashtable org2schema;

  private String sqlLog = System.getProperty("user.home") +
                          System.getProperty("file.separator") + 
                          "art_sql_debug.log";

  private ByteBuffer[] gff_buffer;

  private ByteBuffer gff_buff;

  /** entries to split into */
  private String[] types = { "exon", "gene", "CDS", "transcript" };

  /** true if splitting the GFF into entries */
  private boolean splitGFFEntry;

  private boolean iBatis = false;

  private JPasswordField pfield;

  /**
   * 
   * Create a new Document from a database.
   * 
   * @param location
   *          This should be a URL string giving:
   *          jdbc:postgresql://host:port/datbase_name?user=username
   * 
   */
  public DatabaseDocument(String location, JPasswordField pfield)
  {
    super(location);
    this.pfield = pfield;

    if(System.getProperty("ibatis") != null)
    {
      iBatis = true;
      System.setProperty("chado", location);
    }
  }

  /**
   * 
   * Create a new Document from a database.
   * 
   * @param location
   *          This should be a URL string giving:
   *          jdbc:postgresql://host:port/datbase_name?user=username
   * @param feature_id
   *          ID of a feature to be extracted.
   * 
   */
  public DatabaseDocument(String location, JPasswordField pfield,
                          String feature_id, String schema)
  {
    super(location);
    this.pfield = pfield;

    this.feature_id = feature_id;
    this.schema = schema;

    if(System.getProperty("ibatis") != null)
    {
      iBatis = true;
      System.setProperty("chado", location);
    }
  }

  /**
   * 
   * Create a new Document from a database.
   * 
   * @param location
   *          This should be a URL string giving:
   *          jdbc:postgresql://host:port/datbase_name?user=username
   * @param feature_id
   *          ID of a feature to be extracted.
   * @param splitGFFEntry
   *          split into separate entries based on feature types.
   * @param progress_listener
   *          input stream progress listener
   * 
   */
  public DatabaseDocument(String location, JPasswordField pfield,
                          String feature_id, String schema, boolean splitGFFEntry,
                          InputStreamProgressListener progress_listener)
  {
    super(location);
    this.pfield = pfield;
    this.feature_id = feature_id;
    this.schema = schema;
    this.splitGFFEntry = splitGFFEntry;
    this.progress_listener = progress_listener;
    if(System.getProperty("ibatis") != null)
    {
      iBatis = true;
      System.setProperty("chado", location);
    }
  }

  public DatabaseDocument(String location, JPasswordField pfield,
                          String feature_id, String schema,
                          ByteBuffer gff_buff, String name)
  {
    super(location);
    this.pfield = pfield;
    this.feature_id = feature_id;
    this.schema = schema;
    this.gff_buff = gff_buff;
    this.name = name;
    if(System.getProperty("ibatis") != null)
    {
      iBatis = true;
      System.setProperty("chado", location);
    }
  }

  /**
   * 
   * Append a String to the Document location.
   * 
   * @param name
   *          The name to append.
   * 
   */
  public Document append(String name) throws IOException
  {
    return new DatabaseDocument( ((String)getLocation()) + name, pfield);
  }

  /**
   * 
   * Return the name of this Document (the last element of the Document
   * location).
   * 
   */
  public String getName()
  {
    if(name == null)
    {
      int ind     = ((String) getLocation()).indexOf("?");
      String name = ((String) getLocation()).substring(0, ind);
      ind = name.lastIndexOf("/");
      return name.substring(ind + 1);
    }
    return name;
  }

  /**
   * 
   * Return a Document with the last element stripped off.
   * 
   */
  public Document getParent()
  {
    return null;
  }

  /**
   * 
   * Return true if and only if the Document refered to by this object exists
   * and is readable. Always returns true.
   * 
   */
  public boolean readable()
  {
    return true;
  }

  /**
   * 
   * Return true if and only if the Document refered to by this object exists
   * and can be written to. Always returns false.
   * 
   */
  public boolean writable()
  {
    return true;
  }

  /**
   * Create a new InputStream object from this Document. The contents of the
   * Document can be read from the InputStream.
   * 
   * @exception IOException
   *              Thrown if the Document can't be read from (for example if it
   *              doesn't exist).
   */
  public InputStream getInputStream() throws IOException
  {
    ByteArrayInputStream instream;

    if(gff_buff != null)
    {
      instream = new ByteArrayInputStream(gff_buff.getBytes());
      return instream;
    }

    try
    {
      Connection conn = null;

      if(!iBatis)
        conn = getConnection();

      if(iBatis)
        gff_buffer = getGFFiBatis(feature_id, schema);
      else
        gff_buffer = getGFFJdbc(conn, feature_id, schema);

      ByteBuffer entry = new ByteBuffer();
      if(splitGFFEntry)
      {
        if(gff_buffer[0].size() > 0)
          entry.append(gff_buffer[0]);

        if(iBatis)
          getSequenceIbatis(entry, schema);
        else
          getSequence(conn, entry, schema);
      }
      else
      {
        for(int i = 0; i < gff_buffer.length; i++)
        {
          if(gff_buffer[i].size() > 0)
            entry.append(gff_buffer[i]);
        }

        if(iBatis)
          getSequenceIbatis(entry, schema);
        else
          getSequence(conn, entry, schema);
      }

      if(System.getProperty("debug") != null)
        appendToLogFile(new String(entry.getBytes()), sqlLog);

      instream = new ByteArrayInputStream(entry.getBytes());

      if(conn != null)
        conn.close();
      return instream;
    }
    catch(java.sql.SQLException sqlExp)
    {
      System.out.println("Problems connecting...");
      sqlExp.printStackTrace();
    }

    return null;
  }

  /**
   * 
   * Called (by DatabaseEntrySource) to retrieve all the documents for each
   * entry created.
   * 
   */
  public DatabaseDocument[] getGffDocuments(String location, String id,
                                            String schema)
  {
    int nentries = 0;
    for(int i = 1; i < gff_buffer.length; i++)
    {
      if(gff_buffer[i].size() > 0)
        nentries++;
    }

    DatabaseDocument[] new_docs = new DatabaseDocument[nentries];
    nentries = 0;
    for(int i = 1; i < gff_buffer.length; i++)
    {
      if(gff_buffer[i].size() == 0)
        continue;

      String name;
      if(i >= types.length)
        name = "other";
      else
        name = types[i];

      new_docs[nentries] = new DatabaseDocument(location, pfield, id, schema,
                                                gff_buffer[i], name);
      nentries++;
    }

    return new_docs;
  }

  /**
   * 
   * Return a feature name given the feature_id.
   * 
   */
  private String getFeatureNameIbatis(final Feature feature)
      throws java.sql.SQLException
  {
    SqlMapClient sqlMap = DbSqlConfig.getSqlMapInstance();
    return (String) sqlMap.queryForObject("getFeatureName", feature);
  }

  /**
   * 
   * Return a feature name given the feature_id.
   * 
   */
  private String getFeatureNameJdbc(String feature_id, Connection conn,
      String schema) throws java.sql.SQLException
  {
    Statement st = conn.createStatement();

    String sql = "SELECT name FROM " + schema + ".feature WHERE feature_id= " +
                  feature_id;
    appendToLogFile(sql, sqlLog);
    ResultSet rs = st.executeQuery(sql);
    rs.next();
    return rs.getString("name");
  }

  private ByteBuffer[] getGFFiBatis(String parentFeatureID, String schema)
      throws java.sql.SQLException
  {
    SqlMapClient sqlMap = DbSqlConfig.getSqlMapInstance();

    Feature feature = new Feature();
    feature.setId(Integer.parseInt(parentFeatureID));
    feature.setSchema(schema);

    List featList = sqlMap.queryForList("getGffLine", feature);

    ByteBuffer[] buffers = new ByteBuffer[types.length + 1];
    for(int i = 0; i < buffers.length; i++)
      buffers[i] = new ByteBuffer();

    String parentFeature = getFeatureNameIbatis(feature);
    Hashtable hstore = new Hashtable();
    ByteBuffer this_buff;

    int feature_size = featList.size();

    for(int i = 0; i < feature_size; i++)
    {
      Feature feat = (Feature)featList.get(i);
      int fmin                = feat.getFmin() + 1;
      int fmax                = feat.getFmax();
      long type_id            = feat.getType_id();
      long prop_type_id       = feat.getProp_type_id();
      int strand              = feat.getStrand();
      String name             = feat.getUniquename();
      String typeName         = getCvtermName(null, type_id);
      String propTypeName     = getCvtermName(null, prop_type_id);
      String timelastmodified = feat.getTimelastmodified().toString();
      String feature_id       = Integer.toString(feat.getId());
      hstore.put(feature_id, name);

      String parent_id = feat.getObject_id();
      if(parent_id != null && hstore.containsKey(parent_id))
        parent_id = (String)hstore.get(parent_id);

      // make gff format

      // select buffer
      this_buff = buffers[types.length];
      for(int j = 0; j < types.length; j++)
      {
        if(types[j].equals(typeName))
          this_buff = buffers[j];
      }

      this_buff.append(parentFeature + "\t"); // seqid
      this_buff.append("chado\t");            // source
      this_buff.append(typeName + "\t");      // type
      this_buff.append(fmin + "\t");          // start
      this_buff.append(fmax + "\t");          // end
      this_buff.append(".\t");                // score
      if(strand == -1)                       // strand
        this_buff.append("-\t");
      else if(strand == 1)
        this_buff.append("+\t");
      else
        this_buff.append(".\t");

      this_buff.append(".\t");               // phase
      this_buff.append("ID=" + name + ";");

      if(parent_id != null)
        this_buff.append("Parent=" + parent_id + ";");

      this_buff.append("timelastmodified=" + timelastmodified + ";");

      String value = "";
      if(feat.getValue() != null)
        value = GFFStreamFeature.encode(feat.getValue());

      this_buff.append(propTypeName + "=" + value); // attributes

      // is the next line part of the same feature, if so merge
      boolean rewind = false;
      Feature featNext = null;

      if(i < feature_size - 1)
        featNext = (Feature)featList.get(i + 1);

      // merge next line if part of the same feature
      while(featNext != null && featNext.getUniquename().equals(name))
      {
        prop_type_id = featNext.getProp_type_id();
        propTypeName = getCvtermName(null, prop_type_id);
        value = GFFStreamFeature.encode(featNext.getValue());
        this_buff.append(";" + propTypeName + "=" + value);
        i++;
        if(i < feature_size - 1)
          featNext = (Feature) featList.get(i + 1);
        else
          break;
      }

      this_buff.append("\n");

      progress_listener.progressMade("Read from database: " + name);
    }

    return buffers;
  }

  /**
   * 
   * Given a parent (chromosome, contig, supercontig) retrieve the features in
   * the form of a GFF stream.
   * 
   */
  private ByteBuffer[] getGFFJdbc(Connection conn, String parentFeatureID,
                                  String schema) 
          throws java.sql.SQLException
  {
    Statement st = conn.createStatement(ResultSet.TYPE_SCROLL_INSENSITIVE,
        ResultSet.CONCUR_UPDATABLE);

    String sql = "SELECT timelastmodified, feature.feature_id, object_id, strand, fmin, fmax, uniquename, "
        + schema + ".feature.type_id, "
        + schema + ".featureprop.type_id AS prop_type_id, featureprop.value"
        + " FROM  "
        + schema + ".featureloc, "
        + schema + ".feature"
        + " LEFT JOIN "
        + schema + ".feature_relationship ON "
        + schema + ".feature_relationship.subject_id="
        + schema + ".feature.feature_id"
        + " LEFT JOIN "
        + schema + ".featureprop ON "
        + schema + ".featureprop.feature_id="
        + schema + ".feature.feature_id"
        + " WHERE srcfeature_id = "
        + parentFeatureID + " AND "
        + schema + ".featureloc.feature_id="
        + schema + ".feature.feature_id"
        + " AND ("
        + schema + ".featureloc.rank="
        + schema + ".feature_relationship.rank OR "
        + schema + ".feature_relationship.rank IS NULL)"
        + " ORDER BY "
        + schema + ".feature.type_id,  uniquename";

    appendToLogFile(sql, sqlLog);
    ResultSet rs = st.executeQuery(sql);

    ByteBuffer[] buffers = new ByteBuffer[types.length + 1];
    for(int i = 0; i < buffers.length; i++)
      buffers[i] = new ByteBuffer();

    String parentFeature = getFeatureNameJdbc(parentFeatureID, conn, schema);
    Hashtable hstore = new Hashtable();
    ByteBuffer this_buff;

    while(rs.next())
    {
      int fmin                = rs.getInt("fmin") + 1;
      int fmax                = rs.getInt("fmax");
      long type_id            = rs.getLong("type_id");
      long prop_type_id       = rs.getLong("prop_type_id");
      int strand              = rs.getInt("strand");
      String name             = rs.getString("uniquename");
      String typeName         = getCvtermName(conn, type_id);
      String propTypeName     = getCvtermName(conn, prop_type_id);
      String timelastmodified = rs.getString("timelastmodified");
      String feature_id       = rs.getString("feature_id");
      hstore.put(feature_id, name);

      String parent_id = rs.getString("object_id");
      if(parent_id != null && hstore.containsKey(parent_id))
        parent_id = (String)hstore.get(parent_id);

      // make gff format

      // select buffer
      this_buff = buffers[types.length];
      for(int i = 0; i < types.length; i++)
      {
        if(types[i].equals(typeName))
          this_buff = buffers[i];
      }

      this_buff.append(parentFeature + "\t"); // seqid
      this_buff.append("chado\t");            // source
      this_buff.append(typeName + "\t");      // type
      this_buff.append(fmin + "\t");          // start
      this_buff.append(fmax + "\t");          // end
      this_buff.append(".\t");                // score
      if (strand == -1)                       // strand
        this_buff.append("-\t");
      else if (strand == 1)
        this_buff.append("+\t");
      else
        this_buff.append(".\t");

      this_buff.append(".\t");               // phase
      this_buff.append("ID=" + name + ";");

      if(parent_id != null)
        this_buff.append("Parent=" + parent_id + ";");

      this_buff.append("timelastmodified=" + timelastmodified + ";");

      String value = "";
      if(rs.getString("value") != null)
        value = GFFStreamFeature.encode(rs.getString("value"));

      this_buff.append(propTypeName + "=" + value); // attributes

      // is the next line part of the same feature, if so merge
      boolean rewind = false;
      while((rewind = rs.next()) && rs.getString("uniquename").equals(name))
      {
        prop_type_id = rs.getLong("prop_type_id");
        propTypeName = getCvtermName(conn, prop_type_id);
        value = GFFStreamFeature.encode(rs.getString("value"));
        this_buff.append(";" + propTypeName + "=" + value);
      }

      if(rewind)
        rs.previous();

      this_buff.append("\n");

      progress_listener.progressMade("Read from database: " + name);
    }

    return buffers;
  }

  public static Long getCvtermID(String name)
  {
    Enumeration enum_cvterm = cvterm.keys();
    while(enum_cvterm.hasMoreElements())
    {
      Long key = (Long)enum_cvterm.nextElement();
      if(name.equals(cvterm.get(key)))
        return key;
    }
    return null;
    // return new Long("-1.");
  }

  private String getCvtermName(Connection conn, long id)
  {
    if(cvterm == null)
    {
      if(iBatis)
        getCvtermIbatis(null);
      else
        getCvterm(conn, null);
    }

    return (String)cvterm.get(new Long(id));
  }

  /**
   * 
   * Look up cvterms names and id and return in a hashtable.
   * 
   */
  private Hashtable getCvtermIbatis(String cv_name)
  {
    cvterm = new Hashtable();

    try
    {
      SqlMapClient sqlMap = DbSqlConfig.getSqlMapInstance();

      List cvtem_list = sqlMap.queryForList("getCvterm", null);
      Iterator it = cvtem_list.iterator();

      while(it.hasNext())
      {
        Cvterm cv = (Cvterm)it.next();
        cvterm.put(new Long(cv.getId()), cv.getName());
      }
    }
    catch(SQLException sqle)
    {
      System.err.println(this.getClass() + ": SQLException retrieving CvTerms");
      System.err.println(sqle);
    }

    return cvterm;
  }

  /**
   * 
   * Look up cvterms names and id and return in a hashtable.
   * 
   */
  private Hashtable getCvterm(Connection conn, String cv_name)
  {
    String sql = "SELECT cvterm.cvterm_id, cvterm.name " +
                 "FROM cvterm, cv WHERE cv.cv_id = cvterm.cv_id";

    if(cv_name != null)
      sql = sql + " AND cv.name='" + cv_name + "'";

    appendToLogFile(sql, sqlLog);

    cvterm = new Hashtable();

    try
    {
      Statement s = conn.createStatement();
      ResultSet rs = s.executeQuery(sql);

      while(rs.next())
      {
        long id = rs.getLong("cvterm_id");
        String name = rs.getString("name");

        if(cvterm.get(name) != null)
          System.err.println(this.getClass()
              + ": WARNING - read multiple CvTerms with name = '" + name + "'");

        cvterm.put(new Long(id), name);
      }
    }
    catch(SQLException sqle)
    {
      System.err.println(this.getClass() + ": SQLException retrieving CvTerms");
      System.err.println(sqle);
    }

    return cvterm;
  }

  public ByteBuffer getSequenceIbatis(ByteBuffer buff, String schema)
      throws java.sql.SQLException
  {
    SqlMapClient sqlMap = DbSqlConfig.getSqlMapInstance();

    Feature feature = new Feature();
    feature.setId(Integer.parseInt(feature_id));
    feature.setSchema(schema);

    feature = (Feature)sqlMap.queryForObject("getSequence", 
                                              feature);
    buff.append("##FASTA\n>");
    buff.append(feature.getName());
    buff.append("\n");
    buff.append(feature.getResidues());
    return buff;
  }

  public ByteBuffer getSequence(Connection conn, ByteBuffer buff, String schema)
      throws java.sql.SQLException
  {
    Statement st = conn.createStatement();
    String sql = "SELECT name, residues from " + schema +
                 ".feature where feature_id = '" + feature_id + "'";

    appendToLogFile(sql, sqlLog);

    ResultSet rs = st.executeQuery(sql);
    rs.next();

    buff.append("##FASTA\n>");
    buff.append(rs.getBytes("name"));
    buff.append("\n");
    buff.append(rs.getBytes("residues"));
    return buff;

    // return "##FASTA\n>" + name + "\n" + rs.getString("residues");
  }

  public Hashtable getDatabaseEntries()
  {
    if(iBatis)
      return getDatabaseEntriesIbatis();
    else
      return getDatabaseEntriesJdbc();
  }

  public Hashtable getSchemaEntries()
  {
    return org2schema;
  }

  /**
   * 
   * Create a hashtable of the available entries.
   * 
   */
  private Hashtable getDatabaseEntriesJdbc()
  {
    db = new Hashtable();
    organism = new Vector();
    org2schema = new Hashtable();

    try
    {
      Connection conn = getConnection();
      System.out.println("JDBC Connection");

      Statement st = conn.createStatement();

      String query = "SELECT schema_name FROM information_schema.schemata "+ 
                     "WHERE schema_name=schema_owner ORDER BY schema_name";
      appendToLogFile(query, sqlLog);

      ResultSet rs = st.executeQuery(query);
      Vector schemas = new Vector();

      while(rs.next())
        schemas.add(rs.getString("schema_name"));
      
      for(int i = 0; i < schemas.size(); i++)
      {
        String schema = (String)schemas.get(i);
        appendToLogFile(schema, sqlLog);

        String sql = "SELECT DISTINCT type_id FROM " + schema +
                     ".feature WHERE residues notnull";
        appendToLogFile(sql, sqlLog);

        Vector cvterm_id = new Vector();
        rs = st.executeQuery(sql);

        while(rs.next())
          cvterm_id.add(rs.getString("type_id"));

        if(cvterm_id.size() == 0)  // no residues for this organism
          continue;

        sql = new String(
            "SELECT abbreviation, name, feature_id, type_id FROM organism, "+
            schema + ".feature WHERE (");

        for(int j = 0; j < cvterm_id.size(); j++)
        {
          sql = sql + " type_id = " + (String)cvterm_id.get(j);
          if(j < cvterm_id.size() - 1)
            sql = sql + " OR ";
        }

        sql = sql + ") and organism.organism_id=" + schema
            + ".feature.organism_id " + "and residues notnull "
            + "ORDER BY abbreviation";

        appendToLogFile(sql, sqlLog);

        rs = st.executeQuery(sql);
        while(rs.next())
        {
          String org      = rs.getString("abbreviation");
          String typeName = getCvtermName(conn, rs.getLong("type_id"));
          db.put(org + " - " + typeName + " - " + rs.getString("name"), 
                 rs.getString("feature_id"));
          if(!organism.contains(org))
            organism.add(org);
          if(!org2schema.containsKey(org))
          {
            org2schema.put(org, schema);
          }
        }
      }
      conn.close();
    }
    catch(java.sql.SQLException sqlExp)
    {
      JOptionPane.showMessageDialog(null, "SQL Problems...", "SQL Error",
                                    JOptionPane.ERROR_MESSAGE);
      sqlExp.printStackTrace();
    }
    catch (java.net.ConnectException conn)
    {
      JOptionPane.showMessageDialog(null, "Problems connecting...",
                                    "Database Connection Error - Check Server",
                                    JOptionPane.ERROR_MESSAGE);
      conn.printStackTrace();
    }

    return db;
  }

  /**
   * 
   * Create a hashtable of the available entries.
   * 
   */
  private Hashtable getDatabaseEntriesIbatis()
  {
    db = new Hashtable();
    organism = new Vector();
    org2schema = new Hashtable();

    try
    {
      DbSqlConfig.init(pfield);
      SqlMapClient sqlMap = DbSqlConfig.getSqlMapInstance();

      List schema_list = sqlMap.queryForList("getShema", null);
      Iterator it      = schema_list.iterator();

      while(it.hasNext())
      {
        String schema = (String)it.next();
  
        SchemaCVList schema_CVlist = new SchemaCVList();
        schema_CVlist.setSchema(schema);
         
        List list = sqlMap.queryForList("getResidueType", schema);
        if(list.size() == 0)  // no residues for this organism
          continue;

        schema_CVlist.setCvlist(list);

        List list_residue_features = sqlMap.queryForList("getSchemaResidueFeatures",
                                                         schema_CVlist);
        Iterator it_residue_features = list_residue_features.iterator();
        while(it_residue_features.hasNext())
        {
          Feature feature = (Feature)it_residue_features.next();
          String org      = feature.getAbbreviation();
          String typeName = getCvtermName(null, feature.getType_id());

          db.put(org + " - " + typeName + " - " + feature.getName(), 
                 Integer.toString(feature.getId()));
          if(!organism.contains(org))
            organism.add(org);
          if(!org2schema.containsKey(org))
            org2schema.put(org, schema);
        }
      }
    }
    catch(java.sql.SQLException sqlExp)
    {
      JOptionPane.showMessageDialog(null, "SQL Problems...", "SQL Error",
                                    JOptionPane.ERROR_MESSAGE);
      sqlExp.printStackTrace();
    }
    return db;
  }

  public Vector getOrganism()
  {
    return organism;
  }

  /**
   * 
   * Make a connetion with the jdbc
   * jdbc:postgresql://localhost:13001/chadoCVS?user=es2
   * 
   */
  public Connection getConnection() throws java.sql.SQLException,
      java.net.ConnectException
  {
    String location = (String)getLocation();
    if(pfield == null || pfield.getPassword().length == 0)
      return DriverManager.getConnection(location);

    // assume we have a password
    final int index = location.indexOf("?user=");
    return DriverManager.getConnection(location.substring(0, index), 
                                       location.substring(index + 6),
                                       new String(pfield.getPassword()));
  }

  /**
   * Create a new OutputStream object from this Document. The contents of the
   * Document can be written from the stream.
   * 
   * @exception IOException
   *              Thrown if the Document can't be written.
   */
  public OutputStream getOutputStream() throws IOException
  {
    System.out.println("DatabaseDocument - ReadOnlyException");
    throw new ReadOnlyException("this Database Document can not be written to");
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

  public void commit(Vector sql)
  {
    try
    {
      Connection conn = getConnection();
      int row = 0;

      for(int i = 0; i < sql.size(); i++)
      {
        ChadoTransaction tsn = (ChadoTransaction) sql.get(i);
        String[] sql_array = tsn.getSqlQuery(schema);

        for(int j = 0; j < sql_array.length; j++)
        {
          System.out.println(sql_array[j]);

          Statement st = conn.createStatement();
          row += st.executeUpdate(sql_array[j]);
        }
      }

      conn.close();
    }
    catch (java.sql.SQLException sqlExp)
    {
      sqlExp.printStackTrace();
    }
    catch (java.net.ConnectException conn)
    {
      JOptionPane.showMessageDialog(null, "Problems connecting...",
                                    "Database Connection Error - Check Server",
                                    JOptionPane.ERROR_MESSAGE);
      conn.printStackTrace();
    }

  }

  public static void main(String args[])
  {
    try
    {
      DbSqlConfig.init(new JPasswordField());
      SqlMapClient sqlMap = DbSqlConfig.getSqlMapInstance();

      Feature feature = new Feature();
      feature.setId(Integer.parseInt(args[0]));
      feature.setSchema(args[1]);

      List featureList = sqlMap.queryForList("getGffLine", feature);
 
      for(int i = 0; i < featureList.size(); i++)
      {
        feature = (Feature)featureList.get(i);
        int fmin     = feature.getFmin() + 1;
        int fmax     = feature.getFmax();

        System.out.print(fmin+" "+fmax);
        System.out.print(" "+feature.getType_id());
        System.out.print(" "+feature.getProp_type_id());
        System.out.print(" "+feature.getStrand());
        System.out.print(" "+feature.getUniquename());
        System.out.print(" "+feature.getTimelastmodified().toString());
        System.out.println(" "+Integer.toString(feature.getId()));
      }
    }
    catch(SQLException sqle)
    {
      sqle.printStackTrace();
    }
  }
}
