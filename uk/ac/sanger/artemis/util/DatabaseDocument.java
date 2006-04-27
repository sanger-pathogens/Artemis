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
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.chado.*;


import java.sql.*;
import java.text.SimpleDateFormat;
import java.io.*;
import java.net.ConnectException;
import java.util.Hashtable;
import java.util.Vector;
import java.util.Enumeration;
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

  /** JDBC DAO */
  private JdbcDAO jdbcDAO = null;

  /** iBatis DAO */
  private IBatisDAO connIB = null;

  private ByteBuffer[] gff_buffer;

  private ByteBuffer gff_buff;

  /** entries to split into */
  private String[] types = { "exon", "gene", "CDS", "transcript" };

  /** true if splitting the GFF into entries */
  private boolean splitGFFEntry;

  private boolean iBatis = false;

  private JPasswordField pfield;

  private List schema_list;
  
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
   * Append a String to the Document location.
   * @param name  the name to append.
   */
  public Document append(String name) throws IOException
  {
    return new DatabaseDocument( ((String)getLocation()) + name, pfield);
  }

  /**
   * Return the name of this Document (the last element of the Document
   * location).
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
  *  Set the name of this document.
  */
  public void setName(String name)
  {
    this.name = name;
  }

  /**
   * Return a Document with the last element stripped off.
   */
  public Document getParent()
  {
    return null;
  }

  /**
   * Return true if and only if the Document refered to by this object exists
   * and is readable. Always returns true.
   */
  public boolean readable()
  {
    return true;
  }

  /**
   * Return true if and only if the Document refered to by this object exists
   * and can be written to. Always returns false.
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
      ChadoDAO dao = getDAO();
      gff_buffer = getGff(dao, feature_id);

      ByteBuffer entry = new ByteBuffer();
      if(splitGFFEntry)
      {
        if(gff_buffer[0].size() > 0)
          entry.append(gff_buffer[0]);

        getSequence(dao, entry);
      }
      else
      {
        for(int i = 0; i < gff_buffer.length; i++)
        {
          if(gff_buffer[i].size() > 0)
            entry.append(gff_buffer[i]);
        }

        getSequence(dao, entry);
      }

      instream = new ByteArrayInputStream(entry.getBytes());
      return instream;
    }
    catch(java.sql.SQLException sqlExp)
    {
      JOptionPane.showMessageDialog(null, "Problems Reading...\n" +
          sqlExp.getMessage(),
          "Problems Reading From the Database ",
          JOptionPane.ERROR_MESSAGE);
      
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
   * Create an array of GFF lines.
   * @param dao                 the data access object 
   * @param parentFeatureID     the parent identifier for the features to 
   *                            extract
   * @return   the <code>ByteBuffer</code> array of GFF lines
   * @throws java.sql.SQLException
   */
  private ByteBuffer[] getGff(ChadoDAO dao, String parentFeatureID)
                       throws java.sql.SQLException
  {
    final int srcfeature_id = Integer.parseInt(parentFeatureID);
    List featList = dao.getGff(srcfeature_id, schema);

    ByteBuffer[] buffers = new ByteBuffer[types.length + 1];
    for(int i = 0; i < buffers.length; i++)
      buffers[i] = new ByteBuffer();

    String parentFeature = dao.getFeatureName(srcfeature_id, schema);
    ByteBuffer this_buff;

    int feature_size = featList.size();
    Hashtable id_store = new Hashtable(feature_size);

// build feature name store
    for(int i = 0; i < feature_size; i++)
    {
      ChadoFeature feat = (ChadoFeature)featList.get(i);
      String name       = feat.getUniquename();
      String feature_id = Integer.toString(feat.getId());

      id_store.put(feature_id, name);
    }

    // get all dbrefs
    Hashtable dbxrefs = dao.getDbxref(schema, null);
    
    // get all synonyms
    Hashtable synonym = dao.getAlias(schema, null);

    for(int i = 0; i < feature_size; i++)
    {
      ChadoFeature feat = (ChadoFeature)featList.get(i);
      int fmin                = feat.getFmin() + 1;
      int fmax                = feat.getFmax();
      long type_id            = feat.getType_id();
      int strand              = feat.getStrand();
      int phase               = feat.getPhase();
      String name             = feat.getUniquename();
      String typeName         = getCvtermName(type_id);

//      SimpleDateFormat date_format = new SimpleDateFormat("MM.dd.yyyy hh:mm:ss z");
//      String timelastmodified = date_format.format(feat.getTimelastmodified());
      String timelastmodified = Long.toString(feat.getTimelastmodified().getTime());
      String feature_id       = Integer.toString(feat.getId());

      String parent_id = feat.getObject_id();
      if(parent_id != null && id_store.containsKey(parent_id))
        parent_id = (String)id_store.get(parent_id);

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
      if(strand == -1)                        // strand
        this_buff.append("-\t");
      else if(strand == 1)
        this_buff.append("+\t");
      else
        this_buff.append(".\t");

      if(phase > 3)
        this_buff.append(".\t");               // phase
      else
        this_buff.append(phase+"\t"); 

      this_buff.append("ID=" + name + ";");

      if(parent_id != null)
        this_buff.append("Parent=" + parent_id + ";");

      this_buff.append("timelastmodified=" + timelastmodified + ";");


      // attributes
      Hashtable qualifiers     = feat.getQualifiers();
      if(qualifiers != null)
      {
        Enumeration e_qualifiers = qualifiers.keys();
        while(e_qualifiers.hasMoreElements())
        {
          Long qualifier_type_id = (Long)e_qualifiers.nextElement();
          String qualifier_name = getCvtermName(qualifier_type_id.longValue());
          if(qualifier_name == null)
            continue;
          Vector qualifier_value = (Vector)qualifiers.get(qualifier_type_id);
        
          for(int j=0; j<qualifier_value.size(); j++)
          {
            this_buff.append(qualifier_name+ "=" +
                             GFFStreamFeature.encode((String)qualifier_value.get(j))+";");
          }
        }
      } 

      // append dbxrefs
      if(dbxrefs != null &&
         dbxrefs.containsKey(new Integer(feature_id)))
      {
        this_buff.append("Dbxref=");
        Vector dbxref = (Vector)dbxrefs.get(new Integer(feature_id));
        for(int j=0; j<dbxref.size(); j++)
        {
          this_buff.append((String)dbxref.get(j));
          if(j<dbxref.size()-1)
            this_buff.append(",");
        }
      }
      
      if(synonym != null &&
         synonym.containsKey(new Integer(feature_id)))
      {
        this_buff.append(";");
        
        Alias alias;
        Vector v_synonyms = (Vector)synonym.get(new Integer(feature_id));
        for(int j=0; j<v_synonyms.size(); j++)
        {
          alias = (Alias)v_synonyms.get(j);
          this_buff.append(alias.getCvterm_name()+"=");
          this_buff.append(alias.getName());
          
          //this_buff.append((String)v_synonyms.get(j));
          if(j<v_synonyms.size()-1)
            this_buff.append(",");
        }
      }
      
      this_buff.append("\n");

      progress_listener.progressMade("Read from database: " + name);
    }

    return buffers;
  }

  /**
   * Look up the cvterm_id for a controlled vocabulary name.
   * @param name  
   * @return
   */
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
  }

  /**
   * Look up a cvterm name from the collection of cvterms.
   * @param id  a cvterm_id  
   * @return    the cvterm name
   */
  private String getCvtermName(long id)
  {
    if(cvterm == null)
    {
      try
      {
        getCvterm(getDAO());
      }
      catch(ConnectException ce)
      {
        ce.printStackTrace();
      }
      catch(SQLException sqle)
      {
        JOptionPane.showMessageDialog(null,
            "Problems Looking Up cvterm Name (cvterm_id="+
            Long.toString(id)+") ...\n" +
            sqle.getMessage(),
            "Cvterm Name Look Up",
            JOptionPane.ERROR_MESSAGE);
        sqle.printStackTrace();
      }
    }

    return (String)cvterm.get(new Long(id));
  }

  /**
   * Look up cvterms names and id and return in a hashtable.
   * @param dao the data access object
   * @return    the cvterm <code>Hashtable</code>
   */
  private Hashtable getCvterm(ChadoDAO dao)
  {
    cvterm = new Hashtable();

    try
    {
      List cvtem_list = dao.getCvterm();
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
   * Get the sequence for a feature.
   * @param dao   the data access object
   * @param buff  the buffer to add the sequence to
   * @return      the resulting buffer
   * @throws java.sql.SQLException
   */
  private ByteBuffer getSequence(ChadoDAO dao, ByteBuffer buff)
                     throws java.sql.SQLException
  {
    ChadoFeature feature = dao.getSequence(Integer.parseInt(feature_id),
                                         schema);

    buff.append("##FASTA\n>");
    buff.append(feature.getName());
    buff.append("\n");
    buff.append(feature.getResidues());
    return buff;
  }

  /**
   * Get the <code>List</code> of available schemas.
   * @return  the <code>List</code> of available schemas
   */
  public List getSchema()
  {
    return schema_list;
  }

  /**
   * Create a hashtable of the available entries with residues.
   * @return a <code>Hashtable</code> of the <code>String</code>
   *          representation (schema-type-feature_name) and the
   *          corresponding feature_id
   * @throws ConnectException
   * @throws java.sql.SQLException
   */
  public Hashtable getDatabaseEntries()
                   throws ConnectException, java.sql.SQLException
  {
    db = new Hashtable();
 
    ChadoDAO dao = null;
    
    try
    {
      dao = getDAO();
    }
    catch(ConnectException exp)
    {
      JOptionPane.showMessageDialog(null, "Connection Problems...\n"+
            exp.getMessage(), 
            "Connection Error",
            JOptionPane.ERROR_MESSAGE);
      throw exp;
    }
    catch(java.sql.SQLException sqlExp)
    {
      JOptionPane.showMessageDialog(null, "SQL Problems...\n"+
                                    sqlExp.getMessage(), 
                                    "SQL Error",
                                    JOptionPane.ERROR_MESSAGE);
      throw sqlExp;
    }
      
    try
    {
      schema_list = dao.getSchema();
      Iterator it      = schema_list.iterator();

      while(it.hasNext())
      {
        String schema = (String)it.next();
  
        List list = dao.getResidueType(schema);
         
        if(list.size() == 0)  // no residues for this organism
          continue;

        List list_residue_features = dao.getResidueFeatures(list, schema);
        Iterator it_residue_features = list_residue_features.iterator();
        while(it_residue_features.hasNext())
        {
          ChadoFeature feature = (ChadoFeature)it_residue_features.next();
          String typeName = getCvtermName(feature.getType_id());

          db.put(schema + " - " + typeName + " - " + feature.getName(),
                 Integer.toString(feature.getId()));
        }
      }
    }
    catch(java.sql.SQLException sqlExp)
    {
      JOptionPane.showMessageDialog(null, "SQL Problems...\n"+
                                    sqlExp.getMessage(), 
                                    "SQL Error",
                                    JOptionPane.ERROR_MESSAGE);
      sqlExp.printStackTrace();
    }
    return db;
  }
  
  
  /**
   * Get the data access object (DAO).
   * @return data access object
   */
  private ChadoDAO getDAO()
     throws java.net.ConnectException, SQLException
  {
    if(!iBatis)
    {
      if(jdbcDAO == null)
       jdbcDAO = new JdbcDAO((String)getLocation(), pfield); 
      return jdbcDAO;
    }
    else
    {
      if(connIB == null)
        connIB = new IBatisDAO(pfield);
      return connIB;
    }
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
    final File write_file = new File(System.getProperty("user.dir")+
                                     System.getProperty("file.separator")+
                                     getName());

    final FileOutputStream file_output_stream =
      new FileOutputStream(write_file);

    if(write_file.getName().endsWith(".gz")) 
    {
      // assume this file should be gzipped
      return new java.util.zip.GZIPOutputStream (file_output_stream);
    } 
    else 
      return file_output_stream;
  }

  /**
   * Commit the <code>ChadoTransaction</code> SQL back to the
   * database.
   * @param sql the collection of <code>ChadoTransaction</code> objects
   * @return
   */
  public int commit(Vector sql)
  {
    int i = 0;
    try
    {
      ChadoDAO dao = getDAO();
      boolean unchanged;
      
      //
      // check feature timestamps have not changed
      for(i = 0; i < sql.size(); i++)
      {
        ChadoTransaction tsn = (ChadoTransaction) sql.get(i);
 
        if(tsn.getType() != ChadoTransaction.INSERT_FEATURE ||
           tsn.getType() != ChadoTransaction.DELETE_FEATURE)
        {
          final List uniquename = tsn.getUniquename();
          
          for(int j=0; j<uniquename.size(); j++)
          {  
            System.out.println((String)uniquename.get(j)+" "+tsn.getLastModified());
            unchanged = checkFeatureTimestamp(schema, 
                           (String)uniquename.get(j), 
                         tsn.getLastModified(), dao);
            if(!unchanged)
              return 0;
          }
        }
      }  
       
      //
      // commit to database
      for(i = 0; i < sql.size(); i++)
      {
        ChadoTransaction tsn = (ChadoTransaction) sql.get(i);
        
        if(tsn.getType() == ChadoTransaction.UPDATE)
          dao.updateAttributes(schema, tsn);
        else if(tsn.getType() == ChadoTransaction.INSERT)
          dao.insertAttributes(schema, tsn);
        else if(tsn.getType() == ChadoTransaction.DELETE)
          dao.deleteAttributes(schema, tsn);
        else if(tsn.getType() == ChadoTransaction.INSERT_FEATURE)
          dao.insertFeature(schema, tsn, feature_id);
        else if(tsn.getType() == ChadoTransaction.DELETE_FEATURE)
          dao.deleteFeature(schema, tsn);
        else if(tsn.getType() == ChadoTransaction.DELETE_DBXREF)
          dao.deleteFeatureDbxref(schema, tsn);
        else if(tsn.getType() == ChadoTransaction.INSERT_DBXREF)
          dao.insertFeatureDbxref(schema, tsn);
      }
      
      
      //
      // update timelastmodified timestamp
      Timestamp ts = null;
      Timestamp ts2;
      for(int j = 0; j < sql.size(); j++)
      {
        ChadoTransaction tsn = (ChadoTransaction) sql.get(j);

        if(tsn.getType() != ChadoTransaction.INSERT_FEATURE ||
           tsn.getType() != ChadoTransaction.DELETE_FEATURE)
        {
          final List uniquename = tsn.getUniquename();
          
          // update timelastmodified timestamp
          for(int k=0; k<uniquename.size(); k++)
          {
            dao.writeTimeLastModified(schema, (String)uniquename.get(k), ts);
            ts2 = dao.getTimeLastModified(schema, (String)uniquename.get(k));
            if(ts2 == null)
              continue;
            
            if(ts == null)  
              ts = ts2;
            
            GFFStreamFeature gff_feature = (GFFStreamFeature)tsn.getFeatureObject();
            gff_feature.setLastModified(ts);
          }
        }
      }  
    }
    catch (java.sql.SQLException sqlExp)
    {
      JOptionPane.showMessageDialog(null, "Problems Writing...\n" +
                                    sqlExp.getMessage(),
                                    "Problems Writing to Database ",
                                    JOptionPane.ERROR_MESSAGE);
      sqlExp.printStackTrace();
    }
    catch (java.net.ConnectException conn_ex)
    {
      JOptionPane.showMessageDialog(null, "Problems connecting..."+
                                    conn_ex.getMessage(),
                                    "Database Connection Error - Check Server",
                                    JOptionPane.ERROR_MESSAGE);
      conn_ex.printStackTrace();
    }
      
    return i;
  }
  
  /**
   * Check the <code>Timestamp</code> on a feature (for versioning).
   * @param schema      the schema
   * @param uniquename  the feature uniquename
   * @param timestamp   the last read feature timestamp
   * @throws SQLException
   */
  public boolean checkFeatureTimestamp(final String schema,
                                      final String uniquename,
                                      final Timestamp timestamp,
                                      final ChadoDAO dao)
                                      throws SQLException
  {
    Timestamp now = dao.getTimeLastModified(schema, uniquename);
    
    if(now != null)
    {
      now.setNanos(0);
      timestamp.setNanos(0);
      
      if(now.compareTo(timestamp) != 0)
      {
        SimpleDateFormat date_format = 
                   new SimpleDateFormat("dd.MM.yyyy hh:mm:ss z");
        
        System.out.println(date_format.format(now)+"   "+date_format.format(timestamp));
        JOptionPane.showMessageDialog(null, 
                                      uniquename +
                                      " has been altered at :\n"+
                                      date_format.format(now),
                                      "Feature Changed",
                                      JOptionPane.ERROR_MESSAGE);
        return false;
      }
    }
    return true;
  }

  public static void main(String args[])
  {
    try
    {
      DbSqlConfig.init(new JPasswordField());
      SqlMapClient sqlMap = DbSqlConfig.getSqlMapInstance();

      ChadoFeature feature = new ChadoFeature();
      feature.setId(Integer.parseInt(args[0]));
      feature.setSchema(args[1]);

      List featureList = sqlMap.queryForList("getGffLine", feature);
 
      for(int i = 0; i < featureList.size(); i++)
      {
        feature = (ChadoFeature)featureList.get(i);
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
