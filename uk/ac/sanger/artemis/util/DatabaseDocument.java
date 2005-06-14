/* DatabaseDocument.java
 *
 * created: Fri Dec 18 1998
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 1998,1999,2000  Genome Research Limited
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

import java.sql.*;
import java.io.*;
import java.util.Hashtable;
import java.util.Vector;
import java.util.Date;
import javax.swing.JOptionPane;

/**
 *  Objects of this class are Documents created from a relational database.
 *
 **/

public class DatabaseDocument extends Document 
{
  private String name = null;
  private String feature_id = "1";
  private Hashtable cvterm;
  private InputStreamProgressListener progress_listener;
  private Hashtable db;
  private Vector organism;
  private String sqlLog = System.getProperty("user.home")+
                          System.getProperty("file.separator")+"art_sql.log";

  /**
   *
   *  Create a new Document from a database.
   *  @param location This should be a URL string giving:
   *         jdbc:postgresql://host:port/datbase_name?user=username
   *
   **/
  public DatabaseDocument(String location)
  {
    super(location);
  }

  /**
   *
   *  Create a new Document from a database.
   *  @param location This should be a URL string giving:
   *         jdbc:postgresql://host:port/datbase_name?user=username
   *  @param feature_id ID of a feature to be extracted.
   *
   **/
  public DatabaseDocument(String location, String feature_id)
  {
    super(location);
    this.feature_id = feature_id;
  }

  public DatabaseDocument(String location, String feature_id,
                          InputStreamProgressListener progress_listener)
  {
    super(location);
    this.feature_id = feature_id;
    this.progress_listener = progress_listener;
  }

  /**
   *
   *  Append a String to the Document location.
   *  @param name The name to append.
   *
   **/
  public Document append(String name) throws IOException 
  {
    return new DatabaseDocument(((String)getLocation()) + name);
  }

  /**
   *
   *  Return the name of this Document (the last element of the Document
   *  location).
   *
   **/
  public String getName() 
  {
    if(name == null)
    {
      int ind = ((String)getLocation()).indexOf("?");
      String name = ((String)getLocation()).substring(0,ind);
      ind = name.lastIndexOf("/");
      return name.substring(ind);
    }
    return name;
  }

  /**
   *
   *  Return a Document with the last element stripped off.
   *
   **/
  public String getFeatureName(String feature_id, Connection conn)
                  throws java.sql.SQLException
  { 
    Statement st = conn.createStatement();

    String sql = "SELECT name FROM feature WHERE feature_id= "+feature_id;
    appendToLogFile(sql,sqlLog);
    ResultSet rs = st.executeQuery(sql);
    rs.next();
    return rs.getString("name");
  }
  
  /**
   *
   *  Return a Document with the last element stripped off.
   *
   **/
  public Document getParent()
  {
    return null;
  }

  /**
   *
   *  Return true if and only if the Document refered to by this object exists
   *  and is readable.  Always returns true.
   *
   **/
  public boolean readable () 
  {
    return true;
  }


  /**
   *
   *  Return true if and only if the Document refered to by this object exists
   *  and can be written to.  Always returns false.
   *
   **/
  public boolean writable() 
  {
    return true;
  }


  /**
   *  Create a new InputStream object from this Document.  The contents of the
   *  Document can be read from the InputStream.
   *  @exception IOException Thrown if the Document can't be read from
   *    (for example if it doesn't exist).
   **/
  public InputStream getInputStream() throws IOException 
  {
    try
    {
      Connection conn = getConnection();
      System.out.println("Connected");

      String entry = getGFF(conn,feature_id) + getSequence(conn);
//    String entry = getSequence(conn);
      appendToLogFile(entry,sqlLog);
      ByteArrayInputStream instream = new ByteArrayInputStream(entry.getBytes());
      return instream;
    }
    catch(java.sql.SQLException sqlExp)
    {
      System.out.println("Problems connecting...");
      sqlExp.printStackTrace();
    }

    return null;
  }

  private String getGFF(Connection conn, String parentFeatureID) 
          throws java.sql.SQLException
  {
    Statement st = conn.createStatement();
    String sql = "SELECT strand, fmin, fmax, value, uniquename, feature.type_id, featureprop.type_id AS prop_type_id, strand"+
                 " FROM feature, featureloc, featureprop WHERE srcfeature_id = "+parentFeatureID+
                 " and featureloc.feature_id=featureprop.feature_id"+
                 " and featureloc.feature_id=feature.feature_id";
//                 " and feature.type_id=cvterm.cvterm_id"; // and cvterm.name='gene'";

    appendToLogFile(sql,sqlLog);
    ResultSet rs = st.executeQuery(sql);
    StringBuffer cdsBuffer = new StringBuffer();

    String parentFeature = getFeatureName(parentFeatureID,conn);
    int loop = 1;

    while(rs.next())
    {
      int fmin          = rs.getInt("fmin")+1;
      int fmax          = rs.getInt("fmax");
      long type_id      = rs.getLong("type_id");
      long prop_type_id = rs.getLong("prop_type_id");
      int strand      = rs.getInt("strand");
      String name     = rs.getString("uniquename");
      String typeName = getCvtermName(conn,type_id);
      String propTypeName = getCvtermName(conn,prop_type_id);
// make gff format

      cdsBuffer.append(parentFeature+"\t");    // seqid
      cdsBuffer.append("chado\t");             // source
      cdsBuffer.append(typeName+"\t");         // type
      cdsBuffer.append(fmin+"\t");             // start
      cdsBuffer.append(fmax+"\t");             // end
      cdsBuffer.append(".\t");                 // score
      if(strand == -1)                         // strand
        cdsBuffer.append("-\t");
      else
        cdsBuffer.append("+\t");
      cdsBuffer.append(".\t");                 // phase
      cdsBuffer.append("ID="+name+";"+propTypeName+"="+rs.getString("value")+"\n"); // attributes
      progress_listener.progressMade("Read from database: "+name);
    }

    return cdsBuffer.toString();
  }

  private String getCvtermName(Connection conn, long id)
  {
    if(cvterm == null)
      getCvterm(conn,null);
 
    return (String)cvterm.get(new Long(id));
  }

  private Hashtable getCvterm(Connection conn, String cv_name)
  {
    String sql = "SELECT cvterm.cvterm_id, cvterm.name " +
                 "FROM cvterm, cv WHERE cv.cv_id = cvterm.cv_id";

    if(cv_name != null)
      sql = sql + " AND cv.name='"+cv_name+"'";

    appendToLogFile(sql,sqlLog);

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
          System.err.println(this.getClass() + ": WARNING - read multiple CvTerms with name = '" 
                             + name  + "'");
        
        cvterm.put(new Long(id),name);
      }
    } 
    catch (SQLException sqle) 
    {
      System.err.println(this.getClass() + ": SQLException retrieving CvTerms");
      System.err.println(sqle);
    }

    return cvterm;
  }


  public String getSequence(Connection conn) throws java.sql.SQLException
  {
    Statement st = conn.createStatement();
    String sql = "SELECT name, residues from feature where feature_id = '"+
                                     feature_id+"'";

    appendToLogFile(sql,sqlLog);

    ResultSet rs = st.executeQuery(sql);

    rs.next();
    name = rs.getString("name");
    return "##FASTA\n>" + name + "\n" + rs.getString("residues");
  }


  /**
   *
   *  Create a hashtable of the available entries.
   * 
   **/
  public Hashtable getDatabaseEntries()
  {
    db = new Hashtable();
    organism = new Vector();

    try
    {
      Connection conn = getConnection();
      System.out.println("Connected");

      Statement st = conn.createStatement();

      String sql = "select type_id from feature where residues notnull";
      appendToLogFile(sql,sqlLog);

      ResultSet rs = st.executeQuery(sql);
      rs.next();
      String cvterm_id = rs.getString("type_id");

//    String sql = "select cvterm.cvterm_id, cvterm.name FROM cvterm, cv "+
//                 "WHERE cv.cv_id = cvterm.cv_id and cvterm.name = 'chromosome'";
//    ResultSet rs = st.executeQuery(sql);
//    rs.next();
//    String cvterm_id = rs.getString("cvterm_id");

      sql = new String("SELECT abbreviation, name, feature_id FROM organism, feature WHERE type_id = '"+
                        cvterm_id+"' and organism.organism_id=feature.organism_id "+
                        "ORDER BY abbreviation, name");
      appendToLogFile(sql,sqlLog);

      rs = st.executeQuery(sql);
      while(rs.next())
      {
        String org = rs.getString("abbreviation");
        db.put(org+" - "+rs.getString("name"), rs.getString("feature_id"));
        if(!organism.contains(org))
          organism.add(org);
      }
    }
    catch(java.sql.SQLException sqlExp)
    {
      JOptionPane.showMessageDialog(null, "SQL Problems...",
                                "SQL Error",
                                JOptionPane.ERROR_MESSAGE);
      sqlExp.printStackTrace();
    }
    catch(java.net.ConnectException conn)
    {
      JOptionPane.showMessageDialog(null, "Problems connecting...",
                                "Database Connection Error - Check Server",
                                JOptionPane.ERROR_MESSAGE);
      conn.printStackTrace();    
    }

    return db;
  }

  public Vector getOrganism()
  {
    return organism;
  }

  public Connection getConnection() throws java.sql.SQLException,
                                           java.net.ConnectException
  {
    return DriverManager.getConnection((String)getLocation());
  }


  /**
   *  Create a new OutputStream object from this Document.  The contents of the
   *  Document can be written from the stream.
   *  @exception IOException Thrown if the Document can't be written.
   **/
  public OutputStream getOutputStream() throws IOException
  {
    System.out.println("DatabaseDocument - ReadOnlyException");
    throw new ReadOnlyException ("this Database Document can not be written to");
  }

  /**
  *
  * Appends a log entry to the log file
  * @param logEntry     entry to add to log file
  * @param logFileName  log file name
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
      bw.write(dat+":: "+logEntry);
      bw.newLine();
      bw.flush();
    }
    catch (Exception ioe)
    {
      System.out.println("Error writing to log file "+logFileName);
      ioe.printStackTrace();
    }
    finally                     // always close the file
    {
      if(bw != null)
      try
      {
        bw.close();
      }
      catch (IOException ioe2) {}
    }
  }

}

