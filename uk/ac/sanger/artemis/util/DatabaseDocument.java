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

import uk.ac.sanger.artemis.io.ChadoCanonicalGene;
import uk.ac.sanger.artemis.io.DocumentEntry;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.PartialSequence;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.io.ReadFormatException;

import uk.ac.sanger.artemis.chado.ArtemisUtils;
import uk.ac.sanger.artemis.chado.ChadoCvTermView;
import uk.ac.sanger.artemis.chado.IBatisDAO;
import uk.ac.sanger.artemis.chado.JdbcDAO;
import uk.ac.sanger.artemis.chado.GmodDAO;
import uk.ac.sanger.artemis.chado.ChadoTransaction;
import uk.ac.sanger.artemis.components.database.DatabaseEntrySource;
import uk.ac.sanger.artemis.components.Splash;

import org.gmod.schema.sequence.Feature;
import org.gmod.schema.sequence.FeatureProp;
import org.gmod.schema.sequence.FeatureLoc;
import org.gmod.schema.sequence.FeaturePub;
import org.gmod.schema.sequence.FeatureRelationship;
import org.gmod.schema.sequence.FeatureSynonym;
import org.gmod.schema.sequence.FeatureCvTerm;
import org.gmod.schema.sequence.FeatureCvTermProp;
import org.gmod.schema.sequence.FeatureCvTermDbXRef;
import org.gmod.schema.sequence.FeatureCvTermPub;
import org.gmod.schema.cv.Cv;
import org.gmod.schema.cv.CvTerm;
import org.gmod.schema.general.Db;
import org.gmod.schema.general.DbXRef;
import org.gmod.schema.organism.Organism;
import org.gmod.schema.pub.PubDbXRef;
import org.gmod.schema.pub.Pub;

import java.sql.*;
import java.text.SimpleDateFormat;
import java.io.*;
import java.net.ConnectException;
import java.net.InetAddress;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.HashMap;
import java.util.Set;
import java.util.Vector;
import java.util.Enumeration;
import java.util.List;
import java.util.Iterator;
import java.util.Collection;

import javax.swing.JOptionPane;
import javax.swing.JPasswordField;

/**
 * Objects of this class are Documents created from a relational database.
 * 
 */

public class DatabaseDocument extends Document
{
  private String name = null;

  /** source feature_id */
  private String srcFeatureId = "1";

  /** database schema */
  private String schema = "public";

  private static Hashtable cvterms;
  
  private InputStreamProgressListener progress_listener;

  /** JDBC DAO */
  private JdbcDAO jdbcDAO = null;

  /** iBatis DAO */
  private IBatisDAO connIB = null;

  private ByteBuffer[] gff_buffer;

  private ByteBuffer gff_buff;

  /** entries to split into */
  private String[] types = { "repeat_region", "transcript" };

  /** true if splitting the GFF into entries */
  private boolean splitGFFEntry;

  private boolean iBatis = false;

  private JPasswordField pfield;
  
  private boolean singleSchema = true;

  private List schema_list;
  
  private static List organismNames;
  
  private boolean gene_builder;
  
  // include children in reading from the database
  private boolean readChildren = true;
  
  // range to retrieve features for
  private Range range;
  
  private Feature geneFeature;
  
  private Hashtable idFeatureStore;
  
  private boolean lazyFeatureLoad = true;
  
  public static String EXONMODEL  = "exon-model";
  public static String TRANSCRIPT = "mRNA";
  
  /** list of controlled_curation CV names */
  private static Vector cvControledCuratioNames;
  
  // controlled vocabulary
  /** controlled_curation controlled vocabulary */
  public static String CONTROLLED_CURATION_TAG_CVNAME = 
                                 "CC_";
  /** product controlled vocabulary */
  public static String PRODUCTS_TAG_CVNAME = "genedb_products";
  public static String RILEY_TAG_CVNAME = "RILEY";
  
  private static org.apache.log4j.Logger logger4j = 
    org.apache.log4j.Logger.getLogger(DatabaseDocument.class);

  
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

    if(location.indexOf('=') > -1)
      this.schema = location.substring( location.indexOf('=')+ 1);
    
    if(System.getProperty("ibatis") != null)
    {
      iBatis = true;
      System.setProperty("chado", location);
    }
    initMDC(this);
  }

  
  /**
   * 
   * Create a new Document from a database.
   * 
   * @param location
   *          This should be a URL string giving:
   *          jdbc:postgresql://host:port/database_name?user=username
   * @param feature_id
   *          ID of a feature to be extracted.
   * 
   */
  public DatabaseDocument(String location, JPasswordField pfield,
                          String srcFeatureId, String schema)
  {
    super(location);
    this.pfield = pfield;

    this.srcFeatureId = srcFeatureId;
    this.schema = schema;

    if(System.getProperty("ibatis") != null)
    {
      iBatis = true;
      System.setProperty("chado", location);
    }
    initMDC(this);
  }

  /**
   * 
   * Create a new Document from a database.
   * 
   * @param location
   *          This should be a URL string giving:
   *          jdbc:postgresql://host:port/datbase_name?user=username
   * @param srcFeatureId
   *          ID of a feature to be extracted.
   * @param splitGFFEntry
   *          split into separate entries based on feature types.
   * @param progress_listener
   *          input stream progress listener
   * 
   */
  public DatabaseDocument(String location, JPasswordField pfield,
                          String srcFeatureId, String schema, boolean splitGFFEntry,
                          InputStreamProgressListener progress_listener)
  {
    super(location);
    this.pfield = pfield;
    this.srcFeatureId = srcFeatureId;
    this.schema = schema;
    this.splitGFFEntry = splitGFFEntry;
    this.progress_listener = progress_listener;
    if(System.getProperty("ibatis") != null)
    {
      iBatis = true;
      System.setProperty("chado", location);
    }
    
    reset(location, schema);
    initMDC(this);
  }

  /**
   * Used by the gene builder to read a database entry
   * for a single gene.
   * @param location
   * @param pfield
   * @param srcFeatureId
   * @param schema
   * @param gene_builder
   */
  public DatabaseDocument(String location, JPasswordField pfield,
          String srcFeatureId, String schema, boolean gene_builder)
  {
    super(location);
    this.pfield = pfield;
    this.srcFeatureId = srcFeatureId;
    this.schema = schema;
    this.gene_builder = gene_builder;

    if(System.getProperty("ibatis") != null)
    {
      iBatis = true;
      System.setProperty("chado", location);
    }
    initMDC(this);
  }
  
  public DatabaseDocument(String location, JPasswordField pfield,
                          String srcFeatureId, String schema,
                          ByteBuffer gff_buff, String name)
  {
    super(location);
    this.pfield = pfield;
    this.srcFeatureId = srcFeatureId;
    this.schema = schema;
    this.gff_buff = gff_buff;
    this.name = name;
    if(System.getProperty("ibatis") != null)
    {
      iBatis = true;
      System.setProperty("chado", location);
    }
    initMDC(this);
  }
  
  /**
   * Use another DatabaseDocument to make a new document.
   * @param originalDocument
   * @param srcFeatureId
   * @param schema
   * @param gene_builder
   * @param region_grab
   * @param progress_listener
   */
  public DatabaseDocument (final DatabaseDocument originalDocument,
             final String schema, final Feature geneFeature,
             final Range range,
             final InputStreamProgressListener progress_listener)
  {
    this((String)originalDocument.getLocation(), 
         originalDocument.getPfield(),
         "-1", schema, false);
    this.progress_listener = progress_listener;
    this.range = range;
    this.geneFeature = geneFeature;
  }
  
  /**
   * Use another DatabaseDocument to make a new document.
   * @param originalDocument
   * @param srcFeatureId
   * @param schema
   * @param gene_builder
   * @param region_grab
   * @param progress_listener
   */
  public DatabaseDocument (final DatabaseDocument originalDocument,
             final String srcFeatureId, 
             final String schema, 
             final boolean gene_builder,
             final InputStreamProgressListener progress_listener)
  {
    this((String)originalDocument.getLocation(), 
         originalDocument.getPfield(),
         srcFeatureId, schema, gene_builder);
    this.progress_listener = progress_listener;
  }
  
  
  public static void initMDC(final DatabaseDocument doc)
  {
    // add username & host to MDC data for logging
    try
  	{
      org.apache.log4j.MDC.put("username",doc.getUserName());
  	}
  	catch(NullPointerException npe)
	  {
	    org.apache.log4j.MDC.put("username",System.getProperty("user.name"));
	  }
	
	  try 
	  {
	    org.apache.log4j.MDC.put("host",
    		  InetAddress.getLocalHost().getHostAddress());
	  } 
	  catch(Exception e) {}
  }
  
  public void setReadChildren(final boolean readChildren)
  {
    this.readChildren = readChildren; 
  }
  
  /**
   * Reset the schema.
   * @param location
   * @param schema
   */
  private void reset(String location, String schema)
  {
    this.schema = schema;

    if(!location.endsWith("="+schema))
    {
      int index = location.lastIndexOf('=');
      setLocation(location.substring(0,index+1) + schema);
      if(iBatis && connIB != null)
      {
        try
        {
          connIB.close();
        }
        catch(SQLException e)
        {
          logger4j.warn(e.getMessage());
        }
        connIB  = null;
      }
      
      jdbcDAO = null;
      System.setProperty("chado", (String)getLocation());
      logger4j.debug((String)getLocation());
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


  public DatabaseDocument createDatabaseDocument()
  {
    return new DatabaseDocument( (String)getLocation(), pfield,
                                  srcFeatureId, schema );
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
      GmodDAO dao = getDAO();

      if(gene_builder)
      {
        // creating a gene builder
        List schemaList = new Vector();
        schemaList.add(schema);
          
        ByteBuffer bb = getGeneFeature(srcFeatureId,
                                       schemaList, dao, readChildren);          
        return new ByteArrayInputStream(bb.getBytes());
      }
      else if(range != null)
      {
        //
        // Retrieve all features within a range
       // List schemaList = new Vector();
       // schemaList.add(schema);
        final Feature srcFeature;
        if(geneFeature != null)
        {
          Collection featureLocs = geneFeature.getFeatureLocsForFeatureId();
          Iterator it = featureLocs.iterator();
          final FeatureLoc featureLoc = (FeatureLoc)it.next();

          int srcfeatureid = featureLoc.getFeatureBySrcFeatureId().getFeatureId();
          srcFeature = dao.getFeatureById(srcfeatureid);
          setName(srcFeature.getUniqueName());
          this.srcFeatureId = Integer.toString(srcfeatureid);     
        }
        else
        {
          srcFeature = dao.getFeatureById(Integer.parseInt(srcFeatureId));
        }
        
        final ByteBuffer entryBuffer = getFeaturesInRange(srcFeature, range, dao);
        getChadoSequence(srcFeature, entryBuffer);

        return new ByteArrayInputStream(entryBuffer.getBytes());
      }
      
      ByteBuffer entryBuffer = new ByteBuffer();
      try
      {
        ByteBuffer sequenceBuffer = new ByteBuffer();
        if(dao instanceof IBatisDAO)
          ((IBatisDAO) dao).startTransaction();

        logger4j.debug("RETRIEVE SOURCE FEATURE FROM: "+getLocation());
        Feature srcFeature = getChadoSequence(dao, sequenceBuffer);
        gff_buffer = getGff(dao, srcFeature);
        
        if(splitGFFEntry)
        {
          if(gff_buffer[0].size() > 0)
            entryBuffer.append(gff_buffer[0]);
        }
        else
        {
          for(int i = 0; i < gff_buffer.length; i++)
          {
            if(gff_buffer[i].size() > 0)
              entryBuffer.append(gff_buffer[i]);
          }
        }

        if(dao instanceof IBatisDAO)
          ((IBatisDAO) dao).commitTransaction();
        entryBuffer.append(sequenceBuffer);
      }
      finally
      {
        if(dao instanceof IBatisDAO)
          ((IBatisDAO) dao).endTransaction();
      }

      instream = new ByteArrayInputStream(entryBuffer.getBytes());
      return instream;
    }
    catch(RuntimeException re)
    {
      JOptionPane.showMessageDialog(null, "Problems Reading...\n" +
          re.getMessage(),
          "Problems Reading From the Database ",
          JOptionPane.ERROR_MESSAGE);
      
      re.printStackTrace();
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

      String name = types[i-1];

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
   */
  private ByteBuffer[] getGff(final GmodDAO dao, 
                              final Feature srcFeature)
  {
    //final int srcfeature_id = Integer.parseInt(srcFeatureId);
    
    logger4j.debug("BUILD GFF FEATURES");
    
    // build srcfeature object
    FeatureLoc featureloc = new FeatureLoc();
    featureloc.setFeatureBySrcFeatureId(srcFeature);
    Feature child = new Feature();
    child.setFeatureLoc(featureloc);
    
    final List featList = dao.getFeaturesByLocatedOnFeature(child);
    final ByteBuffer[] buffers = new ByteBuffer[types.length + 1];
    for(int i = 0; i < buffers.length; i++)
      buffers[i] = new ByteBuffer();
    
    ByteBuffer this_buff;

    int feature_size = featList.size();
    final Hashtable id_store = new Hashtable(feature_size);

    // build feature store
    for(int i = 0; i < feature_size; i++)
    {
      Feature feat = (Feature)featList.get(i);
      id_store.put(Integer.toString(feat.getFeatureId()), feat);
    }
    
    if(lazyFeatureLoad)
      idFeatureStore = id_store;
    
    // get all dbrefs & synonyms etc
    final Hashtable dbxrefs;
    final Hashtable synonym;
    final Hashtable featureCvTerms;
    final Hashtable featureCvTermDbXRefs;
    Hashtable featureCvTermPubs = null;
    final Hashtable featurePubs;
    final List pubDbXRefs;
    
    if(lazyFeatureLoad)
    {
      dbxrefs = null;
      synonym = null;
      featureCvTerms = null;
      featureCvTermDbXRefs = null;
      featureCvTermPubs = null;
      featurePubs = null;
      pubDbXRefs = null;
    }
    else
    {
      dbxrefs= IBatisDAO.mergeDbXRef(
        dao.getFeatureDbXRefsBySrcFeature(srcFeature));
      synonym = getAllFeatureSynonyms(
        dao.getFeatureSynonymsBySrcFeature(srcFeature));
      featureCvTerms = getFeatureCvTermsByFeature(dao, 
        dao.getFeatureCvTermsBySrcFeature(srcFeature));
      featureCvTermDbXRefs = getFeatureCvTermDbXRef(dao, 
        dao.getFeatureCvTermDbXRefBySrcFeature(srcFeature));
      
      try
      {
        featureCvTermPubs = getFeatureCvTermPub(dao, 
          dao.getFeatureCvTermPubBySrcFeature(srcFeature));
      } 
      catch(Exception e) 
      { 
        e.printStackTrace();
        if(dao instanceof IBatisDAO)
        {
          try
          {
            ((IBatisDAO) dao).endTransaction();
            ((IBatisDAO) dao).startTransaction();
          }
          catch(SQLException sqle){}  
        }
      }
      featurePubs = getFeaturePubs(dao,
          dao.getFeaturePubsBySrcFeature(srcFeature));
      pubDbXRefs= dao.getPubDbXRef();
    }
    
    // create gff byte stream
    for(int i = 0; i < feature_size; i++)
    { 
      // select buffer based on feature type
      Feature feat = (Feature)featList.get(i);
      int type_id = feat.getCvTerm().getCvTermId();
      String typeName = getCvtermName(type_id, dao, gene_builder);
      this_buff = buffers[0];
      for(int j = 0; j < types.length; j++)
      {
        if(types[j].equals(typeName))
          this_buff = buffers[j+1];
      }

      chadoToGFF(feat, srcFeature.getUniqueName(),
                 dbxrefs, synonym, featureCvTerms,
                 pubDbXRefs, featureCvTermDbXRefs, featureCvTermPubs,
                 featurePubs,
                 id_store, dao, 
                 feat.getFeatureLoc(), this_buff, gene_builder);
       
      if( i%10 == 0 || i == feature_size-1)
        progress_listener.progressMade("Read from database: " + 
                                       feat.getUniqueName());
    }

    return buffers;
  }

  /**
   * Get a <code>Hashtable</code> of feature_id keys and their corresponding 
   * feature_synonym
   * 
   */
  private Hashtable getAllFeatureSynonyms(final List list) 
  {   
    Hashtable synonym = new Hashtable();
    Integer featureId;
    List value;
    FeatureSynonym alias;
    
    for(int i=0; i<list.size(); i++)
    {
      alias = (FeatureSynonym)list.get(i);
      featureId = new Integer(alias.getFeature().getFeatureId());
      if(synonym.containsKey(featureId))
        value = (Vector)synonym.get(featureId);
      else
        value = new Vector();
      
      value.add(alias);
      synonym.put(featureId, value);
    }
    
    return synonym;
  }
  
  /**
   * Get FeaturePub's (i.e. /literature qualifiers).
   * @param dao
   * @param list
   * @return
   */
  private Hashtable getFeaturePubs(final GmodDAO dao,
                                   final List list)
  {
    final Hashtable featurePubs = new Hashtable();
    Integer featureId;
    List value;
    FeaturePub featurePub;
    
    for(int i=0; i<list.size(); i++)
    {
      featurePub = (FeaturePub)list.get(i);
      featureId = new Integer(featurePub.getFeature().getFeatureId());
      if(featurePubs.containsKey(featureId))
        value = (Vector)featurePubs.get(featureId);
      else
        value = new Vector();
      
      value.add(featurePub);
      featurePubs.put(featureId, value);
    }
    
    return featurePubs;
  }
  
  /**
   * 
   * @param dao
   * @param chadoFeature null if we want them all
   * @return
   */
  private Hashtable getFeatureCvTermsByFeature(final GmodDAO dao, 
                                               final List list)
  {
    Hashtable featureCvTerms = new Hashtable();
    Integer featureId;
    List value;
    FeatureCvTerm feature_cvterm;
    
    for(int i=0; i<list.size(); i++)
    {
      feature_cvterm = (FeatureCvTerm)list.get(i);
      featureId = new Integer(feature_cvterm.getFeature().getFeatureId());
      if(featureCvTerms.containsKey(featureId))
        value = (Vector)featureCvTerms.get(featureId);
      else
        value = new Vector();
      
      value.add(feature_cvterm);
      featureCvTerms.put(featureId, value);
    }
    
    return featureCvTerms;
  }
  
  /**
   * 
   * @param dao
   * @param chadoFeature null if we want all
   * @return
   */
  private Hashtable getFeatureCvTermDbXRef(final GmodDAO dao, final List list)
  {
    if(list == null || list.size() == 0)
      return null;
    
    Integer featureCvTermDbXRefId;
    List value;
    
    Hashtable featureCvTermDbXRefs = new Hashtable(list.size());
    for(int i=0; i<list.size(); i++)
    {
      FeatureCvTermDbXRef featureCvTermDbXRef =
        (FeatureCvTermDbXRef)list.get(i);
      
      featureCvTermDbXRefId = new Integer(
          featureCvTermDbXRef.getFeatureCvTerm().getFeatureCvTermId());
      
      if(featureCvTermDbXRefs.containsKey(featureCvTermDbXRefId))
        value = (Vector)featureCvTermDbXRefs.get(featureCvTermDbXRefId);
      else
        value = new Vector();
      
      value.add(featureCvTermDbXRef);
      featureCvTermDbXRefs.put(featureCvTermDbXRefId, value);
    }
     
    return featureCvTermDbXRefs;
  }
  
  private Hashtable getFeatureCvTermPub(final GmodDAO dao,
                                        final List list)
  {
    if(list == null || list.size() == 0)
      return null;
    
    Integer featureCvTermId;
    List value;
    
    Hashtable featureCvTermPubs = new Hashtable(list.size());
    for(int i=0; i<list.size(); i++)
    {
      FeatureCvTermPub featureCvTermPub =
        (FeatureCvTermPub)list.get(i);
      
      featureCvTermId = new Integer(
          featureCvTermPub.getFeatureCvTerm().getFeatureCvTermId());
      
      if(featureCvTermPubs.containsKey(featureCvTermId))
        value = (Vector)featureCvTermPubs.get(featureCvTermId);
      else
        value = new Vector();
      
      value.add(featureCvTermPub);
      featureCvTermPubs.put(featureCvTermId, value);
    }
     
    return featureCvTermPubs;
  }
  
  /**
   * Retrieve the features in a given range
   * @param srcFeature
   * @param range
   * @param dao
   * @return
   */
  private ByteBuffer getFeaturesInRange(final Feature srcFeature,
                                        final Range range, 
                                        final GmodDAO dao)
  { 
    ByteBuffer buff = new ByteBuffer();

    logger4j.debug("GET FEATURES IN RANGE:: "+range.toString());
    List featuresInRange = dao.getFeaturesByRange(range.getStart()-1, 
                                range.getEnd(), 0, srcFeature, null);
    
    List featureIds = new Vector(featuresInRange.size());
    for(int i=0; i<featuresInRange.size(); i++)
    {
      Feature thisFeature = (Feature)featuresInRange.get(i);
      featureIds.add(new Integer(thisFeature.getFeatureId()));
    }
    
    FeatureLoc featureLoc = new FeatureLoc();
    featureLoc.setFmin(new Integer(range.getStart()));
    featureLoc.setFmax(new Integer(range.getEnd()));
    srcFeature.setFeatureLoc(featureLoc);
    
    Hashtable dbxrefs = IBatisDAO.mergeDbXRef(
        dao.getFeatureDbXRefsBySrcFeature(srcFeature));
    Hashtable synonym = getAllFeatureSynonyms(
        dao.getFeatureSynonymsBySrcFeature(srcFeature));
    Hashtable featureCvTerms = getFeatureCvTermsByFeature(dao, 
        dao.getFeatureCvTermsBySrcFeature(srcFeature));
    Hashtable featureCvTermDbXRefs = getFeatureCvTermDbXRef(dao, 
        dao.getFeatureCvTermDbXRefBySrcFeature(srcFeature));
    Hashtable featureCvTermPubs = getFeatureCvTermPub(dao, 
        dao.getFeatureCvTermPubBySrcFeature(srcFeature));
    Hashtable featurePubs = getFeaturePubs(dao,
        dao.getFeaturePubsBySrcFeature(srcFeature));

    List pubDbXRefs = dao.getPubDbXRef();
    
    Hashtable id_store = new Hashtable(featuresInRange.size());

    // build feature name store
    for(int i = 0; i < featuresInRange.size(); i++)
    {
      Feature chadoFeature = (Feature)featuresInRange.get(i);
      String featureId = Integer.toString(chadoFeature.getFeatureId());
      id_store.put(featureId, chadoFeature);
    }
    
    for(int i=0; i<featuresInRange.size(); i++)
    {
      Feature chadoFeature = (Feature)featuresInRange.get(i);
      id_store.put(Integer.toString(chadoFeature.getFeatureId()), chadoFeature);

      chadoToGFF(chadoFeature, srcFeature.getUniqueName(), dbxrefs, synonym, featureCvTerms,
          pubDbXRefs, featureCvTermDbXRefs, featureCvTermPubs, featurePubs,
          id_store, dao, chadoFeature.getFeatureLoc(), buff, gene_builder);
      if( i%10 == 0 || i == featuresInRange.size()-1)
        progress_listener.progressMade("Read from database: " + 
                                       chadoFeature.getUniqueName());
    }
    return buff;
  }
  
  /**
   * Use by the gene editor to retrieve the gene and related
   * features
   * @param search_gene     gene uniquename
   * @param schema_search   schema list to search
   * @param dao             data access method
   * @return  GFF byte buffer
   * @throws SQLException
   * @throws ReadFormatException
   * @throws ConnectException 
   */
  private ByteBuffer getGeneFeature(final String search_gene, 
                                    final List schema_search,
                                    GmodDAO dao, 
                                    final boolean readChildren) 
          throws SQLException, ReadFormatException, ConnectException
  {
    CvTermThread cvThread = null;
    if(DatabaseDocument.cvterms == null)
    {
      cvThread = new CvTermThread(dao);
      cvThread.start();
    }
    
    final Hashtable id_store = new Hashtable();

    boolean singleSchema = true;
    final List pg_schemas = dao.getSchema(); 
    Iterator schemasIt = pg_schemas.iterator();
    while(schemasIt.hasNext())
    {
      String thisSchema = (String)schemasIt.next();
      
      if( thisSchema.equalsIgnoreCase(schema) )
      {
        singleSchema = false;
        break;
      }
    }
    if(singleSchema)
      logger4j.debug("SINGLE SCHEMA");
    else
      reset((String)getLocation(), (String)schema_search.get(0));
    dao = getDAO();
    Feature chadoFeature = 
      (Feature)(dao.getFeaturesByUniqueName(search_gene).get(0));
    
    ChadoCanonicalGene chado_gene = new ChadoCanonicalGene();
    id_store.put(Integer.toString(chadoFeature.getFeatureId()), 
                 chadoFeature);

    List featurelocs = new Vector(chadoFeature.getFeatureLocsForFeatureId());
    FeatureLoc featureloc = (FeatureLoc) featurelocs.get(0);
    int src_id = featureloc.getSrcFeatureId();
    srcFeatureId = Integer.toString(src_id);

    Feature parent = new Feature();
    parent.setFeatureId(src_id);

    logger4j.debug("GET PARENT FEATURE");
    parent = dao.getLazyFeatureNoResiduesById(new Integer(src_id));
    
    chado_gene.setSeqlen(parent.getSeqLen());
    chado_gene.setSrcfeature_id(src_id);

    final ByteBuffer buff = new ByteBuffer();
    
    logger4j.debug("BUILD GENE GFF LINE");
    buildGffLineFromId(dao, chadoFeature.getFeatureId(), 
        id_store, parent.getUniqueName(), src_id, buff, chadoFeature);
    

    if(!readChildren)
    {
      logger4j.debug( new String(buff.getBytes()) );
      return buff;
    }
    
    // get children of gene
    List relations = new Vector(chadoFeature.getFeatureRelationshipsForObjectId());
    Set idsSeen = new HashSet();
    for(int i = 0; i < relations.size(); i++)
    {
      //Feature transcript = new Feature();
      int id = ((FeatureRelationship) relations.get(i)).getFeatureBySubjectId().getFeatureId();
      Integer idInt = new Integer(id);
      if(idsSeen.contains(idInt))
        continue;
      idsSeen.add(idInt);
      Feature transcript = buildGffLineFromId(dao, id, id_store, parent.getUniqueName(), 
                                              src_id, buff, null);

      if( transcript == null || transcript.getCvTerm() == null ||
          transcript.getCvTerm().getName() == null || 
         (transcript.getCvTerm().getName().indexOf("RNA") < 0 &&
          transcript.getCvTerm().getName().indexOf("transcript") < 0 ) )
        continue;
      // get children of transcript - exons and pp
      logger4j.debug("GET CHILDREN OF "+transcript.getName());
      List transcipt_relations = new Vector(
          transcript.getFeatureRelationshipsForObjectId());

      for(int j = 0; j < transcipt_relations.size(); j++)
      {
        id = ((FeatureRelationship) transcipt_relations.get(j)).getFeatureBySubjectId().getFeatureId();

        buildGffLineFromId(dao, id, id_store, parent.getUniqueName(), 
                           src_id, buff, null);
      }
    }

    logger4j.debug( "GFF:\n"+new String(buff.getBytes()) );

    // now wait for cvterm to be loaded
    if(cvThread != null)
    {
      while(cvThread.isAlive())
        try
        {
          Thread.sleep(10);
        }
        catch(InterruptedException e)
        {
          // TODO Auto-generated catch block
          e.printStackTrace();
        }
    }
    
    return buff;
  }
  
  
  /**
   * 
   * @param dao
   * @param featureId
   * @param id_store
   * @param parentName
   * @param srcFeatureId
   * @param this_buff
   * @param chadoFeature
   * @return
   */
  private Feature buildGffLineFromId(final GmodDAO dao, 
                                  final int featureId, 
                                  final Hashtable id_store,
                                  final String parentName,
                                  final int srcFeatureId,
                                  final ByteBuffer this_buff,
                                  Feature chadoFeature)
  {
    if(chadoFeature == null)
      chadoFeature = (Feature)dao.getFeatureById(featureId); 

    id_store.put(Integer.toString(chadoFeature.getFeatureId()), 
                 chadoFeature);

    final FeatureLoc loc = getFeatureLoc(new Vector(
        chadoFeature.getFeatureLocsForFeatureId()), srcFeatureId);
    
    if(loc == null)
    {
      logger4j.debug("FEATURELOC NOT FOUND :: "+chadoFeature.getUniqueName());
      return null;
    }
    final Hashtable dbxrefs = IBatisDAO.mergeDbXRef(
        dao.getFeatureDbXRefsByFeatureUniquename(chadoFeature.getUniqueName()));
    
    final Hashtable synonym = getAllFeatureSynonyms( 
        dao.getFeatureSynonymsByFeatureUniquename(chadoFeature.getUniqueName()));
    
    final Hashtable featureCvTerms = getFeatureCvTermsByFeature(dao, 
                                  dao.getFeatureCvTermsByFeature(chadoFeature));
    
    final Hashtable featureCvTermDbXRefs = getFeatureCvTermDbXRef(dao, 
                             dao.getFeatureCvTermDbXRefByFeature(chadoFeature));
    
    Hashtable featureCvTermPubs = null;
    
    try
    {
      featureCvTermPubs = getFeatureCvTermPub(dao,
                          dao.getFeatureCvTermPubByFeature(chadoFeature));
    }
    catch(RuntimeException re){re.printStackTrace();}

    final Hashtable featurePubs = getFeaturePubs(dao,
        dao.getFeaturePubsByFeature(chadoFeature));
    List pubDbXRefs= new Vector(); //dao.getPubDbXRef();
    chadoToGFF(chadoFeature, parentName, dbxrefs, synonym, featureCvTerms,
        pubDbXRefs, featureCvTermDbXRefs, featureCvTermPubs, featurePubs, 
        id_store, dao, loc, this_buff, gene_builder);  
    return chadoFeature;
  }
  
  /**
   * Convert the chado feature into a GFF line
   * @param feat           Chado feature
   * @param parentFeature  parent of this feature
   * @param dbxrefs        hashtable containing dbxrefs
   * @param synonym        hashtable containing synonynms
   * @param featureCvTerms
   * @param pubDbXRefs
   * @param featureCvTermDbXRefs
   * @param id_store       id store for looking up parent names
   * @param dao            chado data access
   * @param featureloc     feature location for this chado feature
   * @param this_buff      byte buffer of GFF line 
   */
  private static void chadoToGFF(final Feature feat,
                                 final String parentFeature,
                                 final Hashtable dbxrefs,
                                 final Hashtable synonym,
                                 final Hashtable featureCvTerms,
                                 final List pubDbXRefs,
                                 final Hashtable featureCvTermDbXRefs,
                                 final Hashtable featureCvTermPubs,
                                 final Hashtable featurePubs,
                                 final Hashtable id_store,
                                 final GmodDAO dao,
                                 final FeatureLoc featureloc,
                                 final ByteBuffer this_buff,
                                 final boolean gene_builder)
  {
    String gff_source = null;
    
    final int fmin          = featureloc.getFmin().intValue() + 1;
    final int fmax          = featureloc.getFmax().intValue();
    final int type_id       = feat.getCvTerm().getCvTermId();
    final Short strand      = featureloc.getStrand();
    final Integer phase     = featureloc.getPhase();
    final String name       = feat.getUniqueName();
    final String typeName   = getCvtermName(type_id, dao, gene_builder);
    final Integer featureId = new Integer(feat.getFeatureId());
    final String timelastmodified = Long.toString(feat.getTimeLastModified().getTime());

    String parent_id = null;
    String parent_relationship = null;
    int rank = -1;
/*    if(feat.getFeatureRelationship() != null)
    {
      FeatureRelationship feat_relationship = feat.getFeatureRelationship();
      parent_id = Integer.toString(feat_relationship.getFeatureByObjectId().getFeatureId());
      long parent_type_id = feat_relationship.getCvTerm().getCvTermId();
      
      parent_relationship = feat_relationship.getCvTerm().getName();
      
      rank= feat_relationship.getRank();
      if(parent_relationship == null)
        parent_relationship = getCvtermName(parent_type_id, dao);
    }
    else */
    
    ByteBuffer clusterOrthoParalog = null;
    if(feat.getFeatureRelationshipsForSubjectId() != null)
    {
      Collection relations = feat.getFeatureRelationshipsForSubjectId();
      Iterator it = relations.iterator();
      Set featureRelationshipIds = new HashSet();
      //Set duplicates = new HashSet();
      
      while(it.hasNext())
      {
        final FeatureRelationship fr = 
                            (FeatureRelationship)it.next();
        final Integer featureRelationShipId = new Integer( fr.getFeatureRelationshipId() );

        if(featureRelationshipIds.contains( featureRelationShipId ))
          continue;
        
        featureRelationshipIds.add(featureRelationShipId);
        final String cvTermName;
        if( fr.getCvTerm().getName() == null )
        {
          int parent_type_id = fr.getCvTerm().getCvTermId();
          cvTermName = getCvtermName(parent_type_id, dao, gene_builder);
        }
        else
          cvTermName = fr.getCvTerm().getName();
      
        if(cvTermName.equals("derives_from") || cvTermName.equals("part_of") ||
           cvTermName.equals("proper_part_of") || 
           cvTermName.equals("partof") || cvTermName.equals("producedby")) // flybase
        {
          parent_relationship = cvTermName;
          parent_id = Integer.toString(fr.getFeatureByObjectId().getFeatureId());
          rank      = fr.getRank();
        }
        else
        {
          if(clusterOrthoParalog == null)
            clusterOrthoParalog = new ByteBuffer();
          // ortholog/paralog/cluster data
          int orthologueFeature = fr.getFeatureByObjectId().getFeatureId();
          clusterOrthoParalog.append(cvTermName+"="+
              GFFStreamFeature.encode("object_id="+orthologueFeature+"; rank="+fr.getRank())+";");
        }
      }
    }

    // look up parent name
    if(parent_id != null && id_store != null &&  id_store.containsKey(parent_id))
      parent_id = ((Feature)id_store.get(parent_id)).getUniqueName();
 
    // make gff format
    
    Vector dbxref = null;
    // append dbxrefs
    if(dbxrefs != null &&
       dbxrefs.containsKey(featureId))
    {
      dbxref = (Vector)dbxrefs.get(featureId);
      for(int j=0; j<dbxref.size(); j++)
      {
        if(((String)dbxref.get(j)).startsWith("GFF_source:"))
        {
          gff_source = ((String)dbxref.get(j)).substring(11);
          dbxref.removeElementAt(j);
        }
      }
    }

    this_buff.append(parentFeature + "\t"); // seqid
    
    if(gff_source != null)
      this_buff.append(gff_source+"\t");    // source
    else
      this_buff.append("chado\t");
    
    if(typeName.equals("exon"))
      this_buff.append(EXONMODEL + "\t");   // type
    else
      this_buff.append(typeName + "\t");    // type
    this_buff.append(fmin + "\t");          // start
    this_buff.append(fmax + "\t");          // end
    this_buff.append(".\t");                // score
    if(strand.equals( new Short((short)-1)) )                        // strand
      this_buff.append("-\t");
    else if(strand.equals( new Short((short)1)) )
      this_buff.append("+\t");
    else
      this_buff.append(".\t");

    if(phase == null)
      this_buff.append(".\t");               // phase
    else
      this_buff.append(phase+"\t"); 

    this_buff.append("ID=" + name + ";");
    this_buff.append("feature_id=" + featureId.toString() + ";");
    
    if(feat.getName() != null)
      this_buff.append("Name=" + feat.getName() + ";");
   
    if(parent_id != null && !parent_id.equals("0"))
    {
      if(parent_relationship.equals("derives_from"))
        this_buff.append("Derives_from=" + parent_id + ";");
      else
        this_buff.append("Parent=" + parent_id + ";");
    }

    this_buff.append("timelastmodified=" + timelastmodified + ";");
    this_buff.append("isObsolete=" + Boolean.toString(feat.isObsolete()) + ";");
    
    // this is the chado feature_relationship.rank used
    // to order joined features e.g. exons
    if(rank > -1)
      this_buff.append("feature_relationship_rank="+rank+";"); 

    //this_buff.append("feature_id="+feature_id+";");
    
    // attributes
    if(feat.getFeatureProps() != null &&
       feat.getFeatureProps().size() > 0)
    {
      Collection featureprops = feat.getFeatureProps();
      Iterator it = featureprops.iterator();
      
      while(it.hasNext())
      {
        FeatureProp featprop = (FeatureProp)it.next();
        String qualifier_name = getCvtermName(featprop.getCvTerm().getCvTermId(), dao, gene_builder);
        if(qualifier_name == null)
          continue;
        if(featprop.getValue() != null)
          this_buff.append(GFFStreamFeature.encode(qualifier_name)+ "=" +
                           GFFStreamFeature.encode(featprop.getValue())+";");
      }
    }

    if(clusterOrthoParalog != null)
      this_buff.append(clusterOrthoParalog);
    
    // append dbxrefs
    boolean foundPrimaryDbXRef = false;
    if(feat.getDbXRef() != null)
    {
      this_buff.append("Dbxref=");
      this_buff.append(GFFStreamFeature.encode(
          feat.getDbXRef().getDb().getName()+":"+feat.getDbXRef().getAccession()));
      foundPrimaryDbXRef = true;
      if(dbxref == null || dbxref.size() == 0)
        this_buff.append(";");
    }
    
    if(dbxref != null && dbxref.size() > 0)
    {
      if(foundPrimaryDbXRef)
        this_buff.append(",");
      else
        this_buff.append("Dbxref=");
      for(int j=0; j<dbxref.size(); j++)
      {
        this_buff.append(GFFStreamFeature.encode((String)dbxref.get(j)));
        if(j<dbxref.size()-1)
          this_buff.append(",");
      }
      this_buff.append(";");
    }
    
    // append synonyms
    if(synonym != null &&
       synonym.containsKey(featureId))
    {   
      FeatureSynonym alias;
      Vector v_synonyms = (Vector)synonym.get(featureId);
      for(int j=0; j<v_synonyms.size(); j++)
      {
        alias = (FeatureSynonym)v_synonyms.get(j);
        
        this_buff.append( getCvtermName(alias.getSynonym().getCvTerm().getCvTermId(), dao, gene_builder) + "=" );
        //this_buff.append(alias.getSynonym().getCvterm().getName()+"=");
        this_buff.append(alias.getSynonym().getName());
        
        if(!alias.isCurrent())
          this_buff.append(GFFStreamFeature.encode(";current=false"));
        
        //if(j<v_synonyms.size()-1)
        this_buff.append(";");
      }
    }
    
    // /literature
    if(featurePubs != null &&
       featurePubs.containsKey(featureId))
    {
      FeaturePub featurePub;
      Vector v_featurePubs = (Vector)featurePubs.get(featureId);
      for(int j=0; j<v_featurePubs.size(); j++)
      {
        featurePub = (FeaturePub)v_featurePubs.get(j);
        this_buff.append( "literature=" );
        this_buff.append(featurePub.getPub().getUniqueName());
        this_buff.append(";");
      }
    }
    
    // GO, controlled_curation, product
    if(featureCvTerms != null && 
       featureCvTerms.containsKey(featureId))
    {
      FeatureCvTerm feature_cvterm;
      Vector v_feature_cvterms = (Vector)featureCvTerms.get(featureId);
      for(int j=0; j<v_feature_cvterms.size(); j++)
      {
        feature_cvterm = (FeatureCvTerm)v_feature_cvterms.get(j);
        
        Integer featureCvTermId = new Integer( feature_cvterm.getFeatureCvTermId() );
        
        List featureCvTermDbXRefList = null;
        
        if(featureCvTermDbXRefs != null)
          featureCvTermDbXRefList = (List)featureCvTermDbXRefs.get(featureCvTermId);
        
        List featureCvTermPubList = null;
        
        if(featureCvTermPubs != null)
          featureCvTermPubList = (List)featureCvTermPubs.get(featureCvTermId);
          
        appendControlledVocabulary(this_buff, dao, feature_cvterm,
                                   featureCvTermDbXRefList,featureCvTermPubList, pubDbXRefs, gene_builder);
      }
      //System.out.println(new String(this_buff.getBytes()));
    }
      
    this_buff.append("\n");
  }
  
  /**
   * Appends controlled vocabulary terms to the buffer
   * @param attr_buff
   * @param dao
   * @param feature_cvterm
   * @param featureCvTermDbXRef
   */
  public static void appendControlledVocabulary(final ByteBuffer attr_buff,
                                          final GmodDAO dao,
                                          final FeatureCvTerm feature_cvterm,
                                          final List featureCvTermDbXRefs,
                                          final List featureCvTermPubs,
                                          final List pubDbXRefs,
                                          final boolean gene_builder)
  {
    CvTerm cvterm =  getCvTerm( feature_cvterm.getCvTerm().getCvTermId(), dao, gene_builder);
    DbXRef dbXRef = feature_cvterm.getCvTerm().getDbXRef();
    
    if(cvterm.getCv().getName().startsWith(DatabaseDocument.CONTROLLED_CURATION_TAG_CVNAME))
    {
      attr_buff.append("controlled_curation=");
      
      attr_buff.append("term="+
          GFFStreamFeature.encode(feature_cvterm.getCvTerm().getName())+"%3B");
      attr_buff.append("cv="+
          GFFStreamFeature.encode(feature_cvterm.getCvTerm().getCv().getName())+"%3B");   
      
      // N.B. the db_xref may be a FeatureCvTermDbXRef or a Pub for /controlled_curation
      int nfound_dbxref = 0;
      if(feature_cvterm.getPub().getUniqueName() != null &&
         !feature_cvterm.getPub().getUniqueName().equalsIgnoreCase("NULL"))
      {
        // PMID
        Pub pub = feature_cvterm.getPub();
        
        // internal check
        checkPubDbXRef(pubDbXRefs, pub.getPubId(), pub, feature_cvterm);
        
        attr_buff.append("db_xref="+ pub.getUniqueName());
        nfound_dbxref++;
      }
      
      if(featureCvTermDbXRefs != null &&
              featureCvTermDbXRefs.size() > 0)
      {
        for(int i=0; i<featureCvTermDbXRefs.size(); i++)
        {
          FeatureCvTermDbXRef featureCvTermDbXRef =
            (FeatureCvTermDbXRef)featureCvTermDbXRefs.get(i);
    
          if(feature_cvterm.getFeatureCvTermId() != 
            featureCvTermDbXRef.getFeatureCvTerm().getFeatureCvTermId())
            continue;
      
          if(nfound_dbxref == 0)
            attr_buff.append("db_xref=");
          else if(nfound_dbxref > 0)
            attr_buff.append("|");
          
          DbXRef fc_dbXRef = featureCvTermDbXRef.getDbXRef();
          attr_buff.append(fc_dbXRef.getDb().getName()+":");
          attr_buff.append(fc_dbXRef.getAccession());
          nfound_dbxref++;
        }
      }
      
      if(nfound_dbxref > 0)
        attr_buff.append("%3B");
      
      List feature_cvtermprops = (List) feature_cvterm.getFeatureCvTermProps();
      for(int i = 0; i < feature_cvtermprops.size(); i++)
      {
        FeatureCvTermProp feature_cvtermprop = 
          (FeatureCvTermProp)feature_cvtermprops.get(i);
        attr_buff.append(getCvtermName(feature_cvtermprop.getCvTerm()
            .getCvTermId(), dao, gene_builder));
        attr_buff.append("=");
        attr_buff.append(GFFStreamFeature.encode(feature_cvtermprop.getValue()));
        if(i < feature_cvtermprops.size()-1)
          attr_buff.append("%3B");
      }
      
      attr_buff.append(";");
    }
    else if(cvterm.getCv().getName().equals(DatabaseDocument.PRODUCTS_TAG_CVNAME))
    {
      attr_buff.append("product=");
      attr_buff.append(GFFStreamFeature.encode(feature_cvterm.getCvTerm().getName())+";");
    }
    else if(cvterm.getCv().getName().equals(DatabaseDocument.RILEY_TAG_CVNAME))
    {
      // class include the cvTermId as a convenience for looking up the term
      attr_buff.append("class=");
      attr_buff.append(dbXRef.getAccession()+"::"+
                      feature_cvterm.getCvTerm().getCvTermId()+";");
    }
    else
    {
      attr_buff.append("GO=");

      if(cvterm.getCv().getName().equals("molecular_function"))
        attr_buff.append("aspect=F%3B");
      else if(cvterm.getCv().getName().equals("cellular_component"))
        attr_buff.append("aspect=C%3B");
      else if(cvterm.getCv().getName().equals("biological_process"))
        attr_buff.append("aspect=P%3B");

      if(feature_cvterm.isNot())
        attr_buff.append("qualifier=NOT%3B");

      attr_buff.append("GOid="+dbXRef.getDb().getName() + ":"
          + dbXRef.getAccession() + "%3B");
      
      attr_buff.append("term="+
          GFFStreamFeature.encode(feature_cvterm.getCvTerm().getName())+"%3B");
      
      // PMID
      int nfound_pub = 0;
      if(feature_cvterm.getPub() != null &&
         feature_cvterm.getPub().getUniqueName() != null &&
         !feature_cvterm.getPub().getUniqueName().equalsIgnoreCase("NULL"))
      {
        Pub pub = feature_cvterm.getPub();
        attr_buff.append("db_xref="+
            pub.getUniqueName());
        nfound_pub++;
      }
      
      if(featureCvTermPubs != null &&
          featureCvTermPubs.size() > 0)
      {
        for(int i=0; i<featureCvTermPubs.size(); i++)
        {
          FeatureCvTermPub featureCvTermPub =
            (FeatureCvTermPub)featureCvTermPubs.get(i);
          
          if(feature_cvterm.getFeatureCvTermId() != 
            featureCvTermPub.getFeatureCvTerm().getFeatureCvTermId())
            continue;
          
          if(nfound_pub == 0)
            attr_buff.append("db_xref=");
          else if(nfound_pub > 0)
            attr_buff.append("|");

          attr_buff.append(featureCvTermPub.getPub().getUniqueName());
          nfound_pub++;
        }
      }
      
      if(nfound_pub > 0)
        attr_buff.append("%3B");
      
      if(featureCvTermDbXRefs != null &&
          featureCvTermDbXRefs.size() > 0 )
      {  
        int nfound = 0;
        for(int i=0; i<featureCvTermDbXRefs.size(); i++)
        {
          FeatureCvTermDbXRef featureCvTermDbXRef =
            (FeatureCvTermDbXRef)featureCvTermDbXRefs.get(i);
          
          
          if(feature_cvterm.getFeatureCvTermId() != 
            featureCvTermDbXRef.getFeatureCvTerm().getFeatureCvTermId())
          {
            continue;
          }
          
          if(nfound == 0)
            attr_buff.append("with=");
          else if(nfound > 1)
            attr_buff.append("|");
          
          DbXRef fc_dbXRef = featureCvTermDbXRef.getDbXRef();
          attr_buff.append(fc_dbXRef.getDb().getName()+":");
          attr_buff.append(fc_dbXRef.getAccession());
          nfound++;
        }
        
        if(nfound > 0)
          attr_buff.append("%3B");
      }
      
      List feature_cvtermprops = (List)feature_cvterm
          .getFeatureCvTermProps();
      for(int i = 0; i < feature_cvtermprops.size(); i++)
      {
        FeatureCvTermProp feature_cvtermprop = 
          (FeatureCvTermProp)feature_cvtermprops.get(i);
        
        if(feature_cvtermprop.getValue() == null)
          continue;
        
        attr_buff.append(getCvtermName(feature_cvtermprop.getCvTerm()
            .getCvTermId(), dao, gene_builder));
        attr_buff.append("=");
        attr_buff.append(GFFStreamFeature.encode(feature_cvtermprop.getValue()));
        if(i < feature_cvtermprops.size()-1)
          attr_buff.append("%3B");
      }
      
      attr_buff.append(";");
    }
  }
  
  /**
   * Check the PubDbXref contains the Pub in FeatureCvTerm
   * @param pubDbXRefs
   * @param pubId
   * @param pub
   * @param feature_cvterm
   */
  private static void checkPubDbXRef(final List pubDbXRefs, final int pubId,
                                     final Pub pub, final FeatureCvTerm feature_cvterm)
  {
    PubDbXRef pubDbXRef = null;
    for(int i = 0; i < pubDbXRefs.size(); i++)
    {
      pubDbXRef = (PubDbXRef) pubDbXRefs.get(i);
      if(pubDbXRef.getPub().getPubId() == pubId)
      {
        DbXRef dbxref = pubDbXRef.getDbXRef();
        Splash.logger4j.debug("Checking PubDbXRef and found Pub "+dbxref.getDb().getName()+
                              ":"+dbxref.getAccession());
        break;
      }
    }
    
    if(pubDbXRef == null || 
        !pub.getUniqueName().endsWith(pubDbXRef.getDbXRef().getAccession()))
     {
       Splash.logger4j.debug("Checking PubDbXRef and not found Pub "+
                           feature_cvterm.getPub().getUniqueName());
      
       JOptionPane.showMessageDialog(null, "Cannot find pub_dbxref for:\n"+
           feature_cvterm.getPub().getUniqueName(), 
           "Database Error",
           JOptionPane.ERROR_MESSAGE);
     }
  }
  
  /**
   * Look up the cvterm_id for a controlled vocabulary name.
   * @param name  
   * @return
   */
  public static Integer getCvtermID(final String name)
  {
    Enumeration enum_cvterm = cvterms.keys();
    while(enum_cvterm.hasMoreElements())
    {
      Integer key = (Integer)enum_cvterm.nextElement();
      if(name.equalsIgnoreCase( ((CvTerm)cvterms.get(key)).getName() ))
        return key;
    }
    return null;
  }

  /**
   * Look up a cvterm name from the collection of cvterms.
   * @param id  a cvterm_id  
   * @return    the cvterm name
   */
  private static String getCvtermName(final int id, 
                                      final GmodDAO dao,
                                      final boolean gene_builder)
  {
    if(gene_builder)
      return dao.getCvTermById(id).getName();
    
    return getCvTerm(id, dao, gene_builder).getName();
  }
  
  private static CvTerm getCvTerm(final int id, 
                                  final GmodDAO dao,
                                  final boolean gene_builder)
  {
    if(gene_builder)
      return dao.getCvTermById(id);
    if(cvterms == null)
      getCvterms(dao);

    return (CvTerm)cvterms.get(new Integer(id));
  }
  
  /**
   * Use the CvTerm name to return a CvTerm.
   * @param cvTermId
   * @return
   */
  public static CvTerm getCvTermByCvTermName(final String cvterm_name)
  {
    Enumeration enum_cvterm = cvterms.elements();
    while(enum_cvterm.hasMoreElements())
    {
      CvTerm cvterm = (CvTerm)enum_cvterm.nextElement();
      if(cvterm_name.equalsIgnoreCase( cvterm.getName() ))
        return cvterm;
    }
    
    return null;
  }
  
  /**
   * Use the CvTermId to return a CvTerm.
   * @param cvTermId
   * @return
   */
  public static CvTerm getCvTermByCvTermId(final int cvTermId,
                                           final uk.ac.sanger.artemis.io.Feature feature)
  {
    if(cvterms == null)
    {
      try
      {
        DatabaseDocument doc =
          (DatabaseDocument) ((DocumentEntry)feature.getEntry()).getDocument();
        return doc.getDAO().getCvTermById(cvTermId);
      }
      catch(ConnectException e)
      { 
        logger4j.warn(e.getMessage()); 
      }
      catch(SQLException e) 
      { 
        logger4j.warn(e.getMessage()); 
      }
    }
    
    Enumeration enum_cvterm = cvterms.elements();
    while(enum_cvterm.hasMoreElements())
    {
      CvTerm cvterm = (CvTerm)enum_cvterm.nextElement();
      if(cvterm.getCvTermId() == cvTermId)
        return cvterm;
    }
    
    return null;
  }
  
  /**
   * Use the Cv name and CvTerm name to return a CvTerm.
   * @param cvterm_name
   * @param cvName
   * @return
   */
  public static CvTerm getCvTermByCvAndCvTerm(final String cvterm_name,
                                              final String cvName)
  {
    final Enumeration enum_cvterm = cvterms.elements();
    while(enum_cvterm.hasMoreElements())
    {
      final CvTerm cvterm = (CvTerm)enum_cvterm.nextElement();
      if(cvName.equals( cvterm.getCv().getName() ) &&
         cvterm_name.equals( cvterm.getName() ))
        return cvterm;
    }
    return null;
  }

  /**
   * This is different from getCvTermByCvAndCvTerm as it makes sure
   * the Cv name part matches the name supplied to the function, i.e.
   * by matching just the start.
   * @param cvterm_name
   * @param cvName
   * @return
   */
  public static CvTerm getCvTermByCvPartAndCvTerm(final String cvterm_name,
                                              final String cvName)
  {
    Enumeration enum_cvterm = cvterms.elements();
    while(enum_cvterm.hasMoreElements())
    {
      CvTerm cvterm = (CvTerm)enum_cvterm.nextElement();
      if(cvterm.getCv().getName().startsWith( cvName ) &&
         cvterm_name.equals( cvterm.getName() ))
        return cvterm;
    }
    return null;
  }
  /**
   * Look up cvterms names and id and return in a hashtable.
   * @param dao the data access object
   * @return    the cvterm <code>Hashtable</code>
   */
  private static Hashtable getCvterms(final GmodDAO dao)
  {
    try
    {
      final List cvterm_list = dao.getCvTerms();
      final Iterator it = cvterm_list.iterator();
      cvterms = new Hashtable(cvterm_list.size());

      while(it.hasNext())
      {
        final CvTerm cvterm = (CvTerm)it.next();
        cvterms.put(new Integer(cvterm.getCvTermId()), cvterm);
      }
    }
    catch(RuntimeException sqle)
    {
      System.err.println("SQLException retrieving CvTerms");
      System.err.println(sqle);
    }
    logger4j.debug("LOADED CvTerms");
    
    return cvterms;
  }
  
  /**
   * Get CvTerm's in a given CV
   * @param cvName
   * @return
   */
  public Vector getCvTermsByCvName(final String cvName)
  {
    if(cvterms == null)
    {
      logger4j.debug("getCvTermsByCvName LOADING CVTERMS");
      getCvterms(getDAOOnly());
    }
    
    return getCvterms("", cvName);
  }
  
  public List getDatabaseNames()
  {
    GmodDAO dao = getDAOOnly();
    List dbs = dao.getDbs();
    List names = new Vector();
    Iterator it = dbs.iterator();
    while(it.hasNext())
    {
      Db db = (Db)it.next();
      names.add(db.getName());
    }
    return names;
  }
  
  public List getOrganismNames()
  {
    if(organismNames != null && organismNames.size() > 0)
      return organismNames;
    
    GmodDAO dao = getDAOOnly();
    List organisms = dao.getOrganisms();
    organismNames = new Vector();
    Iterator it = organisms.iterator();
    while(it.hasNext())
    {
      Organism organism = (Organism)it.next();
      organismNames.add(organism.getCommonName());
    }
    return organismNames;
  }
  
  public static Vector getCvterms(final String search_str, final String cv_name)
  {
    final Vector cvterm_match = new Vector();
    
    Enumeration enum_cvterm = cvterms.keys();
    while(enum_cvterm.hasMoreElements())
    {
      Integer key = (Integer)enum_cvterm.nextElement();
      CvTerm cvterm = (CvTerm)cvterms.get(key);
      
      if(cvterm.getCv().getName().startsWith(cv_name))
      {
        if(cvterm.getName().indexOf(search_str) > -1)
          cvterm_match.add(cvterm);
      }
    }
    
    return cvterm_match;
  }
  
  /**
   * Get a list of the CV names
   * @return
   */
  public static List getCvControledCurationNames()
  {
    if(cvControledCuratioNames != null)
      return cvControledCuratioNames;
    cvControledCuratioNames = new Vector();
    final Enumeration enum_cvterm = cvterms.elements();
    while(enum_cvterm.hasMoreElements())
    {
      final CvTerm cvTerm = (CvTerm)enum_cvterm.nextElement();
      final String cvNameStr = cvTerm.getCv().getName();
      
      if(cvNameStr.startsWith(DatabaseDocument.CONTROLLED_CURATION_TAG_CVNAME) && 
         !cvControledCuratioNames.contains(cvNameStr))
        cvControledCuratioNames.add(cvNameStr);
    }
    
    return cvControledCuratioNames;
  }

  /**
   * Look up synonym type names e.g. synonym, systematic_id.
   * @return    the synonym tag names
   */
  public static String[] getSynonymTypeNames(final String cv_name, 
                                             final GFFStreamFeature feature)
  {
    if(cvterms == null)
    {
      DatabaseDocument doc = (DatabaseDocument)feature.getDocumentEntry().getDocument();
      try
      {
        Cv cv = new Cv();
        cv.setName(cv_name);
        List synonymCvTerms = doc.getDAO().getCvTermByNameInCv(null, cv);
        String synonymNames[] = new String[synonymCvTerms.size()];
        for(int i=0; i<synonymCvTerms.size(); i++)
          synonymNames[i] = ((CvTerm) synonymCvTerms.get(i)).getName();

        return synonymNames;
      }
      catch(ConnectException e){}
      catch(SQLException e){}   
    }
    
    Vector synonym_names = new Vector();
    Enumeration cvterm_enum = cvterms.elements();
    while(cvterm_enum.hasMoreElements())
    {
      CvTerm cvterm = (CvTerm)cvterm_enum.nextElement();
      if(cvterm.getCv().getName().equals(cv_name))
        synonym_names.add(cvterm.getName());
    }
    
    return (String[])synonym_names.toArray(
                       new String[synonym_names.size()]);
  }
  
  public void insertCvTerm(CvTerm cvTerm)
  {
    final GmodDAO dao = getDAOOnly();
    dao.persist(cvTerm);
    cvTerm = dao.getCvTermByNameAndCvName(cvTerm.getName(), cvTerm.getCv().getName());
    cvterms.put(new Integer(cvTerm.getCvTermId()), cvTerm);
  }
  
  /**
   * Get the sequence for a feature.
   * @param dao   the data access object
   * @param buff  the buffer to add the sequence to
   * @return      the resulting buffer
   * @throws java.sql.SQLException
   */
  private Feature getChadoSequence(GmodDAO dao, ByteBuffer buff)
  {
    Feature feature = dao.getFeatureById(Integer.parseInt(srcFeatureId));
    getChadoSequence(feature, buff);
    return feature;
  }
  
  
  /**
   * Get the sequence for a feature.
   * @param dao   the data access object
   * @param buff  the buffer to add the sequence to
   * @return      the resulting buffer
   * @throws java.sql.SQLException
   */
  private void getChadoSequence(final Feature feature, ByteBuffer buff)
  {
    buff.append("##FASTA\n>");
    buff.append(feature.getUniqueName());
    buff.append("\n");
    buff.append(feature.getResidues());
  }

  /**
   * Get the CDS FeatureLoc's associated with a give protein
   * @param peptideName
   * @return
   */
  public List getCdsFeatureLocsByPeptideName(final String peptideName)
  {
    Feature peptideFeature = getFeatureByUniquename(peptideName);
    Collection frs = peptideFeature.getFeatureRelationshipsForSubjectId();
    Iterator it = frs.iterator();
    Feature transcriptFeature = null;
    while(it.hasNext())
    {
      FeatureRelationship fr = (FeatureRelationship)it.next();
      if(fr.getCvTerm().getName().equalsIgnoreCase("derives_from"))
      {
        transcriptFeature = fr.getFeatureByObjectId();
        logger4j.debug("TRANSCRIPT :: "+transcriptFeature.getUniqueName());
        break;
      }
    }
    
    if(transcriptFeature == null)
      return null;
    return getCdsFeatureLocsByTranscriptName(transcriptFeature.getUniqueName());
  }
  
  /**
   * Get the CDS FeatureLoc's associated with a given transcript
   * @param transcriptName
   * @return
   */
  public List getCdsFeatureLocsByTranscriptName(final String transcriptName)
  {
    Feature transcriptFeature = getFeatureByUniquename(transcriptName);
    if(transcriptFeature == null)
      return null;
    
    Collection frs = transcriptFeature.getFeatureRelationshipsForObjectId();
    Iterator it = frs.iterator();
    List cdsFeatureLocs = new Vector();
    while(it.hasNext())
    {
      FeatureRelationship fr = (FeatureRelationship)it.next();
      org.gmod.schema.sequence.Feature child = fr.getFeatureBySubjectId();
      if(child.getCvTerm().getName().equals("exon") || 
         child.getCvTerm().getName().equals("pseudogenic_exon"))
      {
        Collection featureLocs = child.getFeatureLocsForFeatureId();

        Iterator it2 = featureLocs.iterator();
        while(it2.hasNext())
        {
          FeatureLoc featureLoc = (FeatureLoc) it2.next();
          cdsFeatureLocs.add(featureLoc);
        }
      }
    }

    Collections.sort(cdsFeatureLocs, new LocationComarator());
    return cdsFeatureLocs;
  }
  
  /**
   * Get the sequence for a feature.
   * @param uniqueName   the feature
   * @return      the resulting buffer
   */
  public PartialSequence getChadoSequence(final String uniqueName)
  {
    Feature feature = (Feature) getDAOOnly().getResiduesByUniqueName(uniqueName).get(0);
    char[] c = getChars(feature.getResidues());
    
    PartialSequence ps = new PartialSequence(c, feature.getSeqLen(),
                                feature.getFeatureLoc().getFmin().intValue()+1,
                                feature.getFeatureLoc().getStrand(),
                                feature.getFeatureLoc().getPhase());
    return ps;
  }
  
  /**
   * Convert byte array to char array
   * @param b byte array
   * @return  char array
   */
  private char[] getChars(final byte b[])
  {
    char[] c = new char[b.length];

    for(int i = 0; i < b.length; i++)
      c[i] = (char)b[i];
    return c;
  }
  
  /**
   * Get the <code>List</code> of available schemas.
   * @return  the <code>List</code> of available schemas
   */
  public List getSchema()
  {
    return schema_list;
  }
  
  public Feature getFeatureByUniquename(final String uniqueName) 
  {
    GmodDAO dao = getDAOOnly();
    List features = dao.getFeaturesByUniqueName(uniqueName);
    if(features == null || features.size() < 1)
      return null;
      
    return (Feature)(dao.getFeaturesByUniqueName(uniqueName).get(0));
  }
  
  /**
   * Given a gene unique name return the poplypeptide chado features that belong
   * to that gene
   * @param geneName
   * @return
   */
  public Vector getPolypeptideFeatures(final String geneName)
  {
    Feature geneFeature =  getFeatureByUniquename(geneName);
    if(geneFeature == null)
      return null;
    
    Collection frs = geneFeature.getFeatureRelationshipsForObjectId();
    Iterator it = frs.iterator();
    List transcripts = new Vector(frs.size());
    while(it.hasNext())
    {
      FeatureRelationship fr = (FeatureRelationship)it.next();
      transcripts.add(fr.getFeatureBySubjectId());
    }
    
    Vector polypep = new Vector();
    for(int i=0; i<transcripts.size(); i++)
    {
      org.gmod.schema.sequence.Feature transcript = 
        (org.gmod.schema.sequence.Feature) transcripts.get(i);
      frs = transcript.getFeatureRelationshipsForObjectId();
      it = frs.iterator();
      while(it.hasNext())
      {
        FeatureRelationship fr = (FeatureRelationship)it.next();
        if(fr.getCvTerm().getName().equalsIgnoreCase("derives_from"))
          if(fr.getFeatureBySubjectId().getCvTerm().getName().equalsIgnoreCase("polypeptide"))
            polypep.add(fr.getFeatureBySubjectId());
      }
    }
    return polypep;
  }
  
  /**
   * Given a gene unique name return the poplypeptides that belong
   * to that gene
   * @param geneName
   * @return
   */
  /*public Vector getPolypeptideNames(final String geneName)
  {
    Vector polypeptides = getPolypeptideFeatures(geneName);
    Vector polypeptideNames = new Vector(polypeptides.size());
    for(int i=0; i<polypeptides.size(); i++)
    {
      Feature feature = (Feature)polypeptides.get(i);
      polypeptideNames.add(feature.getUniqueName());
    }
    return polypeptideNames;
  }*/
  
  public List getClustersByFeatureIds(final List featureIds)
  {
    GmodDAO dao = getDAOOnly();
    return dao.getClustersByFeatureIds(featureIds);
  }
  
  public List getParentFeaturesByChildFeatureIds(final List subjectIds)
  {
    GmodDAO dao = getDAOOnly();
    return dao.getParentFeaturesByChildFeatureIds(subjectIds);
  }
  
  public List getFeatureDbXRefsByFeatureId(final List featureIds)
  {
    GmodDAO dao = getDAOOnly();
    return dao.getFeatureDbXRefsByFeatureId(featureIds);
  }
  
  
  /**
   * Used by SimilarityLazyQualifierValue.bulkRetrieve() to get the match features
   * @param featureIds the <code>List</code> of feature_id's
   * @return  the corresponding features
   */
  public List getFeaturesByListOfIds(final List featureIds)
  {
    GmodDAO dao = getDAOOnly();
    return dao.getFeaturesByListOfIds(featureIds);
  }
  
  public List getFeaturePropByFeatureIds(final List featureIds)
  {
    GmodDAO dao = getDAOOnly();
    return dao.getFeaturePropByFeatureIds(featureIds);
  }
  
  public List getSimilarityMatches(List featureIds)
  {
    GmodDAO dao = getDAOOnly();
    if(featureIds == null)
      return dao.getSimilarityMatches(new Integer(srcFeatureId));
    else
      return dao.getSimilarityMatchesByFeatureIds(featureIds);
  }
  
  public List getFeatureLocsByListOfIds(List featureIds)
  {
    GmodDAO dao = getDAOOnly();
    return dao.getFeatureLocsByListOfIds(featureIds);
  }
  
  /**
   * Create a hashtable of the available entries with residues.
   * @return a <code>Hashtable</code> of the <code>String</code>
   *          representation (schema-type-feature_name) and the
   *          corresponding feature_id
   * @throws ConnectException
   * @throws java.sql.SQLException
   */
  public HashMap getDatabaseEntriesMultiSchema()
                   throws ConnectException, java.sql.SQLException
  {
    String schema = null;
    HashMap db    = null;
    try
    { 
      GmodDAO dao = getDAO();
      schema_list = dao.getSchema(); 
      
      Iterator it = schema_list.iterator();
        
      while(it.hasNext())
      {
        schema = (String)it.next();
        
        reset((String)getLocation(),  schema);

        try
        {
          dao = getDAO();
          List list_residue_features = dao.getResidueFeatures();
          
          Iterator it_residue_features = list_residue_features.iterator();
          while(it_residue_features.hasNext())
          {
            Feature feature = (Feature)it_residue_features.next();
            String typeName = getCvtermName(feature.getCvTerm().getCvTermId(), getDAO(), gene_builder); 
          
            if(db == null)
              db = new HashMap();
            db.put(schema + " - " + typeName + " - " + feature.getUniqueName(),
                   Integer.toString(feature.getFeatureId()));
          }
        }
        catch(RuntimeException e){}
        catch(java.sql.SQLException sqlExp){}
      }
      
    }
    catch(RuntimeException sqlExp)
    {
      JOptionPane.showMessageDialog(null, "SQL Problems...\n"+
                                    sqlExp.getMessage(), 
                                    "SQL Error",
                                    JOptionPane.ERROR_MESSAGE);
      
      logger4j.debug(sqlExp.getMessage());
      //sqlExp.printStackTrace();
    }
    catch(ConnectException exp)
    {
      JOptionPane.showMessageDialog(null, "Connection Problems...\n"+
            exp.getMessage(), 
            "Connection Error",
            JOptionPane.ERROR_MESSAGE);
      logger4j.debug(exp.getMessage());
      throw exp;
    }
    catch(java.sql.SQLException sqlExp)
    {
      JOptionPane.showMessageDialog(null, "SQL Problems....\n"+
                                    sqlExp.getMessage(), 
                                    "SQL Error",
                                    JOptionPane.ERROR_MESSAGE);
      logger4j.debug(sqlExp.getMessage());
      throw sqlExp;
    }
    
    return db;
  }
  
  
  /**
   * Create a hashtable of the available entries with residues.
   * 
   * @return a <code>Hashtable</code> of the <code>String</code>
   *          representation (schema-type-feature_name) and the
   *          corresponding feature_id
   * @throws ConnectException
   * @throws java.sql.SQLException
   */
  public HashMap getDatabaseEntries()
                   throws ConnectException, java.sql.SQLException
  {

    HashMap db = null;
    try
    { 
      GmodDAO dao = getDAO();
      CvTermThread cvThread = new CvTermThread(dao);
      cvThread.start();
      
      schema_list = dao.getOrganisms();
      
      /*Organism org = new Organism();
      org.setCommonName("web");
      schema_list.add(org);*/
      
      final List pg_schemas = dao.getSchema();
      
      // build a lookup hash of residue features in the
      // main schema
      Hashtable residueFeaturesLookup = new Hashtable();
      List list_residue_features = dao.getResidueFeatures();
      for(int i=0; i<list_residue_features.size(); i++)
      {
        Feature feature = (Feature)list_residue_features.get(i);
        Integer organismId = new Integer(feature.getOrganism().getOrganismId());
        List features;
        if(residueFeaturesLookup.containsKey(organismId))
          features= (List)residueFeaturesLookup.get(organismId);
        else
          features = new Vector();
        features.add(feature);
        residueFeaturesLookup.put(organismId, features);
      }
      
      // loop over organisms to identify those with features 
      // containing residues
      Iterator it = schema_list.iterator();
      while(it.hasNext())
      {
        final Organism organism = (Organism)it.next();
        String orgName = organism.getCommonName();
        
        if(orgName == null || orgName.equals(""))
        {
          orgName = organism.getGenus() + "." + organism.getSpecies();
          organism.setCommonName(orgName);
        }
        
        // search to see if this is in its own schema
        Iterator schemasIt = pg_schemas.iterator();
        while(schemasIt.hasNext())
        {
          schema = (String)schemasIt.next();         
          if( schema.equalsIgnoreCase(organism.getCommonName()) )
          {
            reset((String)getLocation(),  schema);
            dao = getDAO();
            orgName = schema;

            singleSchema = false;
            break;
          }
        }
        
        try
        {
          if(organism.getOrganismId() > 0)
          {
            Integer organismId = new Integer(organism.getOrganismId());
            
            if(residueFeaturesLookup.containsKey(organismId))
              list_residue_features = (List)residueFeaturesLookup.get(organismId);
            else if(singleSchema && residueFeaturesLookup.size()>0)
              continue;
            else
              list_residue_features = dao.getResidueFeatures(organismId);
          }
          else
            list_residue_features = 
              dao.getResidueFeaturesByOrganismCommonName(organism.getCommonName());
          
          Iterator it_residue_features = list_residue_features.iterator();
          while(it_residue_features.hasNext())
          {
            final Feature feature = (Feature)it_residue_features.next();
            final String typeName = feature.getCvTerm().getName();
                  
            if(db == null)
              db = new HashMap();
            if(organismNames == null)
              organismNames = new Vector();
            
            if(!organismNames.contains(orgName))
              organismNames.add(orgName);         
            
            db.put(orgName + " - " + typeName + " - " + feature.getUniqueName(),
                   Integer.toString(feature.getFeatureId())); 
          }
        }
        catch(RuntimeException e)
        {
          e.printStackTrace();
        }
      }
      
      residueFeaturesLookup.clear();
      list_residue_features.clear();
      
      // now wait for cvterm to be loaded
      while(cvThread.isAlive())
        Thread.sleep(10);
    }
    catch(RuntimeException sqlExp)
    {
      JOptionPane.showMessageDialog(null, "SQL Problems...\n"+
                                    getLocation()+"\n"+
                                    sqlExp.getMessage(), 
                                    "SQL Error",
                                    JOptionPane.ERROR_MESSAGE);
      
      logger4j.debug(sqlExp.getMessage());
      //sqlExp.printStackTrace();
    }
    catch(ConnectException exp)
    {
      JOptionPane.showMessageDialog(null, "Connection Problems...\n"+
            exp.getMessage(), 
            "Connection Error",
            JOptionPane.ERROR_MESSAGE);
      logger4j.debug(exp.getMessage());
      throw exp;
    }
    catch(java.sql.SQLException sqlExp)
    {
      JOptionPane.showMessageDialog(null, "SQL Problems....\n"+
                                    getLocation()+"\n"+
                                    sqlExp.getMessage(), 
                                    "SQL Error",
                                    JOptionPane.ERROR_MESSAGE);
      logger4j.debug(sqlExp.getMessage());
      throw sqlExp;
    }
    catch(InterruptedException e)
    {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    
    return db;
  }
  
  /**
   * 
   */
  public void showCvTermLookUp()
  {
    try
    {
      new ChadoCvTermView( getDAO() );
    }
    catch(ConnectException e)
    {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    catch(SQLException e)
    {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
  }
  
  /**
   * Get the data access object (DAO).
   * @return data access object
   */
  private GmodDAO getDAO()
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
      {
        System.setProperty("chado", (String)getLocation());
        connIB = new IBatisDAO(pfield);
      }
      return connIB;
    }
  }
  
  /**
   * Get the username for this connection
   * @return
   */
  public String getUserName()
  {
    // "localhost:10001/backup?chado"
    
    String url = (String)getLocation();
    int index  = url.indexOf("?");
    
    String userName = url.substring(index+1).trim();
    if(userName.startsWith("user="))
      userName = userName.substring(5);
    
    return userName;
  }
  
  private GmodDAO getDAOOnly()
  {
    GmodDAO dao = null;
    try
    {
      dao = getDAO();
    }
    catch(RuntimeException sqlExp)
    {
      JOptionPane.showMessageDialog(null, "SQL Problems...\n"+
                                    sqlExp.getMessage(), 
                                    "SQL Error",
                                    JOptionPane.ERROR_MESSAGE);
    }
    catch(ConnectException exp)
    {
      JOptionPane.showMessageDialog(null, "Connection Problems...\n"+
            exp.getMessage(), 
            "Connection Error",
            JOptionPane.ERROR_MESSAGE);
    }
    catch(java.sql.SQLException sqlExp)
    {
      JOptionPane.showMessageDialog(null, "SQL Problems....\n"+
                                    sqlExp.getMessage(), 
                                    "SQL Error",
                                    JOptionPane.ERROR_MESSAGE);
    }
    return dao;
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
  public int commit(final Vector sql,
                    final boolean force)
  {
    GmodDAO dao = null;
    int ncommit = -1;
    final Hashtable featureIdStore = new Hashtable();
    boolean useTransactions = false;

    try
    {
      dao = getDAO();
      
      if(!force && dao instanceof IBatisDAO)
        useTransactions = true;
      
      if(useTransactions)
      {
        ((IBatisDAO) dao).startTransaction();
        logger4j.debug("START TRANSACTION");
      }
      boolean unchanged;
      
      //
      // check feature timestamps have not changed
      Vector names_checked = new Vector();
      
      for(int i = 0; i < sql.size(); i++)
      {
        final ChadoTransaction tsn = (ChadoTransaction)sql.get(i);
        final Object uniquenames[] = getUniqueNames(tsn);
        
        if(uniquenames == null)
          continue;
        
        for(int j=0; j<uniquenames.length; j++)
        {
          final String uniquename = (String) uniquenames[j];
          
          if(uniquename == null || names_checked.contains(uniquename))
            continue;

          names_checked.add(uniquename);
          final String keyName = tsn.getFeatureKey();

          unchanged = checkFeatureTimestamp(schema, uniquename, 
              dao, keyName, featureIdStore, tsn);
          if(!unchanged)
          {
            if(useTransactions)
              ((IBatisDAO) dao).endTransaction();
            return 0;
          }
        }
      }  
      
      final Timestamp ts = new Timestamp(new java.util.Date().getTime());
      //
      // commit to database
      for(ncommit = 0; ncommit < sql.size(); ncommit++)
      {
        try
        {
          ChadoTransaction tsn = (ChadoTransaction) sql.get(ncommit);
          commitChadoTransaction(tsn, dao, ts);
        }
        catch (RuntimeException re)
        {
          if(!force)
            throw re;
          logger4j.warn(constructExceptionMessage(re, sql, ncommit));
          logger4j.warn("NOW TRYING TO CONTINUE TO COMMIT");
        }
      }

      //
      // update timelastmodified timestamp
      names_checked = new Vector();
      for(int i = 0; i < sql.size(); i++)
      {
        final ChadoTransaction tsn = (ChadoTransaction) sql.get(i);
        final Object uniquenames[] = getUniqueNames(tsn);
             
        if(uniquenames == null)
          continue;
        
        if(tsn.getType() == ChadoTransaction.UPDATE &&
           tsn.getFeatureObject() instanceof Feature)
        {
          for(int j=0; j<uniquenames.length; j++)
            names_checked.add((String) uniquenames[j]);
          continue;  
        }
        
        for(int j=0; j<uniquenames.length; j++)
        {
          final String uniquename = (String) uniquenames[j];
          if(uniquename == null || names_checked.contains(uniquename))
            continue;

          names_checked.add(uniquename);

          final Feature feature;

          // retieve from featureId store
          if(featureIdStore != null && featureIdStore.containsKey(uniquename))
          {
            Feature f = (Feature) featureIdStore.get(uniquename);

            feature = new Feature();
            feature.setFeatureId(f.getFeatureId());
            feature.setUniqueName(uniquename);
            feature.setObsolete(f.isObsolete());
          }
          else
            feature = dao.getFeatureByUniqueName(uniquename, 
                                       tsn.getFeatureKey());

          if(feature != null)
          {
            feature.setTimeLastModified(ts);
            feature.setName("0");  // do not change name
            dao.merge(feature);
          }
        }

        GFFStreamFeature gff_feature = (GFFStreamFeature) tsn.getGff_feature();
        gff_feature.setLastModified(ts);
      }

      final String nocommit = System.getProperty("nocommit");
      if( useTransactions && 
          (nocommit == null || nocommit.equals("false")))
      {
        ((IBatisDAO) dao).commitTransaction();
        logger4j.debug("TRANSACTION COMPLETE");
      }
      else if(useTransactions && 
          (nocommit != null && nocommit.equals("true")))
        logger4j.debug("TRANSACTION NOT COMMITTED : nocommit property set to true");
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
    catch (RuntimeException re)
    {
      final String msg = constructExceptionMessage(re, sql, ncommit);
      JOptionPane.showMessageDialog(null, msg,
          "Problems Writing to Database ",
          JOptionPane.ERROR_MESSAGE);
      logger4j.error(msg);
      //re.printStackTrace();
    }
    finally
    {
      if(useTransactions)
        try
        {
          ((IBatisDAO) dao).endTransaction();
          logger4j.debug("END TRANSACTION");
        }
        catch(SQLException e){ e.printStackTrace(); }
    }
    
    if(featureIdStore != null)
      featureIdStore.clear();
        
    return ncommit;
  }
  
  /**
   * Get the uniquenames involved in a transaction
   * @param tsn
   * @return
   */
  private Object[] getUniqueNames(final ChadoTransaction tsn)
  {
    if(tsn.getGff_feature() == null)
      return null;
    if(tsn.getGff_feature().getSegmentRangeStore() == null ||
       tsn.getGff_feature().getSegmentRangeStore().size() < 2 ||
       tsn.getFeatureObject() instanceof FeatureProp)
      return new Object[]{ tsn.getUniquename() };
    else 
      return tsn.getGff_feature().getSegmentRangeStore().keySet().toArray();
  }
  
  /**
   * Construct an exeption message from the ChadoTransaction
   * @param re
   * @param sql
   * @param ncommit
   * @return
   */
  private String constructExceptionMessage(final RuntimeException re,
                                           final Vector sql,
                                           final int ncommit)
  {
    String msg = "";
    if(ncommit > -1 && ncommit < sql.size())
    {
      final ChadoTransaction t_failed = (ChadoTransaction)sql.get(ncommit);
      
      if(t_failed.getType() == ChadoTransaction.DELETE)
        msg = "DELETE failed ";
      else if(t_failed.getType() == ChadoTransaction.INSERT)
        msg = "INSERT failed ";
      else if(t_failed.getType() == ChadoTransaction.UPDATE)
        msg = "UPDATE failed ";
        
      if(t_failed.getUniquename() != null)
        msg = msg + "for " + t_failed.getUniquename()+":";
      else if(t_failed.getFeatureObject() != null &&
              t_failed.getFeatureObject() instanceof Feature)
      {
        final Feature chadoFeature = (Feature)t_failed.getFeatureObject();
        if(chadoFeature.getUniqueName() != null)
          msg = msg + "for " + chadoFeature.getUniqueName() +":";
      }
        
      msg = msg+"\n";
    }
    
    return msg + re.getMessage();  
  }
  
  /**
   * Commit a single chado transaction
   * @param tsn
   * @param dao
   */
  private void commitChadoTransaction(final ChadoTransaction tsn,
                                      final GmodDAO dao,
                                      final Timestamp ts)
  {
    if(tsn.getType() == ChadoTransaction.UPDATE)
    {
      if(tsn.getFeatureObject() instanceof Feature)
      {
        Feature feature = (Feature)tsn.getFeatureObject();

        if(feature.getUniqueName() != null)
        {
          final String uniquename;
          if(tsn.getOldUniquename() != null)
            uniquename = (String)tsn.getOldUniquename();
          else
            uniquename = feature.getUniqueName();
          
          Feature old_feature
              = dao.getFeatureByUniqueName(uniquename, tsn.getFeatureKey());
          
          if(old_feature != null)
            feature.setFeatureId( old_feature.getFeatureId() );
          
          tsn.setOldUniquename(feature.getUniqueName());
        }
        feature.setTimeLastModified(ts);
      }
      dao.merge(tsn.getFeatureObject());
      //dao.updateAttributes(tsn);
    }
    else if(tsn.getType() == ChadoTransaction.INSERT)
    {
      if(tsn.getFeatureObject() instanceof FeatureCvTerm)
        ArtemisUtils.inserFeatureCvTerm(dao, (FeatureCvTerm)tsn.getFeatureObject());
      else
      {
        // set srcfeature_id
        if(tsn.getFeatureObject() instanceof Feature &&
            ((Feature) tsn.getFeatureObject()).getFeatureLoc() != null)
        {
          FeatureLoc featureloc = ((Feature) tsn.getFeatureObject()).getFeatureLoc();
          Feature featureBySrcFeatureId = new Feature();
          featureBySrcFeatureId.setFeatureId(Integer.parseInt(srcFeatureId));
          featureloc.setFeatureBySrcFeatureId(featureBySrcFeatureId);
        }
        dao.persist(tsn.getFeatureObject());
      }
      
    }
    else if(tsn.getType() == ChadoTransaction.DELETE)
    {
      if(tsn.getFeatureObject() instanceof FeatureCvTerm)
        ArtemisUtils.deleteFeatureCvTerm(dao, (FeatureCvTerm)tsn.getFeatureObject());
      else
        dao.delete(tsn.getFeatureObject());
    }
  }
  
  /**
   * Check the <code>Timestamp</code> on a feature (for versioning).
   * @param schema      the schema
   * @param uniquename  the feature uniquename
   * @param timestamp   the last read feature timestamp
   */
  public boolean checkFeatureTimestamp(final String schema,
                                       final String uniquename,
                                       final GmodDAO dao,
                                       final String keyName,
                                       final Hashtable featureIdStore,
                                       final ChadoTransaction tsn)
  {
    final Timestamp timestamp  = tsn.getLastModified();
    final Object featureObject = tsn.getFeatureObject();
    
    final Feature feature = dao.getFeatureByUniqueName(uniquename, keyName);
    if(feature == null)
      return true;
    
    featureIdStore.put(uniquename, feature);
    
    if(featureObject instanceof FeatureProp)
      ((FeatureProp)featureObject).setFeature(feature);
    else if(featureObject instanceof FeatureLoc)
    {
      if(((FeatureLoc)featureObject).getFeatureByFeatureId().getUniqueName().equals(uniquename))
      {
        logger4j.debug("Setting featureId for:"  + uniquename );
        ((FeatureLoc)featureObject).setFeatureByFeatureId(feature);
      }
    }
    
    final Timestamp now = feature.getTimeLastModified();
    
    if(now != null && timestamp != null)
    {
      now.setNanos(0);
      timestamp.setNanos(0);
      
      if(now.compareTo(timestamp) != 0)
      {
        final SimpleDateFormat date_format = 
                   new SimpleDateFormat("dd.MM.yyyy hh:mm:ss z");
        
        //System.out.println(date_format.format(now)+"   "+
        //                   date_format.format(timestamp));
        int select = JOptionPane.showConfirmDialog(null, uniquename +
                                      " has been altered at :\n"+
                                      date_format.format(now)+"\nOverwite?", 
                                      "Feature Changed", 
                                      JOptionPane.OK_CANCEL_OPTION);
        if(select == JOptionPane.OK_OPTION)
          return true;
        else
          return false;
      }
    }
    return true;
  }

  
  public static void main(String args[])
  {
    try
    {
      GmodDAO dao;
      DatabaseEntrySource src = new DatabaseEntrySource();
      src.setLocation(true);
      
      if(System.getProperty("ibatis") == null)
        dao = new JdbcDAO(src.getLocation(), src.getPfield()); 
      else
        dao = new IBatisDAO(src.getPfield());
      
      Feature feature = new Feature();
      feature.setUniqueName(args[0]);
      List schemas = new Vector();
      schemas.add(args[1]);
      List featureList = new Vector();
      featureList.add(dao.getFeatureByUniqueName(args[0], "polypeptide")); 
      System.out.println("FINISHED getFeature()");
      for(int i = 0; i < featureList.size(); i++)
      {
        feature = (Feature)featureList.get(i);
        
        String abb  = feature.getOrganism().getAbbreviation();
        String type = feature.getCvTerm().getName();
        int fmin    = feature.getFeatureLoc().getFmin().intValue() + 1;
        int fmax    = feature.getFeatureLoc().getFmax().intValue();
        String featprop = 
          ((FeatureProp)(new Vector(feature.getFeatureProps()).get(0))).getCvTerm().getName();
        
        System.out.print(fmin+".."+fmax);
        System.out.print(" "+type);
        System.out.print(" "+featprop);
        System.out.print(" "+feature.getFeatureLoc().getStrand());
        System.out.print(" "+feature.getUniqueName());
        System.out.print(" "+abb);
        System.out.println(" "+Integer.toString(feature.getFeatureId()));
        
/*      Hashtable synonyms = getAllFeatureSynonyms(dao, null);
        Vector syns = (Vector)synonyms.get(new Integer(feature.getId()));
        for(int j=0; j<syns.size(); j++)
        {
          FeatureSynonym alias = (FeatureSynonym)syns.get(j);
          System.out.print(" "+alias.getSynonym().getCvterm().getName()+
                           "="+alias.getSynonym().getName());
        }*/
        
        System.out.println(" "); 
      }
    }
    catch(SQLException sqle)
    {
      sqle.printStackTrace();
    }
    catch(RuntimeException re)
    {
      re.printStackTrace();
    }
    catch(ConnectException e)
    {
      e.printStackTrace();
    }
  }

  public Document getParent()
  {
    return null;
  }
  
  /**
   * Find from a list the FeatureLoc with a given srcFeature
   * @param locs
   * @param srcfeature_id
   * @return
   */
  public static FeatureLoc getFeatureLoc(List locs, int srcfeature_id)
  {
    for(int i=0; i<locs.size(); i++)
    {
      FeatureLoc loc = (FeatureLoc)locs.get(i);
      if(loc.getFeatureBySrcFeatureId().getFeatureId() == srcfeature_id)
        return loc;
    }
    return null;
  }


  public String getSrcFeatureId()
  {
    return srcFeatureId;
  }


  private JPasswordField getPfield()
  {
    return pfield;
  }


  /**
   * Return true if this looks like a single schema postgres
   * database
   * @return
   */
  public boolean isSingleSchema()
  {
    return singleSchema;
  }
  
  /**
   * Ensure exon featurelocs are in the correct order
   */
  class LocationComarator implements Comparator
  {

    public int compare(Object o1, Object o2)
    {
      int loc1 = ((FeatureLoc)o1).getFmin().intValue();
      int loc2 = ((FeatureLoc)o2).getFmin().intValue();
      
      if(loc2 == loc1)
        return 0;
      int strand = ((FeatureLoc)o1).getStrand().intValue();
      
      if(strand < 0)
      {
        if(loc2 > loc1)
          return 1;
        else
          return -1;
      }
      else
      {
        if(loc2 > loc1)
          return -1;
        else
          return 1;
      }
    } 
  }
  
  class CvTermThread extends Thread 
  {
    private GmodDAO dao;
    CvTermThread(final GmodDAO dao) 
    {
      this.dao = dao;
    }

    public void run() 
    {
      getCvterms(dao);
    }
  }

  public void setRange(Range range)
  {
    this.range = range;
  }


  public Hashtable getIdFeatureStore()
  {
    return idFeatureStore;
  }


  public boolean isLazyFeatureLoad()
  {
    return lazyFeatureLoad;
  }


  public void setLazyFeatureLoad(boolean lazyFeatureLoad)
  {
    this.lazyFeatureLoad = lazyFeatureLoad;
  }
}
