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
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.ReadFormatException;

import uk.ac.sanger.artemis.chado.ArtemisUtils;
import uk.ac.sanger.artemis.chado.IBatisDAO;
import uk.ac.sanger.artemis.chado.JdbcDAO;
import uk.ac.sanger.artemis.chado.GmodDAO;
import uk.ac.sanger.artemis.chado.ChadoTransaction;
import uk.ac.sanger.artemis.components.database.DatabaseEntrySource;
import uk.ac.sanger.artemis.components.Splash;

import org.gmod.schema.sequence.Feature;
import org.gmod.schema.sequence.FeatureProp;
import org.gmod.schema.sequence.FeatureLoc;
import org.gmod.schema.sequence.FeatureRelationship;
import org.gmod.schema.sequence.FeatureSynonym;
import org.gmod.schema.sequence.FeatureCvTerm;
import org.gmod.schema.sequence.FeatureCvTermProp;
import org.gmod.schema.sequence.FeatureCvTermDbXRef;
import org.gmod.schema.sequence.FeatureCvTermPub;
import org.gmod.schema.cv.CvTerm;
import org.gmod.schema.general.DbXRef;
import org.gmod.schema.organism.Organism;
import org.gmod.schema.pub.PubDbXRef;
import org.gmod.schema.pub.Pub;

import java.sql.*;
import java.text.SimpleDateFormat;
import java.io.*;
import java.net.ConnectException;
import java.util.Hashtable;
import java.util.HashMap;
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

  private String feature_id = "1";

  /** database schema */
  private String schema = "public";

  private static Hashtable cvterms;

  private InputStreamProgressListener progress_listener;

  private HashMap db;

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
  
  private boolean gene_builder;
  
  // controlled vocabulary
  /** controlled_curation controlled vocabulary */
  public static String CONTROLLED_CURATION_TAG_CVNAME = 
                                 "CC_";
  /** product controlled vocabulary */
  public static String PRODUCTS_TAG_CVNAME = "genedb_products";
  
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
    
    reset(location, schema);
  }

  /**
   * Used by the gene builder to read a database entry
   * for a single gene.
   * @param location
   * @param pfield
   * @param feature_id
   * @param schema
   * @param gene_builder
   */
  public DatabaseDocument(String location, JPasswordField pfield,
          String feature_id, String schema, boolean gene_builder)
  {
    super(location);
    this.pfield = pfield;
    this.feature_id = feature_id;
    this.schema = schema;
    this.gene_builder = gene_builder;

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
      connIB  = null;
      jdbcDAO = null;
      System.setProperty("chado", (String)getLocation());
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
                                 feature_id, schema );
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
      ByteBuffer entry = new ByteBuffer();
      
      try
      {
        if(dao instanceof IBatisDAO)
          ((IBatisDAO) dao).startTransaction();

        // if creating a gene builder
        if(gene_builder)
        {
          List schemaList = new Vector();
          schemaList.add(schema);
          return new ByteArrayInputStream(getGeneFeature(feature_id,
              schemaList, dao).getBytes());
        }

        gff_buffer = getGff(dao, feature_id);

        
        if(splitGFFEntry)
        {
          if(gff_buffer[0].size() > 0)
            entry.append(gff_buffer[0]);

          getChadoSequence(dao, entry);
        }
        else
        {
          for(int i = 0; i < gff_buffer.length; i++)
          {
            if(gff_buffer[i].size() > 0)
              entry.append(gff_buffer[i]);
          }

          getChadoSequence(dao, entry);
        }

        if(dao instanceof IBatisDAO)
          ((IBatisDAO) dao).commitTransaction();
      }
      finally
      {
        if(dao instanceof IBatisDAO)
          ((IBatisDAO) dao).endTransaction();
      }

      instream = new ByteArrayInputStream(entry.getBytes());
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
   */
  private ByteBuffer[] getGff(GmodDAO dao, String parentFeatureID)
  {
    final int srcfeature_id = Integer.parseInt(parentFeatureID);
    
    Splash.logger4j.debug("Build GFF");
    
    // build srcfeature object
    FeatureLoc featureloc = new FeatureLoc();
    Feature srcFeature = new Feature();
    srcFeature.setFeatureId(srcfeature_id);
    featureloc.setFeatureBySrcFeatureId(srcFeature);
    
    //featureloc.setSrcfeature_id(srcfeature_id);
    Feature parent = new Feature();
    parent.setFeatureLoc(featureloc);
    
    List featList = dao.getFeaturesByLocatedOnFeature(parent);

    ByteBuffer[] buffers = new ByteBuffer[types.length + 1];
    for(int i = 0; i < buffers.length; i++)
      buffers[i] = new ByteBuffer();

    final Feature parentFeature = dao.getFeatureById(srcfeature_id);
    
    ByteBuffer this_buff;

    int feature_size = featList.size();
    Hashtable id_store = new Hashtable(feature_size);

    // build feature name store
    for(int i = 0; i < feature_size; i++)
    {
      Feature feat = (Feature)featList.get(i);
      String name       = feat.getUniqueName();
      String feature_id = Integer.toString(feat.getFeatureId());

      id_store.put(feature_id, name);
    }
    
    // get all dbrefs & synonyms
    Hashtable dbxrefs = IBatisDAO.mergeDbXRef(
                 dao.getFeatureDbXRefsByFeatureUniquename(null));
    Hashtable synonym = getAllFeatureSynonyms(dao, null);

    Hashtable featureCvTerms = getFeatureCvTermsByFeature(dao);
    Hashtable featureCvTermDbXRefs = getFeatureCvTermDbXRef(dao);
    Hashtable featureCvTermPubs = getFeatureCvTermPub(dao);
    
    List pubDbXRefs= dao.getPubDbXRef();

    // create gff byte stream
    for(int i = 0; i < feature_size; i++)
    { 
      // select buffer based on feature type
      Feature feat = (Feature)featList.get(i);
      int type_id = feat.getCvTerm().getCvTermId();
      String typeName = getCvtermName(type_id, dao);
      this_buff = buffers[types.length];
      for(int j = 0; j < types.length; j++)
      {
        if(types[j].equals(typeName))
          this_buff = buffers[j];
      }
      
      
      chadoToGFF(feat, parentFeature.getUniqueName(),
                 dbxrefs, synonym, featureCvTerms,
                 pubDbXRefs, featureCvTermDbXRefs, featureCvTermPubs,
                 id_store, dao, 
                 feat.getFeatureLoc(), this_buff);
       
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
  private Hashtable getAllFeatureSynonyms(final GmodDAO dao, 
          final String uniquename) 
  {
    List list = dao.getFeatureSynonymsByFeatureUniquename(uniquename);  
    
    Hashtable synonym = new Hashtable();
    Integer feature_id;
    List value;
    FeatureSynonym alias;
    
    for(int i=0; i<list.size(); i++)
    {
      alias = (FeatureSynonym)list.get(i);
      feature_id = new Integer(alias.getFeature().getFeatureId());
      if(synonym.containsKey(feature_id))
        value = (Vector)synonym.get(feature_id);
      else
        value = new Vector();
      
      value.add(alias);
      synonym.put(feature_id, value);
    }
    
    return synonym;
  }
  
  private Hashtable getFeatureCvTermsByFeature(final GmodDAO dao)
  {
    List list = dao.getFeatureCvTermsByFeature(null);
    
    Hashtable featureCvTerms = new Hashtable();
    Integer feature_id;
    List value;
    FeatureCvTerm feature_cvterm;
    
    for(int i=0; i<list.size(); i++)
    {
      feature_cvterm = (FeatureCvTerm)list.get(i);
      feature_id = new Integer(feature_cvterm.getFeature().getFeatureId());
      if(featureCvTerms.containsKey(feature_id))
        value = (Vector)featureCvTerms.get(feature_id);
      else
        value = new Vector();
      
      value.add(feature_cvterm);
      featureCvTerms.put(feature_id, value);
    }
    
    return featureCvTerms;
  }
  
  
  private Hashtable getFeatureCvTermDbXRef(final GmodDAO dao)
  {
    List list = dao.getFeatureCvTermDbXRefByFeature(null);
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
  
  private Hashtable getFeatureCvTermPub(final GmodDAO dao)
  {
    List list = dao.getFeatureCvTermPubByFeature(null);
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
                                    GmodDAO dao) 
          throws SQLException, ReadFormatException, ConnectException
  {
    Hashtable id_store = new Hashtable();

    //Feature feature = new Feature();
    //feature.setUniquename(search_gene);

    reset((String)getLocation(), (String)schema_search.get(0));
    dao = getDAO();
    Feature feature = 
      (Feature)dao.getFeatureByUniqueName(search_gene);
    
    ChadoCanonicalGene chado_gene = new ChadoCanonicalGene();
    id_store.put(Integer.toString(feature.getFeatureId()), feature.getUniqueName());

    List featurelocs = new Vector(feature.getFeatureLocsForFeatureId());
    FeatureLoc featureloc = (FeatureLoc) featurelocs.get(0);
    int src_id = featureloc.getFeatureBySrcFeatureId().getFeatureId();

    Feature parent = new Feature();
    parent.setFeatureId(src_id);

    parent = dao.getFeatureById(src_id);  //.getLazyFeature(parent);
    
    chado_gene.setSeqlen(parent.getSeqLen());
    chado_gene.setSrcfeature_id(src_id);

    ByteBuffer buff = new ByteBuffer();
    
    chadoToGFF(feature, null, null, null, null, null, null, null, null, dao,
               featureloc, buff);

    // get children of gene
    List relations = new Vector(feature.getFeatureRelationshipsForObjectId());
    
    for(int i = 0; i < relations.size(); i++)
    {
      //Feature transcript = new Feature();
      
      int id = ((FeatureRelationship) relations.get(i)).getFeatureBySubjectId().getFeatureId();
      //transcript.setId(id);
      Feature transcript = 
        (Feature)dao.getFeatureById(id); //.getLazyFeature(transcript);

      id_store.put(Integer.toString(transcript.getFeatureId()), 
          transcript.getUniqueName());

      FeatureLoc loc = getFeatureLoc(new Vector(
          transcript.getFeatureLocsForFeatureId()), src_id);
      chadoToGFF(transcript, feature.getUniqueName(), null, null, null,
          null, null, null, id_store, dao, loc, buff);

      // get children of transcript - exons and pp
      List transcipt_relations = new Vector(
          transcript.getFeatureRelationshipsForObjectId());

      for(int j = 0; j < transcipt_relations.size(); j++)
      {
        id = ((FeatureRelationship) transcipt_relations.get(j)).getFeatureBySubjectId().getFeatureId();
        //Feature child = new Feature();
        //child.setId(((FeatureRelationship) transcipt_relations.get(j))
        //    .getSubject_id());
        Feature child = 
          (Feature)dao.getFeatureById(id); //dao.getLazyFeature(child);

        id_store.put(Integer.toString(child.getFeatureId()), child.getUniqueName());

        loc = getFeatureLoc(
            new Vector(child.getFeatureLocsForFeatureId()),src_id);
        chadoToGFF(child, transcript.getUniqueName(), null,null, null,
                   null, null, null, id_store, dao, loc, buff);
      }
    }

    return buff;
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
                                 final Hashtable id_store,
                                 final GmodDAO dao,
                                 final FeatureLoc featureloc,
                                 final ByteBuffer this_buff)
  {
    String gff_source = null;
    
    final int fmin           = featureloc.getFmin().intValue() + 1;
    final int fmax           = featureloc.getFmax().intValue();
    final int type_id        = feat.getCvTerm().getCvTermId();
    final Short strand       = featureloc.getStrand();
    final Integer phase      = featureloc.getPhase();
    final String name        = feat.getUniqueName();
    final String typeName    = getCvtermName(type_id, dao);
    final Integer feature_id = new Integer(feat.getFeatureId());
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
      
    if(feat.getFeatureRelationshipsForSubjectId() != null)
    {
      List relations = new Vector(feat.getFeatureRelationshipsForSubjectId());
      for(int i=0; i<relations.size(); i++)
      {
        FeatureRelationship feat_relationship = 
                            (FeatureRelationship)relations.get(i);
        parent_id = Integer.toString(feat_relationship.getFeatureByObjectId().getFeatureId());
        
        if( feat_relationship.getCvTerm().getName() == null )
        {
          int parent_type_id = feat_relationship.getCvTerm().getCvTermId();
          parent_relationship = getCvtermName(parent_type_id, dao);
        }
        else
          parent_relationship = feat_relationship.getCvTerm().getName();
      }
    }
          
    if(parent_id != null && id_store != null &&  id_store.containsKey(parent_id))
      parent_id = (String)id_store.get(parent_id);
 
    // make gff format
    
    Vector dbxref = null;
    // append dbxrefs
    if(dbxrefs != null &&
       dbxrefs.containsKey(feature_id))
    {
      dbxref = (Vector)dbxrefs.get(feature_id);
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
    this_buff.append(typeName + "\t");      // type
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

   
    if(parent_id != null && !parent_id.equals("0"))
    {
      if(parent_relationship.equals("derives_from"))
        this_buff.append("Derives_from=" + parent_id + ";");
      else
        this_buff.append("Parent=" + parent_id + ";");
    }

    this_buff.append("timelastmodified=" + timelastmodified + ";");
    
    // this is the chado feature_relationship.rank used
    // to order features e.g. exons
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
        String qualifier_name = getCvtermName(featprop.getCvTerm().getCvTermId(), dao);
        if(qualifier_name == null)
          continue;
        if(featprop.getValue() != null)
          this_buff.append(qualifier_name+ "=" +
                           GFFStreamFeature.encode(featprop.getValue())+";");
      }
    }

    // append dbxrefs
    if(dbxref != null && dbxref.size() > 0)
    {
      this_buff.append("Dbxref=");
      for(int j=0; j<dbxref.size(); j++)
      {
        this_buff.append((String)dbxref.get(j));
        if(j<dbxref.size()-1)
          this_buff.append(",");
      }
      this_buff.append(";");
    }
    
    // append synonyms
    if(synonym != null &&
       synonym.containsKey(feature_id))
    {   
      FeatureSynonym alias;
      Vector v_synonyms = (Vector)synonym.get(feature_id);
      for(int j=0; j<v_synonyms.size(); j++)
      {
        alias = (FeatureSynonym)v_synonyms.get(j);
        
        this_buff.append( getCvtermName(alias.getSynonym().getCvTerm().getCvTermId(), dao) + "=" );
        //this_buff.append(alias.getSynonym().getCvterm().getName()+"=");
        this_buff.append(alias.getSynonym().getName());
        
        if(j<v_synonyms.size()-1)
          this_buff.append(";");
      }
    }
    
    // GO, controlled_curation, product
    if(featureCvTerms != null && 
       featureCvTerms.containsKey(feature_id))
    {
      FeatureCvTerm feature_cvterm;
      Vector v_feature_cvterms = (Vector)featureCvTerms.get(feature_id);
      for(int j=0; j<v_feature_cvterms.size(); j++)
      {
        feature_cvterm = (FeatureCvTerm)v_feature_cvterms.get(j);
        
        Integer featureCvTermId = new Integer( feature_cvterm.getFeatureCvTermId() );
        
        List featureCvTermDbXRefList = 
            (List)featureCvTermDbXRefs.get(featureCvTermId);
        
        List featureCvTermPubList = (List)featureCvTermPubs.get(featureCvTermId);
          
        appendControlledVocabulary(this_buff, dao, feature_cvterm,
                                   featureCvTermDbXRefList,featureCvTermPubList, pubDbXRefs);
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
                                          final List pubDbXRefs)
  {
    CvTerm cvterm =  getCvTerm( feature_cvterm.getCvTerm().getCvTermId(), dao);
    DbXRef dbXRef = feature_cvterm.getCvTerm().getDbXRef();
    
    if(cvterm.getCv().getName().startsWith(DatabaseDocument.CONTROLLED_CURATION_TAG_CVNAME))
    {
      attr_buff.append("controlled_curation=");
      
      attr_buff.append("term="+feature_cvterm.getCvTerm().getName()+"%3B");
      attr_buff.append("cv="+feature_cvterm.getCvTerm().getCv().getName()+"%3B");   
      
      // N.B. the db_xref may be a FeatureCvTermDbXRef or a Pub for /controlled_curation
      int nfound_dbxref = 0;
      if(feature_cvterm.getPub().getUniqueName() != null &&
         !feature_cvterm.getPub().getUniqueName().equals("NULL"))
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
            .getCvTermId(), dao));
        attr_buff.append("=");
        attr_buff.append(feature_cvtermprop.getValue());
        if(i < feature_cvtermprops.size()-1)
          attr_buff.append("%3B");
      }
      
      attr_buff.append(";");
    }
    else if(cvterm.getCv().getName().equals(DatabaseDocument.PRODUCTS_TAG_CVNAME))
    {
      attr_buff.append("product=");
      attr_buff.append(feature_cvterm.getCvTerm().getName()+";");
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
      
      attr_buff.append("term="+feature_cvterm.getCvTerm().getName()+"%3B");
      
      // PMID
      int nfound_pub = 0;
      if(feature_cvterm.getPub().getUniqueName() != null &&
         !feature_cvterm.getPub().getUniqueName().equals("NULL"))
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
        attr_buff.append(getCvtermName(feature_cvtermprop.getCvTerm()
            .getCvTermId(), dao));
        attr_buff.append("=");
        attr_buff.append(feature_cvtermprop.getValue());
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
  public static Integer getCvtermID(String name)
  {
    Enumeration enum_cvterm = cvterms.keys();
    while(enum_cvterm.hasMoreElements())
    {
      Integer key = (Integer)enum_cvterm.nextElement();
      if(name.equals( ((CvTerm)cvterms.get(key)).getName() ))
        return key;
    }
    return null;
  }

  /**
   * Look up a cvterm name from the collection of cvterms.
   * @param id  a cvterm_id  
   * @return    the cvterm name
   */
  private static String getCvtermName(int id, GmodDAO dao)
  {
    return getCvTerm(id, dao).getName();
  }
  
  private static CvTerm getCvTerm(int id, GmodDAO dao)
  {
    if(cvterms == null)
      getCvterms(dao);

    return (CvTerm)cvterms.get(new Integer(id));
  }
  
  public static CvTerm getCvTermByCvTermName(String cvterm_name)
  {
    Enumeration enum_cvterm = cvterms.elements();
    while(enum_cvterm.hasMoreElements())
    {
      CvTerm cvterm = (CvTerm)enum_cvterm.nextElement();
      if(cvterm_name.equals( cvterm.getName() ))
        return cvterm;
    }
    
    return null;
  }
  
  public static CvTerm getCvTermByCvAndCvTerm(final String cvterm_name,
                                              final String cvName)
  {
    Enumeration enum_cvterm = cvterms.elements();
    while(enum_cvterm.hasMoreElements())
    {
      CvTerm cvterm = (CvTerm)enum_cvterm.nextElement();
      if(cvName.equals( cvterm.getCv().getName() ) &&
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
  private static Hashtable getCvterms(GmodDAO dao)
  {
    cvterms = new Hashtable();

    try
    {
      List cvterm_list = dao.getCvTerms();
      Iterator it = cvterm_list.iterator();

      while(it.hasNext())
      {
        CvTerm cvterm = (CvTerm)it.next();
        cvterms.put(new Integer(cvterm.getCvTermId()), cvterm);
      }
    }
    catch(RuntimeException sqle)
    {
      System.err.println("SQLException retrieving CvTerms");
      System.err.println(sqle);
    }

    return cvterms;
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
   * Look up synonym type names e.g. synonym, systematic_id.
   * @return    the synonym tag names
   */
  public static String[] getSynonymTypeNames(String cv_name)
  {
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
  
  
  /**
   * Get the sequence for a feature.
   * @param dao   the data access object
   * @param buff  the buffer to add the sequence to
   * @return      the resulting buffer
   * @throws java.sql.SQLException
   */
  private ByteBuffer getChadoSequence(GmodDAO dao, ByteBuffer buff)
  {
    Feature feature = dao.getFeatureById(Integer.parseInt(feature_id));
 
    buff.append("##FASTA\n>");
    buff.append(feature.getUniqueName());
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
  public HashMap getDatabaseEntries2()
                   throws ConnectException, java.sql.SQLException
  {
    db = new HashMap();
 
    GmodDAO dao = null;
    
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
      if(dao instanceof IBatisDAO)
        ((IBatisDAO) dao).startTransaction();
      
      schema_list = dao.getSchema();
      Iterator it = schema_list.iterator();

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
          Feature feature = (Feature)it_residue_features.next();
          String typeName = getCvtermName(feature.getCvTerm().getCvTermId(), getDAO()); 
          
          db.put(schema + " - " + typeName + " - " + feature.getUniqueName(),
                 Integer.toString(feature.getFeatureId()));
        }
      }
      
      if(dao instanceof IBatisDAO)
        ((IBatisDAO) dao).commitTransaction();
      
    }
    catch(RuntimeException sqlExp)
    {
      JOptionPane.showMessageDialog(null, "SQL Problems...\n"+
                                    sqlExp.getMessage(), 
                                    "SQL Error",
                                    JOptionPane.ERROR_MESSAGE);
      sqlExp.printStackTrace();
    }
    finally
    {
      if(dao instanceof IBatisDAO)
        ((IBatisDAO) dao).endTransaction();
    }
    
    return db;
  }
  
  
  /**
   * Create a hashtable of the available entries with residues.
   * @return a <code>Hashtable</code> of the <code>String</code>
   *          representation (schema-type-feature_name) and the
   *          corresponding feature_id
   * @throws ConnectException
   * @throws java.sql.SQLException
   */
  public HashMap getDatabaseEntries()
                   throws ConnectException, java.sql.SQLException
  {
    String schema = null;

    try
    { 
      GmodDAO dao = null;
      dao = getDAO();
      schema_list = dao.getOrganisms();
      
      
/*      Organism o = new Organism();
      o.setCommonName("web");
      schema_list.add(o);*/
      
      
      Iterator it = schema_list.iterator();
      db = new HashMap();
      
      while(it.hasNext())
      {
        Organism organism = (Organism)it.next();
        schema = organism.getCommonName();
        
        reset((String)getLocation(),  schema);

        try
        {
          dao = getDAO();
          final List list = dao.getResidueType(schema);
          
          if(list.size() == 0)  // no residues for this organism
            continue;

          List list_residue_features = dao.getResidueFeatures(list, schema);
          Iterator it_residue_features = list_residue_features.iterator();
          while(it_residue_features.hasNext())
          {
            Feature feature = (Feature)it_residue_features.next();
            String typeName = getCvtermName(feature.getCvTerm().getCvTermId(), getDAO()); 
          
            db.put(schema + " - " + typeName + " - " + feature.getUniqueName(),
                   Integer.toString(feature.getFeatureId()));
          }
        }
        catch(RuntimeException e)
        { 
        }
        catch(java.sql.SQLException sqlExp)
        {
        }
      }
      
    }
    catch(RuntimeException sqlExp)
    {
      JOptionPane.showMessageDialog(null, "SQL Problems...\n"+
                                    sqlExp.getMessage(), 
                                    "SQL Error",
                                    JOptionPane.ERROR_MESSAGE);
      sqlExp.printStackTrace();
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
      JOptionPane.showMessageDialog(null, "SQL Problems....\n"+
                                    sqlExp.getMessage(), 
                                    "SQL Error",
                                    JOptionPane.ERROR_MESSAGE);
      throw sqlExp;
    }
    
    return db;
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
      
      GmodDAO dao = getDAO();
      boolean unchanged;
      
      //
      // check feature timestamps have not changed
      Vector names_checked = new Vector();
      for(i = 0; i < sql.size(); i++)
      {
        ChadoTransaction tsn = (ChadoTransaction)sql.get(i);
        if( (tsn.getType() == ChadoTransaction.INSERT ||
             tsn.getType() == ChadoTransaction.DELETE) && 
             tsn.getFeatureObject() instanceof Feature )
          continue;
            
        final String uniquename = tsn.getUniquename();
        
        if(uniquename == null)
          continue;
        if(names_checked.contains(uniquename))
          continue;
            
        names_checked.add(uniquename);
        
        unchanged = checkFeatureTimestamp(schema, 
                         uniquename, 
                         tsn.getLastModified(), dao);
        if(!unchanged)
          return 0;
      }  
      
      try
      {
        if(dao instanceof IBatisDAO)
          ((IBatisDAO) dao).startTransaction();

        //
        // commit to database
        for(i = 0; i < sql.size(); i++)
        {
          ChadoTransaction tsn = (ChadoTransaction) sql.get(i);

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
                    = dao.getFeatureByUniqueName(uniquename);
                feature.setFeatureId( old_feature.getFeatureId() );
                
                tsn.setOldUniquename(feature.getUniqueName());
              }
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
              if(tsn.getFeatureObject() instanceof Feature)
              {
                FeatureLoc featureloc = ((Feature) tsn.getFeatureObject()).getFeatureLoc();
                Feature featureBySrcFeatureId = new Feature();
                featureBySrcFeatureId.setFeatureId(Integer.parseInt(feature_id));
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

        //
        // update timelastmodified timestamp
        Timestamp ts = new Timestamp(new java.util.Date().getTime());
        names_checked = new Vector();
        for(int j = 0; j < sql.size(); j++)
        {
          ChadoTransaction tsn = (ChadoTransaction)sql.get(j);
          
          if( (tsn.getType() == ChadoTransaction.INSERT ||
              tsn.getType() == ChadoTransaction.DELETE) && 
              tsn.getFeatureObject() instanceof Feature )
           continue;
          
          final String uniquename = tsn.getUniquename();
          if(uniquename == null)
            continue;
            
          if(names_checked.contains(uniquename))
            continue;

          names_checked.add(uniquename);

          Feature feature = dao.getFeatureByUniqueName(uniquename);
          if(feature != null)
          {
            feature.setTimeLastModified(ts);
            dao.merge(feature);
          }

          GFFStreamFeature gff_feature = (GFFStreamFeature) tsn
              .getGff_feature();
          gff_feature.setLastModified(ts);
        }
        
        if(dao instanceof IBatisDAO && 
           System.getProperty("nocommit") == null)
           ((IBatisDAO) dao).commitTransaction();
      }
      finally
      {
        if(dao instanceof IBatisDAO)
          ((IBatisDAO) dao).endTransaction();
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
   */
  public boolean checkFeatureTimestamp(final String schema,
                                       final String uniquename,
                                       final Timestamp timestamp,
                                       final GmodDAO dao)
  {
    Feature feature = dao.getFeatureByUniqueName(uniquename);
    if(feature == null)
      return true;
    Timestamp now = feature.getTimeLastModified();
    
    if(now != null && timestamp != null)
    {
      now.setNanos(0);
      timestamp.setNanos(0);
      
      if(now.compareTo(timestamp) != 0)
      {
        SimpleDateFormat date_format = 
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
      featureList.add(dao.getFeatureByUniqueName(args[0])); 
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
}
