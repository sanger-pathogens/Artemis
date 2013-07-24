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

import org.gmod.schema.cv.CvTerm;
import org.gmod.schema.general.DbXRef;
import org.gmod.schema.organism.Organism;
import org.gmod.schema.pub.Pub;
import org.gmod.schema.sequence.Feature;
import org.gmod.schema.sequence.FeatureCvTerm;
import org.gmod.schema.sequence.FeatureCvTermDbXRef;
import org.gmod.schema.sequence.FeatureCvTermProp;
import org.gmod.schema.sequence.FeatureCvTermPub;
import org.gmod.schema.sequence.FeatureDbXRef;
import org.gmod.schema.sequence.FeaturePub;
import org.gmod.schema.sequence.FeatureSynonym;
import org.gmod.schema.sequence.FeatureProp;
import org.gmod.schema.sequence.FeatureLoc;

//import uk.ac.sanger.artemis.io.GFFStreamFeature;
//import uk.ac.sanger.artemis.util.DatabaseDocument;

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.Cursor;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.Dimension;
import java.net.ConnectException;
import java.sql.SQLException;
import java.util.Iterator;
import java.util.List;
import java.util.Collection;
import java.util.Vector;

import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JButton;
import javax.swing.JTabbedPane;
import javax.swing.ListSelectionModel;
import javax.swing.JTable;
import javax.swing.JPasswordField;
import javax.swing.JTextField;
import javax.swing.JTextArea;
import javax.swing.JList;
import javax.swing.JScrollPane;
import javax.swing.JPanel;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.Box;
import javax.swing.event.ListSelectionListener;
import javax.swing.event.ListSelectionEvent;

import uk.ac.sanger.artemis.util.ByteBuffer;
import uk.ac.sanger.artemis.util.DatabaseLocationParser;

/**
 * Chado data access example code. This searches for features by their
 * uniquename and returns their properties and attributes.
 * 
 * @author tjc
 * 
 */
public class ChadoDemo
{
  /** JDBC DAO */
  private JdbcDAO jdbcDAO = null;

  /** iBatis DAO */
  private IBatisDAO connIB = null;

  /** database URL */
  private String location;

  /** password fields */
  private JPasswordField pfield;

  /** results table */
  private JTable result_table;

  /** feature attributes */
  private JTextArea attr_text;

  /** <code>List</code> of <code>Feature</code> objects */
  private List featureList;
  
  /** row data containing results */
  private String rowData[][];
  
  private List pubDbXRefs[];
  
  private JTabbedPane tabbedPane;
  
  private static org.apache.log4j.Logger logger4j = 
	    org.apache.log4j.Logger.getLogger(ChadoDemo.class);

  /** controlled_curation controlled vocabulary */
  public static String CONTROLLED_CURATION_TAG_CVNAME = 
                                 "CC_";
  /** product controlled vocabulary */
  public static String PRODUCTS_TAG_CVNAME = "genedb_products";
  public static String RILEY_TAG_CVNAME = "RILEY";
  
  /** 
   * Chado demo
   */
  public ChadoDemo()
  {
    //uk.ac.sanger.artemis.components.Splash.initLogger();
    
    try
    {
      setLocation();
      final GmodDAO dao = getDAO();
      
      //
      // TESTING - updating feature.residues
      //
      /*((IBatisDAO)dao).startTransaction();
      FeatureForUpdatingResidues chadoFeature = new FeatureForUpdatingResidues();
      chadoFeature.setStartBase(10);
      chadoFeature.setEndBase(20);
      dao.merge(chadoFeature);
      ((IBatisDAO)dao).endTransaction();*/
      //
      //
      
      showFeatureSearchPanel(dao);
      //getCvterm(dao);
    }
    catch(java.net.ConnectException ce)
    {
      ce.printStackTrace();
    }
    catch(SQLException sqlExp)
    {
      JOptionPane.showMessageDialog(null, "SQL Problems...\n"
          + sqlExp.getMessage(), "SQL Error", JOptionPane.ERROR_MESSAGE);
      sqlExp.printStackTrace();
    }

  }

  /**
   * Display a window for searching for features.
   * 
   * @throws java.net.ConnectException
   * @throws SQLException
   */
  private void showFeatureSearchPanel(final GmodDAO dao)
      throws java.net.ConnectException, SQLException
  {
    int index = location.indexOf('=') + 1;
    String schema = location.substring(index);
    
    final List schemas = dao.getOrganismsContainingSrcFeatures();

    Vector v_schemas = new Vector();
    v_schemas.add(0, "All");
    
    for(int i=0; i<schemas.size(); i++)
    {
      Organism o = (Organism)schemas.get(i);
      v_schemas.add(o.getCommonName());
    }

    final JFrame frame = new JFrame("Feature Search");
    final JPanel panel = new JPanel(new BorderLayout());
    final JList schema_list = new JList(v_schemas);
    schema_list.setSelectedValue(schema, true);
    if(schema_list.getSelectedIndex() == -1)
      schema_list.setSelectedValue("All", true);

    Box xbox2 = Box.createHorizontalBox();
    JScrollPane jsp = new JScrollPane(schema_list);

    Box xbox = Box.createHorizontalBox();
    final JTextField gene_text = new JTextField("PFA0005w*",20);
    xbox.add(gene_text);
    gene_text.selectAll();

    result_table = new JTable();

    final JScrollPane scrollpane = new JScrollPane(result_table);
    scrollpane.setPreferredSize(new Dimension(600, 250));

    //panel.add(scrollpane, BorderLayout.CENTER);
    pubDbXRefs = new List[schemas.size()];
    
    JButton findButt = new JButton("FIND");
    findButt.addActionListener(new ActionListener()
    {
      private String columnNames[] = { "schema", "name", "type",
          "feature ID", "loc", "strand", "time modified" };

      public void actionPerformed(ActionEvent event)
      {
        frame.setCursor(new Cursor(Cursor.WAIT_CURSOR));
        String search_gene = gene_text.getText();
        String schema = (String)schema_list.getSelectedValue();
        List schema_search;
        if(schema.equalsIgnoreCase("All"))
          schema_search = schemas;
        else
        {
          schema_search = new Vector();
          schema_search.add(schema);
        }

        try
        {
          rowData = search(search_gene, schema_search, dao);

          result_table = new JTable(rowData, columnNames);
          result_table.getSelectionModel().addListSelectionListener(
                                         new SelectionListener(dao, frame));
          result_table.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);

          result_table.addMouseListener(new MouseAdapter()
          {
            public void mouseClicked(MouseEvent e)
            {
              frame.setCursor(new Cursor(Cursor.WAIT_CURSOR));
              int row = result_table.getSelectedRow();
 
              if(pubDbXRefs[row] == null)
                pubDbXRefs[row] = dao.getPubDbXRef();
              showAttributes(dao, pubDbXRefs[row]);
              frame.setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
            }
          });

          scrollpane.setViewportView(result_table);
        }
        catch(SQLException sqlExp)
        {
          JOptionPane.showMessageDialog(null, "SQL Problems...\n"
              + sqlExp.getMessage(), "SQL Error", JOptionPane.ERROR_MESSAGE);
          sqlExp.printStackTrace();
        }
        catch(ConnectException e)
        {
          // TODO Auto-generated catch block
          e.printStackTrace();
        }
        frame.setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
      }
    });
    xbox.add(findButt);
    xbox.add(Box.createHorizontalGlue());
    
    Box ybox = Box.createVerticalBox();
    xbox2.add(scrollpane);
    xbox2.add(jsp);
    xbox2.add(Box.createHorizontalGlue());
    
    ybox.add(xbox);
    ybox.add(xbox2);
    
    panel.add(ybox, BorderLayout.NORTH);
    
    attr_text = new JTextArea();
    JScrollPane jsp_attr = new JScrollPane(attr_text);
    jsp_attr.setPreferredSize(new Dimension(600, 150));
    
    
    tabbedPane = new JTabbedPane();
    tabbedPane.add("Core", jsp_attr);
    
    panel.add(tabbedPane, BorderLayout.CENTER);

    
    frame.getContentPane().add(panel);
    frame.setJMenuBar(getJMenuBar(dao));
    frame.pack();
    frame.setVisible(true);
  }

  /**
   * Build a <code>JMenuBar</code>.
   * 
   * @return a <code>JMenuBar</code>
   */
  public JMenuBar getJMenuBar(final GmodDAO dao)
  {
    JMenuBar mbar = new JMenuBar();
    JMenu file = new JMenu("File");
    mbar.add(file);

    JMenuItem exit = new JMenuItem("Exit");
    exit.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)
      {
        System.exit(0);
      }
    });
    file.add(exit);

    return mbar;
  }

 

  
  /**
   * Search for a feature/gene name from a <code>List</code> of schemas.
   * 
   * @param search_gene
   *          the feature name
   * @param schema_search
   *          the <code>List</code> to search
   * @param dao
   *          the data access object
   * @return string array of results
   * @throws SQLException
   * @throws ConnectException
   */
  public String[][] search(final String search_gene, final List schema_search,
      GmodDAO dao) throws SQLException, ConnectException
  {
    final String search_name = search_gene.replaceAll("[*]","%");
    Feature feature = new Feature();
    //feature.setUniquename(search_gene.replaceAll("[*]","%"));
    featureList = new Vector();
    
    for(int i=0; i<schema_search.size(); i++)
    {
      featureList.addAll(dao.getFeaturesByAnyCurrentName(search_name));
    }
    
    String rowData[][] = new String[featureList.size()][7];

    for(int i = 0; i < featureList.size(); i++)
    {
      feature = (Feature) featureList.get(i);
      
      // assume only one featureloc
      Collection locs = feature.getFeatureLocsForFeatureId();
      
      if(locs != null && locs.size() > 0)
      {
        Iterator it = locs.iterator();
        FeatureLoc loc = (FeatureLoc)it.next();
        int fmin = loc.getFmin().intValue() + 1;
        int fmax = loc.getFmax().intValue();
        rowData[i][4] = fmin + "..." + fmax;
        rowData[i][5] = Integer.toString(loc.getStrand().shortValue());
      }
      else if(feature.getFeatureLoc() != null)
      {
        FeatureLoc loc = feature.getFeatureLoc();
        int fmin = loc.getFmin().intValue() + 1;
        int fmax = loc.getFmax().intValue();
        rowData[i][4] = fmin + "..." + fmax;
        rowData[i][5] = Integer.toString(loc.getStrand().shortValue());
      }
      
      String schema = feature.getOrganism().getAbbreviation();
      int ind = schema.indexOf('.');
      if(ind > 0)
        schema = schema.substring(0,ind)+schema.substring(ind+1);
      
      System.out.println("\n\nNow get feature type_id.......\n\n");
      rowData[i][0] = schema;
      rowData[i][1] = feature.getUniqueName();
      rowData[i][2] = feature.getCvTerm().getName();
      //rowData[i][2] = (String)cvterm.get(new Long(feature.getType_id()));
      rowData[i][3] = Integer.toString(feature.getFeatureId());
      rowData[i][6] = feature.getTimeLastModified().toString();
    
    }
    return rowData;
  }

 
  /**
   * Get the data access object (DAO).
   * 
   * @return data access object
   */
  private GmodDAO getDAO() throws java.net.ConnectException, SQLException
  {
    if(System.getProperty("ibatis") == null)
    {
      logger4j.debug("Using JDBC");
      if(jdbcDAO == null)
        jdbcDAO = new JdbcDAO(location, pfield);
      return jdbcDAO;
    }
    else
    {
      logger4j.debug("Using iBatis");
      if(connIB == null)
        connIB = new IBatisDAO(pfield);
      return connIB;
    }
  }

  /**
   * Set the database location as:
   * jdbc:postgresql://localhost:13001/chadoCVS?user=es2
   * 
   * @return true if location chosen
   */
  private boolean setLocation()
  {
    Container bacross = new Container();
    bacross.setLayout(new GridLayout(6, 2, 5, 5));

    JLabel lServer = new JLabel("Server : ");
    bacross.add(lServer);
    JTextField inServer = new JTextField("localhost");
    bacross.add(inServer);

    JLabel lPort = new JLabel("Port : ");
    bacross.add(lPort);
    JTextField inPort = new JTextField("5432");
    bacross.add(inPort);

    JLabel lDB = new JLabel("Database : ");
    bacross.add(lDB);
    JTextField inDB = new JTextField("chado");
    bacross.add(inDB);

    JLabel lUser = new JLabel("User : ");
    bacross.add(lUser);
    JTextField inUser = new JTextField("afumigatus");
    bacross.add(inUser);

    JLabel lpasswd = new JLabel("Password : ");
    bacross.add(lpasswd);
    pfield = new JPasswordField(16);
    bacross.add(pfield);

    
    DatabaseLocationParser dlp = new DatabaseLocationParser();
    // given -Dchado=localhost:port/dbname?username
    if(System.getProperty("chado") != null)
    {
      String db_url = System.getProperty("chado").trim();
      dlp.setFromURLString(db_url);
      inServer.setText(dlp.getHost());
      inPort.setText("" + dlp.getPort());
      inDB.setText(dlp.getDatabase());
      inUser.setText(dlp.getUsername());
      
    }

    int n = JOptionPane.showConfirmDialog(null, bacross,
        "Enter Database Address", JOptionPane.OK_CANCEL_OPTION,
        JOptionPane.QUESTION_MESSAGE);
    if(n == JOptionPane.CANCEL_OPTION)
      return false;

    dlp.setHost(inServer.getText());
    dlp.setPort(inPort.getText());
    dlp.setDatabase(inDB.getText());
    dlp.setUsername(inUser.getText());
  
    location = dlp.getCompleteURL();

    System.setProperty("chado", location);
    
    return true;
  }
  
  
  public class SelectionListener implements ListSelectionListener
  {
	  private GmodDAO dao;
    private JFrame frame;
    
    public SelectionListener(final GmodDAO dao, final JFrame frame)
    {
      super();
      this.dao = dao;
      this.frame = frame;
    }
    
    public void valueChanged(ListSelectionEvent e)
    {
      frame.setCursor(new Cursor(Cursor.WAIT_CURSOR));
      int row = result_table.getSelectedRow();
      //reset(location, rowData[row][0]);

      if(pubDbXRefs[row] == null)
        pubDbXRefs[row] = dao.getPubDbXRef();
      showAttributes(dao, pubDbXRefs[row]);
      frame.setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
    }
  }
  
  /**
   * Show the attributes of a selected feature.
   */
  private void showAttributes(final GmodDAO dao,
                              final List pubDbXRefs)
  {
    final int row = result_table.getSelectedRow();
    final ByteBuffer attr_buff = new ByteBuffer();
    final Feature chado_feature = (Feature)featureList.get(row);
    
    Collection locs = chado_feature.getFeatureLocsForFeatureId();
    FeatureLoc loc;
    
    if(locs != null && locs.size() > 0)
    {
      Iterator it = locs.iterator();
      loc = (FeatureLoc)it.next();
    }
    else
      loc = chado_feature.getFeatureLoc();
    /*
    int fmin = loc.getFmin().intValue() + 1;
    int fmax = loc.getFmax().intValue();
    */
    String uniquename = chado_feature.getUniqueName();
    
    attr_buff.append("/ID="+uniquename+"\n");
    List attributes = (List)chado_feature.getFeatureProps();
    List dbxrefs = dao.getFeatureDbXRefsByFeatureUniquename(uniquename);
    List featureCvTerms       = dao.getFeatureCvTermsByFeature(chado_feature);
    List featureCvTermDbXRefs = dao.getFeatureCvTermDbXRefByFeature(chado_feature);
    List featureCvTermPubs    = dao.getFeatureCvTermPubByFeature(chado_feature);
    Feature srcFeature = new Feature();
    srcFeature.setFeatureId(loc.getSrcFeatureId());
    List featurePubs = dao.getFeaturePubsBySrcFeature(srcFeature);

    if(dbxrefs.size() > 0)
    {
      attr_buff.append("/Dbxref=");
      for(int i = 0; i < dbxrefs.size(); i++)
      {
        FeatureDbXRef dbxref = (FeatureDbXRef) dbxrefs.get(i);
        attr_buff.append(dbxref.getDbXRef().getDb().getName() + ":"
            + dbxref.getDbXRef().getAccession() + "; ");
      }
      attr_buff.append("\n");
    }

    Collection synonyms = chado_feature.getFeatureSynonyms();

    // append synonyms
    if(synonyms != null && synonyms.size() > 0)
    {
      FeatureSynonym alias;

      System.out.println("\n\nNow get synonym & type_id.......\n\n");
      Iterator it = synonyms.iterator();
      while(it.hasNext())
      {
        alias = (FeatureSynonym)it.next();
        attr_buff.append("/");
        attr_buff.append(alias.getSynonym().getCvTerm().getName() + "=");
        attr_buff.append(alias.getSynonym().getName());
        attr_buff.append(";");
        attr_buff.append("\n");
      }
    }

    if(attributes != null)
    {
      Iterator it = attributes.iterator();
      while(it.hasNext())
      {
        FeatureProp featprop = (FeatureProp)it.next();

        attr_buff.append("/" + featprop.getCvTerm().getName() + "="
            + decode(featprop.getValue()) + "\n");
      }
    }
    
    if(featureCvTerms != null)
    {
      for(int j=0; j<featureCvTerms.size(); j++)
      {
        attr_buff.append("/");
        FeatureCvTerm feature_cvterm = (FeatureCvTerm)featureCvTerms.get(j);

        appendControlledVocabulary(attr_buff, dao, feature_cvterm,
            featureCvTermDbXRefs, featureCvTermPubs, pubDbXRefs, false);
        
        attr_buff.append("\n");
      }
    }
    
    // /literature
    if(featurePubs != null)
    {
      for(int i=0; i<featurePubs.size(); i++)
      {
        FeaturePub featurePub = (FeaturePub)featurePubs.get(i);
        if(chado_feature.getFeatureId() == featurePub.getFeature().getFeatureId())
        {
          attr_buff.append( "/literature=" );
          attr_buff.append(featurePub.getPub().getUniqueName());
          attr_buff.append("\n");
        }
      }
    }
    
    attr_text.setText(decode((new String(attr_buff.getBytes()))));
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
    CvTerm cvterm = dao.getCvTermById( feature_cvterm.getCvTerm().getCvTermId() );
    DbXRef dbXRef = feature_cvterm.getCvTerm().getDbXRef();
    
    if(cvterm.getCv().getName().startsWith(CONTROLLED_CURATION_TAG_CVNAME))
    {
      attr_buff.append("controlled_curation=");
      
      attr_buff.append("term="+
          feature_cvterm.getCvTerm().getName()+";");
      attr_buff.append("cv="+
          feature_cvterm.getCvTerm().getCv().getName()+";");   
      
      // N.B. the db_xref may be a FeatureCvTermDbXRef or a Pub for /controlled_curation
      int nfound_dbxref = 0;
      if(feature_cvterm.getPub().getUniqueName() != null &&
         !feature_cvterm.getPub().getUniqueName().equalsIgnoreCase("NULL"))
      {
        // PMID
        Pub pub = feature_cvterm.getPub();
        
        // internal check
        //checkPubDbXRef(pubDbXRefs, pub.getPubId(), pub, feature_cvterm);
        
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
        attr_buff.append( dao.getCvTermById( feature_cvtermprop.getCvTerm().getCvTermId() ).getName());
        attr_buff.append("=");
        attr_buff.append(feature_cvtermprop.getValue());
        if(i < feature_cvtermprops.size()-1)
          attr_buff.append(";");
      }
      
      attr_buff.append(";");
    }
    else if(cvterm.getCv().getName().equals(PRODUCTS_TAG_CVNAME))
    {
      attr_buff.append("product=");
      attr_buff.append(feature_cvterm.getCvTerm().getName()+";");
    }
    else if(cvterm.getCv().getName().equals(RILEY_TAG_CVNAME))
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
          feature_cvterm.getCvTerm().getName()+";");
      
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
        
        attr_buff.append(dao.getCvTermById(feature_cvtermprop.getCvTerm().getCvTermId()).getName());
        attr_buff.append("=");
        attr_buff.append(feature_cvtermprop.getValue());
        if(i < feature_cvtermprops.size()-1)
          attr_buff.append(";");
      }
      
      attr_buff.append(";");
    }
  }
  
  
  
  /**
  *
  * For gff-version 3:
  * http://song.sourceforge.net/gff3-jan04.shtml
  *
  * Remove URL escaping rule (e.g. space="%20" or "+")
  *
  */
  private String decode(String s)
  {
    final String map[][] = {
                             { " ",  "%20" },  // white space
                             { ",",  "%2C" },  // comma
                             { ";",  "%3B" },  // semi-colon
                             { "=",  "%3D" },  // equals
                             { "\t", "%09" },  // tab
                             { " ",  "+"   },  // white space
                             { "+",  "%2B" },
                             { "(",  "%28" },  // left bracket
                             { ")",  "%29" }   // right bracket
                           };

    int ind;
    String enc;
    String dec;

    for(int i=0; i<map.length; i++)
    {
      enc = map[i][1];
      dec = map[i][0];
      while( (ind = s.indexOf(enc)) > -1)
        s = s.substring(0,ind) + dec + s.substring(ind+enc.length());
    }

    return s;
  }
  
  public static void main(String args[])
  {
    new ChadoDemo();
  }
}
