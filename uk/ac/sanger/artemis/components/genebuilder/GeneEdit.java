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

package uk.ac.sanger.artemis.components.genebuilder;

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;


import java.sql.SQLException;
import java.util.Hashtable;
import java.util.List;
import java.util.Vector;


import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JOptionPane;

import javax.swing.JPasswordField;
import javax.swing.JTextField;

import javax.swing.JComboBox;
import javax.swing.JScrollPane;
import javax.swing.JPanel;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.Box;

import uk.ac.sanger.artemis.chado.*;
import uk.ac.sanger.artemis.io.ChadoCanonicalGene;
import uk.ac.sanger.artemis.io.InvalidRelationException;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.ReadFormatException;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.ByteBuffer;


/**
 * Chado data access example code. This searches for features by their
 * uniquename and returns their properties and attributes.
 * 
 * @author tjc
 * 
 */
public class GeneEdit
{
  /** JDBC DAO */
  private JdbcDAO jdbcDAO = null;

  /** iBatis DAO */
  private IBatisDAO connIB = null;

  /** database URL */
  private String location;

  /** password fields */
  private JPasswordField pfield;

  /** <code>List</code> of <code>ChadoFeature</code> objects */
  private List featureList;

  /**
   * 
   * 
   */
  public GeneEdit()
  {
    try
    {
      setLocation();
      final ChadoDAO dao = getDAO();
      showFeatureSearchPanel(dao);
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
  private void showFeatureSearchPanel(final ChadoDAO dao)
      throws java.net.ConnectException, SQLException
  {
    int index = location.indexOf('=') + 1;
    String schema = location.substring(index);

    final List schemas = dao.getSchema();

    Vector v_schemas = new Vector(schemas);
    v_schemas.add(0, "All");

    final JPanel panel = new JPanel(new BorderLayout());
    final JComboBox schema_list = new JComboBox(v_schemas);
    schema_list.setSelectedItem(schema);
    

    JScrollPane jsp = new JScrollPane(schema_list);
    panel.add(jsp, BorderLayout.EAST);

    Box xbox = Box.createHorizontalBox();
    final JTextField gene_text = new JTextField(20);
    gene_text.setText("AfA24A6.005"); //"SPAC212.04c");
    xbox.add(gene_text);
    gene_text.selectAll();
    
    
    JButton findButt = new JButton("RETRIEVE");
    findButt.addActionListener(new ActionListener()
    {

      public void actionPerformed(ActionEvent event)
      {
        String search_gene = gene_text.getText();
        String schema = (String)schema_list.getSelectedItem();
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
          search(search_gene, schema_search, dao);
        }
        catch(SQLException sqlExp)
        {
          JOptionPane.showMessageDialog(null, "SQL Problems...\n"
              + sqlExp.getMessage(), "SQL Error", JOptionPane.ERROR_MESSAGE);
          sqlExp.printStackTrace();
        }
        catch(ReadFormatException e)
        {
          // TODO Auto-generated catch block
          e.printStackTrace();
        }
      }
    });
    xbox.add(findButt);
    xbox.add(Box.createHorizontalGlue());
    panel.add(xbox, BorderLayout.NORTH);
    
    JFrame frame = new JFrame("Feature Search");
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
  public JMenuBar getJMenuBar(final ChadoDAO dao)
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
   * @throws ReadFormatException 
   */
  public void search(final String search_gene, final List schema_search,
                     final ChadoDAO dao) throws SQLException, ReadFormatException
  {
    Hashtable id_store = new Hashtable();
    
    ChadoFeature feature = new ChadoFeature();
    feature.setUniquename(search_gene);
    featureList = dao.getLazyFeature(feature,
                                     schema_search);  
    
    ChadoCanonicalGene chado_gene = new ChadoCanonicalGene();
    
    if(featureList.size() > 1)
      System.err.println("More than one feature found!");
    
    feature = (ChadoFeature)featureList.get(0);
    id_store.put(Integer.toString(feature.getId()), feature.getUniquename());
    
    int src = feature.getFeatureloc().getSrcfeature_id();
    
    ChadoFeature parent = new ChadoFeature();
    parent.setId(src);
    List parentList = dao.getLazyFeature(parent, schema_search);
    parent = (ChadoFeature)parentList.get(0);
    chado_gene.setSeqlen( parent.getLength() );
       
    
    ByteBuffer buff = new ByteBuffer();
    DatabaseDocument.chadoToGFF(feature, null, null, null, null, dao, buff);
    
    //System.out.println(new String(buff.getBytes())); 
    
    GFFStreamFeature gff_gene_feature = new GFFStreamFeature(new String(buff.getBytes()));
    chado_gene.setGene(gff_gene_feature);
    
    // get children of gene
    List relations = feature.getFeatureRelationshipsForObjectId();
    
    for(int i=0;i<relations.size(); i++)
    {
      ChadoFeature transcript = new ChadoFeature();
      transcript.setId( ((ChadoFeatureRelationship)relations.get(i)).getSubject_id() );
      featureList = dao.getLazyFeature(transcript,
                                       schema_search);
      
      transcript = (ChadoFeature)featureList.get(0);
      id_store.put(Integer.toString(transcript.getId()), transcript.getUniquename());
      buff = new ByteBuffer();
      DatabaseDocument.chadoToGFF(transcript, feature.getUniquename(), 
                                  null, null, id_store, dao, buff);
      GFFStreamFeature gff_feature = new GFFStreamFeature(new String(buff.getBytes()));
      
      new uk.ac.sanger.artemis.Feature(gff_feature);
      chado_gene.addTranscript(gff_feature);
      
   
      // get children of transcript - exons and pp
      List transcipt_relations = transcript.getFeatureRelationshipsForObjectId();

      for(int j=0; j<transcipt_relations.size(); j++)
      {
        ChadoFeature child = new ChadoFeature();
        child.setId( ((ChadoFeatureRelationship)transcipt_relations.get(j)).getSubject_id() );
        featureList = dao.getLazyFeature(child,
                                         schema_search);
      
        child = (ChadoFeature)featureList.get(0);
        id_store.put(Integer.toString(child.getId()), child.getUniquename());
        buff = new ByteBuffer();
        DatabaseDocument.chadoToGFF(child, transcript.getUniquename(), 
                                    null, null, id_store, dao, buff);
        
        gff_feature = new GFFStreamFeature(new String(buff.getBytes()));
        new uk.ac.sanger.artemis.Feature(gff_feature);
        
        try
        {
          if(child.getCvterm().getName().equalsIgnoreCase("polypeptide"))
            chado_gene.addProtein(transcript.getUniquename(), child);
          else if(child.getCvterm().getName().equalsIgnoreCase("exon"))
            chado_gene.addExon(transcript.getUniquename(), child);
        }
        catch(InvalidRelationException e)
        {
          // TODO Auto-generated catch block
          e.printStackTrace();
        }

      }
    }
    
    gff_gene_feature.setChadoGene(chado_gene);
    
    GeneBuilderFrame frame = 
      new GeneBuilderFrame(new uk.ac.sanger.artemis.Feature(gff_gene_feature), 
                           null, null, null);
    
  }


  
  /**
   * Get the data access object (DAO).
   * 
   * @return data access object
   */
  private ChadoDAO getDAO() throws java.net.ConnectException, SQLException
  {
    if(System.getProperty("ibatis") == null)
    {
      if(jdbcDAO == null)
        jdbcDAO = new JdbcDAO(location, pfield);
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
   * Set the database location as:
   * jdbc:postgresql://localhost:13001/chadoCVS?user=es2
   * 
   * @return true if location chosen
   */
  protected boolean setLocation()
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

    // given -Dchado=localhost:port/dbname?username
    if(System.getProperty("chado") != null)
    {
      String db_url = System.getProperty("chado").trim();
      int index;
      if((index = db_url.indexOf(":")) > -1)
      {
        inServer.setText(db_url.substring(0, index));
        int index2;
        if((index2 = db_url.indexOf("/")) > -1)
        {
          inPort.setText(db_url.substring(index + 1, index2));
          int index3;
          if((index3 = db_url.indexOf("?")) > -1)
          {
            inDB.setText(db_url.substring(index2 + 1, index3));
            inUser.setText(db_url.substring(index3 + 1));

            /*
             * if(!prompt_user) { location = "jdbc:postgresql://"
             * +inServer.getText().trim()+ ":" +inPort.getText().trim()+ "/"
             * +inDB.getText().trim()+ "?user=" +inUser.getText().trim(); return
             * true; }
             */
          }
        }
      }
    }

    int n = JOptionPane.showConfirmDialog(null, bacross,
        "Enter Database Address", JOptionPane.OK_CANCEL_OPTION,
        JOptionPane.QUESTION_MESSAGE);
    if(n == JOptionPane.CANCEL_OPTION)
      return false;

    location = "jdbc:postgresql://" + inServer.getText().trim() + ":"
        + inPort.getText().trim() + "/" + inDB.getText().trim() + "?user="
        + inUser.getText().trim();

    return true;
  }


  public static void main(String args[])
  {
    new GeneEdit();
  }
}
