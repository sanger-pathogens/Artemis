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

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.Dimension;
import java.net.ConnectException;
import java.sql.SQLException;
import java.util.List;
import java.util.Vector;

import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JButton;
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
import uk.ac.sanger.artemis.io.GFFStreamFeature;

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

  /**
   * 
   * 
   */
  public ChadoDemo()
  {
    try
    {
      setLocation();
      final ChadoDAO dao = getDAO();
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
  private void showFeatureSearchPanel(final ChadoDAO dao)
      throws java.net.ConnectException, SQLException
  {
    int index = location.indexOf('=') + 1;
    String schema = location.substring(index);

    final List schemas = dao.getSchema();

    Vector v_schemas = new Vector(schemas);
    v_schemas.add(0, "All");

    final JPanel panel = new JPanel(new BorderLayout());
    final JList schema_list = new JList(v_schemas);
    schema_list.setSelectedValue(schema, true);
    if(schema_list.getSelectedIndex() == -1)
      schema_list.setSelectedValue("All", true);

    JScrollPane jsp = new JScrollPane(schema_list);
    panel.add(jsp, BorderLayout.EAST);

    Box xbox = Box.createHorizontalBox();
    final JTextField gene_text = new JTextField(20);
    gene_text.setText("MAL8P1.204");
    xbox.add(gene_text);
    gene_text.selectAll();

    result_table = new JTable();

    final JScrollPane scrollpane = new JScrollPane(result_table);
    scrollpane.setPreferredSize(new Dimension(600, 200));

    panel.add(scrollpane, BorderLayout.CENTER);

    JButton findButt = new JButton("FIND");
    findButt.addActionListener(new ActionListener()
    {
      private String columnNames[] = { "schema", "name", "type",
          "feature ID", "location", "strand", "time modified" };

      public void actionPerformed(ActionEvent event)
      {
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
          String rowData[][] = search(search_gene, schema_search, dao);

          result_table = new JTable(rowData, columnNames);
          result_table.getSelectionModel().addListSelectionListener(new SelectionListener());
          result_table.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);

          result_table.addMouseListener(new MouseAdapter()
          {
            public void mouseClicked(MouseEvent e)
            {
              int row = result_table.getSelectedRow();
              showAttributes();
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
      }
    });
    xbox.add(findButt);
    xbox.add(Box.createHorizontalGlue());
    panel.add(xbox, BorderLayout.NORTH);
    
    attr_text = new JTextArea();
    JScrollPane jsp_attr = new JScrollPane(attr_text);
    jsp_attr.setPreferredSize(new Dimension(600, 100));
    panel.add(jsp_attr, BorderLayout.SOUTH);

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
   * Show the attributes of a selected feature.
   */
  private void showAttributes()
  {
    int row = result_table.getSelectedRow();
    StringBuffer attr_buff = new StringBuffer();
    Feature chado_feature = (Feature)featureList.get(row);
    List attributes = chado_feature.getFeaturepropList();
    List dbxrefs    = chado_feature.getFeatureDbxrefs();

    if(dbxrefs.size() > 0)
    {
      attr_buff.append("/Dbxref=");
      for(int i = 0; i < dbxrefs.size(); i++)
      {
        FeatureDbxref dbxref = (FeatureDbxref) dbxrefs.get(i);
        attr_buff.append(dbxref.getDbxref().getDb().getName() + ":"
            + dbxref.getDbxref().getAccession() + "; ");
      }
      attr_buff.append("\n");
    }

    List synonyms = chado_feature.getFeatureSynonymsForFeatureId();

    // append synonyms
    if(synonyms.size() > 0)
    {
      FeatureSynonym alias;

      for(int i = 0; i < synonyms.size(); i++)
      {
        alias = (FeatureSynonym) synonyms.get(i);
        attr_buff.append("/");
        attr_buff.append(alias.getSynonym().getCvterm().getName() + "=");
        attr_buff.append(alias.getSynonym().getName());

        attr_buff.append(";");
        attr_buff.append("\n");
      }
    }

    if(attributes != null)
      for(int i = 0; i < attributes.size(); i++)
      {
        FeatureProp featprop = (FeatureProp) attributes.get(i);

        attr_buff.append("/" + featprop.getCvterm().getName() + "="
            + GFFStreamFeature.decode(featprop.getValue()) + "\n");
      }

    attr_text.setText(new String(attr_buff));
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
      ChadoDAO dao) throws SQLException, ConnectException
  {
    final String search_name = search_gene.replaceAll("[*]","%");
    Feature feature = new Feature();
    //feature.setUniquename(search_gene.replaceAll("[*]","%"));
    featureList = new Vector();
    
    for(int i=0; i<schema_search.size(); i++)
    {
      reset(location, (String)schema_search.get(i));
      dao = getDAO();
      featureList.addAll(dao.getFeaturesByAnyCurrentName(search_name));
    }
    
    String rowData[][] = new String[featureList.size()][7];

    for(int i = 0; i < featureList.size(); i++)
    {
      feature = (Feature) featureList.get(i);
      
      // assume only one featureloc
      List locs = feature.getFeaturelocsForFeatureId();
      
      if(locs.size() > 0)
      {
        FeatureLoc loc = (FeatureLoc)locs.get(0);
        int fmin = loc.getFmin() + 1;
        int fmax = loc.getFmax();
        rowData[i][4] = fmin + "..." + fmax;
        rowData[i][5] = Integer.toString(loc.getStrand());
      }
      
      String schema = feature.getOrganism().getAbbreviation().toLowerCase();
      int ind = schema.indexOf('.');
      if(ind > 0)
        schema = schema.substring(0,ind)+schema.substring(ind+1);
      
      System.out.println("get featur type_id.......");
      rowData[i][0] = schema;
      rowData[i][1] = feature.getUniquename();
      rowData[i][2] = feature.getCvterm().getName();
      //rowData[i][2] = (String)cvterm.get(new Long(feature.getType_id()));
      rowData[i][3] = Integer.toString(feature.getId());
      rowData[i][6] = feature.getTimelastmodified().toString();
    
    }
    return rowData;
  }

  /**
   * Reset the schema.
   * @param location
   * @param schema
   */
  private void reset(String location, String schema)
  {
    if(!location.endsWith("="+schema))
    {
      int index = location.lastIndexOf('=');
      location = location.substring(0,index+1) + schema;
      connIB  = null;
      jdbcDAO = null;
      System.setProperty("chado", location);
    }
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


  public class SelectionListener implements ListSelectionListener
  {
    public void valueChanged(ListSelectionEvent e)
    {
      showAttributes();
    }
  }

  public static void main(String args[])
  {
    new ChadoDemo();
  }
}
