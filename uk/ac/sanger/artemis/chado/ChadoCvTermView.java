/* 
 *
 * created: 2007
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2007  Genome Research Limited
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

import org.gmod.schema.cv.Cv;
import org.gmod.schema.cv.CvTerm;

import uk.ac.sanger.artemis.components.genebuilder.JExtendedComboBox;
import uk.ac.sanger.artemis.components.genebuilder.cv.CvTermsComparator;
import uk.ac.sanger.artemis.util.DatabaseLocationParser;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.sql.SQLException;
import java.util.Collections;
import java.util.List;
import java.util.Vector;

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPasswordField;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.ListCellRenderer;
import javax.swing.SwingConstants;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumn;

/**
 * Chado data access example code. This searches for features by their
 * uniquename and returns their properties and attributes.
 * 
 * @author tjc
 * 
 */
public class ChadoCvTermView extends JFrame
{
  /** */
  private static final long serialVersionUID = 1L;

  /** database URL */
  private String location;

  /** password fields */
  private JPasswordField pfield;
  
  private GmodDAO dao;
  
  /** 
   * Chado demo
   */
  public ChadoCvTermView(GmodDAO dao)
  { 
    super("Controlled Vocabulary");
    super.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
    
    try
    {
      if(dao == null)
      {
        setLocation();
        dao = getDAO();
      }

      this.dao = dao;
      
      final List cvs = dao.getAllCvs();
      
      Cv cvSequence = null;
      for(int i=0; i<cvs.size(); i++)
      {
        Cv cv = (Cv)cvs.get(i);
        if(cv.getName().equalsIgnoreCase("sequence"))
        {
          cvSequence = cv;
          break;
        }
      }
      
      final JComboBox cvList = new JComboBox(new Vector(cvs));
      cvList.setRenderer(new CvRenderer());
      
      if(cvSequence != null)
        cvList.setSelectedItem(cvSequence);
      
      final JPanel panel = new JPanel(new BorderLayout());
      final Box searchBox = Box.createVerticalBox();
      final JPanel searchPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
      searchBox.add(searchPanel);
      panel.add(searchBox, BorderLayout.NORTH);
      
      searchPanel.add(new JLabel("CV:"));
      searchPanel.add(cvList);
      
      final JTextField searchTerm = new JTextField(15);
      searchTerm.setText("gene");

      searchPanel.add(new JLabel("Term:"));
      searchPanel.add(searchTerm);
      
      final JButton searchButt = new JButton("GO");
      
      searchButt.addActionListener(new ActionListener()
      {
        private JExtendedComboBox cvTermsList;
        private JTable detailTable;
        
        public void actionPerformed(ActionEvent e)
        {
          final GmodDAO dao = ChadoCvTermView.this.dao;
          
          String search = searchTerm.getText().trim();
          search = search.replaceAll("\\*", "%");

          List cvTerms = dao.getCvTermByNameInCv(search,
                                               (Cv)cvList.getSelectedItem());
          
          Collections.sort(cvTerms, new CvTermsComparator());
          
          if(cvTermsList != null)
            searchBox.remove(cvTermsList);
          
          if(detailTable != null)
            panel.remove(detailTable);
          
          cvTermsList = new JExtendedComboBox(new Vector(cvTerms));
          cvTermsList.setPreferredSize(new Dimension(
              searchPanel.getPreferredSize().width,cvTermsList.getPreferredSize().height));
          cvTermsList.addActionListener(new ActionListener()
          {
            public void actionPerformed(ActionEvent e)
            {
              detailTable = showCvTable(dao, searchTerm, 
                  cvTermsList, cvList, searchBox, panel, 
                  searchPanel, detailTable);
            }
          });
          
          cvTermsList.setRenderer(new CvRenderer());
          searchBox.add(cvTermsList);
          panel.revalidate();
          pack();
        } 
      });
      searchPanel.add(searchButt);

      getContentPane().add(panel);
      pack();
      setVisible(true);
      
      searchTerm.requestFocus();
      searchTerm.selectAll();
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
   * 
   * @param dao
   * @param searchTerm
   * @param cvTermsList
   * @param cvList
   * @param searchBox
   * @param panel
   * @param searchPanel
   * @param detailTable
   * @param frame
   * @return
   */
  private JTable showCvTable(final GmodDAO dao, final JTextField searchTerm,
        final JExtendedComboBox cvTermsList, final JComboBox cvList,
        final Box searchBox, final JPanel panel, final JPanel searchPanel,
        JTable detailTable)
  {
    CvTerm cvTerm = (CvTerm) cvTermsList.getSelectedItem();

    final String columnData[] = { "", "" };
    final String rowData[][] = 
    {
        { "CV Name", cvTerm.getCv().getName() },
        { "Term Name", cvTerm.getName() },
        { "Definition", cvTerm.getDefinition() },
        {
            "DbXRef",
            cvTerm.getDbXRef().getDb().getName() + ":"
                + cvTerm.getDbXRef().getAccession() } 
    };

    if(detailTable != null)
      panel.remove(detailTable);
    detailTable = new JTable(rowData, columnData);
    TableColumn col = detailTable.getColumnModel().getColumn(0);
    col.setMinWidth(35);
    col.setMaxWidth(80);
    col.setPreferredWidth(80);

    TextCellRenderer renderer = new TextCellRenderer();
    detailTable.getColumnModel().getColumn(1).setCellRenderer(renderer);
    detailTable.getColumnModel().getColumn(0).setCellRenderer(
        new LabelTableRender());

    panel.add(detailTable, BorderLayout.CENTER);
    panel.revalidate();
    pack();

    int rowHeight = 0;
    for(int i = 0; i < detailTable.getRowCount(); i++)
    {
      renderer.getTableCellRendererComponent(detailTable, 
          rowData[i][1], false,
          false, i, 1);
      rowHeight += detailTable.getRowHeight(i);
    }

    detailTable.setPreferredSize(new Dimension(
        detailTable.getPreferredSize().width, rowHeight));
    panel.revalidate();
    pack();
    return detailTable;
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
   * Get the data access object (DAO).
   * 
   * @return data access object
   */
  private GmodDAO getDAO() throws java.net.ConnectException, SQLException
  {
    if(System.getProperty("ibatis") == null)
      return new JdbcDAO(location, pfield);
    else
      return new IBatisDAO(pfield);
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

    
    DatabaseLocationParser dlp; 
    // given -Dchado=localhost:port/dbname?username
    if(System.getProperty("chado") != null)
    {
      String db_url = System.getProperty("chado").trim();
      dlp = new DatabaseLocationParser(db_url);
      inServer.setText(dlp.getHost());
      inPort.setText("" + dlp.getPort());
      inDB.setText(dlp.getDatabase());
      inUser.setText(dlp.getUsername());
   
    }
    else
    {
        // Still need to initialise the object
        dlp = new DatabaseLocationParser();
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

  class CvRenderer extends JLabel implements ListCellRenderer
  {
    /** */
    private static final long serialVersionUID = 1L;

    public CvRenderer() 
    {
      setOpaque(true);
    }
    
    public Component getListCellRendererComponent(
        JList list,
        Object value,
        int index,
        boolean isSelected,
        boolean cellHasFocus)
    {
      if(value instanceof Cv)
      {
        Cv cv = (Cv)value;
        setText(cv.getName());
      }
      else if(value instanceof CvTerm)
      {
        CvTerm cvTerm = (CvTerm)value;
        setText(cvTerm.getName());
      }
      else
        setText((String)value);
      
      setBackground(isSelected ? Color.red : Color.white);
      setForeground(isSelected ? Color.white : Color.black);
      return this;
    }
  }
  
  
  class TextCellRenderer extends JTextArea implements TableCellRenderer 
  {
    /** */
    private static final long serialVersionUID = 1L;

    public TextCellRenderer() 
    {
      setLineWrap(true);
      setWrapStyleWord(true);
    }

    public Component getTableCellRendererComponent(JTable table, Object
          value, boolean isSelected, boolean hasFocus, int row, int column) 
    {
      if(value == null)
        setText("");
      else
        setText(value.toString());
      setSize(table.getColumnModel().getColumn(column).getWidth(),
              getPreferredSize().height);

      if(table.getRowHeight(row) != getPreferredSize().height) 
        table.setRowHeight(row, getPreferredSize().height);
      
      return this;
    }

  }
  
  public class LabelTableRender extends DefaultTableCellRenderer 
  {
    /** */
    private static final long serialVersionUID = 1L;

    public LabelTableRender() 
    {
      super();
      this.setVerticalAlignment(SwingConstants.TOP);
    }
  }

  
  public static void main(String args[])
  {
    new ChadoCvTermView(null);
  }
}
