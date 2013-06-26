/********************************************************************
*
*  This library is free software; you can redistribute it and/or
*  modify it under the terms of the GNU Library General Public
*  License as published by the Free Software Foundation; either
*  version 2 of the License, or (at your option) any later version.
*
*  This library is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
*  Library General Public License for more details.
*
*  You should have received a copy of the GNU Library General Public
*  License along with this library; if not, write to the
*  Free Software Foundation, Inc., 59 Temple Place - Suite 330,
*  Boston, MA  02111-1307, USA.
*
*  Copyright (C) Genome Research Limited
*
********************************************************************/

package uk.ac.sanger.artemis.components.filetree;

import uk.ac.sanger.artemis.components.SwingWorker;
import uk.ac.sanger.artemis.components.database.DatabaseEntrySource;
import uk.ac.sanger.artemis.components.database.DatabaseJPanel;
import uk.ac.sanger.artemis.j2ssh.SshLogin;
import uk.ac.sanger.artemis.j2ssh.SshFileManager;
import uk.ac.sanger.artemis.util.StringVector;
import uk.ac.sanger.artemis.Options;

import javax.swing.table.TableColumn;
import javax.swing.*;

import java.io.File;
import java.io.FileFilter;
import java.awt.event.*;
import java.awt.*;
import java.util.Properties;
import java.util.Enumeration;
import java.util.Vector;
import javax.swing.border.Border;

public class LocalAndRemoteFileManager extends JFrame
{

  /***/
  private static final long serialVersionUID = 1L;
  private JScrollPane remoteTree;
  private SshJTreeTable sshtree;
  private JSplitPane treePane = null;
  private DatabaseEntrySource entry_source;
  public static JCheckBoxMenuItem lazyLoad = 
    new JCheckBoxMenuItem("Lazy load feature data", false);
  private static JCheckBoxMenuItem automaticHistory =
    new JCheckBoxMenuItem("Automatic History Annotation", false);
  
  public static JCheckBoxMenuItem domainLoad = 
    new JCheckBoxMenuItem("Display protein domains", false);

  public LocalAndRemoteFileManager(JFrame frame)
  {
    this(frame,getArtemisFilter());
  }

  /**
  *
  * File Manager Frame
  * @param frame  parent frame
  * @param filter file name filter
  *
  */
  public LocalAndRemoteFileManager(JFrame frame, FileFilter filter)
  {
    super();

    final JPanel localPanel = new JPanel(new BorderLayout());
    
    final SshLogin ssh_login = new SshLogin();
    JTreeTable ftree = new JTreeTable(new FileSystemModel(getLocalDirectories(), 
                                      filter, this));
    JScrollPane localTree = new JScrollPane(ftree);
    localTree.getViewport().setBackground(Color.white);
    localPanel.add(localTree,BorderLayout.CENTER);

    final JLabel local_status_line = getStatusLabel("LOCAL");
    localPanel.add(local_status_line,BorderLayout.NORTH);

    final JPanel remotePanel = new JPanel(new BorderLayout());

    //
    final Dimension screen    = Toolkit.getDefaultToolkit().getScreenSize();
    final Dimension panelSize = new Dimension((int)(screen.getWidth()/3),
                                          (int)(screen.getHeight()/4));
    String remote_name = "";
    final JLabel remote_status_line = getStatusLabel("");

    if(FileList.ssh_client == null)  // if no connection etablished yet
    {
      final Box bdown = Box.createVerticalBox();
      JButton connect = new JButton("Connect");

      connect.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          login(remotePanel, bdown, ssh_login, panelSize, 
                local_status_line, remote_status_line);
        }
      });
 
      bdown.add(ssh_login.getLogin());
      // listen to passwd field for return press
      JPasswordField pwf = ssh_login.getJPasswordField();
      pwf.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          login(remotePanel, bdown, ssh_login, panelSize, 
                local_status_line, remote_status_line); 
        }
      });

      bdown.add(connect);
      /*int ypos = panelSize.height-connect.getPreferredSize().height;
      if(ypos>0)
        bdown.add(Box.createVerticalStrut(ypos/2));*/
      bdown.add(Box.createVerticalGlue());
      
      remotePanel.add(bdown, BorderLayout.SOUTH);
      remotePanel.setPreferredSize(panelSize);
    }
    else
    {
      FileList flist = new FileList();
      setRemoteTree(flist, sshtree, remoteTree, remotePanel,
                    panelSize, remote_status_line);
    }

    remote_status_line.setText("REMOTE "+remote_name);
    remotePanel.add(remote_status_line,BorderLayout.NORTH);

    treePane = new JSplitPane(JSplitPane.VERTICAL_SPLIT,
                                 localPanel,remotePanel);
    treePane.setOneTouchExpandable(true);
    
    JPanel pane = (JPanel)getContentPane();
    pane.setLayout(new BorderLayout());
    
    DbConnectionThread dbthread = null;
    if(System.getProperty("chado") != null)
    { 
      setTitle("Database and File Manager");
      entry_source = new DatabaseEntrySource();
      
      boolean promptUser = true;
      if(System.getProperty("read_only") != null)
      {
        promptUser = false;
        entry_source.setReadOnly(true);
      }
      
      if(!entry_source.setLocation(promptUser))
        return;
    
      JLabel label  = new JLabel(" Database Loading...");
      JPanel dbPane = new JPanel();
      dbPane.add(label);
      dbPane.setBackground(Color.white);
      dbPane.setPreferredSize(panelSize);
      
      JSplitPane mainSplit = new JSplitPane(JSplitPane.VERTICAL_SPLIT,
                                            dbPane, treePane);
      
      dbthread = new DbConnectionThread(mainSplit, panelSize, 
                                        entry_source, dbPane);
      dbthread.start();

      treePane.setDividerLocation((int)(screen.getHeight()/4));
      mainSplit.setOneTouchExpandable(true);
      mainSplit.setDividerLocation((int)(screen.getHeight()/4));
      pane.add(mainSplit, BorderLayout.CENTER);
    }
    else
    {
      setTitle("File Manager");
      pane.add(treePane, BorderLayout.CENTER);
      treePane.setDividerLocation((int)(screen.getHeight()/4));
    }
    setJMenuBar(makeMenuBar(pane,ftree,sshtree,localPanel,
                            remotePanel,treePane,panelSize,dbthread));
    localPanel.add(getFileFileterComboBox(ftree), BorderLayout.SOUTH);

    localTree.setPreferredSize(panelSize);

    // Set the column width
    int width = panelSize.width;
    setColumnWidth(ftree, width);

    pack();
  
    int yloc = (int)((screen.getHeight()-getHeight())/2);
    setLocation(0,yloc);  
    setVisible(true);
  }

  private void login(JPanel remotePanel, Box bdown, SshLogin ssh_login,
                     Dimension panelSize, JLabel local_status_line,
                     JLabel remote_status_line)
  {
    setCursor(new Cursor(Cursor.WAIT_CURSOR));

    final SshFileManager ssh_fm;
    try
    {
      ssh_fm = new SshFileManager(ssh_login);
    }
    catch(NullPointerException npe)
    {
      setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
      JOptionPane.showMessageDialog(LocalAndRemoteFileManager.this,
                                    "Check login details and try again.",
                                    "Failed Login", JOptionPane.ERROR_MESSAGE);
      return;
    }
    FileList flist = new FileList(ssh_fm);
    remotePanel.remove(bdown);
    int divider_loc = treePane.getDividerLocation();
    setRemoteTree(flist, sshtree, remoteTree, remotePanel,
                  panelSize, remote_status_line);

    if(treePane.getOrientation() == JSplitPane.HORIZONTAL_SPLIT)
      treePane.setBottomComponent(remotePanel);
    else
      treePane.setRightComponent(remotePanel);

    treePane.setDividerLocation(divider_loc);

    setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
  }

  private void setRemoteTree(final FileList flist, SshJTreeTable sshtree, 
                          JScrollPane remoteTree, JPanel remotePanel,
                          final Dimension panelSize, final JLabel remote_status_line)
  {
    sshtree = new SshJTreeTable(new FileSystemModel( 
                      getRemoteDirectories(flist.pwd()), LocalAndRemoteFileManager.this),
                      LocalAndRemoteFileManager.this);
    remoteTree = new JScrollPane(sshtree);
    remoteTree.setPreferredSize(panelSize);
    remoteTree.getViewport().setBackground(Color.white);
    remotePanel.add(remoteTree,BorderLayout.CENTER);

    String remote_name = SshLogin.getHostname();
    if(!SshLogin.getPort().equals(""))
      remote_name = remote_name + ":" + SshLogin.getPort();
    remote_status_line.setText("REMOTE "+remote_name);
    setColumnWidth(sshtree, panelSize.width);
  }


  private void setColumnWidth(JTable table, int width)
  {
    TableColumn col0 = table.getColumnModel().getColumn(0);
    col0.setPreferredWidth( (int)(width*0.60) );

    TableColumn col1  = table.getColumnModel().getColumn(1);
    col1.setPreferredWidth( (int)(width*0.12) );

    TableColumn col2 = table.getColumnModel().getColumn(2);
    col2.setPreferredWidth( (int)(width*0.28) );
  }

  /**
  *
  * Create a status JLabel with bevelled border
  *
  */
  private JLabel getStatusLabel(String status)
  {
    final JLabel status_line = new JLabel(status);
    Border loweredbevel = BorderFactory.createLoweredBevelBorder();
    Border raisedbevel = BorderFactory.createRaisedBevelBorder();
    Border compound = BorderFactory.createCompoundBorder(raisedbevel,loweredbevel);
    status_line.setBorder(compound);

    final FontMetrics fm =
      this.getFontMetrics(status_line.getFont());
    final int font_height = fm.getHeight()+10;

    status_line.setMinimumSize(new Dimension(100, font_height));
    status_line.setPreferredSize(new Dimension(100, font_height));
    return status_line;
  }

  /**
  *
  * Look in j2ssh.properties for local directories.
  *
  */
  private File[] getLocalDirectories()
  {
    final Properties settings = SshLogin.getProperties();
    Enumeration enum_prop = settings.propertyNames();
    Vector dirs = new Vector();

    dirs.add(new File(System.getProperty("user.home")));
    dirs.add(new File(System.getProperty("user.dir")));

    while(enum_prop.hasMoreElements())
    {
      final String property = (String)enum_prop.nextElement();
      File f = new File(settings.getProperty(property));
      if(property.startsWith("localdir") && f.exists())
        dirs.add(f);
    }

    File fdirs[] = new File[dirs.size()];
    for(int i=0; i<dirs.size(); i++)
      fdirs[i] = (File)dirs.get(i);

    return fdirs;
  }

  /**
  *
  * Look in j2ssh.properties for remote directories.
  *
  */
  private String[] getRemoteDirectories(String pwd)
  {
    final Properties settings = SshLogin.getProperties();
    Enumeration enum_prop = settings.propertyNames();
    Vector dirs = new Vector();
    dirs.add(pwd);
    while(enum_prop.hasMoreElements())
    {
      final String property = (String)enum_prop.nextElement();
      if(property.startsWith("remotedir"))
        dirs.add(settings.getProperty(property));
    }

    String sdirs[] = new String[dirs.size()];
    for(int i=0; i<dirs.size(); i++)
      sdirs[i] = (String)dirs.get(i);

    return sdirs;
  }

  protected JComboBox getFileFileterComboBox(final JTreeTable ftree)
  {
    String[] filters = { "Artemis Files", "Sequence Files", 
                         "Feature Files", "All Files" };
    final JComboBox comboFilter = new JComboBox(filters);
    comboFilter.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        FileSystemModel model = (FileSystemModel)(ftree.getTree().getModel());
        String select = (String)comboFilter.getSelectedItem(); 
        if(select.equals("Artemis Files"))
          model.setFilter(getArtemisFilter());
        else if(select.equals("Sequence Files"))
          model.setFilter(getSequenceFilter());
        else if(select.equals("Feature Files"))
          model.setFilter(getFeatureFilter());
        else if(select.equals("All Files"))
        {
          model.setFilter(new FileFilter()
          {
            public boolean accept(File pathname)
            {
              if(pathname.getName().startsWith("."))
                return false;
              return true;
            }
          });
        }
        ftree.refreshAll();
        ftree.revalidate();
      }
    });
    return comboFilter;
  }

  /**
  *
  * Get a file filter for sequence and feature suffixes.
  * @return file filter
  */
  protected static FileFilter getArtemisFilter()
  {
    final StringVector sequence_suffixes =
      Options.getOptions().getOptionValues("sequence_file_suffixes");

    final StringVector feature_suffixes =
      Options.getOptions().getOptionValues("feature_file_suffixes");

    final FileFilter artemis_filter = new FileFilter()
    {
      public boolean accept(File pathname)
      {
        if(pathname.isDirectory() &&
           !pathname.getName().startsWith("."))
          return true;
          
        for(int i = 0; i<sequence_suffixes.size(); ++i)
        {
          final String suffix = (String)sequence_suffixes.elementAt(i);

          if(pathname.getName().endsWith("." + suffix) ||
             pathname.getName().endsWith("." + suffix + ".gz"))
            return true;
        }

        for(int i = 0; i<feature_suffixes.size(); ++i)
        {
          final String suffix = (String)feature_suffixes.elementAt(i);

          if(pathname.getName().endsWith("." + suffix) ||
             pathname.getName().endsWith("." + suffix + ".gz"))
            return true;
        }
        return false;
      }
    };
    return artemis_filter;
  }


  /**
  *
  * Get a file filter for feature suffixes.
  * @return file filter
  */
  protected static FileFilter getFeatureFilter()
  {
    final StringVector feature_suffixes =
      Options.getOptions().getOptionValues("feature_file_suffixes");

    final FileFilter feature_filter = new FileFilter()
    {
      public boolean accept(File pathname)
      {
        if(pathname.isDirectory() &&
           !pathname.getName().startsWith("."))
          return true;

        for(int i = 0; i<feature_suffixes.size(); ++i)
        {
          final String suffix = (String)feature_suffixes.elementAt(i);

          if(pathname.getName().endsWith("." + suffix) ||
             pathname.getName().endsWith("." + suffix + ".gz"))
            return true;
        }
        return false;
      }
    };
    return feature_filter;
  }

  /**
  *
  * Get a file filter for sequence suffixes.
  * @return file filter
  */
  protected static FileFilter getSequenceFilter()
  {
    final StringVector sequence_suffixes =
      Options.getOptions().getOptionValues("sequence_file_suffixes");

    final FileFilter seq_filter = new FileFilter()
    {
      public boolean accept(File pathname)
      {
        if(pathname.isDirectory() &&
           !pathname.getName().startsWith("."))
          return true;
         
        for(int i = 0; i<sequence_suffixes.size(); ++i)
        {
          final String suffix = (String)sequence_suffixes.elementAt(i);

          if(pathname.getName().endsWith("." + suffix) ||
             pathname.getName().endsWith("." + suffix + ".gz"))
            return true;
        }

        return false;
      }
    };
    return seq_filter;
  }

  /**
   * @return the automatic history
   */
  public static boolean isAutomaticHistory()
  {
    return automaticHistory.isSelected();
  }
  
  public static void setAutomaticHistory(boolean flag)
  {
    automaticHistory.setSelected(flag);
  }

  /**
  *
  * Set up a menu and tool bar
  * @param pane   panel to add toolbar to
  * @param ftree  file tree display
  *
  */
  private JMenuBar makeMenuBar(JPanel pane, 
                               final JTreeTable ftree, final SshJTreeTable sshtree,
                               final JPanel localPanel, final JPanel remotePanel,
                               final JSplitPane treePane, final Dimension panelSize,
                               final DbConnectionThread dbthread)
  {
    JMenuBar mBar = new JMenuBar();
    JMenu fileMenu = new JMenu("File");
    mBar.add(fileMenu);
    
    JRadioButtonMenuItem prefV = new JRadioButtonMenuItem("Vertical Split");
    fileMenu.add(prefV);
    prefV.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        treePane.remove(remotePanel);
        treePane.remove(localPanel);
        treePane.setOrientation(JSplitPane.VERTICAL_SPLIT);
        treePane.setTopComponent(localPanel);
        treePane.setBottomComponent(remotePanel);
        remotePanel.setPreferredSize(panelSize);
        localPanel.setPreferredSize(panelSize);

        pack();
        treePane.setDividerLocation(0.5);
      }
    });
    prefV.setSelected(true);
    ButtonGroup group = new ButtonGroup();
    group.add(prefV);

    JRadioButtonMenuItem prefH = new JRadioButtonMenuItem("Horizontal Split");
    fileMenu.add(prefH);
    prefH.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        treePane.remove(remotePanel);
        treePane.remove(localPanel);
        treePane.setOrientation(JSplitPane.HORIZONTAL_SPLIT);
        treePane.setLeftComponent(localPanel);
        treePane.setRightComponent(remotePanel);

        remotePanel.setPreferredSize(panelSize);
        localPanel.setPreferredSize(panelSize);

        pack();
        treePane.setDividerLocation(0.5);
      }
    });
    group.add(prefH);
//  prefH.setSelected(true);

    if(System.getProperty("chado") != null)
    {
      fileMenu.add(new JSeparator());
      final JMenuItem fileShow = new JMenuItem("Open Selected Database Sequence ...");
      fileShow.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          if(dbthread.getDatabaseJPanel() != null)
            dbthread.getDatabaseJPanel().showSelected(entry_source, null);
        }
      });
      fileMenu.add(fileShow);

      final JCheckBoxMenuItem splitGFF = new JCheckBoxMenuItem(
          "Split into entries ...", false);
      splitGFF.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          if(splitGFF.isSelected())
          {
            DatabaseEntryFilterPanel messagePanel = new DatabaseEntryFilterPanel();
            int val = JOptionPane.showConfirmDialog(
                LocalAndRemoteFileManager.this, messagePanel, "Define Entry",
                JOptionPane.OK_CANCEL_OPTION);

            if(val == JOptionPane.OK_OPTION)
              messagePanel.setTypesForEntries();
            else
              splitGFF.setSelected(false);
          }
          
          dbthread.setSplitGFFEntry(splitGFF.isSelected());
        }
      });
      fileMenu.add(splitGFF);
      
      if(Options.getOptions().getPropertyTruthValue("show_polypeptide_domains"))
        domainLoad.setSelected(true);
      fileMenu.add(domainLoad);
      fileMenu.add(lazyLoad);
      
      if(Options.getOptions().getPropertyTruthValue("automatic_history_annotation"))
        automaticHistory.setSelected(true);
      fileMenu.add(automaticHistory);
      
      JMenuItem clearCache = new JMenuItem("Clear database manager cache");
      clearCache.addActionListener(new ActionListener()
      {

        public void actionPerformed(ActionEvent e)
        {
          File cacheDir = new File(Options.CACHE_PATH);
          
          if(cacheDir.exists())
          {
            File cacheFiles[] = cacheDir.listFiles();
            if(cacheFiles != null)
            {
              for(int i=0; i<cacheFiles.length; i++)
                cacheFiles[i].delete();
            }
          }
        }
      });
      fileMenu.add(clearCache);
    }
    
    
    final JMenuItem validate = new JMenuItem("Validate Selected Sequence / Organism ...");
    validate.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        SwingWorker entryWorker = new SwingWorker()
        {
          public Object construct()
          {
            if(dbthread.getDatabaseJPanel() != null)
            {
              dbthread.getDatabaseJPanel().setCursor(new Cursor(Cursor.WAIT_CURSOR));
              try
              {
                dbthread.getDatabaseJPanel().validate(entry_source);
              }
              finally
              {
                dbthread.getDatabaseJPanel().setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
              }
            }
            return null;
          }
        };
        entryWorker.start();
      }
    });
    fileMenu.add(new JSeparator());
    fileMenu.add(validate);
    
    JMenuItem fileMenuClose = new JMenuItem("Close");
    fileMenuClose.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        setVisible(false);
      }
    });
    fileMenu.add(new JSeparator());
    fileMenu.add(fileMenuClose);

    // remote tool bar set up
//  JToolBar remoteToolBar  = new JToolBar();
//  remotePanel.add(remoteToolBar, BorderLayout.NORTH);

    // local tool bar set up
//  JToolBar toolBar  = new JToolBar();
//  localPanel.add(toolBar, BorderLayout.NORTH);

    return mBar;
  }
  
  
  private class DbConnectionThread extends Thread
  {
    private JSplitPane dbSplitPane;
    private Dimension panelSize;
    private DatabaseEntrySource entry_source;
    private JPanel topPanel;
    private DatabaseJPanel dbPane;
    private boolean splitGFFEntry = false;
    
    public DbConnectionThread(final JSplitPane dbSplitPane,
                              final Dimension panelSize,
                              final DatabaseEntrySource entry_source,
                              final JPanel topPanel)
    {
      this.dbSplitPane = dbSplitPane;
      this.panelSize = panelSize;
      this.entry_source = entry_source;
      this.topPanel = topPanel;
    }

    public void run()
    {
      topPanel.setCursor(new Cursor(Cursor.WAIT_CURSOR));
      dbPane = new DatabaseJPanel(entry_source, null);
      dbPane.setPreferredSize(panelSize);
      dbPane.setSplitGFFEntry(splitGFFEntry);
      dbSplitPane.setTopComponent(dbPane);
      topPanel.setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
    }
    
    protected void setSplitGFFEntry(final boolean splitGFFEntry)
    {
      if(dbPane != null)
        dbPane.setSplitGFFEntry(splitGFFEntry);
      else
        this.splitGFFEntry = splitGFFEntry;
    }
    
    protected DatabaseJPanel getDatabaseJPanel()
    {
      return dbPane;
    }
  }


  public static void main(String args[])
  {
    //final javax.swing.LookAndFeel look_and_feel =
    //  javax.swing.UIManager.getLookAndFeel();

    final javax.swing.plaf.FontUIResource font_ui_resource =
      Options.getOptions().getFontUIResource();

    java.util.Enumeration keys = UIManager.getDefaults().keys();
    while(keys.hasMoreElements())
    {
      Object key = keys.nextElement();
      Object value = UIManager.get(key);
      if(value instanceof javax.swing.plaf.FontUIResource)
        UIManager.put(key, font_ui_resource);
    }

    JFrame frame = new LocalAndRemoteFileManager(null);
    frame.pack();
    frame.setVisible(true);
  }
  
  public DatabaseEntrySource getDatabaseEntrySource()
  {
    return entry_source;
  }

}

