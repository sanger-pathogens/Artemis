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

import uk.ac.sanger.artemis.j2ssh.SshLogin;
import uk.ac.sanger.artemis.util.StringVector;
import uk.ac.sanger.artemis.Options;

import javax.swing.*;
import java.io.File;
import java.io.FileFilter;
import java.awt.event.*;
import java.awt.geom.*;
import java.awt.*;
import java.util.Properties;
import java.util.Enumeration;
import java.util.Vector;
import javax.swing.border.Border;

public class LocalAndRemoteFileManager extends JFrame
{

  /** busy cursor */
  private Cursor cbusy = new Cursor(Cursor.WAIT_CURSOR);
  /** done cursor */
  private Cursor cdone = new Cursor(Cursor.DEFAULT_CURSOR);

  public LocalAndRemoteFileManager(JFrame frame)
  {
    this(frame,getArtemisFilter());
  }

  /**
  *
  * File Manager Frame
  * @param frame  parent frame
  *
  */
  public LocalAndRemoteFileManager(JFrame frame, FileFilter filter)
  {
    super("File Manager");

    FileList flist = new FileList();
    
    if(flist.isConnected())
    {
      final JPanel localPanel = new JPanel(new BorderLayout());
    
      JTreeTable ftree = new JTreeTable(new FileSystemModel(getLocalDirectories(), filter));
      JScrollPane localTree = new JScrollPane(ftree);
      localPanel.add(localTree,BorderLayout.CENTER);

      final JLabel local_status_line = getStatusLabel("LOCAL");
      localPanel.add(local_status_line,BorderLayout.NORTH);

      final JPanel remotePanel = new JPanel(new BorderLayout());

      SshJTreeTable sshtree = new SshJTreeTable(
                       new FileSystemModel( getRemoteDirectories(flist.pwd()) ));
      JScrollPane remoteTree = new JScrollPane(sshtree);
      remotePanel.add(remoteTree,BorderLayout.CENTER);
    
      String remote_name = SshLogin.getHostname();
      if(!SshLogin.getPort().equals(""))
        remote_name = remote_name + ":" + SshLogin.getPort();

      final JLabel remote_status_line = getStatusLabel("REMOTE "+remote_name);
      remotePanel.add(remote_status_line,BorderLayout.NORTH);

      final JSplitPane treePane = new JSplitPane(JSplitPane.VERTICAL_SPLIT,
                                                 localPanel,remotePanel);

      JPanel pane = (JPanel)getContentPane();
      pane.setLayout(new BorderLayout());
      pane.add(treePane, BorderLayout.CENTER);

      Dimension screen    = Toolkit.getDefaultToolkit().getScreenSize();
      Dimension panelSize = new Dimension((int)(screen.getWidth()/3), 
                                          (int)(screen.getHeight()/3));
      treePane.setDividerLocation((int)(screen.getHeight()/3));
      setJMenuBar(makeMenuBar(pane,ftree,sshtree,localPanel,remotePanel,treePane,panelSize));
      localPanel.add(getFileFileterComboBox(ftree), BorderLayout.SOUTH);
  
      remoteTree.setPreferredSize(panelSize);
      localTree.setPreferredSize(panelSize);

      pack();
    
      int yloc = (int)((screen.getHeight()-getHeight())/2);
      setLocation(0,yloc);  
      setVisible(true);
    }
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
  *
  * Set up a menu and tool bar
  * @param pane   panel to add toolbar to
  * @param ftree  file tree display
  *
  */
  private JMenuBar makeMenuBar(JPanel pane, final JTreeTable ftree, final SshJTreeTable sshtree,
                               final JPanel localPanel, final JPanel remotePanel,
                               final JSplitPane treePane, final Dimension panelSize)
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
//  prefV.setSelected(true);
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
    prefH.setSelected(true);

    JMenuItem fileMenuClose = new JMenuItem("Close");
    fileMenuClose.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        setVisible(false);
      }
    });
    fileMenu.add(fileMenuClose);

    // remote tool bar set up
//  JToolBar remoteToolBar  = new JToolBar();
//  remotePanel.add(remoteToolBar, BorderLayout.NORTH);

    // local tool bar set up
//  JToolBar toolBar  = new JToolBar();
//  localPanel.add(toolBar, BorderLayout.NORTH);

    return mBar;
  }

  /**
  *
  * Used to draw a Shape.
  *
  */
  public static GeneralPath makeShape(float loc[][]) 
  {
    GeneralPath path = new GeneralPath(GeneralPath.WIND_NON_ZERO);

    path.moveTo(loc[0][0],loc[0][1]);

    for(int i=1; i<loc.length; i++)
      path.lineTo(loc[i][0],loc[i][1]);
    
    return path;
  }

  class UpButton extends JButton
  {
    public void paintComponent(Graphics g)
    {
      super.paintComponent(g);
      Graphics2D g2 = (Graphics2D)g; 

      g2.setColor(new Color(0,128,0));
      float loc1[][] = { {11,18}, {7,18}, {7,14},
                         {3,14},  {11,4} };
                
      g2.fill(makeShape(loc1));
      g2.setColor(Color.green);

      float loc2[][] = { {11,18}, {15,18}, {15,14},
                         {19,14},  {11,4} };
      g2.fill(makeShape(loc2));

      setSize(22,24);
    }
  }

  class TextButton extends JButton
  {
    private String s1;
    private String s2;
    private Color c;

    public TextButton(String text)
    {
      this(text, Color.black);
    }

    public TextButton(String text, Color c)
    {
      super();
      this.s1 = text.substring(0);
      this.s2 = text.substring(1);
    }

    public void paintComponent(Graphics g)
    {
      super.paintComponent(g);
      Graphics2D g2 = (Graphics2D)g;
      Font font = new Font("Monospaced", Font.BOLD, 14);
      g2.setFont(font);

      g2.setColor(c);
      g2.drawString(s1,4,18);
      g2.setColor(Color.red);
      g2.drawString(s2,10,15);
      setSize(22,24);
    }
  }


  public static void main(String args[])
  {
    final javax.swing.LookAndFeel look_and_feel =
      javax.swing.UIManager.getLookAndFeel();

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
}

