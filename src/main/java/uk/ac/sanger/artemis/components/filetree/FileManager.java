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

import uk.ac.sanger.artemis.util.StringVector;
import uk.ac.sanger.artemis.Options;

import javax.swing.table.TableColumn;
import javax.swing.*;
import java.io.File;
import java.io.FileFilter;
import java.awt.event.*;
import java.awt.geom.*;
import java.awt.*;
import java.util.Properties;
import java.util.Vector;
import java.util.Enumeration;

public class FileManager extends JFrame
{

  /** */
  private static final long serialVersionUID = 1L;

  public FileManager(JFrame frame)
  {
    this(frame,getArtemisFilter());
  }

  /**
  *
  * File Manager Frame
  * @param frame  parent frame
  *
  */
  public FileManager(JFrame frame, FileFilter filter)
  {
    super("File Manager");

    FileSystemModel model = new FileSystemModel(getLocalDirectories(), filter, this);
    JTreeTable ftree = new JTreeTable(model);
    JScrollPane jsp = new JScrollPane(ftree);
    jsp.getViewport().setBackground(Color.white);

    JPanel pane = (JPanel)getContentPane();
    pane.setLayout(new BorderLayout());
    pane.add(jsp, BorderLayout.CENTER);
    setJMenuBar(makeMenuBar(pane,ftree));
    pane.add(getFileFileterComboBox(model, ftree), BorderLayout.SOUTH);

    Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();

    Dimension panelSize = new Dimension((int)(screen.getWidth()/3),
                                        (int)(screen.getHeight()/2));
    jsp.setPreferredSize(panelSize);

    int width = panelSize.width;
    TableColumn col0 = ftree.getColumnModel().getColumn(0);
    col0.setPreferredWidth( (int)(width*0.60) );

    TableColumn col1  = ftree.getColumnModel().getColumn(1);
    col1.setPreferredWidth( (int)(width*0.12) );

    TableColumn col2 = ftree.getColumnModel().getColumn(2);
    col2.setPreferredWidth( (int)(width*0.28) );

    pack();
    
    int yloc = (int)((screen.getHeight()-getHeight())/2);
    setLocation(0,yloc);  
    setVisible(true);
  }

  /**
  *
  * Look in j2ssh.properties for local directories.
  *
  */
  private File[] getLocalDirectories()
  {
    final Properties settings = new Properties();
    ClassLoader cl = getClass().getClassLoader();
    // try out of the classpath
    try
    {
      settings.load(cl.getResourceAsStream("j2ssh.properties"));
    }
    catch (Exception e)
    {
    }

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

  protected JComboBox getFileFileterComboBox(final FileSystemModel model,
                                             final JTreeTable ftree)
  {
    String[] filters = { "Artemis Files", "Sequence Files", 
                         "Feature Files", "All Files" };
    final JComboBox comboFilter = new JComboBox(filters);
    comboFilter.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
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
  *
  * Set up a menu and tool bar
  * @param pane   panel to add toolbar to
  * @param ftree  file tree display
  *
  */
  private JMenuBar makeMenuBar(JPanel pane, final JTreeTable ftree)
  {
    JMenuBar mBar = new JMenuBar();
    JMenu fileMenu = new JMenu("File");
    mBar.add(fileMenu);
    
    JMenuItem fileMenuClose = new JMenuItem("Close");
    fileMenuClose.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        setVisible(false);
      }
    });
    fileMenu.add(fileMenuClose);

    return mBar;
  }

}

