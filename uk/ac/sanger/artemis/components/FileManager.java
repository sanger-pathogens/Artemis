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

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.util.StringVector;
import uk.ac.sanger.artemis.Options;

import javax.swing.*;
import java.io.File;
import java.io.FileFilter;
import java.awt.event.*;
import java.awt.geom.*;
import java.awt.*;

public class FileManager extends JFrame
{

  /** busy cursor */
  private Cursor cbusy = new Cursor(Cursor.WAIT_CURSOR);
  /** done cursor */
  private Cursor cdone = new Cursor(Cursor.DEFAULT_CURSOR);

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

    FileTree ftree  = new FileTree(new File(System.getProperty("user.dir")),
                                   frame, filter);
    JScrollPane jsp = new JScrollPane(ftree);
    JPanel pane = (JPanel)getContentPane();
    pane.setLayout(new BorderLayout());
    pane.add(jsp, BorderLayout.CENTER);
    setJMenuBar(makeMenuBar(pane,ftree));
    pane.add(getFileFileterComboBox(ftree), BorderLayout.SOUTH);

    Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
    jsp.setPreferredSize(new Dimension(210,
                         (int)(screen.getHeight()/2)));
    pack();
    
    int yloc = (int)((screen.getHeight()-getHeight())/2);
    setLocation(0,yloc);  
    setVisible(true);
  }

  protected JComboBox getFileFileterComboBox(final FileTree ftree)
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
          ftree.setFilter(getArtemisFilter());
        else if(select.equals("Sequence Files"))
          ftree.setFilter(getSequenceFilter());
        else if(select.equals("Feature Files"))
          ftree.setFilter(getFeatureFilter());
        else if(select.equals("All Files"))
        {
          ftree.setFilter(new FileFilter()
          {
            public boolean accept(File pathname)
            {
              if(pathname.getName().startsWith("."))
                return false;
              return true;
            }
          });
        }
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
          final String suffix = sequence_suffixes.elementAt(i);

          if(pathname.getName().endsWith("." + suffix) ||
             pathname.getName().endsWith("." + suffix + ".gz"))
            return true;
        }

        for(int i = 0; i<feature_suffixes.size(); ++i)
        {
          final String suffix = feature_suffixes.elementAt(i);

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
          final String suffix = feature_suffixes.elementAt(i);

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
          final String suffix = sequence_suffixes.elementAt(i);

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
  private JMenuBar makeMenuBar(JPanel pane, final FileTree ftree)
  {
    JMenuBar mBar = new JMenuBar();
    JMenu fileMenu = new JMenu("File");
    mBar.add(fileMenu);
    
    JMenuItem fileMenuGoto = new JMenuItem("Go to Directory ...");
    fileMenuGoto.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        String dir = ftree.getRoot().getAbsolutePath();
        String newDir = JOptionPane.showInputDialog(FileManager.this,
                                             "Go to Directory:",dir);      

        if(newDir == null)
          return;
       
        newDir = newDir.trim();
        File newDirFile = new File(newDir);
        
        if(newDirFile.exists() &&
           newDirFile.canRead() &&
           !newDir.equals(dir))
          ftree.newRoot(newDir);
        else
        {
          String error = null;
          if(!newDirFile.exists())
            error = new String(newDir+" doesn't exist!");
          else if(!newDirFile.canRead())
            error = new String(newDir+" cannot be read!");
          else if(newDir.equals(dir))
            error = new String("Same directory!");

          if(error != null)
            JOptionPane.showMessageDialog(FileManager.this,
                                        error, "Warning",
                                        JOptionPane.WARNING_MESSAGE);
        }
      }
    });
    fileMenu.add(fileMenuGoto);
    fileMenu.add(new JSeparator());
    
    JMenuItem fileMenuClose = new JMenuItem("Close");
    fileMenuClose.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        setVisible(false);
      }
    });
    fileMenu.add(fileMenuClose);

    // tool bar set up
    JToolBar toolBar  = new JToolBar();
    Dimension buttonSize = new Dimension(22,24);

    JButton upBt = new JButton()
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
    };
    upBt.setPreferredSize(buttonSize);
    upBt.setMinimumSize(buttonSize);

    upBt.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        FileManager.this.setCursor(cbusy);
        File root = ftree.getRoot();
        String parent = root.getParent();
        if(parent != null)
          ftree.newRoot(parent);
        FileManager.this.setCursor(cdone);
      }
    });
    toolBar.add(upBt);

// yeastpub
    JButton shortCut1 = new JButton()
    {
      public void paintComponent(Graphics g)
      {
        super.paintComponent(g);
        Graphics2D g2 = (Graphics2D)g;
        Font font = new Font("Monospaced", Font.BOLD, 14);
        g2.setFont(font);

        g2.setColor(Color.black);
        g2.drawString("Y",4,18);
        g2.setColor(Color.red);
        g2.drawString("P",10,15);
        setSize(22,24);
      }
    };
    shortCut1.setPreferredSize(buttonSize);
    shortCut1.setMinimumSize(buttonSize);
    shortCut1.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        ftree.newRoot("/nfs/disk222/yeastpub");
      }
    });

    if((new File("/nfs/disk222/yeastpub")).exists())
      toolBar.add(shortCut1);

// pathdata
   JButton shortCut2 = new JButton()
    {
      public void paintComponent(Graphics g)
      {
        super.paintComponent(g);
        Graphics2D g2 = (Graphics2D)g;
        Font font = new Font("Monospaced", Font.BOLD, 14);
        g2.setFont(font);

        g2.setColor(Color.black);
        g2.drawString("P",4,18);
        g2.setColor(Color.red);
        g2.drawString("D",10,15);
        setSize(22,24);
      }
    };
    shortCut2.setPreferredSize(buttonSize);
    shortCut2.setMinimumSize(buttonSize);
    shortCut2.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        ftree.newRoot("/nfs/pathdata/");
      }
    });

    if((new File("/nfs/pathdata/")).exists())
      toolBar.add(shortCut2);

// home button
    JButton homeBt = new JButton()
    {
      public void paintComponent(Graphics g)
      {
        super.paintComponent(g);
        Graphics2D g2 = (Graphics2D)g;
                                                                                
        g2.setColor(Color.blue);
        float loc1[][] = { {3,14}, {11,3}, {19,14},
                           {17,14}, {17,18}, {5,18}, {5,14} };
        g2.fill(makeShape(loc1));
                                                                                
        setSize(22,24);
      }
    };
    homeBt.setPreferredSize(buttonSize);
    homeBt.setMinimumSize(buttonSize);
    homeBt.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        ftree.newRoot(System.getProperty("user.home"));
      }
    });
    toolBar.add(homeBt);

    toolBar.add(Box.createVerticalStrut(35));
    pane.add(toolBar, BorderLayout.NORTH);

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


}

