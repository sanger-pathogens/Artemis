/*
 *
 * created: Wed Aug 3 2004
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2000  Genome Research Limited
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or(at your option) any later version.
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

package uk.ac.sanger.artemis.editor;

import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Insets;
import java.awt.Toolkit;

import java.util.Vector;

import javax.swing.*;

public class BigPane extends JFrame
{

  protected static Font font = new Font("Monospaced",Font.PLAIN,11);
  protected static Font font_sm = new Font("Monospaced",Font.PLAIN,10);
  protected static JCheckBoxMenuItem srsBrowser;
  protected static JCheckBoxMenuItem srsTabPane;
  protected static JCheckBoxMenuItem srsWin;
  protected static JInternalFrame srsFrame;
  private JTextArea qualifier;
  private DataViewInternalFrame dataView;
  
  public BigPane(Object dataFile[], JTextArea qualifier)
  {
    this(dataFile,qualifier.getText());
    this.qualifier = qualifier;
  }

  public BigPane(Object dataFile[], String qualifier_txt)
  {
    super("Object Editor");

    MultiLineToolTipUI.initialize();
    setFont(font);
    JDesktopPane desktop = new JDesktopPane();
    desktop.setDragMode(JDesktopPane.LIVE_DRAG_MODE);
    getContentPane().add(desktop);

    //Make the big window be indented 80 pixels from each edge
    //of the screen.
    int inset = 80;
    Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
    setBounds(inset, inset,
              screenSize.width  - inset*2,
              screenSize.height - inset*2);

    addWindowListener(new winExit());

    // data set
    int hgt = getHeight()-85;
    int wid = getWidth()-100;
    dataView = new DataViewInternalFrame(dataFile,desktop,
                                       wid,qualifier_txt);
    dataView.setLocation(50,0);
    dataView.setSize(wid,hgt);
    dataView.setVisible(true);
    desktop.add(dataView);

    JMenuBar menuBar = createMenuBar(desktop);
    setJMenuBar(menuBar);

// toolbar
    JToolBar toolBar = createToolbar();
    getContentPane().add(toolBar,BorderLayout.NORTH);

    setVisible(true);
  }


  /**
  *
  * Create a toolbar
  * @return toolbar.
  *
  */
  private JToolBar createToolbar()
  {
    JToolBar toolBar = new JToolBar();
    
    JButton applyButt = new JButton("APPLY");
    applyButt.setToolTipText("Apply annotation changed to feature editor");
    applyButt.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        qualifier.setText(dataView.getFeatureText());
      }
    });
    applyButt.setBackground(new Color(0,0,81));
    applyButt.setForeground(Color.white);
    applyButt.setBorderPainted(false);
    applyButt.setMargin(new Insets(0,0,0,0));
    applyButt.setFont(font);
    toolBar.add(applyButt);
    
    return toolBar;
  }

  /**
  *
  * Create a menu bar. 
  * @param desktop pane.
  * @return menu bar.
  *
  */
  private JMenuBar createMenuBar(final JDesktopPane desktop)
  {
    JMenuBar menuBar = new JMenuBar();
    JMenu fileMenu = new JMenu("File");
    menuBar.add(fileMenu);

    JMenuItem reReadMenu = new JMenuItem("Re-read selected results");
    reReadMenu.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        dataView.reReadSelectedResults();
      }
    });
    fileMenu.add(reReadMenu);
    fileMenu.add(new JSeparator());
 
    JMenuItem applyMenu = new JMenuItem("Apply to Feature Editor");
    applyMenu.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        qualifier.setText(dataView.getFeatureText());
      }
    });
    fileMenu.add(applyMenu);
    fileMenu.add(new JSeparator());
//
    JMenuItem exitMenu = new JMenuItem("Close");
    exitMenu.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        if(qualifier == null)
          System.exit(0);

        String oldTxt = qualifier.getText().trim();
        String newTxt = dataView.getFeatureText().trim();

        // changes have been made to feature annotation
        if(!oldTxt.equals(newTxt))
        {
          int ok = JOptionPane.showConfirmDialog(BigPane.this, 
                              "Apply changes now?",
                              "Apply Changes",
                              JOptionPane.YES_NO_CANCEL_OPTION,
                              JOptionPane.QUESTION_MESSAGE);

          if(ok == JOptionPane.CANCEL_OPTION)
            return;

          if(ok == JOptionPane.OK_OPTION)
            qualifier.setText(newTxt);
        }

        BigPane.srsFrame = null;
        dispose();
      }
    });
        
    fileMenu.add(exitMenu);

   //srs menu items
    JMenu optionMenu = new JMenu("Options");
    menuBar.add(optionMenu);

    JMenu srsMenu = new JMenu("Show SRS in");
    optionMenu.add(srsMenu);
    
    srsBrowser = new JCheckBoxMenuItem("Browser",false);
    srsTabPane = new JCheckBoxMenuItem("Tab Pane",true);
    srsWin     = new JCheckBoxMenuItem("New Window",false);

    srsMenu.add(srsBrowser);
    srsMenu.add(srsTabPane);
    srsMenu.add(srsWin);

   //drag mode
    JMenu dragMenu = new JMenu("Drag Mode");
    optionMenu.add(dragMenu);

    JRadioButtonMenuItem liveDrag = new JRadioButtonMenuItem("Live",true);
    liveDrag.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        desktop.setDragMode(JDesktopPane.LIVE_DRAG_MODE);
      }
    });

    dragMenu.add(liveDrag);

    JRadioButtonMenuItem outlineDrag = new JRadioButtonMenuItem("Outline",false);
    outlineDrag.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        desktop.setDragMode(JDesktopPane.OUTLINE_DRAG_MODE);
      }
    });
    dragMenu.add(outlineDrag);
    ButtonGroup buttGroup = new ButtonGroup();
    buttGroup.add(liveDrag);
    buttGroup.add(outlineDrag);

    return menuBar;
  }

  /**
  *
  * Set up the tabbed SRS frame
  *
  */ 
  protected static void setUpSRSFrame(int hgt, JDesktopPane desktop)
  {
    BigPane.srsFrame = new JInternalFrame("SRS",
                                           true, //resizable
                                           true, //closable
                                           true, //maximizable
                                           true);//iconifiable
    BigPane.srsFrame.setDefaultCloseOperation(JInternalFrame.HIDE_ON_CLOSE);
    
    BigPane.srsFrame.setLocation(0,0);
    BigPane.srsFrame.setSize(500,hgt);
    final JTabbedPane jtab = new JTabbedPane();
    BigPane.srsFrame.getContentPane().add(jtab);

    JMenuBar menuBar = new JMenuBar();
    CommonMenu cmen = new CommonMenu(BigPane.srsFrame);
    menuBar.add(cmen);
    JMenuItem closeTabMenu = new JMenuItem("Close tab");
    cmen.add(closeTabMenu);
    closeTabMenu.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        int select = jtab.getSelectedIndex(); 
        if(select > -1)
          jtab.removeTabAt(select);
      }
    });

    BigPane.srsFrame.setJMenuBar(menuBar);
    desktop.add(BigPane.srsFrame);
  }

  /**
  *
  * Extends WindowAdapter to close window.
  *
  */
  class winExit extends WindowAdapter
  {
     public void windowClosing(WindowEvent we)
     {
        if(qualifier == null)
          System.exit(0);

        qualifier.setText(dataView.getFeatureText());
        BigPane.srsFrame = null;
        dispose();

     }
  }


  public static void main(String args[])
  {
    if(args.length < 1)
    {
      System.out.println("Usage:: java BigPane data_file");
      System.exit(0);
    }
    new BigPane(args,"");
  }
}
