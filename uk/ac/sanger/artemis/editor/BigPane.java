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

import java.util.StringTokenizer;
import java.util.Vector;
import javax.swing.*;

import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.io.Location;

public class BigPane extends JFrame
{
  protected static Font font    = new Font("Monospaced",Font.PLAIN,11);
  protected static Font font_sm = new Font("Monospaced",Font.PLAIN,10);
  protected static JCheckBoxMenuItem srsBrowser;
  protected static JCheckBoxMenuItem srsTabPane;
  protected static JCheckBoxMenuItem srsWin;
  protected static JInternalFrame srsFrame;
  protected static JCheckBox addNote = new JCheckBox("Add Note");

  private JTextArea qualifier;
  private DataViewInternalFrame dataView;
  private FeatureVector overlapFeature;
  private Feature edit_feature;
  private JDesktopPane desktop = null;

  public BigPane()
  {
    super("Object Editor");
  }

  public void set(Object dataFile[], JTextArea qualifier,
             FeatureVector overlapFeature, 
             final Feature edit_feature) 
  {
    set(dataFile,qualifier.getText(),overlapFeature,edit_feature);
    this.qualifier      = qualifier;
  }

  public void set(Object dataFile[], String qualifier_txt,
             FeatureVector overlapFeature,
             final Feature edit_feature)
  {
    this.overlapFeature = overlapFeature;
    this.edit_feature   = edit_feature;

    addNote.setSelected(false);
    setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
    MultiLineToolTipUI.initialize();
    setFont(font);

    if(desktop == null)
    {
      desktop = new JDesktopPane();
      desktop.setDragMode(JDesktopPane.LIVE_DRAG_MODE);
      getContentPane().add(desktop);
    }

    //Make the big window be indented 80 pixels from each edge
    //of the screen.
    final int inset = 80;
    final Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
    setBounds(inset, inset,
              screenSize.width  - inset*2,
              screenSize.height - inset*2);

    addWindowListener(new winExit());

    final JScrollPane scrollEvidence = new JScrollPane();
    // data set
    final int hgt = getHeight()-85;
    final int wid = getWidth()/2-10;
    
    dataView = new DataViewInternalFrame(dataFile,desktop, scrollEvidence,
                                         wid,hgt,qualifier_txt);
    dataView.setLocation(5,0);
    dataView.setSize(wid,hgt);
    dataView.setVisible(true);
    desktop.add(dataView);

    // evidence
    final JInternalFrame evidence = new JInternalFrame("Evidence", true,
                                                      true, true, true);

    Box evidenceBox = dataView.getEvidenceBox();
    if(overlapFeature != null)
      evidenceBox.add(getOverlapFeatures(overlapFeature,desktop),0);
    evidenceBox.add(Box.createVerticalGlue());

    ScrollPanel scroller = new ScrollPanel();
    scroller.add(evidenceBox);
    scrollEvidence.setViewportView(scroller);
    evidence.getContentPane().add(scrollEvidence);
    evidence.setLocation(wid+10,0);
    evidence.setSize(wid,hgt);
    evidence.setVisible(true);
    desktop.add(evidence);   

    final JMenuBar menuBar = createMenuBar(desktop);
    setJMenuBar(menuBar);

    // toolbar
    final JToolBar toolBar = createToolbar();
    getContentPane().add(toolBar,BorderLayout.NORTH);

    setVisible(true);
  }


  /**
  *
  * Display for overlapping Pfam features.
  *
  */
  private Box getOverlapFeatures(FeatureVector overlapFeature,
                                 JDesktopPane desktop)
  {
    Box bdown = Box.createVerticalBox();
    bdown.add(new EvidenceViewer(edit_feature,overlapFeature,desktop));
    return bdown;
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

    addNote.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        if(addNote.isSelected())
          dataView.updateNote();
        else
          dataView.deleteNote();
      }
    });
    toolBar.add(addNote);
    
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

    final JMenuItem reReadMenu = new JMenuItem("Re-read selected results");
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
    final JMenuItem exitMenu = new JMenuItem("Close");
    exitMenu.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        onClose();
      }
    });
        
    fileMenu.add(exitMenu);

   //srs menu items
    final JMenu optionMenu = new JMenu("Options");
    menuBar.add(optionMenu);

    final JMenu srsMenu = new JMenu("Show SRS in");
    optionMenu.add(srsMenu);
    
    srsBrowser = new JCheckBoxMenuItem("Browser",false);
    srsTabPane = new JCheckBoxMenuItem("Tab Pane",true);
    srsWin     = new JCheckBoxMenuItem("New Window",false);

    srsMenu.add(srsBrowser);
    srsMenu.add(srsTabPane);
    srsMenu.add(srsWin);

   //drag mode
    final JMenu dragMenu = new JMenu("Drag Mode");
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

    final JRadioButtonMenuItem outlineDrag = new JRadioButtonMenuItem("Outline",false);
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
  * Routine to call when the editor is closed.
  *
  */
  private void onClose()
  {
    if(qualifier == null)
      System.exit(0);

    // remember the splitpane divider locations
    dataView.setDataDividerLocation();
    dataView.setAnnotationDividerLocation();

    // update feature text
    final String oldTxt = qualifier.getText().trim();
    final String newTxt = dataView.getFeatureText().trim();

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

    // stop getz processes
    setVisible(false);
    dataView.stopGetz();
    dataView.dispose();
    BigPane.srsFrame = null;
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

    final JMenuBar menuBar = new JMenuBar();
    final CommonMenu cmen = new CommonMenu(BigPane.srsFrame);
    menuBar.add(cmen);
    final JMenuItem closeTabMenu = new JMenuItem("Close tab");
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
       onClose();
     }
  }


  public static void main(String args[])
  {
    if(args.length < 1)
    {
      System.out.println("Usage:: java BigPane data_file");
      System.exit(0);
    }
    BigPane bp = new BigPane();
    bp.set(args,"",null,null);
  }
}
