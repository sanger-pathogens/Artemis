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

import javax.swing.*;
import java.awt.*;
import javax.swing.event.ChangeListener;
import javax.swing.event.ChangeEvent;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.util.Vector;
import java.io.File;

public class DataViewInternalFrame extends JInternalFrame
{
  private JTabbedPane tabPane = new JTabbedPane();
  private Annotation ann;
  private Box evidenceBox;

  public DataViewInternalFrame(Object dataFile[], JDesktopPane desktop,
                               final JScrollPane scrollEvidence,
                               int wid, int hgt, String qualifier_txt)
  {
    super("Document " + dataFile[0], 
              true, //resizable
              true, //closable
              true, //maximizable
              true);//iconifiable

// graphical evidence display
    JInternalFrame evidence = new JInternalFrame("Evidence", true,
                                                 true, true, true);
    JPanel evidencePanel = (JPanel)evidence.getContentPane();
    evidenceBox = Box.createVerticalBox();

//
    ann = new Annotation(desktop);

    StringBuffer annFormat = new StringBuffer();
    annFormat.append(htmlBreaks(qualifier_txt.trim()));

    for(int i=0; i<dataFile.length; i++)
    {
      //ensure results file exists
      File fdata = new File((String)dataFile[i]);
      if(!fdata.exists())
      {
        fdata = new File((String)dataFile[i]+".gz");
        
        if(!fdata.exists())
        {
          JOptionPane.showMessageDialog(desktop, "Results file: \n"+
                                      dataFile[i] + "\ndoes not exist!",
                                      "File Not Found",
                                      JOptionPane.WARNING_MESSAGE);
          continue;
        }
      }
  
      String tabName = (String)dataFile[i];
      int ind = tabName.lastIndexOf("/");
      if(ind > -1)
      {
        String go = "";

        if(tabName.indexOf("blastp+go") > -1)
          go = ":: GO :: ";
        tabName = go + tabName.substring(ind+1);
      }

      // add fasta results internal frame
      FastaTextPane fastaPane = new FastaTextPane((String)dataFile[i]);

      if(qualifier_txt.indexOf("/"+fastaPane.getFormat()+"_file=\"") == -1)
      {
        if(i > 0)
          annFormat.append("\n<br>");
        annFormat.append("/"+fastaPane.getFormat()+"_file=\""+
                                     dataFile[i]+"\"");
      }

      // graphical view
      final JScrollPane dbviewScroll = new JScrollPane();
      final DBViewer dbview = new DBViewer(fastaPane,dbviewScroll);
      dbviewScroll.setViewportView(dbview);

      final Dimension d = new Dimension((int)dbviewScroll.getPreferredSize().getWidth(), 
                                       hgt/3);
      final Box xBox = Box.createHorizontalBox();
      final MouseOverButton hide = new MouseOverButton("X");
      hide.setForeground(Color.blue);
      hide.setBackground(Color.white);
      hide.setFont(BigPane.font);
      hide.setMargin(new Insets(0,0,0,0));
      hide.setBorderPainted(false);
      hide.setActionCommand("HIDE");

      final Box bacross = Box.createHorizontalBox();
      bacross.add(dbviewScroll);

      hide.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent event)
        {
          if(hide.getActionCommand().equals("HIDE"))
          {
            bacross.remove(dbviewScroll);
            bacross.add(xBox);
            hide.setText("+");
            scrollEvidence.setViewportView(evidenceBox);
            hide.setActionCommand("SHOW");
          }
          else
          {
            bacross.remove(xBox);
            dbviewScroll.setColumnHeaderView(xBox);
            bacross.add(dbviewScroll);
            hide.setText("X");
            scrollEvidence.setViewportView(evidenceBox);
            hide.setActionCommand("HIDE");
          }
        }
      });

      xBox.add(hide);
      JLabel tabLabel = new JLabel(fastaPane.getFormat()+" "+tabName);
      tabLabel.setFont(BigPane.font);

      tabLabel.setOpaque(true);
      xBox.add(tabLabel);
      xBox.add(Box.createHorizontalGlue());

      dbviewScroll.setPreferredSize(d);
      dbviewScroll.setColumnHeaderView(xBox);
      fastaPane.addFastaListener(dbview);
      evidenceBox.add(bacross);
    
      // add data pane
      DataCollectionPane dataPane =
         new DataCollectionPane(fastaPane,ann,desktop);
      fastaPane.addFastaListener(dataPane);

      ActiveJSplitPane split = new ActiveJSplitPane(JSplitPane.VERTICAL_SPLIT,
                                                    fastaPane,dataPane);
      split.setLabel(tabLabel);
      split.setDividerLocation(150);
      split.setOneTouchExpandable(true);
      if(i == 0)
        split.setActive(true);

      tabPane.add(fastaPane.getFormat()+" "+tabName,split);
    }

//  evidenceBox.add(Box.createVerticalGlue());
  
    // add tab pane listener
    tabPane.addChangeListener(new TabChangeListener());

    // add setActivennotator text pane
    ann.setAnnotation(annFormat.toString().trim());

    JScrollPane annotationScroll = new JScrollPane(ann);   
    annotationScroll.setPreferredSize(new Dimension(500,150));

    JSplitPane split = new JSplitPane(JSplitPane.VERTICAL_SPLIT,
                                      annotationScroll,tabPane);

    split.setDividerLocation(150);
    getContentPane().add(split);
     
    setVisible(true);
    evidence.setVisible(true);
    desktop.add(evidence);
  }


  protected Box getEvidenceBox()
  {
    return evidenceBox;
  }

  protected String getFeatureText()
  {
    return ann.getFeatureText();
  }

  protected void reReadSelectedResults()
  {
    ActiveJSplitPane split = (ActiveJSplitPane)tabPane.getSelectedComponent();
    Component comps[] = split.getComponents();
    for(int i=0; i<comps.length;i++)
    {
      if(comps[i] instanceof FastaTextPane)
      {
        ((FastaTextPane)comps[i]).reRead();
        return;
      }
    }
  }

  private String htmlBreaks(String t)
  {
    int ind = 0;
    while((ind = t.indexOf("\n",ind+5)) > -1)
      t = t.substring(0,ind) + "<br>" +
          t.substring(ind);

    return t;
  }

  public class TabChangeListener implements ChangeListener
  {
    ActiveJSplitPane lastSelected = (ActiveJSplitPane)tabPane.getSelectedComponent();
    public void stateChanged(ChangeEvent e)
    {
      ActiveJSplitPane split = (ActiveJSplitPane)tabPane.getSelectedComponent();
      lastSelected.setActive(false);
      split.setActive(true);
      lastSelected = split;
    }
  }


  public class ActiveJSplitPane extends JSplitPane
  {
    private JLabel tabLabel;
    private Color bg;

    public ActiveJSplitPane(int newOrientation,
                  Component newLeftComponent,
                  Component newRightComponent)
    {
      super(newOrientation,newLeftComponent,newRightComponent);
    }
    
    public void setLabel(JLabel tabLabel)
    {
      this.tabLabel = tabLabel;
      this.bg = tabLabel.getBackground();
    }

    public void setActive(boolean active)
    {
      if(active)
        tabLabel.setBackground(Color.yellow);
      else
        tabLabel.setBackground(bg);
    }
  }
}

