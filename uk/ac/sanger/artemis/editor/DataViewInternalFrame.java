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
import java.awt.BorderLayout;
import java.awt.Dimension;
import java.util.Vector;
import java.io.File;
import java.awt.Component;

public class DataViewInternalFrame extends JInternalFrame
{
  private JTabbedPane tabPane = new JTabbedPane();
  private Annotation ann;

  public DataViewInternalFrame(Object dataFile[], JDesktopPane desktop,
                               int wid, String qualifier_txt)
  {
    super("Document " + dataFile[0], 
              true, //resizable
              true, //closable
              true, //maximizable
              true);//iconifiable

    ann   = new Annotation(desktop);

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
      JScrollPane dbviewScroll = new JScrollPane();
      DBViewer dbview = new DBViewer(fastaPane,dbviewScroll);
      dbviewScroll.setViewportView(dbview);
      fastaPane.addFastaListener(dbview);
     
      // add data pane
      DataCollectionPane dataPane =
         new DataCollectionPane(fastaPane,ann,desktop);
      fastaPane.addFastaListener(dataPane);

      JSplitPane split = new JSplitPane(JSplitPane.VERTICAL_SPLIT,
                                        fastaPane,dataPane);
      split.setDividerLocation(150);
      split.setOneTouchExpandable(true);

      JSplitPane tabPanel = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT,
                                           split,dbviewScroll);
      tabPanel.setDividerLocation(wid/2);
      tabPanel.setOneTouchExpandable(true);

      String tabName = (String)dataFile[i];
      int ind = tabName.lastIndexOf("/");
      if(ind > -1)
      {
        String go = "";

        if(tabName.indexOf("blastp+go") > -1)
          go = ":: GO :: ";         
        tabName = go + tabName.substring(ind+1);
      }

      tabPane.add(fastaPane.getFormat()+" "+tabName,tabPanel);

    }
  
    // add annotator text pane
    ann.setAnnotation(annFormat.toString().trim());
    JScrollPane annotationScroll = new JScrollPane(ann);   
    annotationScroll.setPreferredSize(new Dimension(500,300));

    JSplitPane split = new JSplitPane(JSplitPane.VERTICAL_SPLIT,
                                      annotationScroll,tabPane);
    split.setDividerLocation(250);
    getContentPane().add(split);
     
    setVisible(true);
  }

  protected String getFeatureText()
  {
    return ann.getFeatureText();
  }

  protected void reReadSelectedResults()
  {
    JSplitPane split = (JSplitPane)tabPane.getSelectedComponent();
    Component comps[] = split.getComponents();
    for(int i=0; i<comps.length;i++)
    {
      if(comps[i] instanceof JSplitPane) 
      {
        Component comps2[] = ((JSplitPane)comps[i]).getComponents();
        for(int j=0; j<comps2.length;j++)
        {
          if(comps2[j] instanceof FastaTextPane)
          {
            ((FastaTextPane)comps2[j]).reRead();
            return;
          }
        }
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

}

