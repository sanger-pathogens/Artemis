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

public class DataViewInternalFrame extends JInternalFrame
{
  private JTabbedPane tabPane = new JTabbedPane();

  public DataViewInternalFrame(String dataFile[],
                               JDesktopPane desktop)
  {
    super("Document " + dataFile[0], 
              true, //resizable
              true, //closable
              true, //maximizable
              true);//iconifiable

    Annotation ann   = new Annotation(desktop);
    JPanel dataPanel = (JPanel)getContentPane();
    dataPanel.setLayout(new BorderLayout());

    StringBuffer annFormat = new StringBuffer();

    for(int i=0; i<dataFile.length; i++)
    {
      // add fasta results internal frame
      FastaTextPane fastaPane = new FastaTextPane(dataFile[i]);

      annFormat.append("/"+fastaPane.getFormat()+"_file=\""+
                                    dataFile[i]+"\"\n");
      // add data pane
      Vector hits = fastaPane.getHitCollection();
      DataCollectionPane dataPane =
         new DataCollectionPane(hits, dataFile[i], fastaPane, ann, desktop);

      JSplitPane split = new JSplitPane(JSplitPane.VERTICAL_SPLIT,
                                        fastaPane,dataPane);
      split.setDividerLocation(150);

      tabPane.add(fastaPane.getFormat(),split);
    }
  
    // add annotator text pane
    ann.setAnnotation(annFormat.toString().trim());
    JScrollPane annotationScroll = new JScrollPane(ann);   
    annotationScroll.setPreferredSize(new Dimension(500,300));

    JSplitPane split = new JSplitPane(JSplitPane.VERTICAL_SPLIT,
                                      annotationScroll,tabPane);
    split.setDividerLocation(250);
    dataPanel.add(split);
     
    setVisible(true);
  }

}

