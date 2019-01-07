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

import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.components.ViewMenu;
import uk.ac.sanger.artemis.util.Document;

import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.util.Hashtable;
import java.util.Vector;
import java.util.StringTokenizer;
import java.util.Enumeration;
import java.io.File;
import java.io.BufferedReader;
import java.io.StringReader;
import java.io.IOException;

public class DataViewInternalFrame extends JInternalFrame
{
  /** */
  private static final long serialVersionUID = 1L;
  protected static int dataDividerLocation = 250;
  protected static int annotationDividerLocation = 150;

  private JTabbedPane tabPane = new JTabbedPane();
  private Annotation ann;
  private Box evidenceBox;
  private Vector fastaCollection = new Vector();
  private JSplitPane split;

  public DataViewInternalFrame(Hashtable dataFile, JDesktopPane desktop,
                               final JScrollPane scrollEvidence,
                               int wid, int hgt, String qualifier_txt,
                               final Feature edit_feature)
  {
    super("Document " , 
              true, //resizable
              true, //closable
              true, //maximizable
              true);//iconifiable

// graphical evidence display
    JInternalFrame evidence = new JInternalFrame("Evidence", true,
                                                 true, true, true);
    //JPanel evidencePanel = (JPanel)evidence.getContentPane();
    evidenceBox = Box.createVerticalBox();

//
    ann = new Annotation(desktop);

    StringBuffer annFormat = new StringBuffer();
    annFormat.append(htmlBreaks(qualifier_txt.trim()));

    int icount = 0;
    Enumeration enumPrograms = dataFile.keys();
    while(enumPrograms.hasMoreElements())
    {
      String programName = (String) enumPrograms.nextElement();
      Vector files = (Vector) dataFile.get(programName);

      for(int i = 0; i < files.size(); i++)
      {
        String fileName = (String) files.get(i);
        Document document = null;
        
        try
        {
          String thisFileName = fileName;
          int ind = fileName.indexOf(programName + File.separatorChar);
          if(ind > -1) 
            thisFileName = fileName.substring (ind + programName.length () + 1);

          document = ViewMenu.getSearchDocument(edit_feature, programName, thisFileName);
        }
        catch (IOException e)
        {
          e.printStackTrace();
        }

        if(document == null)
        {
          JOptionPane.showMessageDialog(desktop, "Results file: \n"
              + fileName + "\ndoes not exist!", "File Not Found",
              JOptionPane.WARNING_MESSAGE);
          continue;
        }

        String tabName = (String) fileName;
        int ind = tabName.lastIndexOf("/");
        if(ind > -1)
        {
          String go = "";

          if(tabName.indexOf("blastp+go") > -1)
            go = ":: GO :: ";
          tabName = go + tabName.substring(ind + 1);
        }

        // add fasta results internal frame
        FastaTextPane fastaPane = new FastaTextPane(document);
        if(fastaPane.getFormat() != null)
          fastaCollection.add(fastaPane);
        else
          continue;
        
        /*
        if(qualifier_txt.indexOf("/" + fastaPane.getFormat() + "_file=\"") == -1)
        {
          if(icount > 0)
            annFormat.append("\n<br>");
          annFormat.append("/" + fastaPane.getFormat() + "_file=\"" + fileName
              + "\"");
        }
        */

        // graphical view
        final JScrollPane dbviewScroll = new JScrollPane();
        final DBViewer dbview = new DBViewer(fastaPane, dbviewScroll);
        dbviewScroll.setViewportView(dbview);

        final Dimension d = new Dimension((int) dbviewScroll.getPreferredSize()
            .getWidth(), hgt / 3);

        final Box yBox = Box.createVerticalBox();
        final Box xBox = Box.createHorizontalBox();
        final MouseOverButton hide = new MouseOverButton("X");
        hide.setForeground(Color.blue);
        hide.setBackground(Color.white);
        hide.setFont(BigPane.font);
        hide.setMargin(new Insets(0, 0, 0, 0));
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
              bacross.add(yBox);
              hide.setText("+");
              scrollEvidence.setViewportView(evidenceBox);
              hide.setActionCommand("SHOW");
            }
            else
            {
              bacross.remove(yBox);
              dbviewScroll.setColumnHeaderView(yBox);
              bacross.add(dbviewScroll);
              hide.setText("X");
              scrollEvidence.setViewportView(evidenceBox);
              hide.setActionCommand("HIDE");
            }
          }
        });

        xBox.add(hide);
        JLabel tabLabel = new JLabel(fastaPane.getFormat() + " " + tabName);
        tabLabel.setFont(BigPane.font);

        tabLabel.setOpaque(true);
        xBox.add(tabLabel);
        xBox.add(Box.createHorizontalGlue());

        yBox.add(xBox);
        yBox.add(dbview.getRuler());

        dbviewScroll.setPreferredSize(d);
        dbviewScroll.setColumnHeaderView(yBox);
        fastaPane.addFastaListener(dbview);

        evidenceBox.add(bacross);

        // add data pane
        DataCollectionPane dataPane = new DataCollectionPane(fastaPane, ann,
            desktop, this);
        fastaPane.addFastaListener(dataPane);

        ActiveJSplitPane split = new ActiveJSplitPane(
            JSplitPane.VERTICAL_SPLIT, fastaPane, dataPane);
        split.setLabel(tabLabel);
        split.setDividerLocation(DataViewInternalFrame.dataDividerLocation);
        split.setOneTouchExpandable(true);
        if(icount == 0)
          split.setActive(true);

        tabPane.add(fastaPane.getFormat() + " " + tabName, split);
        icount++;
      }
    }

// evidenceBox.add(Box.createVerticalGlue());
  
    // add tab pane listener
    tabPane.addChangeListener(new TabChangeListener());

    // add setActive annotator text pane
    ann.setAnnotation(annFormat.toString().trim());

    JScrollPane annotationScroll = new JScrollPane(ann);   
    annotationScroll.setPreferredSize(new Dimension(500,150));

    split = new JSplitPane(JSplitPane.VERTICAL_SPLIT,
                           annotationScroll,tabPane);

    split.setDividerLocation(annotationDividerLocation);
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


  protected void stopGetz()
  {
    Enumeration fastaEnum = fastaCollection.elements();

    while(fastaEnum.hasMoreElements())
      ((FastaTextPane)fastaEnum.nextElement()).stopGetz();
  }


  /**
  *
  * Delete note field in for all similar hits.
  *
  */
  protected void deleteNote()
  {
    ann.deleteNote();
  }


  /**
  * 
  * Add a note field in for all similar hits.
  *
  */
  protected void updateNote()
  {
    StringReader in     = new StringReader(getFeatureText());
    BufferedReader buff = new BufferedReader(in);
    String line       = null;
    StringBuffer note = null;

    try
    {
      while((line = buff.readLine()) != null)
      {
        if(line.startsWith("/similarity="))
        {
          if(note == null)
            note = new StringBuffer("\n/note=\"Similar to ");
          else
            note.append(", and to ");

          StringTokenizer tok = new StringTokenizer(line,";");
          String type = tok.nextToken();
          int ind1 = type.indexOf("\"");
          type = type.substring(ind1+1);

          String id   = tok.nextToken();
          ind1 = id.indexOf("with=")+4;
          if(ind1 == 3)
            ind1 = 0;         
            
//        ind1 = id.indexOf(":");
          id = id.substring(ind1+1).trim();

          String next = tok.nextToken().trim();
          if(next.endsWith("."))
            next = next.substring(0,next.length()-1);

          note.append(next);
          note.append(tok.nextToken().toLowerCase());

          String length = tok.nextToken().trim();
          if(!length.startsWith("length"))
            note.append(" "+length.toLowerCase());
          note.append(" "+id);

          while(!length.startsWith("length"))
            length = tok.nextToken().trim();

          ind1 = length.indexOf("=");
          if(ind1 == -1)
            ind1 = length.indexOf(" ");

          length = length.substring(ind1+1);

          if(!length.endsWith(" aa"))
            length = length + " aa";

          note.append(" ("+length+") ");
          note.append(type+" scores: ");
          
          String pid  = tok.nextToken().trim();  // percent id

          if(pid.endsWith("%"))  // fasta
          {
            ind1 = pid.indexOf(" ");
            pid  = pid.substring(ind1+1);

            tok.nextToken();                       // ungapped id
            String eval = tok.nextToken().trim();
            String len  = tok.nextToken().trim();

            while(!len.endsWith("aa overlap"))
              len  = tok.nextToken().trim();

            len = len.substring(0,len.length()-8);

            note.append(eval);
            note.append(", "+pid+" id in ");
            note.append(len);
          }
          else                  // blastp
            note.append(pid);
        }
      }
    }
    catch(IOException ioe){}

    if(note != null)
      note.append("\"");
    else
      return;

    ann.insert(note.toString(),true);
//  System.out.println(note.toString());
//  Enumeration enumHits = hits.elements();
//  while(enumHits.hasMoreElements())
//  {
//    String id = (String)enumHits.nextElement();
//    System.out.println(id);
//    findHitInfo(String
//  }
    return;
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

  protected void setDataDividerLocation()
  {
    ActiveJSplitPane lastSelected = (ActiveJSplitPane)tabPane.getSelectedComponent();
    if(lastSelected != null)
      dataDividerLocation = lastSelected.getDividerLocation();
  }

  protected void setAnnotationDividerLocation()
  {
    annotationDividerLocation = split.getDividerLocation();
  }


  public class ActiveJSplitPane extends JSplitPane
  {
    /** */
    private static final long serialVersionUID = 1L;
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

