/*
 * Copyright (C) 2008  Genome Research Limited
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
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
 *  @author: Tim Carver
 */

package uk.ac.sanger.artemis.circular;


import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.util.Vector;
import java.util.Hashtable;

import javax.swing.Box;
import javax.swing.ButtonGroup;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;

import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.SimpleEntryGroup;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.io.RangeVector;
import uk.ac.sanger.artemis.sequence.NoSequenceException;
import uk.ac.sanger.artemis.util.OutOfRangeException;


/**
* DNA draw wizard 
*/
public class Wizard
{
  private DNADraw dna = null;
     
  public static Track TRACK_1 = new Track(0.95d, "CDS", "pseudo", false, true, false, null);
  public static Track TRACK_2 = new Track(0.9d,  "CDS", "pseudo", false, false, true, null);
  public static Track TRACK_3 = new Track(0.85d, "CDS", "pseudo", true, true, true, null);
  public static Track TRACK_4 = new Track(0.8d,  "misc_feature", true, true, null);
  public static Track TRACK_5 = new Track(0.75d, null, true, true, null);

  
  public static Track[] tracks = { TRACK_1, TRACK_2, TRACK_3, TRACK_4, TRACK_5 };
  
  public Wizard(DNADraw dna_current)
  {
    int n = getOption(dna_current);  // option 0 - read data file
                                     // option 1 - create dna display
                                     // option 2 - edit existing dna
    if(n == 0)
      dna = getFeaturesFromFile(dna_current);
    else if(n == 1 || n == 2)
    {
      Vector block = new Vector();
      Vector restrictionEnzyme = new Vector();
      if(dna_current == null)
        dna = new DNADraw();
      else
      {
        dna = dna_current;
        block = dna_current.getGeneticMarker();
        restrictionEnzyme = dna_current.getRestrictionEnzyme();
      }

      LineAttribute la = new LineAttribute(dna);

      GeneticMarker gm;
      if(dna_current != null)
        gm = new GeneticMarker(dna_current,block);
      else
        gm = new GeneticMarker(dna,block);

      RestrictionEnzyme re;
      if(dna_current != null)
        re = new RestrictionEnzyme(dna_current,restrictionEnzyme);
      else
        re = new RestrictionEnzyme(dna,restrictionEnzyme);

      Ticks tk = new Ticks(dna_current,false);

      la.setMinimumSize(la.getPreferredSize());
      la.setMaximumSize(la.getPreferredSize());

      re.setMinimumSize(re.getPreferredSize());
      re.setMaximumSize(re.getPreferredSize());

      ScrollPanel pane = new ScrollPanel(new BorderLayout());
      Box bdown = Box.createVerticalBox();
      bdown.add(new JLabel("Properties"));
      Box bacross = Box.createHorizontalBox();
      bacross.add(la);
      bacross.add(tk);
      bacross.add(Box.createHorizontalGlue());
      bdown.add(bacross);

      bdown.add(new JSeparator());
      bdown.add(Box.createVerticalStrut(10));
      bdown.add(new JLabel("Features"));
      bacross = Box.createHorizontalBox();
      bacross.add(gm);
      bacross.add(Box.createHorizontalGlue());
      bdown.add(bacross);

      bdown.add(new JSeparator());
      bdown.add(Box.createVerticalStrut(10));
      bdown.add(new JLabel("Restriction Enzymes"));
      bacross = Box.createHorizontalBox();
      bacross.add(re);
      bacross.add(Box.createHorizontalGlue());
      bdown.add(bacross);
      pane.add(bdown,BorderLayout.CENTER);
    
      JScrollPane createWizScroll = new JScrollPane(pane);
    
      Dimension dscreen = createWizScroll.getToolkit().getScreenSize();
      int wid = (int)dscreen.getWidth();
      if(wid > 700)
        wid = 700;

      int hgt = (int)dscreen.getHeight();
      if(hgt > 750)
        hgt = 700;
      hgt-=50;

      Dimension d = new Dimension(wid,hgt); 
      createWizScroll.setPreferredSize(d);

      JOptionPane.showMessageDialog(null,
                      createWizScroll, "DNA Wizard",
                      JOptionPane.PLAIN_MESSAGE);

      dna.setGeneticMarker(block);
      dna.setRestrictionEnzyme(restrictionEnzyme);
      dna.setLineAttributes(la.getLineAttr());
      dna.setStartTick(tk.getStartTick());
      dna.setMinorTickInterval(tk.getMinorTickInterval());
      dna.setTickInterval(tk.getTickInterval());

      int s = la.getStart();
      dna.setStart(s);
  
      s = la.getEnd();
      dna.setEnd(s);

    }
  }

  protected static DNADraw getFeaturesFromFile(DNADraw dna_current)
  {
    Options.getOptions();
    uk.ac.sanger.artemis.components.FileDialogEntrySource entrySource = 
      new uk.ac.sanger.artemis.components.FileDialogEntrySource(null, null);

    try
    {
      final EntryGroup entryGroup = new SimpleEntryGroup();
      final Entry entry = entrySource.getEntry(true);
      entryGroup.add(entry);
      
      for(int i=0; i<tracks.length; i++)
        tracks[i].setEntry(entry);
      
      
      FeatureVector features = entry.getAllFeatures();
      Vector block = new Vector();

      if(dna_current == null)
        dna_current = new DNADraw();
      
      dna_current.setArtemisEntryGroup(entryGroup);
      
      Hashtable lineAttr = new Hashtable();
      lineAttr.put("lsize",new Integer(1));
      lineAttr.put("circular",new Boolean(true));
      lineAttr.put("start",new Integer(0));
      lineAttr.put("end",new Integer(entry.getBases().getLength()));

      dna_current.setLineAttributes(lineAttr);       

      for(int i=0; i<features.size(); i++)
      {
        Feature f = features.elementAt(i);
        
        RangeVector ranges = f.getLocation().getRanges();
        
        for(int j=0; j<ranges.size(); j++)
        {
          Range range = (Range) ranges.get(j);

          Color col = f.getColour();
          if(col == null || col.equals(Color.white))
            col = Color.lightGray;

          Track track;
          
          if(TRACK_1.isOnTrack(f))
            track = TRACK_1;
          else if(TRACK_2.isOnTrack(f))
            track = TRACK_2;
          else if(TRACK_3.isOnTrack(f))
            track = TRACK_3;
          else if(TRACK_4.isOnTrack(f))
            track = TRACK_4;
          else
            track = TRACK_5;
          
          Block drawBlock = new Block(f.getIDString(), 
              range.getStart(),
              range.getEnd(), 
              col, 
              10.f, 
              track, dna_current);
          
          drawBlock.setFeature(f);
          block.add(drawBlock);
        }
      }
      
      int div;
      if(entry.getBases().getLength() < 1000)
        div = 100;
      else if(entry.getBases().getLength() < 10000)
        div = 1000;
      else if(entry.getBases().getLength() < 100000)
        div = 10000;
      else
        div = 100000;
      int tick = entry.getBases().getLength()/div;
      tick = tick*(div/10);
  
      dna_current.setGeneticMarker(block);
      dna_current.setRestrictionEnzyme(new Vector());
      dna_current.setMinorTickInterval(tick);
      dna_current.setTickInterval(tick);
    }
    catch(OutOfRangeException e)
    {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    catch(NoSequenceException e)
    {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    return dna_current;
  }
  
  
  
  
  protected static DNADraw readEntry(final DNADraw dna_current)
  {
    Options.getOptions();
    uk.ac.sanger.artemis.components.FileDialogEntrySource entrySource = 
      new uk.ac.sanger.artemis.components.FileDialogEntrySource(null, null);

    try
    {
      final Entry entry = entrySource.getEntry(true);
      dna_current.getArtemisEntryGroup().add(entry);

      FeatureVector features = entry.getAllFeatures();

      Track track = addTrack(entry);
      
      for(int i=0; i<features.size(); i++)
      {
        Feature f = features.elementAt(i);
        RangeVector ranges = f.getLocation().getRanges();
        
        for(int j=0; j<ranges.size(); j++)
        {
          Range range = (Range) ranges.get(j);

          Color col = f.getColour();
          if(col == null || col.equals(Color.white))
            col = Color.lightGray;
          
          Block drawBlock = new Block(f.getIDString(), 
              range.getStart(),
              range.getEnd(), 
              col, 
              10.f, 
              track, dna_current);
          
          drawBlock.setFeature(f);
          dna_current.getGeneticMarker().add(drawBlock);
        }
      }
    }
    catch(OutOfRangeException e)
    {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    catch(NoSequenceException e)
    {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    return dna_current;
  }
  
  
  
  
  

  public DNADraw getDNADraw()
  {
    return dna;
  }


  private int getOption(DNADraw dna_current)
  {
    Box bdown = Box.createVerticalBox();

    JRadioButton[] radioButtons;

    if(dna_current !=  null)
      radioButtons = new JRadioButton[3];
    else
      radioButtons = new JRadioButton[2];
    final ButtonGroup group = new ButtonGroup();
    radioButtons[0] = new JRadioButton("Read in data file");
    group.add(radioButtons[0]);
    radioButtons[1] = new JRadioButton("Create new dna display");
    group.add(radioButtons[1]);
    radioButtons[1].setSelected(true);
    bdown.add(radioButtons[0]);
    bdown.add(radioButtons[1]);

    if(dna_current !=  null)
    {
      radioButtons[2] = new JRadioButton("Edit current dna display");
      group.add(radioButtons[2]);
      radioButtons[2].setSelected(true);
      bdown.add(radioButtons[2]);
    }

    JPanel pane = new JPanel(new BorderLayout());
    pane.add(bdown);
    JOptionPane.showMessageDialog(null, 
                      pane, "DNA Viewer Wizard", 
                      JOptionPane.QUESTION_MESSAGE);

    if(radioButtons[0].isSelected())
      return 0;
    else if(radioButtons[1].isSelected())
      return 1;
    else if(radioButtons[2].isSelected())
      return 2;
     
    return 1;
  }

  protected static Track[] getTracks()
  {
    return tracks;
  }
  
  
  protected static Track addTrack(Entry entry)
  {
    Track[] tracks = getTracks();
    Track[] newTracks = new Track[tracks.length+1];
    for(int i=0; i<tracks.length; i++)
      newTracks[i] = tracks[i];
    
    Track newTrack = new Track(
        tracks[tracks.length-1].getPosition()-0.05, entry);
    newTracks[tracks.length] = newTrack;
    Wizard.tracks = newTracks;
    return newTrack;
  }
}

