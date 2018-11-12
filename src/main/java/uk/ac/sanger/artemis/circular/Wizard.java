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
import java.awt.Desktop;
import java.awt.Dimension;
import java.awt.desktop.AboutEvent;
import java.awt.desktop.AboutHandler;
import java.awt.desktop.QuitEvent;
import java.awt.desktop.QuitHandler;
import java.awt.desktop.QuitResponse;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.Vector;
import java.util.Hashtable;

import javax.swing.Box;
import javax.swing.ButtonGroup;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;

import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.EntryVector;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.SimpleEntryGroup;
import uk.ac.sanger.artemis.components.EntryFileDialog;
import uk.ac.sanger.artemis.components.MessageDialog;
import uk.ac.sanger.artemis.components.StickyFileChooser;
import uk.ac.sanger.artemis.components.SwingWorker;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.EntryInformationException;
import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.io.Location;
import uk.ac.sanger.artemis.io.MSPcrunchDocumentEntry;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.io.RangeVector;
import uk.ac.sanger.artemis.sequence.Bases;
import uk.ac.sanger.artemis.sequence.NoSequenceException;
import uk.ac.sanger.artemis.util.Document;
import uk.ac.sanger.artemis.util.DocumentFactory;
import uk.ac.sanger.artemis.util.IconManager;
import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.util.ReadOnlyException;
import uk.ac.sanger.artemis.util.WorkingGZIPInputStream;


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

  private SwingWorker workerGraph;
  public static Track[] tracks = { TRACK_1, TRACK_2, TRACK_3, TRACK_4, TRACK_5 };
  
  /*
   * Default constructor.
   */
  protected Wizard()
  {
	  IconManager.setApplicationIcon(IconManager.DNAPLOTTER_NAME);
		
	  if (System.getProperty("mrj.version") != null ||
	            System.getProperty("os.name").toLowerCase().indexOf("mac") >= 0)
	  {
		  initMac();
	  }
  }
  
  public Wizard(DNADraw dna_current)
  {
	this();
	    
    int n = getOption(dna_current);  // option 0 - read data file
                                     // option 1 - edit existing dna
                                     // option 2 - read template
    if(n == 0)
      dna = getDNADrawFromFile(dna_current);
    else if(n == 2)
    {
      StickyFileChooser chooser = new StickyFileChooser();
      chooser.showOpenDialog(null);

      File fileTemplate = chooser.getSelectedFile();
      if(!fileTemplate.exists())
        JOptionPane.showMessageDialog(null, 
            fileTemplate.getName()+" cannot be found!", 
            "Missing File", JOptionPane.WARNING_MESSAGE);
      loadTemplate(chooser.getSelectedFile());
    }
    else if(n == 1)
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
  
  /**
   * Open a DNA plot based on a template file
   * @param template
   */
  public Wizard(final String templateName)
  {
	this();
    loadTemplate(templateName);
  }
  
  /**
   * Load from a template file
   * @param template
   */
  private void loadTemplate(final String templateName)
  {
    final ProgressFrame progress = new ProgressFrame();
    progress.setString("Reading from "+templateName+"   ");
    progress.setValue(2);
    
    if(dna == null)
      dna = new DNADraw();
    Options.getOptions();
    final BufferedReader inputStream = getReader(templateName);
    loadTemplate(inputStream, templateName, progress);
  }
  
  /**
   * Load from a template file
   * @param template
   */
  private void loadTemplate(final File templateFile)
  {
    final ProgressFrame progress = new ProgressFrame();
    progress.setString("Reading from "+templateFile.getName()+"   ");
    progress.setValue(2);
    
    if(dna == null)
      dna = new DNADraw();
    Options.getOptions();
    FileReader reader;
    try
    {
      reader = new FileReader(templateFile);
      final BufferedReader inputStream = new BufferedReader(reader);
      loadTemplate(inputStream, templateFile.getName(), progress);
    }
    catch(FileNotFoundException e)
    {
      e.printStackTrace();
    }
  }
  
  private void loadTemplate(final BufferedReader inputStream,
                            final String templateName,
                            final ProgressFrame progress)
  {
    try
    {
      final EntryGroup entryGroup = new SimpleEntryGroup();
      final Hashtable fileEntrys = new Hashtable();
      Vector v_tracks = new Vector();
      String inLine = null;
      String lineAttrStr[]  = null;
      String tickMarksStr[] = null;
      String gcGraphStr[]      = null;
      String gcSkewGraphStr[]  = null;
      String userGraphStr[]    = null;
      
      String lineAttrStart    = "# line attributes:";
      String tickMarksStart   = "# tick marks:";
      String gcGraphStart     = "# GC Graph:";
      String gcSkewGraphStart = "# GC Skew Graph:";
      String userGraphStart   = "# User Graph:";
      String mergeBlastFeatures = "# merge.blast";
      
      while((inLine = inputStream.readLine()) != null)
      {
        if(inLine.startsWith("#") || inLine.trim().equals(""))
        {
          if(inLine.startsWith(lineAttrStart))
            lineAttrStr = inLine.substring(lineAttrStart.length()).trim().split("[=\\s]");
          else if(inLine.startsWith(tickMarksStart))
            tickMarksStr = inLine.substring(tickMarksStart.length()).trim().split("[=\\s]");
          else if(inLine.startsWith(gcGraphStart))
            gcGraphStr = inLine.substring(gcGraphStart.length()).trim().split("[=\\s]");
          else if(inLine.startsWith(gcSkewGraphStart))
            gcSkewGraphStr = inLine.substring(gcSkewGraphStart.length()).trim().split("[=\\s]");
          else if(inLine.startsWith(userGraphStart))
            userGraphStr = inLine.substring(userGraphStart.length()).trim().split("[=\\s]");
          else if(inLine.startsWith(mergeBlastFeatures))
            System.setProperty("merge.blast", "true");
          continue;
        }

        String properties[] = inLine.split("\t");
        
        final String separator;
        if (properties[11].indexOf ("://") != -1)
          separator = "/";
        else
          separator = File.separator;
        String fileName = properties[11] + separator  + properties[10];
        Entry entry;
        if(!fileEntrys.containsKey(fileName))
        {
          progress.setString("Reading "+properties[10]);
          progress.setValue(4);
          entry = getEntry(fileName, entryGroup);
          if(entry == null)
            continue;
          fileEntrys.put(fileName, entry);
        }
        else
        {
          entry = (Entry)fileEntrys.get(fileName);
        }
        
        Track track = new Track(.9,entry);
        track.setPropertiesFromTemplate(inLine);
        v_tracks.add(track);
      }
      inputStream.close();
      
      progress.setString("Read template "+templateName);
      progress.setValue(7);
      Track[] newTracks = new Track[v_tracks.size()];
      for(int i=0; i<v_tracks.size(); i++)
        newTracks[i] = (Track) v_tracks.get(i);
      
      Wizard.tracks = newTracks;

      dna.setArtemisEntryGroup(entryGroup);

      int sequenceLength = entryGroup.getSequenceEntry().getBases().getLength();
      
      Hashtable lineAttr = new Hashtable();
      lineAttr.put("lsize", new Integer(1));
      lineAttr.put("circular", new Boolean(true));
      lineAttr.put("start", new Integer(0));
      lineAttr.put("end", new Integer(sequenceLength));
      if(lineAttrStr != null)
      {
        for(int i=0; i<lineAttrStr.length; i++)
        {
          if(lineAttrStr[i].startsWith("line_size"))
            lineAttr.put("lsize", new Integer(lineAttrStr[i+1]));
          else if(lineAttrStr[i].startsWith("circular"))
            lineAttr.put("circular", new Boolean(lineAttrStr[i+1]));
          else if(lineAttrStr[i].startsWith("line_height"))
            dna.setLineHeight(Float.parseFloat(lineAttrStr[i+1]));
          else if(lineAttrStr[i].startsWith("bases_per_line"))
            dna.setBasesPerLine(Integer.parseInt(lineAttrStr[i+1]));
        }
      }
      dna.setLineAttributes(lineAttr);

      final int div;
      if(sequenceLength < 1000)
        div = 100;
      else if(sequenceLength < 10000)
        div = 1000;
      else if(sequenceLength < 100000)
        div = 10000;
      else
        div = 100000;
      int tick = sequenceLength / div;
      tick = tick * (div / 10);
      int tick2 = tick / 2;
      tick = tick2 * 2;
      
      if(tickMarksStr != null)
      {
        for(int i=0; i<tickMarksStr.length; i++)
        {
          if(tickMarksStr[i].startsWith("major"))
            tick = Integer.parseInt(tickMarksStr[i+1]);
          else if(tickMarksStr[i].startsWith("minor"))
            tick2 = Integer.parseInt(tickMarksStr[i+1]);
        }
      }

      dna.setGeneticMarker(new Vector());
      dna.setRestrictionEnzyme(new Vector());
      dna.setMinorTickInterval(tick2);
      dna.setTickInterval(tick);
      
      TrackManager trackManager = dna.getTrackManager();
      if(trackManager == null)
      {
        trackManager = new TrackManager(dna);
        dna.setTrackManager(trackManager);
      }
      trackManager.update(tracks);
      
      final String[] this_gcGraphStr = gcGraphStr;
      final String[] this_gcSkewGraphStr = gcSkewGraphStr;
      final String[] this_userGraphStr = userGraphStr;
      workerGraph = new SwingWorker()
      {
        public Object construct()
        {
          loadGraphs(this_gcGraphStr, this_gcSkewGraphStr, this_userGraphStr, dna, progress);
          return null;
        } 
      };
    }
    catch(FileNotFoundException e)
    {
      e.printStackTrace();
    }
    catch(IOException e)
    {
      e.printStackTrace();
    }
    catch(NoSequenceException e)
    {
      e.printStackTrace();
    }
  }
  
  /**
   * Load graphs using template file details
   * @param gcGraphStr
   * @param gcSkewGraphStr
   * @param userGraphStr
   * @param dna
   */
  private void loadGraphs(final String gcGraphStr[], 
                          final String gcSkewGraphStr[], 
                          final String userGraphStr[],
                          final DNADraw dna,
                          final ProgressFrame progress)
  {
    
    if(gcGraphStr != null)
    {
      if(progress != null)
      {
        progress.setString("Calculating GC graph points");
      }
      GCGraph gcGraph = new GCGraph(dna);
      gcGraph.setOptionsStr(gcGraphStr);
      dna.setGcGraph(gcGraph);
      gcGraph.calcGraphValues();
      dna.add(gcGraph);
      dna.repaint();
      dna.revalidate();
    }
    if(gcSkewGraphStr != null)
    {
      if(progress != null)
      {
        progress.setString("Calculating GC Skew graph points");
        progress.setValue(8);
      }
      GCSkewGraph gcSkewGraph = new GCSkewGraph(dna);
      gcSkewGraph.setOptionsStr(gcSkewGraphStr);
      dna.setGcSkewGraph(gcSkewGraph);
      gcSkewGraph.calcGraphValues();
      dna.add(gcSkewGraph);
      dna.repaint();
      dna.revalidate();
    }
    if(userGraphStr != null)
    {
      String fileName = null;
      for(int i=0;i<userGraphStr.length; i++)
      {
        if(userGraphStr[i].startsWith("file_name"))
          fileName = userGraphStr[i+1];
      }
      //final uk.ac.sanger.artemis.util.Document document =
      //  new uk.ac.sanger.artemis.util.FileDocument(new File(fileName));
      try
      {
        if(progress != null)
        {
          progress.setString("Calculating user graph points");
          progress.setValue(9);
        }
        UserGraph userGraph = new UserGraph(dna, fileName);
        userGraph.setOptionsStr(userGraphStr);
        dna.setUserGraph(userGraph);
        userGraph.calcGraphValues();
        dna.add(userGraph);
        dna.repaint();
        dna.revalidate();
      }
      catch(IOException e)
      {
        e.printStackTrace();
        return;
      }
    }
    progress.dispose();
  }
  
  /**
   * Return an Artemis entry from a file 
   * @param entryFileName
   * @param entryGroup
   * @return
   * @throws NoSequenceException
   */
  private Entry getEntry(final String entryFileName, final EntryGroup entryGroup) 
                   throws NoSequenceException
  {
    final Document entry_document = DocumentFactory.makeDocument(entryFileName);
    final EntryInformation artemis_entry_information =
      Options.getArtemisEntryInformation();
    
    uk.ac.sanger.artemis.io.Entry new_embl_entry =
      EntryFileDialog.getEntryFromFile(null, entry_document,
                                       artemis_entry_information,
                                       false);

    if(new_embl_entry == null)  // the read failed
      return null;

    new_embl_entry = mergeOption(new_embl_entry);
    
    Entry entry = null;
    try
    {
      Bases bases = null;
      if(entryGroup.getSequenceEntry() != null)
        bases = entryGroup.getSequenceEntry().getBases();
      if(bases == null)
        entry = new Entry(new_embl_entry);
      else
        entry = new Entry(bases,new_embl_entry);
      
      entryGroup.add(entry);
    } 
    catch(OutOfRangeException e) 
    {
      new MessageDialog(null, "read failed: one of the features in " +
          entryFileName + " has an out of range " +
                        "location: " + e.getMessage());
    }
    return entry;
  }
  
  
  private static uk.ac.sanger.artemis.io.Entry mergeOption(uk.ac.sanger.artemis.io.Entry embl_entry)
  {
    if(embl_entry instanceof uk.ac.sanger.artemis.io.MSPcrunchDocumentEntry &&
        System.getProperty("merge.blast") == null)
    {
      int status = JOptionPane.showConfirmDialog(null, 
          "This looks like a BLAST file. Do you want to merge\n"+ 
          "overlapping BLAST hits to improve performance?", "Read BLAST", 
          JOptionPane.OK_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE);
      if(status == JOptionPane.OK_OPTION)
        System.setProperty("merge.blast", "");
      else
        System.setProperty("merge.blast", "false");
    }
    
    if( embl_entry instanceof uk.ac.sanger.artemis.io.MSPcrunchDocumentEntry &&
        System.getProperty("merge.blast") != null &&
       !System.getProperty("merge.blast").equals("false") )
      embl_entry = mergeOverlappingFeatures(embl_entry);
    return embl_entry;
  }
  
  /**
   * Merge overlapping BLAST hit features to improve the performance of
   * DNAPlotter.
   * @param embl_entry
   * @return
   */
  private static uk.ac.sanger.artemis.io.Entry mergeOverlappingFeatures(uk.ac.sanger.artemis.io.Entry embl_entry)
  {
    final String name = embl_entry.getName();
    final RangeVector ranges = new RangeVector();
    uk.ac.sanger.artemis.io.FeatureVector features = embl_entry.getAllFeatures();
    System.out.print("Number of features: before merge = "+features.size());
    
    for(int i=0; i<features.size(); i++)
    {
      Range ri = ((uk.ac.sanger.artemis.io.Feature) features.elementAt(i)).
          getLocation().getTotalRange();
      boolean overlaps = false;
      Range rjOld = null;
      Range rj = null;
      int j;
      for(j=0; j<ranges.size(); j++)
      {
        rjOld = (Range) ranges.get(j);
        if(ri.overlaps(rjOld))
        {
          int start = rjOld.getStart();
          int end = rjOld.getEnd();
          if(start > ri.getStart())
            start = ri.getStart();
          if(end < ri.getEnd())
            end = ri.getEnd();
          
          try
          {
            rj = new Range(start, end);
            overlaps = true;
            break;
          }
          catch (OutOfRangeException e){}
        }
      }
      if(!overlaps)
        ranges.add(ri);
      else
      {
        ranges.remove(rjOld);
        ranges.add(rj);
      }
    }
    
    embl_entry.dispose();
    embl_entry = new MSPcrunchDocumentEntry(
        Options.getOptions().getArtemisEntryInformation())
    {
      public boolean isReadOnly () {
        return false;
      }
    };
    embl_entry.setName(name);
  
    final Key key = new Key("CRUNCH_D");
    for(int i=0; i<ranges.size(); i++)
    {
      Range r = (Range)ranges.get(i);
      try
      {
        Feature f = new uk.ac.sanger.artemis.Feature(
            new uk.ac.sanger.artemis.io.EmblStreamFeature(key, new Location(r), null));
        embl_entry.add(f.getEmblFeature());
      }
      catch (ReadOnlyException e){}
      catch (EntryInformationException e){}
      catch (OutOfRangeException e)
      {
        e.printStackTrace();
      }
    }
    
    System.out.println(" after merge = "+ranges.size());
    return embl_entry;
  }

  /**
   * Create a DNADraw panel from a file
   * @param dna_current
   * @return
   */
  protected static DNADraw getDNADrawFromFile(DNADraw dna_current)
  {
    Options.getOptions();
    uk.ac.sanger.artemis.components.FileDialogEntrySource entrySource = 
      new uk.ac.sanger.artemis.components.FileDialogEntrySource(null, null);
    
    final EntryGroup entryGroup = new SimpleEntryGroup();
    Entry entry;
    try
    {
      entry = entrySource.getEntry(true);
      entryGroup.add(entry);
      return getDNADrawFromArtemisEntry(dna_current, entryGroup, entry);
    }
    catch(OutOfRangeException e)
    {
      e.printStackTrace();
    }
    catch(NoSequenceException e)
    {
      JOptionPane.showMessageDialog(null, "No sequence found!", 
          "Sequence Missing", JOptionPane.WARNING_MESSAGE);
    }
    return null;
  }
  
  /**
   * Create a DNADraw panel from an entry
   * @param dna_current
   * @param entryGroup
   * @param entry
   * @return
   */
  public static DNADraw getDNADrawFromArtemisEntry(DNADraw dna_current,
                                             final EntryGroup entryGroup,
                                             final Entry entry)
  {
    for(int i=0; i<tracks.length; i++)
      tracks[i].setEntry(entry);
          
    FeatureVector features = entry.getAllFeatures();
    Vector block = new Vector();

    if(dna_current == null)
      dna_current = new DNADraw();

    dna_current.setArtemisEntryGroup(entryGroup);

    Hashtable lineAttr = new Hashtable();
    lineAttr.put("lsize", new Integer(1));
    lineAttr.put("circular", new Boolean(true));
    lineAttr.put("start", new Integer(0));
    lineAttr.put("end", new Integer(entry.getBases().getLength()));

    dna_current.setLineAttributes(lineAttr);

    for(int i = 0; i < features.size(); i++)
    {
      Feature f = features.elementAt(i);

      RangeVector ranges = f.getLocation().getRanges();

      for(int j = 0; j < ranges.size(); j++)
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

        Block drawBlock = new Block(f.getIDString(), range.getStart(), range
            .getEnd(), col, 10.f, track, dna_current);

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
    int tick = entry.getBases().getLength() / div;
    tick = tick * (div / 10);
    int tick2 = tick / 2;
    tick = tick2 * 2;

    dna_current.setGeneticMarker(block);
    dna_current.setRestrictionEnzyme(new Vector());
    dna_current.setMinorTickInterval(tick2);
    dna_current.setTickInterval(tick);

    EntryVector entries = entryGroup.getActiveEntries();
    for(int i=0; i<entries.size(); i++)
    {
      Entry this_entry = entries.elementAt(i);
      if(!this_entry.getName().equals(entry.getName()))
        addFeaturesFromEntry(this_entry, dna_current);
    }
    return dna_current;
  }
  
  /**
   * Read a new entry from a file
   * @param dna_current
   * @param bases
   * @return
   */ 
  public static DNADraw readEntry(final DNADraw dna_current,
                                  final Bases bases)
  {
    Options.getOptions();
    uk.ac.sanger.artemis.components.FileDialogEntrySource entrySource = 
      new uk.ac.sanger.artemis.components.FileDialogEntrySource(null, null);

    try
    {
      Entry entry = entrySource.getEntry(bases,true);
      
      if(entry.getEMBLEntry() instanceof uk.ac.sanger.artemis.io.MSPcrunchDocumentEntry)
      {
        uk.ac.sanger.artemis.io.Entry new_embl_entry = 
          mergeOption(entry.getEMBLEntry());
        entry = new Entry(bases, new_embl_entry);
      }
      
      dna_current.getArtemisEntryGroup().add(entry);

      addFeaturesFromEntry(entry, dna_current);
    }
    catch(OutOfRangeException e)
    {
      JOptionPane.showMessageDialog(null, 
          "Feature found out of range:\n"+
          e.getMessage(),"Out of Range", 
          JOptionPane.WARNING_MESSAGE);
    }

    return dna_current;
  }
  
  /**
   * Add features from an entry to a new track
   * @param entry
   * @param dna_current
   */
  private static void addFeaturesFromEntry(final Entry entry, 
                                   final DNADraw dna_current)
  {
    FeatureVector features = entry.getAllFeatures();

    Track track = addTrack(entry);
    track.setAny(true);
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
  
  

  public DNADraw getDNADraw()
  {
    return dna;
  }


  private int getOption(DNADraw dna_current)
  {
    Box bdown = Box.createVerticalBox();

    JRadioButton[] radioButtons;

    radioButtons = new JRadioButton[2];

    final ButtonGroup group = new ButtonGroup();
    radioButtons[0] = new JRadioButton("Read in sequence file");
    group.add(radioButtons[0]);

    radioButtons[0].setSelected(true);
    bdown.add(radioButtons[0]);


    radioButtons[0].setSelected(true);
    if(dna_current !=  null)
    {
      radioButtons[1] = new JRadioButton("Edit current dna display");
      group.add(radioButtons[1]);
      radioButtons[1].setSelected(true);
    }
    else
    {
      radioButtons[1] = new JRadioButton("Read template file");
      group.add(radioButtons[1]);  
    }
    bdown.add(radioButtons[1]);
    
    JPanel pane = new JPanel(new BorderLayout());
    pane.add(bdown);
    JOptionPane.showMessageDialog(null, 
                      pane, "DNA Viewer Wizard", 
                      JOptionPane.QUESTION_MESSAGE);

    if(radioButtons[0].isSelected())
      return 0;
    else if(radioButtons[1].isSelected() && dna_current !=  null)
      return 1;
    else if(radioButtons[1].isSelected())
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
    
    final double position;
    if(tracks.length > 1)
      position = tracks[tracks.length-1].getPosition()-0.05;
    else
      position = 0.95;
    Track newTrack = new Track(position, entry);
    newTracks[tracks.length] = newTrack;
    Wizard.tracks = newTracks;
    return newTrack;
  }
  
  protected static void deleteTrack(int trackIndex)
  {
    Track[] tracks = getTracks();
    Track[] newTracks = new Track[tracks.length-1];
    int count = 0;
    for(int i=0; i<tracks.length; i++)
    {
      if(i == trackIndex)
        continue;
      
      newTracks[count] = tracks[i];
      count++;
    }

    Wizard.tracks = newTracks;
  }
  
  /**
   * Return reader for a file or a URL
   * @param templateName
   * @return
   */
  protected static BufferedReader getReader(final String templateName)
  {
    final File fileTemplate = new File(templateName);
    BufferedReader inputStream = null;
    if(!fileTemplate.exists())
    {
      if (templateName.indexOf ("://") != -1) 
      {
        URL template;
        try
        {
          template = new URL(templateName);
          inputStream = new BufferedReader(
              new InputStreamReader(template.openStream()));
        }
        catch(MalformedURLException e)
        {
          e.printStackTrace();
        }
        catch(IOException e)
        {
          e.printStackTrace();
        }
      } 
    }
    else
    {
      try
      {
        if(templateName.endsWith(".gz"))
          inputStream = new BufferedReader(
              new InputStreamReader(new WorkingGZIPInputStream( 
                  new FileInputStream(fileTemplate) )));
        else
          inputStream = new BufferedReader(new FileReader(fileTemplate));
      }
      catch(FileNotFoundException e)
      {
        e.printStackTrace();
      }
      catch (IOException e)
      {
        e.printStackTrace();
      }
    }
    return inputStream;
  }

  public SwingWorker getWorkerGraph()
  {
    return workerGraph;
  }
  
  /**
   * Set up about and exit Mac menu options.
   */
  protected void initMac()
  {
        try 
        { 
          // Special changes for Java 9...
          Desktop desktop = Desktop.getDesktop();
          desktop.setAboutHandler(new AboutHandler() 
          {
          	@Override
  			public void handleAbout(AboutEvent e)
  			{
  				about();
  			}
          });
          
          desktop.setQuitHandler(new QuitHandler()
          {
  			@Override
  			public void handleQuitRequestWith(QuitEvent e, QuitResponse response)
  			{
  				System.exit(0);
  			}
          });
          
        } 
        catch (Exception e)
        {
      	  e.printStackTrace();
        }
  }
  
  /**
   * Display aboout dialog.
   */
  protected void about()
  {
		final JOptionPane pane = new JOptionPane("Circular-Plot\nthis is free software and is distributed"
				+ "\nunder the terms of the GNU General Public License.", JOptionPane.INFORMATION_MESSAGE);
		final JDialog d = pane.createDialog((JFrame) null, "About");
		d.setLocation(10, 10);
		d.setVisible(true);
  }
}

