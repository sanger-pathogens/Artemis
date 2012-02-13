/* BamView
 *
 * created: 2009
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2009  Genome Research Limited
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
package uk.ac.sanger.artemis.components.alignment;


import java.awt.AlphaComposite;
import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Composite;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Insets;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.Stroke;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.lang.management.ManagementFactory;
import java.lang.management.MemoryMXBean;
import java.net.URL;
import java.util.Collections;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.List;
import java.util.Vector;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JScrollBar;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JTextField;
import javax.swing.UIManager;
import javax.swing.border.Border;
import javax.swing.border.EmptyBorder;
import javax.swing.text.JTextComponent;

import org.apache.log4j.Level;

import net.sf.picard.sam.BuildBamIndex;
import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMException;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.util.CloseableIterator;

import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.Selection;
import uk.ac.sanger.artemis.SelectionChangeEvent;
import uk.ac.sanger.artemis.SelectionChangeListener;
import uk.ac.sanger.artemis.SimpleEntryGroup;
import uk.ac.sanger.artemis.components.DisplayAdjustmentEvent;
import uk.ac.sanger.artemis.components.DisplayAdjustmentListener;
import uk.ac.sanger.artemis.components.EntryEdit;
import uk.ac.sanger.artemis.components.EntryFileDialog;
import uk.ac.sanger.artemis.components.FeatureDisplay;
import uk.ac.sanger.artemis.components.FileViewer;
import uk.ac.sanger.artemis.components.MessageDialog;
import uk.ac.sanger.artemis.components.NonModalDialog;
import uk.ac.sanger.artemis.components.SwingWorker;
import uk.ac.sanger.artemis.components.genebuilder.AutoCompleteComboDocument;
import uk.ac.sanger.artemis.components.variant.FeatureContigPredicate;
import uk.ac.sanger.artemis.editor.MultiLineToolTipUI;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.sequence.Bases;
import uk.ac.sanger.artemis.sequence.MarkerRange;
import uk.ac.sanger.artemis.sequence.NoSequenceException;
import uk.ac.sanger.artemis.util.Document;
import uk.ac.sanger.artemis.util.DocumentFactory;
import uk.ac.sanger.artemis.util.FTPSeekableStream;
import uk.ac.sanger.artemis.util.OutOfRangeException;

public class BamView extends JPanel
                     implements DisplayAdjustmentListener, SelectionChangeListener
{
  private static final long serialVersionUID = 1L;

  private List<SAMRecord> readsInView;
  private Hashtable<String, SAMFileReader> samFileReaderHash = new Hashtable<String, SAMFileReader>();

  private Hashtable<String, Integer> seqLengths = new Hashtable<String, Integer>();
  private Hashtable<String, Integer> offsetLengths;
  private Vector<String> seqNames = new Vector<String>();
  protected List<String> bamList;
  private List<Integer> hideBamList = new Vector<Integer>();

  private SAMRecordPredicate samRecordFlagPredicate;
  private SAMRecordMapQPredicate samRecordMapQPredicate;

  private SAMRecordFilter filterFrame;

  private Bases bases;
  private JScrollPane jspView;
  private JScrollBar scrollBar;
  
  private JComboBox combo;
  private boolean isOrientation = false;
  private boolean isSingle = false;
  private boolean isSNPs = false;
  
  private boolean isCoverage = false;
  private boolean isSNPplot = false;
  
  private EntryEdit entry_edit;
  private FeatureDisplay feature_display;
  private Selection selection;
  private JPanel mainPanel = new JPanel();
  private CoveragePanel coveragePanel;
  private SnpPanel snpPanel;

  private boolean logScale = false;
  private Ruler ruler;
  private int nbasesInView;
  
  private int startBase = -1;
  private int endBase   = -1;
  private int laststart;
  private int lastend;

  private boolean asynchronous = true;
  private boolean showBaseAlignment = false;
  
  private JMenu bamFilesMenu = new JMenu("BAM files");
  private JCheckBoxMenuItem logMenuItem = new JCheckBoxMenuItem("Use Log Scale", logScale);
  
  private JCheckBoxMenuItem cbStackView = new JCheckBoxMenuItem("Stack", true);
  private JCheckBoxMenuItem cbPairedStackView = new JCheckBoxMenuItem("Paired Stack");
  private JCheckBoxMenuItem cbStrandStackView = new JCheckBoxMenuItem("Strand Stack");
  private JCheckBoxMenuItem cbIsizeStackView = new JCheckBoxMenuItem("Inferred Size", false);
  private JCheckBoxMenuItem cbCoverageView = new JCheckBoxMenuItem("Coverage", false);
  private JCheckBoxMenuItem cbLastSelected;
  
  private ButtonGroup buttonGroup = new ButtonGroup();
  
  private JCheckBoxMenuItem colourByCoverageColour = new JCheckBoxMenuItem("Coverage Plot Colours");
  private JCheckBoxMenuItem baseQualityColour = new JCheckBoxMenuItem("Base Quality");
  private JCheckBoxMenuItem markInsertions = new JCheckBoxMenuItem("Mark Insertions", true);
  private AlphaComposite translucent = 
    AlphaComposite.getInstance(AlphaComposite.SRC_OVER, 0.6f);
  
  private CoveragePanel coverageView = new CoveragePanel();
  
  /** Used to colour the frames. */
  private static Color LIGHT_GREY = new Color(200, 200, 200);
  private static Color DARK_GREEN = new Color(0, 150, 0);
  private static Color DARK_ORANGE = new Color(255,140,0);
  private static Color DEEP_PINK   = new Color(139,10,80);
  
  private Point lastMousePoint = null;
  private SAMRecord mouseOverSAMRecord = null;
  private SAMRecord highlightSAMRecord = null;
  private String mouseOverInsertion;
  // record of where a mouse drag starts
  private int dragStart = -1;
  
  private static int MAX_BASES = 26000;
  private int maxHeight = 800;
  
  private boolean concatSequences = false;
  private int ALIGNMENT_PIX_PER_BASE;
  private int BASE_HEIGHT;

  private JPopupMenu popup;
  private PopupMessageFrame popFrame = new PopupMessageFrame();
  private PopupMessageFrame waitingFrame = new PopupMessageFrame("waiting...");
  
  public static org.apache.log4j.Logger logger4j = 
    org.apache.log4j.Logger.getLogger(BamView.class);
  
  public BamView(List<String> bamList, 
                String reference,
                int nbasesInView,
                final EntryEdit entry_edit,
                final FeatureDisplay feature_display,
                final Bases bases,
                final JPanel containerPanel,
                final JFrame frame)
  {
    this(bamList, reference, nbasesInView, feature_display, bases, containerPanel, frame);
    this.entry_edit = entry_edit;
  }
  
  public BamView(List<String> bamList, 
                 String reference,
                 int nbasesInView,
                 final FeatureDisplay feature_display,
                 final Bases bases,
                 final JPanel containerPanel,
                 final JFrame frame)
  {
    super();
    setBackground(Color.white);
    this.bamList = bamList;
    this.nbasesInView = nbasesInView;
    this.feature_display = feature_display;
    this.bases = bases;

    containerPanel.setLayout(new BoxLayout(containerPanel, BoxLayout.Y_AXIS));
    containerPanel.add(mainPanel);
    
    // filter out unmapped reads by default
    setSamRecordFlagPredicate(
        new SAMRecordFlagPredicate(SAMRecordFlagPredicate.READ_UNMAPPED_FLAG));
    
    if(reference != null)
    {
      EntryGroup entryGroup = new SimpleEntryGroup();
      try
      {
        getEntry(reference,entryGroup);
      }
      catch (NoSequenceException e)
      {
        e.printStackTrace();
      }
    }
    
    try
    {
      readHeaderPicard();
    }
    catch(java.lang.UnsupportedClassVersionError err)
    {
      JOptionPane.showMessageDialog(null, 
          "This requires Java 1.6 or higher.", 
          "Check Java Version", JOptionPane.WARNING_MESSAGE);
    }
    catch (IOException e)
    {
      e.printStackTrace();
    }

    final javax.swing.plaf.FontUIResource font_ui_resource =
      Options.getOptions().getFontUIResource();
    
    Enumeration<Object> keys = UIManager.getDefaults().keys();
    while(keys.hasMoreElements()) 
    {
      Object key = keys.nextElement();
      Object value = UIManager.get(key);
      if(value instanceof javax.swing.plaf.FontUIResource) 
        UIManager.put(key, font_ui_resource);
    }

    setFont(Options.getOptions().getFont());
    FontMetrics fm  = getFontMetrics(getFont());
    ALIGNMENT_PIX_PER_BASE = fm.charWidth('M');
    BASE_HEIGHT = fm.getMaxAscent();
    selection = new Selection(null);
    
    MultiLineToolTipUI.initialize();
    setToolTipText("");
    
    buttonGroup.add(cbStackView);
    buttonGroup.add(cbPairedStackView);
    buttonGroup.add(cbStrandStackView);
    buttonGroup.add(cbIsizeStackView);
    buttonGroup.add(cbCoverageView);
    addMouseListener(new PopupListener());

    jspView = new JScrollPane(this, 
        JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
        JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
    
    jspView.setViewportBorder(BorderFactory.createMatteBorder(0, 0, 1, 0, Color.DARK_GRAY));
    Border empty = new EmptyBorder(0,0,0,0);
    jspView.setBorder(empty);
    jspView.getVerticalScrollBar().setUnitIncrement(8);

    addBamToPanel(frame);
  }
  
  public String getToolTipText()
  {
    if(mouseOverSAMRecord == null)
      return null;
    
    String msg = 
        mouseOverSAMRecord.getReadName() + "\n" + 
        mouseOverSAMRecord.getAlignmentStart() + ".." +
        mouseOverSAMRecord.getAlignmentEnd() + "\nisize=" +
        mouseOverSAMRecord.getInferredInsertSize() + "\nmapq=" +
        mouseOverSAMRecord.getMappingQuality()+"\nrname="+
        mouseOverSAMRecord.getReferenceName();

    if( mouseOverSAMRecord.getReadPairedFlag() && 
        mouseOverSAMRecord.getProperPairFlag() && 
       !mouseOverSAMRecord.getMateUnmappedFlag())
    {
      msg = msg +
        "\nstrand (read/mate): "+
       (mouseOverSAMRecord.getReadNegativeStrandFlag() ? "-" : "+")+" / "+
       (mouseOverSAMRecord.getMateNegativeStrandFlag() ? "-" : "+");
    }
    else
      msg = msg +
        "\nstrand (read/mate): "+
       (mouseOverSAMRecord.getReadNegativeStrandFlag() ? "-" : "+");
    
    if(msg != null && mouseOverInsertion != null)
      msg = msg + "\nInsertion at:" +mouseOverInsertion;
    
    return msg;
  }
  
  /**
   * Get the BAM index file from the list
   * @param bam
   * @return
   * @throws IOException
   */
  private File getBamIndexFile(String bam) throws IOException
  {
    File bamIndexFile = null;
    if (bam.startsWith("http") || bam.startsWith("ftp"))
    {
      final URL urlBamIndexFile = new URL(bam + ".bai");
      InputStream is = urlBamIndexFile.openStream();

      // Create temp file.
      bamIndexFile = File.createTempFile(urlBamIndexFile.getFile().replaceAll(
          "[\\/\\s]", "_"), ".bai");
      bamIndexFile.deleteOnExit();

      FileOutputStream out = new FileOutputStream(bamIndexFile);
      int c;
      while ((c = is.read()) != -1)
        out.write(c);
      out.flush();
      out.close();
      is.close();

      logger4j.debug("create... " + bamIndexFile.getAbsolutePath());
    }
    else
      bamIndexFile = new File(bam + ".bai");

    return bamIndexFile;
  }
    
  /**
   * Get the SAM file reader.
   * @param bam
   * @return
   * @throws IOException
   */
  private SAMFileReader getSAMFileReader(final String bam) throws IOException
  {  
    if(samFileReaderHash.containsKey(bam))
      return samFileReaderHash.get(bam);

    File bamIndexFile = getBamIndexFile(bam);
    if(!bamIndexFile.exists())
    {
      try
      {
        logger4j.warn("Index file not found so using picard to index the BAM.");
        // Use Picard to index the file
        // requires reads to be sorted by coordinate
        new BuildBamIndex().instanceMain(
          new String[]{ "I="+bam, "O="+bamIndexFile.getAbsolutePath(), "MAX_RECORDS_IN_RAM=50000", "VALIDATION_STRINGENCY=SILENT" });
      }
      catch(SAMException e)
      {
        String ls = System.getProperty("line.separator");
        String msg = 
            "BAM index file is missing. The BAM file needs to be sorted and indexed"+ls+
            "This can be done using samtools (http://samtools.sf.net/):"+ls+ls+
            "samtools sort <in.bam> <out.prefix>"+ls+
            "samtools index <sorted.bam>";
        
        throw new SAMException(msg);
      }
    }
    
    final SAMFileReader samFileReader;
    
    if(bam.startsWith("ftp"))
    {
      FTPSeekableStream fss = new FTPSeekableStream(new URL(bam));
      samFileReader = new SAMFileReader(fss, bamIndexFile, false);
    }
    else if(!bam.startsWith("http"))
    {
      File bamFile = new File(bam);
      samFileReader = new SAMFileReader(bamFile, bamIndexFile);
    }
    else
    {
      final URL urlBamFile = new URL(bam);
      samFileReader = new SAMFileReader(urlBamFile, bamIndexFile, false);
    }
    samFileReader.setValidationStringency(ValidationStringency.SILENT);
    samFileReaderHash.put(bam, samFileReader);
    
    return samFileReader;
  }

  private void readHeaderPicard() throws IOException
  {
    String bam = bamList.get(0);
    final SAMFileReader inputSam = getSAMFileReader(bam);
    
    //final SAMFileReader inputSam = new SAMFileReader(bamFile, indexFile);
    SAMFileHeader header = inputSam.getFileHeader();
    List<SAMSequenceRecord> readGroups = header.getSequenceDictionary().getSequences();
    
    for(int i=0; i<readGroups.size(); i++)
    {
      seqLengths.put(readGroups.get(i).getSequenceName(),
                     readGroups.get(i).getSequenceLength());
      seqNames.add(readGroups.get(i).getSequenceName());
    }
    //inputSam.close();
  }

  
  /**
   * Read a SAM or BAM file.
   * @throws IOException 
   */
  private void readFromBamPicard(int start, int end, int bamIndex, float pixPerBase) 
          throws IOException
  {
    // Open the input file.  Automatically detects whether input is SAM or BAM
    // and delegates to a reader implementation for the appropriate format.
    String bam = bamList.get(bamIndex);  
    final SAMFileReader inputSam = getSAMFileReader(bam);
    
    //final SAMFileReader inputSam = new SAMFileReader(bamFile, indexFile);


    if(concatSequences)
    {
      int len = 0;
      int lastLen = 1;
      for(int i=0; i<seqNames.size(); i++)
      {
        int thisLength = seqLengths.get(seqNames.get(i));
        len += thisLength;

        if( (lastLen >= start && lastLen < end) ||
            (len >= start && len < end) ||
            (start >= lastLen && start < len) ||
            (end >= lastLen && end < len) )
        {
          int offset = getSequenceOffset(seqNames.get(i)); 
          int thisStart = start - offset;
          if(thisStart < 1)
            thisStart = 1;
          int thisEnd   = end - offset;
          if(thisEnd > thisLength)
            thisEnd = thisLength;
          
          //System.out.println("READ "+seqNames.get(i)+"  "+thisStart+".."+thisEnd);
          iterateOverBam(inputSam, seqNames.get(i), thisStart, thisEnd, bamIndex, pixPerBase, bam);
        }
        lastLen = len;
      }
    }
    else
    {
      String refName = (String) combo.getSelectedItem();
      iterateOverBam(inputSam, refName, start, end, bamIndex, pixPerBase, bam);
    }
    
    //inputSam.close();
    //System.out.println("readFromBamPicard "+start+".."+end);
    //System.out.println("Reads in view ... "+readsInView.size());
  }
  
  /**
   * Iterate over BAM file and load into the <code>List</code> of
   * <code>SAMRecord</code>.
   * @param inputSam
   * @param refName
   * @param start
   * @param end
   */
  private void iterateOverBam(final SAMFileReader inputSam, 
                              String refName, int start, int end,
                              int bamIndex, float pixPerBase,
                              String bam)
  { 
    boolean multipleBAM = false;
    if(bamList.size() > 1)
      multipleBAM = true;

    CloseableIterator<SAMRecord> it = inputSam.queryOverlapping(refName, start, end);
    MemoryMXBean memory = ManagementFactory.getMemoryMXBean();
    int checkMemAfter = 8000;
    int cnt = 0;
    int seqOffset = getSequenceOffset(refName);
    int offset = seqOffset- getBaseAtStartOfView();
    
    while ( it.hasNext() )
    {
      try
      {
        cnt++;
        SAMRecord samRecord = it.next();
        if( samRecordFlagPredicate == null ||
           !samRecordFlagPredicate.testPredicate(samRecord))
        {
          if(samRecordMapQPredicate == null ||
             samRecordMapQPredicate.testPredicate(samRecord))
          {
            if(multipleBAM)
              samRecord.setAttribute("FL", bamIndex);
            
            if(isCoverageView(pixPerBase))
              coverageView.addRecord(samRecord, offset, bam);
            
            if(isCoverage)
              coveragePanel.addRecord(samRecord, offset, bam);
            
            if(isSNPplot)
              snpPanel.addRecord(samRecord, seqOffset);
            
            if(!isCoverageView(pixPerBase))
              readsInView.add(samRecord);
          }
        }
        
        if(cnt > checkMemAfter)
        {
          cnt = 0;
          float heapFraction =
            (float)((float)memory.getHeapMemoryUsage().getUsed()/
                    (float)memory.getHeapMemoryUsage().getMax());
          logger4j.debug("Heap memory usage (used/max): "+heapFraction);
          
          if(readsInView.size() > checkMemAfter*2 && !waitingFrame.isVisible())
            waitingFrame.showWaiting("loading...", mainPanel);
          
          if(heapFraction > 0.90) 
          {
            popFrame.show(
              "Using > 90 % of the maximum memory limit:"+
              (memory.getHeapMemoryUsage().getMax()/1000000.f)+" Mb.\n"+
              "Not all reads in this range have been read in. Zoom in or\n"+
              "consider increasing the memory for this application.",
              mainPanel,
              15000);
            break;
          }
        }
      }
      catch(Exception e)
      {
        System.out.println(e.getMessage());
      }
    }

    it.close();
  }

  private int getSequenceLength()
  {
    if(concatSequences)
    {
      int len = 0;
      for(int i=0; i<seqNames.size(); i++)
        len += seqLengths.get(seqNames.get(i));
      return len;
    }
    else
      return seqLengths.get((String) combo.getSelectedItem());
  }
  
  /**
   * For BAM files with multiple references sequences, calculate
   * the offset from the start of the concatenated sequence for 
   * a given reference.
   * @param refName
   * @return
   */
  protected int getSequenceOffset(String refName)
  {
    if(!concatSequences)
      return 0;
    
    if(offsetLengths == null)
    {   
/*    offsetLengths = new Hashtable<String, Integer>(combo.getItemCount());
      int offset = 0;
      for(int i=0; i<combo.getItemCount(); i++)
      {
        String thisSeqName = (String) combo.getItemAt(i);
        offsetLengths.put(thisSeqName, offset);
        offset += seqLengths.get(combo.getItemAt(i));
      }*/

      FeatureVector features = feature_display.getEntryGroup().getAllFeatures();
      offsetLengths = new Hashtable<String, Integer>(seqNames.size());
      for(int i=0; i<seqNames.size(); i++)
      {
        FeatureContigPredicate predicate = new FeatureContigPredicate(seqNames.get(i).trim());
        for(int j=0; j<features.size(); j++)
        {
          if(predicate.testPredicate(features.elementAt(j)))
          {
            offsetLengths.put(seqNames.get(i), features.elementAt(j).getFirstBase()-1);
            break;
          }
        }
      }
      
      if(offsetLengths.size() != seqNames.size())
      {
        System.err.println("Found: "+offsetLengths.size() +" of "+ seqNames.size());
        JOptionPane.showMessageDialog(this, 
            "There is a problem matching the reference sequences\n"+
            "to the names in the BAM file. This may mean the labels\n"+
            "on the reference features do not match those in the in\n"+
            "the BAM file.", 
            "Problem Found", JOptionPane.WARNING_MESSAGE);
        concatSequences = false;
        return 0;
      }
    }
    return offsetLengths.get(refName);
  }
  
  /**
   * Override
   */
  protected void paintComponent(Graphics g)
  {
	super.paintComponent(g);
	Graphics2D g2 = (Graphics2D)g;

	mouseOverSAMRecord = null;
    int seqLength = getSequenceLength();
	float pixPerBase = getPixPerBaseByWidth();
	
    int start;
    int end;
    
    if(startBase > 0)
      start = startBase;
    else
      start = getBaseAtStartOfView();
    
    if(endBase > 0)
      end = endBase;
    else
    {
      end   = start + nbasesInView - 1;
      if(end > seqLength)
        end = seqLength;
    }

    boolean changeToStackView = false;
    MemoryMXBean memory = ManagementFactory.getMemoryMXBean();
    if(laststart != start ||
       lastend   != end ||
       CoveragePanel.isRedraw())
    {
      if(isCoverageView(pixPerBase))
        coverageView.init(this, pixPerBase, start, end);
      if(isCoverage)
        coveragePanel.init(this, pixPerBase, start, end);
      if(isSNPplot)
        snpPanel.init(this, pixPerBase, start, end);

      synchronized (this)
      {
        try
        {
          float heapFractionUsedBefore = (float) ((float) memory.getHeapMemoryUsage().getUsed() / 
                                                  (float) memory.getHeapMemoryUsage().getMax());

          if(readsInView == null)
            readsInView = new Vector<SAMRecord>();
          else
            readsInView.clear();

          for(int i=0; i<bamList.size(); i++)
          {
            if(!hideBamList.contains(i))
              readFromBamPicard(start, end, i, pixPerBase);
          }
          float heapFractionUsedAfter = (float) ((float) memory.getHeapMemoryUsage().getUsed() / 
                                                 (float) memory.getHeapMemoryUsage().getMax());

          // System.out.println("Heap Max  : "+memory.getHeapMemoryUsage().getMax());
          // System.out.println("Heap Used : "+memory.getHeapMemoryUsage().getUsed());
          // System.out.println("Heap memory used "+heapFractionUsedAfter);

          if ((heapFractionUsedAfter - heapFractionUsedBefore) > 0.06
              && !isStackView() && heapFractionUsedAfter > 0.8)
          {
            cbStackView.setSelected(true);
            changeToStackView = true;
          }

          if((!isStackView() && !isStrandStackView()) || isBaseAlignmentView(pixPerBase))
          {
            Collections.sort(readsInView, new SAMRecordComparator());
          }
          else if( (isStackView() || isStrandStackView()) &&
              bamList.size() > 1)
          {
            // merge multiple BAM files
            Collections.sort(readsInView, new SAMRecordPositionComparator(BamView.this));
          }
        }
        catch (OutOfMemoryError ome)
        {
          JOptionPane.showMessageDialog(this, "Out of Memory");
          readsInView.clear();
          return;
        }
        catch(IOException me)
        {
          me.printStackTrace();
        }
        catch(net.sf.samtools.util.RuntimeIOException re)
        {
          JOptionPane.showMessageDialog(this, re.getMessage());
        }
      }
    }

    laststart = start;
    lastend   = end;
    
    if(showBaseAlignment)
	  drawBaseAlignment(g2, seqLength, pixPerBase, start, end);
	else
	{
	  if(isCoverageView(pixPerBase))
	    drawCoverage(g2,start, end, pixPerBase);
	  else if(isStackView())  
	    drawStackView(g2, seqLength, pixPerBase, start, end);
	  else if(isPairedStackView())
	    drawPairedStackView(g2, seqLength, pixPerBase, start, end);
	  else if(isStrandStackView())
	    drawStrandStackView(g2, seqLength, pixPerBase, start, end);
	  else
	    drawLineView(g2, seqLength, pixPerBase, start, end);
	}
    
    if(isCoverage)
      coveragePanel.repaint();
    if(isSNPplot)
      snpPanel.repaint();

	if(waitingFrame.isVisible())
      waitingFrame.hideFrame();
	if(changeToStackView)
	{
	  popFrame.show(
          "Note :: Changed to the stack view to save memory.\n"+
          "Currently this is using "+ 
          (memory.getHeapMemoryUsage().getUsed()/1000000.f)+" Mb "+
          "and the maximum\nmemory limit is "+
          (memory.getHeapMemoryUsage().getMax()/1000000.f)+" Mb.",
          mainPanel,
          15000);
	}
  }
  
  
  private float getPixPerBaseByWidth()
  {
    return (float)mainPanel.getWidth() / (float)nbasesInView;
  }
  
  
  private int getMaxBasesInPanel(int seqLength)
  {
    if(feature_display == null)
      return seqLength+nbasesInView/3;
    else
      return seqLength+nbasesInView;
  }
  
  /**
   * Draw the zoomed-in base view.
   * @param g2
   * @param seqLength
   * @param pixPerBase
   * @param start
   * @param end
   */
  private void drawBaseAlignment(Graphics2D g2, 
                                 int seqLength, 
                                 float pixPerBase, 
                                 final int start, 
                                 int end)
  {
    ruler.start = start;
    ruler.end = end;
    ruler.repaint();
    
    int ypos = 0;

    String refSeq = null;
    int refSeqStart = start;
    
    end = start + ( mainPanel.getWidth() * ALIGNMENT_PIX_PER_BASE );
    if(bases != null)
    {
      // draw the reference sequence
      ypos+=11;

      try
      {
        int seqEnd = end+2;
        if(seqEnd > bases.getLength())
          seqEnd = bases.getLength();

        if(refSeqStart < 1)
          refSeqStart = 1;
        refSeq = 
          bases.getSubSequence(new Range(refSeqStart, seqEnd), Bases.FORWARD).toUpperCase();
        
        g2.setColor(LIGHT_GREY);
        g2.fillRect(0, ypos-11, mainPanel.getWidth(), 11);
        drawSelectionRange(g2, ALIGNMENT_PIX_PER_BASE, start, end);
        g2.setColor(Color.black);
        g2.drawString(refSeq, 0, ypos);
      }
      catch (OutOfRangeException e)
      {
        e.printStackTrace();
      }
    }
    else
      drawSelectionRange(g2, ALIGNMENT_PIX_PER_BASE, start, end);

    g2.setStroke(new BasicStroke (2.f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_ROUND));
    
    boolean drawn[] = new boolean[readsInView.size()];
    for(int i=0; i<readsInView.size(); i++)
      drawn[i] = false;
    
    Rectangle r = jspView.getViewport().getViewRect();
    int nreads = readsInView.size();
    
    for (int i = 0; i < nreads; i++)
    {
      try
      {
        if (!drawn[i])
        {
          ypos += 11;

          SAMRecord thisRead = readsInView.get(i);
          if (ypos < r.getMaxY() || ypos > r.getMinY())
            drawSequence(g2, thisRead, ypos, refSeq, refSeqStart);
          drawn[i] = true;

          int thisEnd = thisRead.getAlignmentEnd();
          if (thisEnd == 0)
            thisEnd = thisRead.getAlignmentStart() + thisRead.getReadLength();

          for (int j = i + 1; j < nreads; j++)
          {
            if (!drawn[j])
            {
              SAMRecord nextRead = readsInView.get(j);
              int nextStart = nextRead.getAlignmentStart();
              if (nextStart > thisEnd + 1)
              {
                if (ypos < r.getMaxY() || ypos > r.getMinY())
                  drawSequence(g2, nextRead, ypos, refSeq, refSeqStart);

                drawn[j] = true;
                thisEnd = nextRead.getAlignmentEnd();
                if (thisEnd == 0)
                  thisEnd = nextStart + nextRead.getReadLength();
              }
              else if (ypos > r.getMaxY() || ypos < r.getMinY())
                break;
            }
          }
        }
      }
      catch (ArrayIndexOutOfBoundsException ae)
      {
        System.err.println(readsInView.size()+"  "+nreads);
        ae.printStackTrace();
      }
    }
    
    if(ypos > getHeight())
    {
      Dimension d = getPreferredSize();
      d.setSize(getPreferredSize().getWidth(), ypos);
      setPreferredSize(d);
      revalidate();
    }
  }

  
  /**
   * Draw the query sequence
   * @param g2
   * @param read
   * @param pixPerBase
   * @param ypos
   */
  private void drawSequence(Graphics2D g2, SAMRecord samRecord, 
                            int ypos, String refSeq, int refSeqStart)
  {
    if (!samRecord.getReadPairedFlag() ||  // read is not paired in sequencing
        samRecord.getMateUnmappedFlag() )  // mate is unmapped )  // mate is unmapped 
      g2.setColor(Color.black);
    else
      g2.setColor(Color.blue);
    
    Color col = g2.getColor();
    int xpos;
    int len    = 0;
    int refPos = 0;
    String readSeq = samRecord.getReadString();
    int offset = getSequenceOffset(samRecord.getReferenceName());

    byte[] phredQuality = null;
    if(baseQualityColour.isSelected())
      phredQuality = samRecord.getBaseQualities();

    Hashtable<Integer, String> insertions = null;
    List<AlignmentBlock> blocks = samRecord.getAlignmentBlocks();
    for(int i=0; i<blocks.size(); i++)
    {
      AlignmentBlock block = blocks.get(i);
      int blockStart = block.getReadStart();
      len += block.getLength();
      for(int j=0; j<block.getLength(); j++)
      {
        int readPos = blockStart-1+j;
        xpos = block.getReferenceStart() - 1 + j + offset;
        refPos = xpos - refSeqStart + 1;

        if(phredQuality != null)
          setColourByBaseQuality(g2, phredQuality[readPos]);

        if(isSNPs && refSeq != null && refPos > 0 && refPos < refSeq.length())
        { 
          if(readSeq.charAt(readPos) != refSeq.charAt(refPos))
            g2.setColor(Color.red);
          else
            g2.setColor(col);
        }

        g2.drawString(readSeq.substring(readPos, readPos+1), 
                      refPos*ALIGNMENT_PIX_PER_BASE, ypos);
        

      }
          
      // look for insertions
      if(markInsertions.isSelected() && i < blocks.size()-1)
      {
        int blockEnd = blockStart+block.getLength();
        int nextBlockStart = blocks.get(i+1).getReadStart();
        int insertSize = nextBlockStart - blockEnd;
        if(insertSize > 0)
        {
          if(insertions == null)
            insertions = new Hashtable<Integer, String>();

          g2.setColor(DEEP_PINK);

          int xscreen = (refPos+1)*ALIGNMENT_PIX_PER_BASE;
          insertions.put(xscreen, 
              (refPos+refSeqStart+1)+" "+
              readSeq.substring(blockEnd-1, nextBlockStart-1));
          g2.drawLine(xscreen, ypos, xscreen, ypos-BASE_HEIGHT);
          
          // mark on reference sequence as well
          if(bases != null)
            g2.drawLine(xscreen, 11, xscreen, 11-BASE_HEIGHT);
          g2.setColor(col);
        }
      }
      
      // highlight
      if(highlightSAMRecord != null &&
         highlightSAMRecord.getReadName().equals(samRecord.getReadName()))
      {
        refPos =  block.getReferenceStart() + offset - refSeqStart;
        int xstart = refPos*ALIGNMENT_PIX_PER_BASE;
        int width  = block.getLength()*ALIGNMENT_PIX_PER_BASE;
        Color col1 = g2.getColor();
        g2.setColor(Color.red);
        g2.drawRect(xstart, ypos-BASE_HEIGHT, width, BASE_HEIGHT);        
        if(i < blocks.size()-1)
        {
          int nextStart = 
            (blocks.get(i+1).getReferenceStart() + offset - refSeqStart)*ALIGNMENT_PIX_PER_BASE;
          g2.drawLine(xstart+width, ypos-(BASE_HEIGHT/2), nextStart, ypos-(BASE_HEIGHT/2));
        }
        
        g2.setColor(col1);
      }
      else if(i < blocks.size()-1)
      {
        refPos =  block.getReferenceStart() + offset - refSeqStart;
        int xstart = refPos*ALIGNMENT_PIX_PER_BASE;
        int width  = block.getLength()*ALIGNMENT_PIX_PER_BASE;
        int nextStart = 
          (blocks.get(i+1).getReferenceStart() + offset - refSeqStart)*ALIGNMENT_PIX_PER_BASE;
        g2.drawLine(xstart+width, ypos-(BASE_HEIGHT/2), nextStart, ypos-(BASE_HEIGHT/2));
      }
    }

    if(lastMousePoint != null && blocks.size() > 0)
    {
      refPos = blocks.get(0).getReferenceStart()+offset-refSeqStart;
      int xstart = refPos*ALIGNMENT_PIX_PER_BASE;
      
      refPos = blocks.get(blocks.size()-1).getReferenceStart()+
               blocks.get(blocks.size()-1).getLength()+offset-refSeqStart;
      int xend   = (refPos+len)*ALIGNMENT_PIX_PER_BASE;

      if(lastMousePoint.getY() > ypos-11 && lastMousePoint.getY() < ypos)
      if(lastMousePoint.getX() > xstart &&
         lastMousePoint.getX() < xend)
      {
        mouseOverSAMRecord = samRecord;

        if(insertions != null)
          mouseOverInsertion = insertions.get((int)lastMousePoint.getX());
      }
    }
  }
  
  /**
   * Colour bases on their mapping quality.
   * @param g2
   * @param baseQuality
   */
  private void setColourByBaseQuality(Graphics2D g2, byte baseQuality)
  {
    if (baseQuality < 10)
      g2.setColor(Color.blue);
    else if (baseQuality < 20)
      g2.setColor(DARK_GREEN);
    else if (baseQuality < 30)
      g2.setColor(DARK_ORANGE);
    else
      g2.setColor(Color.black);
  }
  
  /**
   * Draw zoomed-out view.
   * @param g2
   * @param seqLength
   * @param pixPerBase
   * @param start
   * @param end
   */
  private void drawLineView(Graphics2D g2, int seqLength, float pixPerBase, int start, int end)
  {
    drawSelectionRange(g2, pixPerBase,start, end);
    if(isShowScale())
      drawScale(g2, start, end, pixPerBase, getHeight());
    
    Stroke stroke =
      new BasicStroke (1.3f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_ROUND);
    g2.setStroke(stroke);
    
    int scaleHeight;
    if(isShowScale())
      scaleHeight = 15;
    else
      scaleHeight = 0;
    
    int baseAtStartOfView = getBaseAtStartOfView();
    Rectangle r = jspView.getViewport().getViewRect();
    
    for(int i=0; i<readsInView.size(); i++)
    {
      SAMRecord samRecord = readsInView.get(i);
      SAMRecord samNextRecord = null;      

      if( !samRecord.getReadPairedFlag() ||  // read is not paired in sequencing
          samRecord.getMateUnmappedFlag() )  // mate is unmapped
      {
        if(isSingle)
        {
          int ypos = getYPos(scaleHeight, samRecord.getReadString().length()); // (getHeight() - scaleHeight) - samRecord.getReadString().length();
          if(ypos > r.getMaxY() || ypos < r.getMinY())
            continue;
          
          g2.setColor(Color.black);
          drawRead(g2, samRecord, pixPerBase, ypos, baseAtStartOfView);
        }
        continue;
      }

      int ypos = getYPos(scaleHeight, Math.abs(samRecord.getInferredInsertSize()));
      if( (ypos > r.getMaxY() || ypos < r.getMinY()) && ypos > 0 )
        continue;
      
      if(i < readsInView.size()-1)
      {
        samNextRecord = readsInView.get(++i);

        if(samRecord.getReadName().equals(samNextRecord.getReadName()))
        { 
          // draw connection between paired reads
          if(samRecord.getAlignmentEnd() < samNextRecord.getAlignmentStart() && 
              (samNextRecord.getAlignmentStart()-samRecord.getAlignmentEnd())*pixPerBase > 2.f)
          {
        	g2.setColor(Color.LIGHT_GRAY);

            int offset1 = getSequenceOffset(samRecord.getReferenceName());
            int end1   = samRecord.getAlignmentEnd()+offset1-baseAtStartOfView;
            
            int offset2 = getSequenceOffset(samNextRecord.getReferenceName());
            int start2  = samNextRecord.getAlignmentStart()+offset2-baseAtStartOfView;
            
            drawTranslucentLine(g2, 
                   (int)(end1*pixPerBase), (int)(start2*pixPerBase), ypos);
          }
          
          if(colourByCoverageColour.isSelected())
            g2.setColor(getColourByCoverageColour(samRecord));
          else if( (samRecord.getReadNegativeStrandFlag() && // strand of the query (1 for reverse)
                    samNextRecord.getReadNegativeStrandFlag()) ||
                   (!samRecord.getReadNegativeStrandFlag() && 
                    !samNextRecord.getReadNegativeStrandFlag()))
            g2.setColor(Color.red);
          else
            g2.setColor(Color.blue);

          drawRead(g2, samRecord, pixPerBase, ypos, baseAtStartOfView);
          drawRead(g2, samNextRecord, pixPerBase, ypos, baseAtStartOfView);
        }
        else
        {
          drawLoneRead(g2, samRecord, ypos, pixPerBase, baseAtStartOfView, scaleHeight);
          i--;
        }
      }
      else
      {
        drawLoneRead(g2, samRecord, ypos, pixPerBase, baseAtStartOfView, scaleHeight);
      }
    }
    
    drawYScale(g2, scaleHeight);
  }
  
  private int getYPos(int scaleHeight, int size)
  {
    int ypos;
    
    if(!logScale)
      ypos = (getHeight() - scaleHeight) - size;
    else
    {
      int logInfSize = (int)( Math.log(size) * 100);
      ypos = (getHeight() - scaleHeight) - logInfSize;
    }
    
    return ypos;
  }
 
  /**
   * Draw the reads as lines in vertical stacks. The reads are colour 
   * coded as follows:
   * 
   * blue  - reads are unique and are paired with a mapped mate
   * black - reads are unique and are not paired or have an unmapped mate
   * green - reads are duplicates
   * 
   * @param g2
   * @param seqLength
   * @param pixPerBase
   * @param start
   * @param end
   */
  private void drawStackView(Graphics2D g2, 
                             int seqLength, 
                             float pixPerBase, 
                             int start, 
                             int end)
  {
    drawSelectionRange(g2, pixPerBase,start, end);
    if(isShowScale())
      drawScale(g2, start, end, pixPerBase, getHeight());

    BasicStroke stroke = new BasicStroke(
        1.3f,
        BasicStroke.CAP_BUTT, 
        BasicStroke.JOIN_MITER);
    g2.setStroke(stroke);
    
    int scaleHeight;
    if(isShowScale())
      scaleHeight = 15;
    else
      scaleHeight = 0;
    
    int ypos = (getHeight() - scaleHeight);
    int ydiff = 2;
    if(isOrientation)
      ydiff= 4;
    int maxEnd = 0;
    int lstStart = 0;
    int lstEnd = 0;
    int baseAtStartOfView = getBaseAtStartOfView();
    g2.setColor(Color.blue);
    Rectangle r = jspView.getViewport().getViewRect();
    
    for(int i=0; i<readsInView.size(); i++)
    {
      SAMRecord samRecord = readsInView.get(i);
      int offset = getSequenceOffset(samRecord.getReferenceName());

      int recordStart = samRecord.getAlignmentStart()+offset;
      int recordEnd = samRecord.getAlignmentEnd()+offset;
      
      if(colourByCoverageColour.isSelected() ||
         lstStart != recordStart || lstEnd != recordEnd)
      { 
        if(colourByCoverageColour.isSelected())
          g2.setColor(getColourByCoverageColour(samRecord));
        else if (!samRecord.getReadPairedFlag() ||   // read is not paired in sequencing
                  samRecord.getMateUnmappedFlag() )  // mate is unmapped )  // mate is unmapped 
          g2.setColor(Color.black);
        else
          g2.setColor(Color.blue);
        
        if(maxEnd < recordStart)
        {
          ypos = (getHeight() - scaleHeight)-ydiff;
          maxEnd = recordEnd+2;
        }
        else
          ypos = ypos-ydiff;
      }
      else
        g2.setColor(DARK_GREEN);

      lstStart = recordStart;
      lstEnd   = recordEnd;
      
      if(ypos > r.getMaxY() || ypos < r.getMinY())
        continue;
      drawRead(g2, samRecord, pixPerBase, ypos, baseAtStartOfView);
    }
  }
  
  /**
   * Draw the reads as lines in vertical stacks. The reads are colour 
   * coded as follows:
   * 
   * blue  - reads are unique and are paired with a mapped mate
   * black - reads are unique and are not paired or have an unmapped mate
   * green - reads are duplicates
   * 
   * @param g2
   * @param seqLength
   * @param pixPerBase
   * @param start
   * @param end
   */
  private void drawStrandStackView(Graphics2D g2, 
                                   int seqLength, 
                                   float pixPerBase, 
                                   int start, 
                                   int end)
  {
    drawSelectionRange(g2, pixPerBase,start, end);   
    BasicStroke stroke = new BasicStroke(
        1.3f,
        BasicStroke.CAP_BUTT, 
        BasicStroke.JOIN_MITER);
    
    int scaleHeight = 15;
    drawScale(g2, start, end, pixPerBase, ((getHeight()+scaleHeight)/2));

    int ymid = (getHeight()/ 2);
    int ydiff = 2;
    if(isOrientation)
      ydiff= 4;
    // positive strand    
    drawStrand(g2, false, scaleHeight, ymid-(scaleHeight/2), -ydiff, pixPerBase, stroke);
    
    // negative strand
    drawStrand(g2, true, scaleHeight, ymid+(scaleHeight/2), ydiff, pixPerBase, stroke);
  }
  
  private void drawStrand(Graphics2D g2, 
                          boolean isStrandNegative, 
                          int scaleHeight,
                          int ymid,
                          int ystep,
                          float pixPerBase,
                          Stroke stroke)
  {
    int ypos = (getHeight() - scaleHeight);
    int maxEnd = 0;
    int lstStart = 0;
    int lstEnd = 0;
    int baseAtStartOfView = getBaseAtStartOfView();
    g2.setColor(Color.blue);
    Rectangle r = jspView.getViewport().getViewRect();
    
    for(int i=0; i<readsInView.size(); i++)
    {
      SAMRecord samRecord = readsInView.get(i);
      
      if( samRecord.getReadNegativeStrandFlag() == isStrandNegative )
      {
        int offset = getSequenceOffset(samRecord.getReferenceName());
        int recordStart = samRecord.getAlignmentStart()+offset;
        int recordEnd   = samRecord.getAlignmentEnd()+offset;
      
        if(colourByCoverageColour.isSelected() ||
            lstStart != recordStart || lstEnd != recordEnd)
        { 
          if(colourByCoverageColour.isSelected())
            g2.setColor(getColourByCoverageColour(samRecord));
          else if (!samRecord.getReadPairedFlag() ||   // read is not paired in sequencing
                    samRecord.getMateUnmappedFlag() )  // mate is unmapped 
            g2.setColor(Color.black);
          else
            g2.setColor(Color.blue);
        
          if(maxEnd < recordStart)
          {
            ypos = ymid + ystep;
            maxEnd = recordEnd+2;
          }
          else
            ypos = ypos + ystep;
        }
        else
          g2.setColor(DARK_GREEN);

        lstStart = recordStart;
        lstEnd   = recordEnd;
        
        if(ypos > r.getMaxY() || ypos < r.getMinY())
          continue;
        drawRead(g2, samRecord, pixPerBase, ypos, baseAtStartOfView);
      }
    }
  }
  
  /**
   * Draw the reads as lines in vertical stacks. The reads are colour 
   * coded as follows:
   * 
   * blue  - reads are unique and are paired with a mapped mate
   * black - reads are unique and are not paired or have an unmapped mate
   * green - reads are duplicates
   * 
   * @param g2
   * @param seqLength
   * @param pixPerBase
   * @param start
   * @param end
   */
  private void drawPairedStackView(Graphics2D g2, 
                                   int seqLength, 
                                   float pixPerBase, 
                                   int start, 
                                   int end)
  {
    drawSelectionRange(g2, pixPerBase,start, end);
    if(isShowScale())
      drawScale(g2, start, end, pixPerBase, getHeight());
    
    Vector<PairedRead> pairedReads = new Vector<PairedRead>();
    for(int i=0; i<readsInView.size(); i++)
    {
      SAMRecord samRecord = readsInView.get(i);

      if( !samRecord.getReadPairedFlag() ||  // read is not paired in sequencing
          samRecord.getMateUnmappedFlag() )  // mate is unmapped
        continue;

      SAMRecord samNextRecord = null;      
      if(i < readsInView.size()-1)
      {
        samNextRecord = readsInView.get(++i);
        PairedRead pr = new PairedRead();
        if(samRecord.getReadName().equals(samNextRecord.getReadName()))
        { 
          if(samRecord.getAlignmentStart() < samNextRecord.getAlignmentStart())
          {
            pr.sam1 = samRecord;
            pr.sam2 = samNextRecord;
          }
          else
          {
            pr.sam2 = samRecord;
            pr.sam1 = samNextRecord;
          }
          
        }
        else
        {
          --i;
          pr.sam1 = samRecord;
          pr.sam2 = null;
        }
        pairedReads.add(pr);
      }
    }
    Collections.sort(pairedReads, new PairedReadComparator());
    
    Stroke originalStroke = new BasicStroke (1.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_ROUND); 
    Stroke stroke =
            new BasicStroke (1.3f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_ROUND);
    g2.setStroke(stroke);
    
    int scaleHeight;
    if(isShowScale())
      scaleHeight = 15;
    else
      scaleHeight = 0;
    
    int ydiff = 3;
    if(isOrientation)
      ydiff= 5;
    int ypos = getHeight() - scaleHeight - ydiff;
    int lastEnd = 0;
    int baseAtStartOfView = getBaseAtStartOfView();
    Rectangle r = jspView.getViewport().getViewRect();
    
    for(int i=0; i<pairedReads.size(); i++)
    {
      PairedRead pr = pairedReads.get(i);
      
      if(pr.sam1.getAlignmentStart() > lastEnd)
      {
        ypos = getHeight() - scaleHeight - ydiff;
        
        if(pr.sam2 != null)
        {  
          lastEnd = pr.sam2.getAlignmentEnd();
        }
        else
          lastEnd = pr.sam1.getAlignmentEnd();
      }
      else
        ypos = ypos - ydiff;
      
      if(ypos > r.getMaxY() || ypos < r.getMinY())
        continue;
      
      g2.setStroke(originalStroke);
      g2.setColor(Color.gray);
      
      if(pr.sam2 != null)
      {
        int offset1 = getSequenceOffset(pr.sam1.getReferenceName());
        int offset2 = getSequenceOffset(pr.sam2.getReferenceName());
        drawTranslucentJointedLine(g2, 
                (int)((pr.sam1.getAlignmentEnd()+offset1-getBaseAtStartOfView())*pixPerBase),
                (int)((pr.sam2.getAlignmentStart()+offset2-getBaseAtStartOfView())*pixPerBase), ypos);
      }
      else
      {
        if(!pr.sam1.getMateUnmappedFlag() &&
            pr.sam1.getMateReferenceName().equals(pr.sam1.getReferenceName()))
        {
          int prStart;
          
          if(pr.sam1.getAlignmentStart() > pr.sam1.getMateAlignmentStart())
            prStart = pr.sam1.getAlignmentEnd();
          else
            prStart = pr.sam1.getAlignmentStart();
          
          int offset = getSequenceOffset(pr.sam1.getReferenceName());
          drawTranslucentJointedLine(g2, 
              (int)( (prStart+offset-getBaseAtStartOfView())*pixPerBase),
              (int)( (pr.sam1.getMateAlignmentStart()+offset-getBaseAtStartOfView())*pixPerBase), ypos);
        }
      }
      
      if(colourByCoverageColour.isSelected())
        g2.setColor(getColourByCoverageColour(pr.sam1));
      else if( pr.sam1.getReadNegativeStrandFlag() && // strand of the query (1 for reverse)
          ( pr.sam2 != null && pr.sam2.getReadNegativeStrandFlag() ) )
        g2.setColor(Color.red);
      else
        g2.setColor(Color.blue);
      
      drawRead(g2, pr.sam1, pixPerBase, ypos, baseAtStartOfView);
      
      if(pr.sam2 != null)
        drawRead(g2, pr.sam2, pixPerBase, ypos, baseAtStartOfView);
    }
  }
  
  /**
   * Draw the read coverage.
   * @param g2
   * @param start
   * @param end
   * @param pixPerBase
   */
  private void drawCoverage(Graphics2D g2, int start, int end, float pixPerBase)
  {
    int scaleHeight = 0;
    if(isShowScale())
    {
      drawScale(g2, start, end, pixPerBase, getHeight());
      scaleHeight = 15;
    }

    int hgt = jspView.getVisibleRect().height-scaleHeight;
    g2.translate(0, getHeight()-hgt-scaleHeight);
    coverageView.draw(g2, getWidth(), hgt);
    coverageView.drawMax(g2);  
  }
  
  /**
   * Draw a read that apparently has a read mate that is not in view.
   * @param g2
   * @param thisRead
   * @param ypos
   * @param pixPerBase
   * @param originalStroke
   * @param stroke
   */
  private void drawLoneRead(Graphics2D g2, SAMRecord samRecord, int ypos, 
      float pixPerBase, int baseAtStartOfView, int scaleHeight)
  {
    boolean offTheTop = false;
    int offset = getSequenceOffset(samRecord.getReferenceName());
    int thisStart = samRecord.getAlignmentStart()+offset;
    int thisEnd   = thisStart + samRecord.getReadString().length() -1;
    
    if(ypos <= 0)
    {
      offTheTop = true;
      ypos = samRecord.getReadString().length();
    }
    
    if(samRecord.getInferredInsertSize() == 0)
    {
      offTheTop = true;
      ypos = getHeight() - scaleHeight - 5;
    }
      
    if(samRecord.getInferredInsertSize() != 0 &&
      Math.abs(samRecord.getMateAlignmentStart()-samRecord.getAlignmentEnd())*pixPerBase > 2.f)
    {
      g2.setColor(Color.LIGHT_GRAY);
      
      if(samRecord.getAlignmentEnd() < samRecord.getMateAlignmentStart())
      {
        int nextStart = 
          (int)((samRecord.getMateAlignmentStart()-getBaseAtStartOfView()+offset)*pixPerBase);
        drawTranslucentLine(g2, 
          (int)((thisEnd-getBaseAtStartOfView())*pixPerBase), nextStart, ypos);
      }
      else
      {
        int nextStart = 
            (int)((samRecord.getMateAlignmentStart()-getBaseAtStartOfView()+offset)*pixPerBase);
        drawTranslucentLine(g2, 
            (int)((thisStart-getBaseAtStartOfView())*pixPerBase), nextStart, ypos);
      }
    }
    
    if(colourByCoverageColour.isSelected())
      g2.setColor(getColourByCoverageColour(samRecord));
    else if(offTheTop)
      g2.setColor(DARK_ORANGE); 
    else if(samRecord.getReadNegativeStrandFlag() &&
            samRecord.getMateNegativeStrandFlag()) // strand of the query (1 for reverse)
      g2.setColor(Color.red);
    else
      g2.setColor(Color.blue);
 
    drawRead(g2, samRecord, pixPerBase, ypos, baseAtStartOfView);
    
    if (isSNPs)
      showSNPsOnReads(g2, samRecord, pixPerBase, ypos, offset);
  }

  
  private void drawScale(Graphics2D g2, int start, int end, float pixPerBase, int ypos)
  {
    g2.setColor(Color.black);
    g2.drawLine( 0, ypos-14,
                 (int)((end - getBaseAtStartOfView())*pixPerBase),   ypos-14);
    int interval = end-start;
    
    if(interval > 256000)
      drawTicks(g2, start, end, pixPerBase, 512000, ypos);
    else if(interval > 64000)
      drawTicks(g2, start, end, pixPerBase, 12800, ypos);
    else if(interval > 16000)
      drawTicks(g2, start, end, pixPerBase, 3200, ypos);
    else if(interval > 4000)
      drawTicks(g2, start, end, pixPerBase, 800, ypos);
    else if(interval > 1000)
      drawTicks(g2, start, end, pixPerBase, 400, ypos);
    else
      drawTicks(g2, start, end, pixPerBase, 100, ypos);
  }
  
  private void drawTicks(Graphics2D g2, int start, int end, float pixPerBase, int division, int ypos)
  {
    int markStart = (Math.round(start/division)*division);
    
    if(markStart < 1)
      markStart = 1;
    
    int sm = markStart-(division/2);
    float x;
    if(sm > start)
    {
      x = (sm-getBaseAtStartOfView())*pixPerBase;
      g2.drawLine((int)x, ypos-14,(int)x, ypos-12);
    }
    
    for(int m=markStart; m<end; m+=division)
    {
      x = (m-getBaseAtStartOfView())*pixPerBase;
      g2.drawString(Integer.toString(m), x, ypos-1);
      g2.drawLine((int)x, ypos-14,(int)x, ypos-11);
      
      sm = m+(division/2);
      
      if(sm < end)
      {
        x = (sm-getBaseAtStartOfView())*pixPerBase;
        g2.drawLine((int)x, ypos-14,(int)x, ypos-12);
      }
      
      if(m == 1)
        m = 0;
    }
  }
  
  /**
   * Draw a y-scale for inferred size (isize) of reads.
   * @param g2
   * @param xScaleHeight
   */
  private void drawYScale(Graphics2D g2, int xScaleHeight)
  {
    g2.setColor(Color.black);
    int maxY = getPreferredSize().height-xScaleHeight;
    
    if(logScale)
    {
      int start = 10;
      int count = 0;
      int ypos = getYPos(xScaleHeight, start);
      
      while(ypos > 0 && count < 15 && start > 1)
      {
        g2.drawLine(0, ypos, 2, ypos);
        g2.drawString(Integer.toString(start), 3, ypos);
        start = start*5;
        ypos = getYPos(xScaleHeight, start);
        count++;
      }
      return;
    }
    
    for(int i=100; i<maxY; i+=100)
    {
      int ypos = getHeight()-i-xScaleHeight;
      g2.drawLine(0, ypos, 2, ypos);
      g2.drawString(Integer.toString(i), 3, ypos);
    }
  }
  
  /**
   * Draw a given read.
   * @param g2
   * @param thisRead
   * @param pixPerBase
   * @param ypos
   * @param baseAtStartOfView
   */
  private void drawRead(Graphics2D g2, SAMRecord thisRead,
		                float pixPerBase,
		                int ypos,
		                int baseAtStartOfView)
  {
    int offset = getSequenceOffset(thisRead.getReferenceName());

    int thisStart = thisRead.getAlignmentStart()+offset-baseAtStartOfView;
    int thisEnd   = thisRead.getAlignmentEnd()+offset-baseAtStartOfView;
    
    if(highlightSAMRecord != null && 
       highlightSAMRecord.getReadName().equals(thisRead.getReadName()))
    {
       Stroke originalStroke = g2.getStroke();
       Stroke stroke =
         new BasicStroke (3.f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_ROUND);
       g2.setStroke(stroke);
       Color c = g2.getColor();
       g2.setColor(Color.black);
       g2.drawLine((int)( thisStart * pixPerBase), ypos,
                   (int)( thisEnd * pixPerBase), ypos);
       g2.setColor(c);
       g2.setStroke(originalStroke);
    }

    if(thisRead.getCigar().getCigarElements().size() == 1)
      g2.drawLine((int)( thisStart * pixPerBase), ypos,
                  (int)( thisEnd * pixPerBase), ypos);
    else
    {
      List<AlignmentBlock> blocks = thisRead.getAlignmentBlocks();
      Color c = g2.getColor();
      int lastEnd = 0;
      for(int i=0; i<blocks.size(); i++)
      {
        AlignmentBlock block = blocks.get(i);
        int blockStart = block.getReferenceStart()+offset-baseAtStartOfView;
        int blockEnd = blockStart + block.getLength() - 1;

        g2.drawLine((int)( blockStart * pixPerBase), ypos,
                    (int)( blockEnd * pixPerBase), ypos);
        if(i > 0 && blockStart != lastEnd)
        {
          g2.setColor(Color.gray);
          g2.drawLine((int)( blockStart * pixPerBase), ypos,
                      (int)( lastEnd * pixPerBase), ypos);
          g2.setColor(c);
        }
        lastEnd = blockEnd;
      }
    }
    
    if(isOrientation)
      drawArrow(g2, thisRead, thisStart, thisEnd, pixPerBase, ypos);

    // test if the mouse is over this read
    if(lastMousePoint != null)
    {
      if(lastMousePoint.getY()+2 > ypos && lastMousePoint.getY()-2 < ypos)
      if(lastMousePoint.getX() > thisStart * pixPerBase &&
         lastMousePoint.getX() < thisEnd * pixPerBase)
      {
        mouseOverSAMRecord = thisRead;
      }
    }
    
    if (isSNPs)
      showSNPsOnReads(g2, thisRead, pixPerBase, ypos, offset);
  }
  
  /**
   * Draw arrow on the read to indicate orientation.
   * @param g2
   * @param thisRead
   * @param thisStart
   * @param thisEnd
   * @param pixPerBase
   * @param ypos
   */
  private void drawArrow(Graphics2D g2,
                         SAMRecord thisRead, 
                         int thisStart, 
                         int thisEnd, 
                         float pixPerBase, 
                         int ypos)
  {
    if(thisRead.getReadNegativeStrandFlag())
    {
      int apos = ypos + 2;
      g2.drawLine((int)( (thisStart+5) * pixPerBase), apos,
                  (int)( thisStart * pixPerBase), ypos);
    }
    else
    {
      int apos = ypos - 2;
      g2.drawLine((int)( (thisEnd-5) * pixPerBase), apos,
                  (int)( thisEnd * pixPerBase), ypos);
    }  
  }
  
  /**
   * Highlight a selected range
   * @param g2
   * @param pixPerBase
   * @param start
   * @param end
   */
  private void drawSelectionRange(Graphics2D g2, float pixPerBase, int start, int end)
  {
    if(getSelection() != null)
    {
      Range selectedRange = getSelection().getSelectionRange();

      if(selectedRange != null)
      {
        int rangeStart = selectedRange.getStart();
        int rangeEnd   = selectedRange.getEnd();
        
        if(end < rangeStart || start > rangeEnd)
          return;
        
        int x = (int) (pixPerBase*(rangeStart-getBaseAtStartOfView()));
        int width = (int) (pixPerBase*(rangeEnd-rangeStart+1));
        
        g2.setColor(Color.pink);
        g2.fillRect(x, 0, width, getHeight());
      }
    }
  }
  
  /**
   * Draw a translucent line
   * @param g2
   * @param start
   * @param end
   * @param ypos
   */
  private void drawTranslucentLine(Graphics2D g2, int start, int end, int ypos)
  {
    Composite origComposite = g2.getComposite();
    g2.setComposite(translucent);
    g2.drawLine(start, ypos, end, ypos);
    g2.setComposite(origComposite);
  }
  
  /**
   * Draw a translucent line
   * @param g2
   * @param start
   * @param end
   * @param ypos
   */
  private void drawTranslucentJointedLine(Graphics2D g2, int start, int end, int ypos)
  {
    Composite origComposite = g2.getComposite();
    g2.setComposite(translucent);
    
    int mid = (int) ((end-start)/2.f)+start;
    //g2.drawLine(start, ypos, end, ypos);
    g2.drawLine(start, ypos, mid, ypos-5);
    g2.drawLine(mid, ypos-5, end, ypos);
    g2.setComposite(origComposite);
  }
  
  /**
   * Display the SNPs for the given read.
   * @param g2
   * @param thisRead
   * @param pixPerBase
   * @param ypos
   */
  private void showSNPsOnReads(Graphics2D g2, SAMRecord thisRead,
                               float pixPerBase, int ypos, int offset)
  {
    int thisStart = thisRead.getAlignmentStart();
    int thisEnd   = thisRead.getAlignmentEnd();
    
    // use alignment blocks of the contiguous alignment of
    // subsets of read bases to a reference sequence
    List<AlignmentBlock> blocks = thisRead.getAlignmentBlocks();
    try
    {
      char[] refSeq = bases.getSubSequenceC(
          new Range(thisStart+offset, thisEnd+offset), Bases.FORWARD);
      byte[] readSeq = thisRead.getReadBases();

      Color col = g2.getColor();
      g2.setColor(Color.red);

      offset = offset - getBaseAtStartOfView();
      for(int i=0; i<blocks.size(); i++)
      {
        AlignmentBlock block = blocks.get(i);
        for(int j=0; j<block.getLength(); j++)
        {
          int readPos = block.getReadStart()-1+j;
          int refPos  = block.getReferenceStart()+j;

          if (Character.toUpperCase(refSeq[refPos-thisStart]) != readSeq[readPos])
          {
            g2.drawLine((int) ((refPos+offset) * pixPerBase), ypos + 2,
                        (int) ((refPos+offset) * pixPerBase), ypos - 2);
          }
        }
        
      }

      g2.setColor(col);
    }
    catch (OutOfRangeException e)
    {
      System.err.println(thisRead.getReadName()+" "+e.getMessage());
    }
  }
  
  /**
   * Add the alignment view to the supplied <code>JPanel</code> in
   * a <code>JScrollPane</code>.
   * @param mainPanel  panel to add the alignment to
   * @param frame
   * @param autohide automatically hide the top panel containing the buttons
   * @param feature_display
   */
  private void addBamToPanel(final JFrame frame)
  {
    final JComponent topPanel = bamTopPanel(frame);
    mainPanel.setPreferredSize(new Dimension(900, 400));
    
    setDisplay(1, nbasesInView, null);
    mainPanel.setLayout(new BorderLayout());

    if(topPanel instanceof JPanel)
      mainPanel.add(topPanel, BorderLayout.NORTH);
    mainPanel.add(jspView, BorderLayout.CENTER);

    JPanel bottomPanel = new JPanel(new BorderLayout());
    coveragePanel = new CoveragePanel(this);
    bottomPanel.add(coveragePanel, BorderLayout.CENTER);

    //
    snpPanel = new SnpPanel(this, bases);
    bottomPanel.add(snpPanel, BorderLayout.NORTH);
    
    if(feature_display == null)
    {
      scrollBar = new JScrollBar(JScrollBar.HORIZONTAL, 1, nbasesInView, 1,
          getMaxBasesInPanel(getSequenceLength()));
      scrollBar.setUnitIncrement(nbasesInView/20);
      scrollBar.addAdjustmentListener(new AdjustmentListener()
      {
        public void adjustmentValueChanged(AdjustmentEvent e)
        {
          repaint();

          if(isSNPplot)
            snpPanel.repaint();
          if(isCoverage)
            coveragePanel.repaint();
        }
      });
      bottomPanel.add(scrollBar, BorderLayout.SOUTH);
    }
    else
    {
      if(!concatSequences)
      {
        int seqLen = seqLengths.get((String) combo.getSelectedItem());
        int artemisSeqLen = feature_display.getSequenceLength();
        if(seqLen != artemisSeqLen)
        {
          int newIndex = -1;
          for(int i=0; i<seqNames.size(); i++)
          {
            if(seqLengths.get(seqNames.get(i)) == artemisSeqLen)
            {
              // this looks like the correct sequence
              combo.setSelectedIndex(i);
              newIndex = i;
            }
          }

          if(!Options.isBlackBeltMode())
          {
            String label[] = {
                "The length of the sequence loaded does not match the length of",
                "the default reference sequence in the BAM ("+seqNames.get(0)+").",
                (newIndex == -1 ? "" : "The length does match the reference "+
                    seqNames.get(newIndex)+" so this has been set as the default.") 
            };
            new NonModalDialog(frame, label);
          }
        }
      }
    }

    mainPanel.add(bottomPanel, BorderLayout.SOUTH);
    coveragePanel.setPreferredSize(new Dimension(900, 100));
    coveragePanel.setVisible(false);
    snpPanel.setPreferredSize(new Dimension(900, 100));
    snpPanel.setVisible(false);

    mainPanel.revalidate();
    jspView.getVerticalScrollBar().setValue(
        jspView.getVerticalScrollBar().getMaximum());
  }
  
  private void addToViewMenu(final int thisBamIndex)
  {
    File f = new File(bamList.get(thisBamIndex));
    final JCheckBoxMenuItem cbBam = new JCheckBoxMenuItem(
                                     f.getName(), true);
    bamFilesMenu.add(cbBam);
    cbBam.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        if(cbBam.isSelected())
          hideBamList.remove(new Integer(thisBamIndex));
        else
          hideBamList.add(new Integer(thisBamIndex));
        laststart = -1;
        repaint();
      } 
    });
  }

  
  private void createMenus(JComponent menu)
  {
    final JMenuItem addBam = new JMenuItem("Add BAM ...");
    menu.add(addBam);
    addBam.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        FileSelectionDialog bamFileSelection = new FileSelectionDialog(
            null, false, "BamView", "BAM");
        List<String> bamFiles = bamFileSelection.getFiles(".*\\.bam$");
        int count = bamList.size();
       
        bamList.addAll(bamFileSelection.getFiles(".*\\.bam$"));
        
        for(int i=0; i<bamFiles.size(); i++)
          addToViewMenu(i+count);
        laststart = -1; 
        repaint();
      } 
    });
    
    bamFilesMenu.setFont(addBam.getFont());
    menu.add(bamFilesMenu);
    
    
    final JMenu analyse = new JMenu("Analyse");
    menu.add(analyse);
    final JMenuItem readCount = new JMenuItem("Read count of selected features ...");
    analyse.add(readCount);
    if(feature_display == null)
      readCount.setEnabled(false);
    readCount.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        FeatureVector features = feature_display.getSelection().getAllFeatures();

        JCheckBox overlap = new JCheckBox("Include all overlapping reads", true);
        overlap.setToolTipText("Include reads that partially overlap the feature");
        JCheckBox spliced = new JCheckBox("Introns included", true);
        Box yBox = Box.createVerticalBox();
        yBox.add(overlap);
        yBox.add(spliced);
        JOptionPane.showMessageDialog(null, yBox, "Read Count Option", JOptionPane.INFORMATION_MESSAGE);
        
        new MappedReads(features, (String)combo.getSelectedItem(), samFileReaderHash, bamList,
            seqNames, offsetLengths, concatSequences, seqLengths, 
            samRecordFlagPredicate, samRecordMapQPredicate,
            !overlap.isSelected(), spliced.isSelected());
      } 
    });
    
    final JMenuItem rpkm = new JMenuItem("RPKM value of selected features ...");
    analyse.add(rpkm);
    if(feature_display == null)
      rpkm.setEnabled(false);
    rpkm.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        FeatureVector features = feature_display.getSelection().getAllFeatures();

        JCheckBox overlap = new JCheckBox("Include all overlapping reads", true);
        overlap.setToolTipText("Include reads that partially overlap the feature");
        JCheckBox spliced = new JCheckBox("Introns included", true);
        JCheckBox allRefSeqs = new JCheckBox("Use reads mapped to all reference sequences", false);
        
        Box yBox = Box.createVerticalBox();
        yBox.add(overlap);
        yBox.add(spliced);
        
        if(seqLengths.size() > 1)
          yBox.add(allRefSeqs);
        
        JOptionPane.showMessageDialog(null, yBox, "Read Count Option", JOptionPane.INFORMATION_MESSAGE);
        
        int seqlen = 0;
        if(feature_display != null)
          seqlen = feature_display.getSequenceLength();
        else if(bases != null)
          seqlen = bases.getLength();
        
        new MappedReads(features, (String)combo.getSelectedItem(),
            samFileReaderHash, bamList, seqNames, offsetLengths, concatSequences, 
            seqLengths, seqlen, samRecordFlagPredicate, samRecordMapQPredicate,
            !overlap.isSelected(), spliced.isSelected(), allRefSeqs.isSelected());
      } 
    });

    for(int i=0; i<bamList.size(); i++)
      addToViewMenu(i);
    
    menu.add(new JSeparator());
    
    JMenu viewMenu = new JMenu("Views");
    cbStackView.setFont(viewMenu.getFont());
    cbIsizeStackView.setFont(viewMenu.getFont());
    cbPairedStackView.setFont(viewMenu.getFont());
    cbStrandStackView.setFont(viewMenu.getFont());
    cbCoverageView.setFont(viewMenu.getFont());
    
    baseQualityColour.setFont(viewMenu.getFont());
    colourByCoverageColour.setFont(viewMenu.getFont());
    markInsertions.setFont(viewMenu.getFont());
    
    cbIsizeStackView.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        laststart = -1;
        logMenuItem.setEnabled(isIsizeStackView());
        repaint();
      }
    });
    viewMenu.add(cbIsizeStackView);
    
    
    cbStackView.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        laststart = -1;
        logMenuItem.setEnabled(isIsizeStackView());
        repaint();
      }
    });
    viewMenu.add(cbStackView);
    

    cbPairedStackView.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        laststart = -1;
        logMenuItem.setEnabled(isIsizeStackView());
        repaint();
      }
    });
    viewMenu.add(cbPairedStackView);
    
    cbStrandStackView.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        laststart = -1;
        if(isStrandStackView())
        {
          setViewportMidPoint();
        }
        logMenuItem.setEnabled(isIsizeStackView());
        repaint();
      }
    });
    viewMenu.add(cbStrandStackView);
    
    cbCoverageView.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        laststart = -1;
        logMenuItem.setEnabled(isIsizeStackView());
        repaint();
      }
    });
    viewMenu.add(cbCoverageView);
    
    menu.add(viewMenu);
 
    final JCheckBoxMenuItem checkBoxSNPs = new JCheckBoxMenuItem("SNP marks");
    // 
    JMenu colourMenu = new JMenu("Colour By");
    colourMenu.add(colourByCoverageColour);
    
    baseQualityColour.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        if(baseQualityColour.isSelected())
        {
          checkBoxSNPs.setSelected(false);
          isSNPs = false;
        }
        repaint();
      }
    });
    colourMenu.add(baseQualityColour);
    menu.add(colourMenu);
    
    //
    JMenu showMenu = new JMenu("Show");
    JCheckBoxMenuItem checkBoxOrientation = new JCheckBoxMenuItem("Orientation");
    checkBoxOrientation.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        isOrientation = !isOrientation;
        repaint();
      }
    });
    showMenu.add(checkBoxOrientation);
    
    JCheckBoxMenuItem checkBoxSingle = new JCheckBoxMenuItem("Single Reads");
    checkBoxSingle.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        repaint();
        isSingle = !isSingle;
      }
    });
    showMenu.add(checkBoxSingle);
    
    checkBoxSNPs.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        if (isSNPs && bases == null)
        {
          JOptionPane.showMessageDialog(null,
              "No reference sequence supplied to identify SNPs.", "SNPs",
              JOptionPane.INFORMATION_MESSAGE);
        }
        isSNPs = !isSNPs;
        
        if(isSNPs)
          baseQualityColour.setSelected(false);
        repaint();
      }
    });
    showMenu.add(checkBoxSNPs);
    
    markInsertions.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        repaint();
      }
    });
    showMenu.add(markInsertions);
    menu.add(showMenu);
    
    //
    JMenu graphMenu = new JMenu("Graph");
    JCheckBoxMenuItem checkBoxCoverage = new JCheckBoxMenuItem("Coverage");
    checkBoxCoverage.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        isCoverage = !isCoverage;
        coveragePanel.setVisible(isCoverage);
        
        if(isCoverage && !cbCoverageView.isSelected())
          laststart = -1;
        repaint();
      }
    });
    graphMenu.add(checkBoxCoverage);
    
    JCheckBoxMenuItem checkBoxSNP = new JCheckBoxMenuItem("SNP");
    checkBoxSNP.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        isSNPplot = !isSNPplot;
        snpPanel.setVisible(isSNPplot);
        laststart = -1;
        repaint();
      }
    });
    graphMenu.add(checkBoxSNP);
    menu.add(graphMenu);
    
    
    if(feature_display != null)
    {
      final JCheckBoxMenuItem checkBoxSync =
        new JCheckBoxMenuItem("Asynchronous", asynchronous);
      checkBoxSync.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          asynchronous = checkBoxSync.isSelected();
        }
      });
      menu.add(checkBoxSync);
    }
    
    menu.add(new JSeparator());

    JMenu maxHeightMenu = new JMenu("BamView Height");
    menu.add(maxHeightMenu);
    
    final String hgts[] =
       {"500", "800", "1000", "1500", "2500", "5000", "50000"};
    
    ButtonGroup bgroup = new ButtonGroup();
    for(int i=0; i<hgts.length; i++)
    {
      final String hgt = hgts[i];
      final JCheckBoxMenuItem maxHeightMenuItem = new JCheckBoxMenuItem(hgt);
      bgroup.add(maxHeightMenuItem);
      maxHeightMenuItem.setSelected(hgts[i].equals(Integer.toString(maxHeight)));
      maxHeightMenu.add(maxHeightMenuItem);
      maxHeightMenuItem.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          if(maxHeightMenuItem.isSelected())
            maxHeight = Integer.parseInt(hgt);
          int start = getBaseAtStartOfView();
          setDisplay(start, nbasesInView+start, null);
        }
      });
    }
    
    menu.add(new JSeparator());
    logMenuItem.setFont(menu.getFont());
    menu.add(logMenuItem);
    logMenuItem.setEnabled(isIsizeStackView());
    
    logMenuItem.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        logScale = logMenuItem.isSelected();
        repaint();
      }
    });
    
    JMenuItem filter = new JMenuItem("Filter Reads ...");
    menu.add(filter);
    filter.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        if(filterFrame == null)
          filterFrame = new SAMRecordFilter(BamView.this);
        else
          filterFrame.setVisible(true);
      } 
    });
    
    JMenuItem readList = new JMenuItem("List Reads ...");
    menu.add(readList);
    readList.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        new SAMRecordList(BamView.this);
      }
    });

    final JMenuItem bamSplitter = new JMenuItem("Clone window");
    bamSplitter.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        BamView bamView = new BamView(new Vector<String>(bamList), 
            null, nbasesInView, entry_edit,
            feature_display, bases, (JPanel) mainPanel.getParent(), null);
        bamView.getJspView().getVerticalScrollBar().setValue(
            bamView.getJspView().getVerticalScrollBar().getMaximum());

        int start = getBaseAtStartOfView();
        setDisplay(start, nbasesInView+start, null);
        if(feature_display != null)
        {
          feature_display.addDisplayAdjustmentListener(bamView);
          feature_display.getSelection().addSelectionChangeListener(bamView);
        }
      } 
    });
    menu.add(new JSeparator());
    menu.add(bamSplitter);

    //
    JMenu coverageMenu = new JMenu("Coverage Options");
    coverageView.init(this, 0.f, 0, 0);
    coverageView.createMenus(coverageMenu);
    viewMenu.add(new JSeparator());
    viewMenu.add(coverageMenu);
  }
  
  private JComponent bamTopPanel(final JFrame frame)
  {
    final JComponent topPanel;
    if(frame == null)
    {
      topPanel = new JPanel(new FlowLayout(FlowLayout.LEADING, 0, 0));
      if(feature_display != null)
        this.selection = feature_display.getSelection(); 
    }
    else
    { 
      topPanel = new JMenuBar();
      frame.setJMenuBar((JMenuBar)topPanel);
      
      JMenu fileMenu = new JMenu("File");
      topPanel.add(fileMenu);

      JMenuItem readBam = new JMenuItem("Open new BamView ...");
      fileMenu.add(readBam);
      readBam.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          String[] s = { "NEW-BAMVIEW" };
          BamView.main(s);
        } 
      });
      
      JMenuItem close = new JMenuItem("Close");
      fileMenu.add(close);
      close.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          BamView.this.setVisible(false);
          Component comp = BamView.this;
          
          while( !(comp instanceof JFrame) )
            comp = comp.getParent();
          ((JFrame)comp).dispose();
        } 
      });
      
      JMenuItem exit = new JMenuItem("Exit");
      fileMenu.add(new JSeparator());
      fileMenu.add(exit);
      exit.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          int status = JOptionPane.showConfirmDialog(BamView.this, 
              "Exit BamView?", "Exit", 
              JOptionPane.OK_CANCEL_OPTION);
          if(status != JOptionPane.OK_OPTION)
            return;
          System.exit(0);
        } 
      });
      
      addKeyListener(new KeyAdapter()
      {
        public void keyPressed(final KeyEvent event)
        {
          switch (event.getKeyCode())
          {
            case KeyEvent.VK_UP:
              setZoomLevel((int) (BamView.this.nbasesInView * 1.1));
              break;
            case KeyEvent.VK_DOWN:
              if (showBaseAlignment)
                break;
              setZoomLevel((int) (BamView.this.nbasesInView * .9));
              break;
            default:
              break;
          }
        }
      });
    }
    
    if(seqNames.size() > 1)
    {
      int len = 0;
      for(int i=0; i<seqNames.size(); i++)
        len += seqLengths.get(seqNames.get(i));
      
      if(feature_display != null &&
         len == feature_display.getSequenceLength())
        concatSequences = true;
      else if(bases != null &&
          len == bases.getLength() )
        concatSequences = true;
    }

    // auto hide top panel
    final JCheckBox buttonAutoHide = new JCheckBox("Hide", (frame == null));
    buttonAutoHide.setToolTipText("Auto-Hide");
    final MouseMotionListener mouseMotionListener = new MouseMotionListener()
    {
      public void mouseDragged(MouseEvent event)
      {
        handleCanvasMouseDrag(event);
      }
      
      public void mouseMoved(MouseEvent e)
      {
        lastMousePoint = e.getPoint();
        
        int thisHgt = HEIGHT;
        if (thisHgt < 5)
          thisHgt = 15;

        int y = (int) (e.getY() - jspView.getViewport().getViewRect().getY());
        if (y < thisHgt)
        {
          topPanel.setVisible(true);
        }
        else
        {
          if (buttonAutoHide.isSelected())
            topPanel.setVisible(false);
        }
        mainPanel.repaint();
        mainPanel.revalidate();
      }
    };
    addMouseMotionListener(mouseMotionListener);

    combo = new JComboBox(seqNames);
    JTextComponent editor = (JTextComponent) combo.getEditor().getEditorComponent();
    editor.setDocument(new AutoCompleteComboDocument(combo));
    combo.setEditable(true);
    combo.setMaximumRowCount(20);
    
    combo.addItemListener(new ItemListener()
    {
      public void itemStateChanged(ItemEvent e)
      {
        laststart = -1;
        if(feature_display != null)
          setZoomLevel(feature_display.getMaxVisibleBases());
        else
          setZoomLevel(BamView.this.nbasesInView);
      }
    });
    topPanel.add(combo);

    if(feature_display == null)
    {
      final JTextField baseText = new JTextField(8);
      JButton goTo = new JButton("GoTo:");
      goTo.setToolTipText("Go to base position");
      goTo.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          try
          {
            int basePosition = Integer.parseInt(baseText.getText());
            scrollBar.setValue(basePosition);
          }
          catch (NumberFormatException nfe)
          {
            JOptionPane.showMessageDialog(BamView.this,
                "Expecting a base number!", "Number Format",
                JOptionPane.WARNING_MESSAGE);
          }
        }
      });
      topPanel.add(goTo);
      topPanel.add(baseText);

      JButton zoomIn = new JButton("-");
      zoomIn.setToolTipText("Zoom out (up arrow)");
      Insets ins = new Insets(1,1,1,1);
      zoomIn.setMargin(ins);
      zoomIn.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          setZoomLevel((int) (BamView.this.nbasesInView * 1.1));
        }
      });
      topPanel.add(zoomIn);

      JButton zoomOut = new JButton("+");
      zoomOut.setToolTipText("Zoom in (down arrow)");
      zoomOut.setMargin(ins);
      zoomOut.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          if (showBaseAlignment)
            return;
          setZoomLevel((int) (BamView.this.nbasesInView * .9));
        }
      });
      topPanel.add(zoomOut);
    }
    
    topPanel.add(buttonAutoHide);
    
    if(feature_display != null)
    {
      JButton close = new JButton("Close");
      topPanel.add(close);
      close.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          final JPanel containerPanel = (JPanel) mainPanel.getParent();
          feature_display.removeDisplayAdjustmentListener(BamView.this);
          feature_display.getSelection().removeSelectionChangeListener(BamView.this);
          containerPanel.remove(mainPanel);
          
          if(containerPanel.getComponentCount() > 0)
            containerPanel.revalidate();
          else
          {
            if(entry_edit != null)
              entry_edit.setNGDivider();
            else
              containerPanel.setVisible(false);
          }
        }
      });
    }
    return topPanel;
  }
  
  public void setVisible(boolean visible)
  {
    super.setVisible(visible);
    mainPanel.setVisible(visible);
  }
  
  private void setViewportMidPoint()
  {
    Point p = jspView.getViewport().getLocation();
    p.y = (getHeight() - jspView.getViewport().getViewRect().height)/2;
    jspView.getViewport().setViewPosition(p);
  }
  
  protected int getBaseAtStartOfView()
  {
    if(feature_display != null)
      return feature_display.getForwardBaseAtLeftEdge();
    else
      return scrollBar.getValue();
  }
  
  /**
   * Set the panel size based on the number of bases visible
   * and repaint.
   * @param nbasesInView
   */
  private void setZoomLevel(final int nbasesInView)
  {
    int startValue = getBaseAtStartOfView();
    this.nbasesInView = nbasesInView;
    float pixPerBase = getPixPerBaseByWidth(); 

    if(isBaseAlignmentView(pixPerBase))
    {
      pixPerBase = ALIGNMENT_PIX_PER_BASE;
      this.nbasesInView = (int)(mainPanel.getWidth()/pixPerBase);
      jspView.getVerticalScrollBar().setValue(0);
      
      if(ruler == null)
        ruler = new Ruler();
      jspView.setColumnHeaderView(ruler);
      showBaseAlignment = true;
      baseQualityColour.setEnabled(true);
      markInsertions.setEnabled(true);
    }
    else if(jspView != null)
    {
      if(!cbCoverageView.isSelected() && nbasesInView >= MAX_BASES)
      {
        cbLastSelected = getSelectedCheckBoxMenuItem();
        cbCoverageView.setSelected(true);
      }
      else if(cbCoverageView.isSelected() && nbasesInView < MAX_BASES && cbLastSelected != null)
      {
        cbLastSelected.setSelected(true);
        cbLastSelected = null;
      }
      
      jspView.setColumnHeaderView(null);
      
      if(!isStrandStackView())
        jspView.getVerticalScrollBar().setValue(
            jspView.getVerticalScrollBar().getMaximum());
      else
        setViewportMidPoint();
      showBaseAlignment = false;
      baseQualityColour.setEnabled(false);
      markInsertions.setEnabled(false);
    }
    
    if(scrollBar != null)
    {
      scrollBar.setValues(startValue, nbasesInView, 1, 
             getMaxBasesInPanel(getSequenceLength()));
      scrollBar.setUnitIncrement(nbasesInView/20);
      scrollBar.setBlockIncrement(nbasesInView);
    }
  }

  
  /**
   * Set the start and end base positions to display.
   * @param start
   * @param end
   * @param event
   */
  public void setDisplay(int start,
                         int end,
                         DisplayAdjustmentEvent event)
  {
    this.startBase = start;
    this.endBase   = end;
    this.nbasesInView = end-start+1;
    lastMousePoint = null;

    float pixPerBase;
    if(jspView.getViewport().getViewRect().width > 0)
      pixPerBase = getPixPerBaseByWidth();
    else
    {
      if(feature_display == null)
        pixPerBase = 1000.f/(float)(end-start+1);
      else
        pixPerBase = feature_display.getWidth()/(float)(end-start+1);
    }
    
    Dimension d = new Dimension();
    d.setSize(nbasesInView*pixPerBase, maxHeight);
    setPreferredSize(d);
    
    if(event == null)
    {
      this.startBase = -1;
      this.endBase   = -1;
    }
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
    
    //System.out.println(entryFileName);
    final uk.ac.sanger.artemis.io.Entry new_embl_entry =
      EntryFileDialog.getEntryFromFile(null, entry_document,
                                       artemis_entry_information,
                                       false);

    if(new_embl_entry == null)  // the read failed
      return null;

    Entry entry = null;
    try
    {
      if(entryGroup.getSequenceEntry() != null)
        bases = entryGroup.getSequenceEntry().getBases();
      
      if(bases == null)
      {
        entry = new Entry(new_embl_entry);
        bases = entry.getBases();
      }
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
  
  private boolean isShowScale()
  {
    return (feature_display == null ? true : false);
  }

  public JScrollPane getJspView()
  {
    return jspView;
  }
  
  /**
   *  Handle a mouse drag event on the drawing canvas.
   **/
  private void handleCanvasMouseDrag(final MouseEvent event)
  {
    if(event.getButton() == MouseEvent.BUTTON3 || bases == null) 
      return;
    
    highlightSAMRecord = null;
    if(event.getClickCount() > 1)
    {
      getSelection().clear();
      repaint();
      return;  
    }
    
    highlightRange(event, 
        MouseEvent.BUTTON1_DOWN_MASK & MouseEvent.BUTTON2_DOWN_MASK);
  }
  
  /**
   * 
   * @param event
   * @param onmask
   */
  private void highlightRange(MouseEvent event, int onmask)
  {
    int seqLength = getSequenceLength();
    float pixPerBase = getPixPerBaseByWidth();
    int start = (int) ( ( (event.getPoint().getX())/pixPerBase) + getBaseAtStartOfView() );
    
    if(start < 1)
      start = 1;
    if(start > seqLength)
      start = seqLength;
    
    if (dragStart < 0 && (event.getModifiersEx() & onmask) == onmask)
      dragStart = start;
    else if((event.getModifiersEx() & onmask) != onmask)
      dragStart = -1;
    
    MarkerRange drag_range;
    try
    {
      if(dragStart < 0)
        drag_range = new MarkerRange (bases.getForwardStrand(), start, start);
      else
        drag_range = new MarkerRange (bases.getForwardStrand(), dragStart, start);
      getSelection().setMarkerRange(drag_range);
      repaint();
    }
    catch (OutOfRangeException e)
    {
      e.printStackTrace();
    }
  }
  
  /**
   * Get the colour for the given read given to it by the coverage plot.
   * @param samRecord
   * @return
   */
  private Color getColourByCoverageColour(SAMRecord samRecord)
  {
    LineAttributes lines[] = CoveragePanel.getLineAttributes(bamList.size());
    int fileIndex = 0;
    if(bamList.size()>1)
      fileIndex = (Integer) samRecord.getAttribute("FL");
    return lines[fileIndex].getLineColour(); 
  }

  protected int getMaxBases()
  {
    return MAX_BASES;
  }
  
  protected void setMaxBases(int max)
  {
    MAX_BASES = max;
  }
  
  private boolean isStackView()
  {
    return cbStackView.isSelected();  
  }
  
  private boolean isPairedStackView()
  {
    return cbPairedStackView.isSelected();
  }
  
  private boolean isStrandStackView()
  {
    return cbStrandStackView.isSelected();
  }
  
  private boolean isCoverageView(float pixPerBase)
  {
    if(isBaseAlignmentView(pixPerBase))
      return false;
    return cbCoverageView.isSelected();
  }
  
  private boolean isIsizeStackView()
  {
    return cbIsizeStackView.isSelected();
  }
  
  private boolean isBaseAlignmentView(float pixPerBase)
  {
    if(pixPerBase*1.08f >= ALIGNMENT_PIX_PER_BASE)
      return true;
    return false;
  }
  
  private JCheckBoxMenuItem getSelectedCheckBoxMenuItem()
  {
    if(isStackView())
      return cbStackView;
    if(isPairedStackView())
      return cbPairedStackView;
    if(isStrandStackView())
      return cbStrandStackView;
    if(isIsizeStackView())
      return cbIsizeStackView;
    return cbCoverageView;
  }
  
  private Selection getSelection()
  {
    return selection;
  }
  
  protected List<SAMRecord> getReadsInView()
  {
    return readsInView;
  }
  
  protected int getBasesInView()
  {
    return nbasesInView;
  }
  
  protected void setHighlightSAMRecord(SAMRecord highlightSAMRecord)
  {
    this.highlightSAMRecord = highlightSAMRecord;
  }
  
  protected SAMRecord getHighlightSAMRecord()
  {
    return highlightSAMRecord;
  }
  
  protected FeatureDisplay getFeatureDisplay()
  {
    return feature_display;
  }
  
  private String getVersion()
  {
    final ClassLoader cl = this.getClass().getClassLoader();
    try
    {
      String line;
      InputStream in = cl.getResourceAsStream("etc/versions");
      BufferedReader reader = new BufferedReader(new InputStreamReader(in));
      while((line = reader.readLine()) != null)
      {
        if(line.startsWith("BamView"))
          return line.substring( "BamView".length() ).trim();
      }
      reader.close();
      in.close();
    }
    catch (Exception ex)
    {
    }
    return null;
  }
  
  /**
   * Artemis event notification
   */
  public void displayAdjustmentValueChanged(final DisplayAdjustmentEvent event)
  {
    if(event.getType() == DisplayAdjustmentEvent.REV_COMP_EVENT &&
       event.isRevCompDisplay())
      JOptionPane.showMessageDialog(this, 
          "Flipping the display is not supported by BamView.", "Warning", 
          JOptionPane.WARNING_MESSAGE);

    if(!asynchronous)
    {
      // if not asynchronous
      displayAdjustmentWork(event);
      return;
    }
    
    SwingWorker worker = new SwingWorker()
    {
      public Object construct()
      {
        try
        {
          Thread.sleep(500);
        }
        catch (InterruptedException e)
        {
          e.printStackTrace();
        }
        
        if(event.getStart() != ((FeatureDisplay)event.getSource()).getForwardBaseAtLeftEdge())
        {
          waitingFrame.showWaiting("waiting...", mainPanel);
          return null;
        }
      
        displayAdjustmentWork(event);
        waitingFrame.setVisible(false);
        return null;
      }
    };
    worker.start();
  }
  
  /**
   * Carry out the display agjustment event action.
   * @param event
   */
  private void displayAdjustmentWork(final DisplayAdjustmentEvent event)
  {
    if(event.getType() == DisplayAdjustmentEvent.SCALE_ADJUST_EVENT)
    {
      laststart = -1;

      BamView.this.startBase = event.getStart();
      BamView.this.endBase   = event.getEnd();

      int width = feature_display.getMaxVisibleBases();
      setZoomLevel(width);
      repaint();
    }
    else
    {
      setDisplay(event.getStart(), 
        event.getStart()+feature_display.getMaxVisibleBases(), event);
      repaint();
    }
  }
  
  public void selectionChanged(SelectionChangeEvent event)
  {
    repaint();
  }
  
  private class Ruler extends JPanel
  {
    private static final long serialVersionUID = 1L;
    int start;
    int end;

    public Ruler()
    {
      super();
      setPreferredSize(new Dimension(mainPanel.getWidth(), 15));
      setBackground(Color.white);
      setFont(getFont().deriveFont(11.f));
    }

    public void paintComponent(Graphics g)
    {
      super.paintComponent(g);
      Graphics2D g2 = (Graphics2D)g;
      drawBaseScale(g2, start, end, 12);
    }

    private void drawBaseScale(Graphics2D g2, int start, int end, int ypos)
    {
      int startMark = (((int)(start/10))*10)+1;

      if(end > getSequenceLength())
        end = getSequenceLength();
      
      for(int i=startMark; i<end; i+=10)
      {
        int xpos = (i-start)*ALIGNMENT_PIX_PER_BASE;
        g2.drawString(Integer.toString(i), xpos, ypos);
        
        xpos+=(ALIGNMENT_PIX_PER_BASE/2);
        g2.drawLine(xpos, ypos+1, xpos, ypos+5);
      }
    }
  }
  
  /**
  * Popup menu listener
  */
  class PopupListener extends MouseAdapter
  {
	JMenuItem gotoMateMenuItem;
	JMenuItem showDetails;
	
    public void mouseClicked(MouseEvent e)
    {
      if(e.isPopupTrigger() ||
         e.getButton() == MouseEvent.BUTTON3)
        return;
      
      BamView.this.requestFocus();
      
      if(e.getClickCount() > 1)
        getSelection().clear(); 
      else if(e.getButton() == MouseEvent.BUTTON1)
        highlightSAMRecord = mouseOverSAMRecord;
      else
        highlightRange(e, MouseEvent.BUTTON2_DOWN_MASK);
      repaint();
    }
    
    public void mousePressed(MouseEvent e)
    {
      maybeShowPopup(e);
    }

    public void mouseReleased(MouseEvent e)
    {
      dragStart = -1;
      maybeShowPopup(e);
    }

    private void maybeShowPopup(MouseEvent e)
    {
      if(e.isPopupTrigger())
      {
        if(popup == null)
        {
          popup = new JPopupMenu();
          createMenus(popup);
        }
        
        if(gotoMateMenuItem != null)
          popup.remove(gotoMateMenuItem);

        if(showDetails != null)
          popup.remove(showDetails);
        
        if( mouseOverSAMRecord != null && 
            mouseOverSAMRecord.getReadPairedFlag() &&
           !mouseOverSAMRecord.getMateUnmappedFlag() )
        {
          final SAMRecord thisSAMRecord = mouseOverSAMRecord;
          gotoMateMenuItem = new JMenuItem("Go to mate of : "+
              thisSAMRecord.getReadName());
          gotoMateMenuItem.addActionListener(new ActionListener()
          {
			public void actionPerformed(ActionEvent e) 
			{
			  String name = thisSAMRecord.getMateReferenceName();
			  if(name.equals("="))
			    name = thisSAMRecord.getReferenceName();
			  int offset = getSequenceOffset(name);
			  if(feature_display != null)
			    feature_display.makeBaseVisible(
			        thisSAMRecord.getMateAlignmentStart()+offset);
			  else
			    scrollBar.setValue(
			        thisSAMRecord.getMateAlignmentStart()+offset-
			        (nbasesInView/2));
			  
			  highlightSAMRecord = thisSAMRecord; 
			}  
          });
          popup.add(gotoMateMenuItem);
        }  
          
        if( mouseOverSAMRecord != null)
        {
          final SAMRecord thisSAMRecord = mouseOverSAMRecord;
          showDetails = new JMenuItem("Show details of : "+
              thisSAMRecord.getReadName());
          showDetails.addActionListener(new ActionListener()
          {
            public void actionPerformed(ActionEvent e) 
            {
              openFileViewer(thisSAMRecord, getMate(thisSAMRecord));
            }
          });
          popup.add(showDetails);
        }
        popup.show(e.getComponent(),
                e.getX(), e.getY());
      }
    }
  }
  
  protected static void openFileViewer(SAMRecord readRecord, SAMRecord mateRecord)
  {
    FileViewer viewDetail = new FileViewer(readRecord.getReadName(), true, false, false);
    appendToDetailView(readRecord, mateRecord, viewDetail);
  }
  
  private static void appendToDetailView(SAMRecord thisSAMRecord, SAMRecord thisSAMRecordMate, FileViewer viewDetail)
  {
    viewDetail.appendString("Read Name             "+thisSAMRecord.getReadName()+"\n", Level.INFO);
    viewDetail.appendString("Coordinates           "+thisSAMRecord.getAlignmentStart()+".."+
                                                     thisSAMRecord.getAlignmentEnd()+"\n", Level.DEBUG);
    viewDetail.appendString("Length                "+thisSAMRecord.getReadLength()+"\n", Level.DEBUG);
    viewDetail.appendString("Reference Name        "+thisSAMRecord.getReferenceName()+"\n", Level.DEBUG);
    viewDetail.appendString("Inferred Size         "+thisSAMRecord.getInferredInsertSize()+"\n", Level.DEBUG);
    viewDetail.appendString("Mapping Quality       "+thisSAMRecord.getMappingQuality()+"\n", Level.DEBUG);
    viewDetail.appendString("Cigar String          "+thisSAMRecord.getCigarString()+"\n", Level.DEBUG);
    viewDetail.appendString("Strand                "+
        (thisSAMRecord.getReadNegativeStrandFlag() ? "-\n\n" : "+\n\n"), Level.DEBUG);
    
    if(!thisSAMRecord.getMateUnmappedFlag())
    {     
      if(thisSAMRecordMate != null)
      {
        viewDetail.appendString("Mate Coordinates      "+thisSAMRecordMate.getAlignmentStart()+".."+
            thisSAMRecordMate.getAlignmentEnd()+"\n", Level.DEBUG);
        viewDetail.appendString("Mate Length           "+thisSAMRecordMate.getReadLength()+"\n", Level.DEBUG);
        viewDetail.appendString("Mate Reference Name   "+thisSAMRecordMate.getReferenceName()+"\n", Level.DEBUG);
        viewDetail.appendString("Mate Inferred Size    "+thisSAMRecordMate.getInferredInsertSize()+"\n", Level.DEBUG);
        viewDetail.appendString("Mate Mapping Quality  "+thisSAMRecordMate.getMappingQuality()+"\n", Level.DEBUG);
        viewDetail.appendString("Mate Cigar String     "+thisSAMRecordMate.getCigarString()+"\n", Level.DEBUG);
      }
      else
      {
        viewDetail.appendString("Mate Start Coordinate "+thisSAMRecord.getMateAlignmentStart()+"\n", Level.DEBUG);
        viewDetail.appendString("Mate Reference Name   "+thisSAMRecord.getMateReferenceName()+"\n", Level.DEBUG);
      }
      viewDetail.appendString("Mate Strand           "+
          (thisSAMRecord.getMateNegativeStrandFlag() ? "-" : "+"), Level.DEBUG);
    }
    else
    {
      viewDetail.appendString("Mate Unmapped ", Level.DEBUG);
    }
    
    viewDetail.appendString("\n\nFlags:", Level.INFO);
    viewDetail.appendString("\nDuplicate Read    "+
        (thisSAMRecord.getDuplicateReadFlag() ? "yes" : "no"), Level.DEBUG);
    
    viewDetail.appendString("\nRead Paired       "+
        (thisSAMRecord.getReadPairedFlag() ? "yes" : "no"), Level.DEBUG);
    if(thisSAMRecord.getReadPairedFlag())
    {
      viewDetail.appendString("\nFirst of Pair     "+
        (thisSAMRecord.getFirstOfPairFlag() ? "yes" : "no"), Level.DEBUG);
      viewDetail.appendString("\nMate Unmapped     "+
        (thisSAMRecord.getMateUnmappedFlag() ? "yes" : "no"), Level.DEBUG);  
      viewDetail.appendString("\nProper Pair       "+
        (thisSAMRecord.getProperPairFlag() ? "yes" : "no"), Level.DEBUG);
    }
    viewDetail.appendString("\nRead Fails Vendor\nQuality Check     "+
        (thisSAMRecord.getReadFailsVendorQualityCheckFlag() ? "yes" : "no"), Level.DEBUG);
    viewDetail.appendString("\nRead Unmapped     "+
        (thisSAMRecord.getReadUnmappedFlag() ? "yes" : "no"), Level.DEBUG);
    
    if(thisSAMRecord.getReadPairedFlag())
      viewDetail.appendString("\nSecond Of Pair    "+
        (thisSAMRecord.getSecondOfPairFlag() ? "yes" : "no"), Level.DEBUG);
    
    viewDetail.appendString("\n\nRead Bases:\n", Level.INFO);
    
    // wrap the read bases
    String seq = new String(thisSAMRecord.getReadBases());
    int MAX_SEQ_LINE_LENGTH = 100;
    for(int i=0; i<=seq.length(); i+=MAX_SEQ_LINE_LENGTH)
    {
      int iend = i+MAX_SEQ_LINE_LENGTH;
      if(iend > seq.length())
    	iend = seq.length();
      viewDetail.appendString(seq.substring(i, iend)+"\n", Level.DEBUG);
    }
  }

  /**
   * Query for the mate of a read
   * @param mate
   * @return
   */
  protected SAMRecord getMate(SAMRecord thisSAMRecord)
  {
    SAMRecord mate = null;
    try
    {
      int fileIndex = 0;
      if(bamList.size()>1)
        fileIndex = (Integer) thisSAMRecord.getAttribute("FL");
      String bam = bamList.get(fileIndex);  
      final SAMFileReader inputSam = getSAMFileReader(bam);
      mate = inputSam.queryMate(thisSAMRecord);
    }
    catch (Exception e)
    {
      e.printStackTrace();
    }
    return mate;
  }
  
  protected SAMRecordPredicate getSamRecordFlagPredicate()
  {
    return samRecordFlagPredicate;
  }

  protected void setSamRecordFlagPredicate(
      SAMRecordPredicate samRecordFlagPredicate)
  {
    laststart = -1;
    lastend = -1;
    this.samRecordFlagPredicate = samRecordFlagPredicate;
  }
  
  protected SAMRecordMapQPredicate getSamRecordMapQPredicate()
  {
    return samRecordMapQPredicate;
  }

  protected void setSamRecordMapQPredicate(
      SAMRecordMapQPredicate samRecordMapQPredicate)
  {
    laststart = -1;
    lastend = -1;
    this.samRecordMapQPredicate = samRecordMapQPredicate;
  }
  
  class PairedRead
  {
    SAMRecord sam1;
    SAMRecord sam2;
  } 
 
  public static void main(String[] args)
  {
    BamFrame frame = new BamFrame();
    if(args.length == 0 && BamFrame.isMac())
    {
      try
      {
        Thread.sleep(1000);
      }
      catch (InterruptedException e1) {}
      if(frame.getBamFile() != null)
        args = new String[]{ frame.getBamFile() };
    }
      
    List<String> bam = new Vector<String>();
    String reference = null;
    if(args.length == 0 || args[0].equals("NEW-BAMVIEW"))
    {
      System.setProperty("default_directory", System.getProperty("user.dir"));
      FileSelectionDialog fileSelection = new FileSelectionDialog(
          null, true, "BamView", "BAM");
      bam = fileSelection.getFiles(".*\\.bam$");
      reference = fileSelection.getReferenceFile();
      if(reference == null || reference.equals(""))
        reference = null;
      
      if(bam == null || bam.size() < 1)
      {
        if(args.length > 0 && args[0].equals("NEW-BAMVIEW"))
          return;
        System.err.println("No files found.");
        System.exit(0);
      }
    }
    else if(!args[0].startsWith("-"))
    {
      for(int i=0; i< args.length; i++)
        bam.add(args[i]);
    }
    int nbasesInView = 2000;
    String chr = null;
    String vw  = null;
    boolean orientation = false;
    boolean covPlot     = false;
    boolean snpPlot     = false;
    int base = 0;
    
    for(int i=0;i<args.length; i++)
    {
      if(args[i].equals("-a"))
      {
        while(i < args.length-1 && !args[++i].startsWith("-"))
        {
          String filename = args[i];
          if(FileSelectionDialog.isListOfFiles(filename))
            bam.addAll(FileSelectionDialog.getListOfFiles(filename));
          else
            bam.add(filename);
        }
        --i;
      }
      else if(args[i].equals("-r"))
        reference = args[++i];
      else if(args[i].equals("-n"))
        nbasesInView = Integer.parseInt(args[++i]);
      else if(args[i].equals("-s"))
        System.setProperty("samtoolDir", args[++i]);
      else if(args[i].equals("-c"))
        chr = args[++i].trim();
      else if(args[i].equals("-b"))
        base = Integer.parseInt(args[++i].trim());
      else if(args[i].equals("-v"))
        vw = args[++i].trim();
      else if(args[i].equals("-o"))
        orientation = true;
      else if(args[i].equals("-pc"))
        covPlot = true;
      else if(args[i].equals("-ps"))
        snpPlot = true;
      else if(args[i].startsWith("-h"))
      { 
        System.out.println("-h\t show help");
        
        System.out.println("-a\t BAM/SAM file to display");
        System.out.println("-r\t reference file (optional)");
        System.out.println("-n\t number of bases to display in the view (optional)");
        System.out.println("-c\t chromosome name (optional)");
        System.out.println("-v\t view (optional - IS (inferred size), S (stack, default), PS (paired stack), ST (strand), C (coverage))");
        System.out.println("-b\t base position (optional)");
        System.out.println("-o\t show orientation (optional)");
        System.out.println("-pc\t plot coverage (optional)");
        System.out.println("-ps\t plot SNP (optional and only with -r)");
        System.exit(0);
      }
    }

    final BamView view = new BamView(bam, reference, nbasesInView, null, null,
        (JPanel)frame.getContentPane(), frame);
    frame.setTitle("BamView v"+view.getVersion());
    
    if(chr != null)
      view.combo.setSelectedItem(chr);
    if(vw != null)
    {
      if(vw.equalsIgnoreCase("IS"))
        view.cbIsizeStackView.setSelected(true);
      if(vw.equalsIgnoreCase("PS"))
        view.cbPairedStackView.setSelected(true);
      if(vw.equalsIgnoreCase("ST"))
        view.cbStrandStackView.setSelected(true);
      if(vw.equalsIgnoreCase("C"))
        view.cbCoverageView.setSelected(true);
    }
    if(base > 0)
      view.scrollBar.setValue(base);
    if(orientation)
      view.isOrientation = true;
    if(covPlot)
    {
      view.isCoverage = true;
      view.coveragePanel.setVisible(true);
    }
    if(snpPlot)
    {
      view.isSNPplot = true;
      view.snpPanel.setVisible(true);
    }

    // translucent
    //frame.getRootPane().putClientProperty("Window.alpha", new Float(0.9f));
    /*frame.addWindowFocusListener(new WindowFocusListener()
    {
      public void windowGainedFocus(WindowEvent e)
      {
        view.requestFocus();
      }
      public void windowLostFocus(WindowEvent e){}
    });*/

    frame.pack();

    view.jspView.getVerticalScrollBar().setValue(
        view.jspView.getVerticalScrollBar().getMaximum());
    frame.setVisible(true);
  }
}
