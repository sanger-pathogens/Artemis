/* JamView
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

import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.Stroke;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionListener;
import java.awt.event.WindowEvent;
import java.awt.event.WindowFocusListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.StringReader;
import java.util.Collections;
import java.util.Comparator;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.List;
import java.util.Vector;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.Scrollable;
import javax.swing.SwingConstants;
import javax.swing.UIManager;
import javax.swing.border.Border;
import javax.swing.border.EmptyBorder;

import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.util.CloseableIterator;

import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.Selection;
import uk.ac.sanger.artemis.SelectionChangeEvent;
import uk.ac.sanger.artemis.SelectionChangeListener;
import uk.ac.sanger.artemis.SimpleEntryGroup;
import uk.ac.sanger.artemis.components.DisplayAdjustmentEvent;
import uk.ac.sanger.artemis.components.DisplayAdjustmentListener;
import uk.ac.sanger.artemis.components.EntryFileDialog;
import uk.ac.sanger.artemis.components.FeatureDisplay;
import uk.ac.sanger.artemis.components.MessageDialog;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.sequence.Bases;
import uk.ac.sanger.artemis.sequence.MarkerRange;
import uk.ac.sanger.artemis.sequence.NoSequenceException;
import uk.ac.sanger.artemis.util.Document;
import uk.ac.sanger.artemis.util.DocumentFactory;
import uk.ac.sanger.artemis.util.OutOfRangeException;

public class JamView extends JPanel
                     implements Scrollable, DisplayAdjustmentListener, SelectionChangeListener
{
  private static final long serialVersionUID = 1L;
  private List<SAMRecord> readsInView;
  private Hashtable<String, Integer> seqLengths = new Hashtable<String, Integer>();
  private Vector<String> seqNames = new Vector<String>();
  private String bam;

  private Bases bases;
  private JScrollPane jspView;
  private JComboBox combo;
  private JCheckBox checkBoxSingle;
  private JCheckBox checkBoxSNPs;
  
  private FeatureDisplay feature_display;
  private Selection selection;
  private JPanel mainPanel;
  private boolean showScale = true;
  private Ruler ruler = new Ruler();
  private int nbasesInView;
  
  private int startBase = -1;
  private int endBase   = -1;
  private int laststart;
  private int lastend;
  private int maxUnitIncrement = 8;
  private Cursor cbusy = new Cursor(Cursor.WAIT_CURSOR);
  private Cursor cdone = new Cursor(Cursor.DEFAULT_CURSOR);
  
  private boolean showBaseAlignment = false;
  /** Used to colour the frames. */
  private Color light_grey = new Color(200, 200, 200);
  private int ALIGNMENT_PIX_PER_BASE;
  
  private JPopupMenu popup;
  // use the picard library otherwise call samtools
  public static boolean PICARD = true;

 
  public JamView(String bam, 
                 String reference,
                 int nbasesInView)
  {
    super();
    setBackground(Color.white);
    this.bam = bam;
    this.nbasesInView = nbasesInView;
    
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
    
    if(PICARD)
      readHeaderPicard();
    else
      readHeader();
    
    // set font size
    //setFont(getFont().deriveFont(12.f));

    //Options.getOptions().getFontUIResource();
    //final javax.swing.plaf.FontUIResource font_ui_resource =
    //  new javax.swing.plaf.FontUIResource(getFont());
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
    selection = new Selection(null);
  }
  
  /**
   * Read the BAM/SAM header
   */
  private void readHeader()
  {
    String samtoolCmd = "";
    if(System.getProperty("samtoolDir") != null)
      samtoolCmd = System.getProperty("samtoolDir");
    String cmd[] = { samtoolCmd+File.separator+"samtools",  
			     "view", "-H", bam };
	
    RunSamTools samtools = new RunSamTools(cmd, null, null, null);
	
    if(samtools.getProcessStderr() != null)
      System.out.println(samtools.getProcessStderr());
 
    String header = samtools.getProcessStdout();
    
    StringReader samReader = new StringReader(header);
	BufferedReader buff = new BufferedReader(samReader);
    
	String line;
	try 
	{
	  while((line = buff.readLine()) != null)
	  {
	    if(line.indexOf("LN:") > -1)
	    {
	      String parts[] = line.split("\t");
	      String name = "";
	      int seqLength = 0;
	      for(int i=0; i<parts.length; i++)
	      {
	        if(parts[i].startsWith("LN:"))
	          seqLength = Integer.parseInt( parts[i].substring(3) );
	        else if(parts[i].startsWith("SN:"))
	          name = parts[i].substring(3);
	      }
	      seqLengths.put(name, seqLength);
	      seqNames.add(name);
	    }
	  }
	} 
	catch (IOException e) 
	{
	  e.printStackTrace();
	}
  }

  private void readHeaderPicard()
  {
    File bamFile = new File(bam);
    File indexFile = new File(bam+".bai");
    final SAMFileReader inputSam = new SAMFileReader(bamFile, indexFile);
    SAMFileHeader header = inputSam.getFileHeader();
    List<SAMSequenceRecord> readGroups = header.getSequenceDictionary().getSequences();
    
    for(int i=0; i<readGroups.size(); i++)
    {
      seqLengths.put(readGroups.get(i).getSequenceName(),
                     readGroups.get(i).getSequenceLength());
      seqNames.add(readGroups.get(i).getSequenceName());
    }
    inputSam.close();
  }

  /**
   * Read data from BAM/SAM file for a region.
   * @param start
   * @param end
   * @param pair_sort
   */
  private void readFromBam(int start, int end)
  {
    String refName = (String) combo.getSelectedItem();
    
    String samtoolCmd = "";
    if(System.getProperty("samtoolDir") != null)
      samtoolCmd = System.getProperty("samtoolDir");
    
	String cmd[] = { samtoolCmd+File.separator+"samtools",  
				     "view", 
				     bam, refName+":"+start+"-"+end };
		
	for(int i=0; i<cmd.length;i++)
	  System.out.print(cmd[i]+" ");
	System.out.println();
	
    if(readsInView == null)
      readsInView = new Vector<SAMRecord>();
    else
      readsInView.clear();
	RunSamTools samtools = new RunSamTools(cmd, null, null, readsInView);
		
	if(samtools.getProcessStderr() != null)
      System.out.println(samtools.getProcessStderr());
	    
	samtools.waitForStdout();
  }
  
  /**
   * Read a SAM or BAM file.
   */
  private void readFromBamPicard(int start, int end)
  {
    // Open the input file.  Automatically detects whether input is SAM or BAM
    // and delegates to a reader implementation for the appropriate format.
    File bamFile = new File(bam);
    File indexFile = new File(bam+".bai");
    final SAMFileReader inputSam = new SAMFileReader(bamFile, indexFile);
    inputSam.setValidationStringency(ValidationStringency.SILENT);
    
    if(readsInView == null)
      readsInView = new Vector<SAMRecord>();
    else
      readsInView.clear();
    String refName = (String) combo.getSelectedItem();
    CloseableIterator<SAMRecord> it = inputSam.queryOverlapping(refName, start, end);
    while ( it.hasNext() )
    {
      try
      {
        SAMRecord samRecord = it.next();
        readsInView.add(samRecord);
      }
      catch(Exception e)
      {
        System.out.println(e.getMessage());
      }
    }
    inputSam.close();
  }
  
  /**
   * Override
   */
  protected void paintComponent(Graphics g)
  {
	super.paintComponent(g);
	Graphics2D g2 = (Graphics2D)g;

	String refName = (String) combo.getSelectedItem();
    int seqLength = seqLengths.get(refName);
	float pixPerBase = ((float)getWidth())/(float)(seqLength);
	
    int start;
    final int end;
    
    if(startBase > 0)
      start = startBase;
    else
    {
      double x = jspView.getViewport().getViewRect().getX();
      start = (int) (seqLength * ( (float)x / (float)getWidth()));
    }
    
    if(endBase > 0)
      end = endBase;
    else
      end   = (int) (start + ((float)jspView.getViewport().getWidth() / 
                              (float)pixPerBase));

    //System.out.println("paintComponent "+start+".."+end);
    
    if(laststart != start ||
       lastend   != end)
    {
      try
      {
        setCursor(cbusy);

        if(PICARD)
          readFromBamPicard(start, end);
        else
          readFromBam(start, end);

        Collections.sort(readsInView, new SAMRecordComparator());
        setCursor(cdone);
      }
      catch(OutOfMemoryError ome)
      {
        JOptionPane.showMessageDialog(this, "Out of Memory");
        return;
      }
    }
    
    laststart = start;
    lastend   = end;
	if(pixPerBase*2 >= ALIGNMENT_PIX_PER_BASE)
	  drawBaseAlignment(g2, seqLength, pixPerBase, start, end);
	else
	  drawLineView(g2, seqLength, pixPerBase, start, end);
  }
  
  /**
   * Draw the zoomed-in base view.
   * @param g2
   * @param seqLength
   * @param pixPerBase
   * @param start
   * @param end
   */
  private void drawBaseAlignment(Graphics2D g2, int seqLength, 
                                 float pixPerBase, final int start, final int end)
  {
    ruler.start = start;
    ruler.end = end;
    ruler.repaint();
    
    int ypos = 0;
    
    
    String refSeq = null;
    int refSeqStart = start;
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
        int xpos = (refSeqStart-1)*ALIGNMENT_PIX_PER_BASE;
        
        g2.setColor(light_grey);
        g2.fillRect(xpos, ypos-11, 
            jspView.getViewport().getWidth()+(ALIGNMENT_PIX_PER_BASE*2), 11);
        drawSelectionRange(g2, pixPerBase, start, end);
        g2.setColor(Color.black);
        g2.drawString(refSeq, xpos, ypos);
        
        //for(int i=0;i<refSeq.length(); i++)
        //{
          //xpos = ((refSeqStart-1) + i)*ALIGNMENT_PIX_PER_BASE;
          //g2.drawString(refSeq.substring(i, i+1), xpos, ypos);
        //}
      }
      catch (OutOfRangeException e)
      {
        e.printStackTrace();
      }
    }
    else
      drawSelectionRange(g2, pixPerBase, start, end);

    
    boolean drawn[] = new boolean[readsInView.size()];
    for(int i=0; i<readsInView.size(); i++)
      drawn[i] = false;
    
    for(int i=0; i<readsInView.size(); i++)
    {
      if (!drawn[i])
      {
        SAMRecord thisRead = readsInView.get(i);
        ypos+=11;

        drawSequence(g2, thisRead, pixPerBase, ypos, refSeq, refSeqStart);
        drawn[i] = true;
        
        int thisEnd = thisRead.getAlignmentEnd();
        if(thisEnd == 0)
          thisEnd = thisRead.getAlignmentStart()+thisRead.getReadLength();
        
        for(int j=i+1; j<readsInView.size(); j++)
        {
          if (!drawn[j])
          {
            SAMRecord nextRead = readsInView.get(j);
            if(nextRead.getAlignmentStart() > thisEnd+1)
            {
              drawSequence(g2, nextRead, pixPerBase, ypos, refSeq, refSeqStart);
              drawn[j] = true;
              thisEnd = nextRead.getAlignmentEnd();
              if(thisEnd == 0)
                thisEnd = nextRead.getAlignmentStart()+nextRead.getReadLength();
            }
          }
        }
      }
    }
    
    if(ypos > getHeight())
    {
      setPreferredSize(new Dimension(getWidth(), ypos));
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
                            float pixPerBase, int ypos, String refSeq, int refSeqStart)
  {
    if (!samRecord.getReadPairedFlag() ||  // read is not paired in sequencing
        samRecord.getMateUnmappedFlag() )  // mate is unmapped )  // mate is unmapped 
      g2.setColor(Color.black);
    else
      g2.setColor(Color.blue);
    
    Color col = g2.getColor();
    int xpos;
    String readSeq = samRecord.getReadString();

    List<AlignmentBlock> blocks = samRecord.getAlignmentBlocks();
    for(int i=0; i<blocks.size(); i++)
    {
      AlignmentBlock block = blocks.get(i);
      for(int j=0; j<block.getLength(); j++)
      {
        int readPos = block.getReadStart()-1+j;
        xpos  = block.getReferenceStart()-1+j;
        int refPos = xpos-refSeqStart+1;
        
        if(checkBoxSNPs.isSelected() && refSeq != null && refPos > 0 && refPos < refSeq.length())
        { 
          if(readSeq.charAt(readPos) != refSeq.charAt(refPos))
            g2.setColor(Color.red);
          else
            g2.setColor(col);
        }
        g2.drawString(readSeq.substring(readPos, readPos+1), xpos*ALIGNMENT_PIX_PER_BASE, ypos);
      }
    }
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
      drawScale(g2, start, end, pixPerBase);
    
    Stroke originalStroke = new BasicStroke (1.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_ROUND); 
    Stroke stroke =
            new BasicStroke (1.3f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_ROUND);
    int scaleHeight;
    if(isShowScale())
      scaleHeight = 15;
    else
      scaleHeight = 0;
    
    for(int i=0; i<readsInView.size(); i++)
    {
      SAMRecord samRecord = readsInView.get(i);
      SAMRecord samNextRecord = null;      

      if( !samRecord.getReadPairedFlag() ||  // read is not paired in sequencing
          samRecord.getMateUnmappedFlag() )  // mate is unmapped
      {
        if(checkBoxSingle.isSelected())
        {
          int ypos = (getHeight() - scaleHeight) - samRecord.getReadString().length();
          g2.setColor(Color.orange);
          drawRead(g2, samRecord, pixPerBase, stroke, ypos);
        }
        continue;
      }

      int ypos = (getHeight() - scaleHeight) - ( Math.abs(samRecord.getInferredInsertSize()) );
      if(i < readsInView.size()-1)
      {
        samNextRecord = readsInView.get(++i);

        if(samRecord.getReadName().equals(samNextRecord.getReadName()))
        { 
          // draw connection between paired reads
          if(samRecord.getAlignmentEnd() < samNextRecord.getAlignmentStart() && 
              (samNextRecord.getAlignmentStart()-samRecord.getAlignmentEnd())*pixPerBase > 2.f)
          {
            g2.setStroke(originalStroke);
            g2.setColor(Color.LIGHT_GRAY);
            g2.drawLine((int)(samRecord.getAlignmentEnd()*pixPerBase), ypos, 
                        (int)(samNextRecord.getAlignmentStart()*pixPerBase), ypos);
          }
          
          if( samRecord.getReadNegativeStrandFlag() && // strand of the query (1 for reverse)
              samNextRecord.getReadNegativeStrandFlag() )
            g2.setColor(Color.red);
          else
            g2.setColor(Color.blue);
          
          drawRead(g2, samRecord, pixPerBase, stroke, ypos);
          drawRead(g2, samNextRecord, pixPerBase, stroke, ypos);
        }
        else
        {
          drawLoneRead(g2, samRecord, ypos, pixPerBase, originalStroke, stroke);
          i--;
        }
      }
      else
      {
        drawLoneRead(g2, samRecord, ypos, pixPerBase, originalStroke, stroke);
      }
    }
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
      float pixPerBase, Stroke originalStroke, Stroke stroke)
  {
    boolean drawLine = true;   
    int thisStart = samRecord.getAlignmentStart()-1;
    int thisEnd   = thisStart + samRecord.getReadString().length();
    
    if(drawLine &&
        (samRecord.getMateAlignmentStart()-samRecord.getAlignmentEnd())*pixPerBase > 2.f)
    {
      g2.setStroke(originalStroke);
      g2.setColor(Color.LIGHT_GRAY);
      
      if(samRecord.getAlignmentEnd() < samRecord.getMateAlignmentStart())
      {
        int nextStart = (int) ((samRecord.getMateAlignmentStart()-1)*pixPerBase);
        g2.drawLine((int)(thisStart*pixPerBase), ypos, nextStart, ypos);
      }
      else
      {
        int nextStart = (int) ((samRecord.getMateAlignmentStart()-1)*pixPerBase);
        g2.drawLine((int)(thisEnd*pixPerBase), ypos, nextStart, ypos);
      }
    }
    
    if(samRecord.getReadNegativeStrandFlag()) // strand of the query (1 for reverse)
      g2.setColor(Color.red);
    else
      g2.setColor(Color.blue);
    
    if(ypos <= 0)
    {
      ypos = samRecord.getReadString().length();
      drawLine = false;
      g2.setColor(Color.orange); 
    }  
    drawRead(g2, samRecord, pixPerBase, stroke, ypos);
  }

  
  private void drawScale(Graphics2D g2, int start, int end, float pixPerBase)
  {
    g2.setColor(Color.black);
    g2.drawLine( (int)(start*pixPerBase), getHeight()-14,
                 (int)(end*pixPerBase),   getHeight()-14);
    int interval = end-start;
    
    if(interval > 256000)
      drawTicks(g2, start, end, pixPerBase, 512000);
    else if(interval > 64000)
      drawTicks(g2, start, end, pixPerBase, 12800);
    else if(interval > 16000)
      drawTicks(g2, start, end, pixPerBase, 3200);
    else if(interval > 4000)
      drawTicks(g2, start, end, pixPerBase, 800);
    else if(interval > 1000)
      drawTicks(g2, start, end, pixPerBase, 400);
    else
      drawTicks(g2, start, end, pixPerBase, 100);
  }
  
  private void drawTicks(Graphics2D g2, int start, int end, float pixPerBase, int division)
  {
    int markStart = (Math.round(start/division)*division);
    
    if(markStart < 1)
      markStart = 1;
    
    int sm = markStart-(division/2);
    
    if(sm > start)
      g2.drawLine((int)(sm*pixPerBase), getHeight()-14,(int)(sm*pixPerBase), getHeight()-12);
    
    for(int m=markStart; m<end; m+=division)
    {
      g2.drawString(Integer.toString(m), m*pixPerBase, getHeight()-1);
      g2.drawLine((int)(m*pixPerBase), getHeight()-14,(int)(m*pixPerBase), getHeight()-11);
      
      sm = m+(division/2);
      
      if(sm < end)
        g2.drawLine((int)(sm*pixPerBase), getHeight()-14,(int)(sm*pixPerBase), getHeight()-12);
      
      if(m == 1)
        m = 0;
    }
  }
  
  private void drawRead(Graphics2D g2, SAMRecord thisRead,
		                float pixPerBase, Stroke stroke, int ypos)
  {
    int thisStart = thisRead.getAlignmentStart()-1;
    int thisEnd   = thisRead.getAlignmentEnd();
    g2.setStroke(stroke);
    g2.drawLine((int) (thisStart * pixPerBase), ypos,
                (int) (thisEnd * pixPerBase), ypos);

    if (checkBoxSNPs.isSelected())
      showSNPsOnReads(g2, thisRead, pixPerBase, ypos);
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
        
        int x = (int) (pixPerBase*(rangeStart-1));
        int width = (int) (pixPerBase*(rangeEnd-rangeStart+1));
        
        g2.setColor(Color.pink);
        g2.fillRect(x, 0, width, getHeight());
      }
    }
  }
  
  /**
   * Display the SNPs for the given read.
   * @param g2
   * @param thisRead
   * @param pixPerBase
   * @param ypos
   */
  private void showSNPsOnReads(Graphics2D g2, SAMRecord thisRead,
                               float pixPerBase, int ypos)
  {
    int thisStart = thisRead.getAlignmentStart();
    int thisEnd   = thisRead.getAlignmentEnd();
    
    // use alignment blocks of the contiguous alignment of
    // subsets of read bases to a reference sequence
    List<AlignmentBlock> blocks = thisRead.getAlignmentBlocks();
    try
    {
      char[] refSeq = bases.getSubSequenceC(
          new Range(thisStart, thisEnd), Bases.FORWARD);
      byte[] readSeq = thisRead.getReadBases();

      Color col = g2.getColor();
      g2.setColor(Color.red);

      for(int i=0; i<blocks.size(); i++)
      {
        AlignmentBlock block = blocks.get(i);
        for(int j=0; j<block.getLength(); j++)
        {
          int readPos = block.getReadStart()-1+j;
          int refPos  = block.getReferenceStart()+j;

          if (Character.toUpperCase(refSeq[refPos-thisStart]) != readSeq[readPos])
          {
            g2.drawLine((int) ((thisStart + i) * pixPerBase), ypos + 2,
                        (int) ((thisStart + i) * pixPerBase), ypos - 2);
          }
        }
        
      }

      g2.setColor(col);
    }
    catch (OutOfRangeException e)
    {
      e.printStackTrace();
    }
  }
  
  /**
   * Add the alignment view to the supplied <code>JPanel</code> in
   * a <code>JScrollPane</code>.
   * @param mainPanel  panel to add the alignment to
   * @param autohide automatically hide the top panel containing the buttons
   */
  public void addJamToPanel(final JPanel mainPanel,
                            final boolean autohide,
                            final FeatureDisplay feature_display)
  {
    this.mainPanel = mainPanel;
    
    if(feature_display != null)
    {
      this.feature_display = feature_display;
      this.selection = feature_display.getSelection();
    }
    
    final JPanel topPanel = new JPanel(new GridBagLayout());
    GridBagConstraints gc = new GridBagConstraints();

    // auto hide top panel
    final JCheckBox buttonAutoHide = new JCheckBox("Auto-hide", autohide);
    final MouseMotionListener mouseMotionListener = new MouseMotionListener()
    {
      public void mouseDragged(MouseEvent event)
      {
        handleCanvasMouseDragOrClick(event);
      }
      
      public void mouseMoved(MouseEvent e)
      {
        int thisHgt = HEIGHT;
        if (thisHgt < 5)
          thisHgt = 15;

        int y = (int) (e.getY() - jspView.getViewport().getViewRect().getY());
        if (y < thisHgt)
        {
          if (!containsComponent(topPanel, mainPanel))
            mainPanel.add(topPanel, BorderLayout.NORTH);
        }
        else
        {
          if (buttonAutoHide.isSelected() && containsComponent(topPanel, mainPanel))
            mainPanel.remove(topPanel);
        }
        mainPanel.repaint();
        mainPanel.revalidate();
      }
    };
    addMouseMotionListener(mouseMotionListener);

    combo = new JComboBox(seqNames);
    combo.setEditable(false);

    combo.addItemListener(new ItemListener()
    {
      public void itemStateChanged(ItemEvent e)
      {
        laststart = -1;
        lastend = -1;
        setZoomLevel(JamView.this.nbasesInView);
      }
    });
    gc.fill = GridBagConstraints.NONE;
    gc.anchor = GridBagConstraints.NORTHWEST;
    topPanel.add(combo, gc);

    checkBoxSingle = new JCheckBox("Single Reads");
    checkBoxSingle.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        repaint();
      }
    });
    topPanel.add(checkBoxSingle, gc);

    checkBoxSNPs = new JCheckBox("SNPs");
    checkBoxSNPs.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        if (checkBoxSNPs.isSelected() && bases == null)
        {
          JOptionPane.showMessageDialog(null,
              "No reference sequence supplied to identify SNPs.", "SNPs",
              JOptionPane.INFORMATION_MESSAGE);
        }
        repaint();
      }
    });
    topPanel.add(checkBoxSNPs, gc);

    if (feature_display == null)
    {
      final JTextField baseText = new JTextField(10);
      JButton goTo = new JButton("GoTo:");
      goTo.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          try
          {
            int basePosition = Integer.parseInt(baseText.getText());
            goToBasePosition(basePosition);
          }
          catch (NumberFormatException nfe)
          {
            JOptionPane.showMessageDialog(JamView.this,
                "Expecting a base number!", "Number Format",
                JOptionPane.WARNING_MESSAGE);
          }
        }
      });
      topPanel.add(goTo, gc);
      topPanel.add(baseText, gc);

      JButton zoomIn = new JButton("+");
      Insets ins = new Insets(0,0,0,0);
      zoomIn.setMargin(ins);
      zoomIn.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          int startBase = getBaseAtStartOfView();
          setZoomLevel((int) (JamView.this.nbasesInView * 1.1));
          goToBasePosition(startBase);
          repaint();
        }
      });
      topPanel.add(zoomIn, gc);

      JButton zoomOut = new JButton("-");
      zoomOut.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          if (showBaseAlignment)
            return;
          int startBase = getBaseAtStartOfView();
          setZoomLevel((int) (JamView.this.nbasesInView * .9));
          goToBasePosition(startBase);
          repaint();
        }
      });
      topPanel.add(zoomOut, gc);
    }
    
    topPanel.add(buttonAutoHide, gc);

    mainPanel.setPreferredSize(new Dimension(800, 400));
    jspView = new JScrollPane(this, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
        JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);

    setDisplay(1, nbasesInView, null);
    mainPanel.setLayout(new BorderLayout());

    mainPanel.add(topPanel, BorderLayout.NORTH);
    mainPanel.add(jspView, BorderLayout.CENTER);

    jspView.getVerticalScrollBar().setValue(
        jspView.getVerticalScrollBar().getMaximum());

    if (feature_display == null)
    {
      addKeyListener(new KeyAdapter()
      {
        public void keyPressed(final KeyEvent event)
        {
          switch (event.getKeyCode())
          {
          case KeyEvent.VK_UP:
            int startBase = getBaseAtStartOfView();
            setZoomLevel((int) (JamView.this.nbasesInView * 1.1));
            goToBasePosition(startBase);
            repaint();
            break;
          case KeyEvent.VK_DOWN:
            if (showBaseAlignment)
              break;
            startBase = getBaseAtStartOfView();
            setZoomLevel((int) (JamView.this.nbasesInView * .9));
            goToBasePosition(startBase);
            repaint();
            break;
          default:
            break;
          }
        }
      });
    }

    addMouseListener(new PopupListener());
    setFocusable(true);
    requestFocusInWindow();
    addFocusListener(new FocusListener()
    {
      public void focusGained(FocusEvent fe){}
      public void focusLost(FocusEvent fe){}
    });
  }
  
  /**
   * Check to see if this component is contained by the display
   * (FeatureDisplay) component.
   * @return
   */
  private boolean containsComponent(JPanel topPanel, JPanel mainPanel)
  {
    Component[] c = mainPanel.getComponents();
    for(int i=0; i<c.length; i++)
    {
      if(c[i].equals(topPanel))
        return true;
    }
    
    return false;
  }
  
  public void setVisible(boolean visible)
  {
    super.setVisible(visible);
    mainPanel.setVisible(visible);
  }
  
  private int getBaseAtStartOfView()
  {
    String refName = (String) combo.getSelectedItem();
    int seqLength = seqLengths.get(refName);
    double x = jspView.getViewport().getViewRect().getX();
    return (int) (seqLength * ( x / getWidth()));
  }
  
  /**
   * Set the panel size based on the number of bases visible
   * and repaint.
   * @param nbasesInView
   */
  private void setZoomLevel(final int nbasesInView)
  {
    this.nbasesInView = nbasesInView;
    int start = getBaseAtStartOfView();
    setDisplay(start,start+nbasesInView, null);
    revalidate();
    repaint();
  }
  
  /**
   * Set the ViewPort so it starts at the given base position.
   * @param base
   */
  private void goToBasePosition(int base)
  {
    Point p = jspView.getViewport().getViewPosition();
    String refName = (String) combo.getSelectedItem();
    int seqLength = seqLengths.get(refName);
    p.x = (int) ((getPreferredSize().width)*(((float)base-1)/(float)seqLength));
    jspView.getViewport().setViewPosition(p);
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
    
    String refName = (String) combo.getSelectedItem();
    int seqLength = seqLengths.get(refName);

    double pixPerBase = getPixPerBase(start, end);
    if(pixPerBase*2 > ALIGNMENT_PIX_PER_BASE)
    {
      pixPerBase = ALIGNMENT_PIX_PER_BASE;
      checkBoxSingle.setVisible(false);
      jspView.getVerticalScrollBar().setValue(0);
      jspView.setColumnHeaderView(ruler);
      showBaseAlignment = true;
    }
    else if(jspView != null)
    {
      checkBoxSingle.setVisible(true);
      jspView.setColumnHeaderView(null);
      showBaseAlignment = false;
    }
    
    Dimension d = new Dimension();
    d.setSize((seqLength*pixPerBase), 800.d);
    setPreferredSize(d);
    goToBasePosition(startBase);
    
    repaint();
    if(event == null || 
       event.getType() == DisplayAdjustmentEvent.SCALE_ADJUST_EVENT)
      revalidate();
    
    this.startBase = -1;
    this.endBase   = -1;
  }
  
  private double getPixPerBase(int start, int end)
  {
    if(feature_display == null)
      return 1000.d/(double)(end-start+1);
    else
      return feature_display.getWidth()/(double)(end-start+1);
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
    
    System.out.println(entryFileName);
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
    return showScale;
  }

  public void setShowScale(boolean showScale)
  {
    this.showScale = showScale;
  }
  
  public JScrollPane getJspView()
  {
    return jspView;
  }
  
  public void setBases(Bases bases)
  {
    this.bases = bases;
  }
  
  /**
   * Remove JScrollPane border
   */
  public void removeBorder()
  {
    Border empty = new EmptyBorder(0,0,0,0);
    jspView.setBorder(empty);
  }
  
  /**
   *  Handle a mouse drag event on the drawing canvas.
   **/
  private void handleCanvasMouseDragOrClick(final MouseEvent event)
  {
    if(!showBaseAlignment || event.isShiftDown()) 
      return;

    String refName = (String) combo.getSelectedItem();
    int seqLength = seqLengths.get(refName);
    float pixPerBase = ((float)getWidth())/(float)(seqLength);    
    int start = (int) Math.round(event.getPoint().getX()/pixPerBase);
    if(start < 1)
      start = 1;
    if(start > seqLength)
      start = seqLength;
    
    MarkerRange drag_range;
    try
    {
      drag_range = new MarkerRange (bases.getForwardStrand(), start, start);
      getSelection().setMarkerRange(drag_range);
      repaint();
    }
    catch (OutOfRangeException e)
    {
      e.printStackTrace();
    }
  }
  
  private Selection getSelection()
  {
    return selection;
  }
  
  private class Ruler extends JPanel
  {
    private static final long serialVersionUID = 1L;
    int start;
    int end;

    public Ruler()
    {
      super();
      setPreferredSize(new Dimension(getPreferredSize().width, 15));
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

      for(int i=startMark; i<end; i+=10)
      {
        int xpos = (i-1-start)*ALIGNMENT_PIX_PER_BASE;
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
    public void mouseClicked(MouseEvent e)
    {
      JamView.this.requestFocus();
      handleCanvasMouseDragOrClick(e);
    }
    
    public void mousePressed(MouseEvent e)
    {
      maybeShowPopup(e);
    }

    public void mouseReleased(MouseEvent e)
    {
      maybeShowPopup(e);
    }

    private void maybeShowPopup(MouseEvent e)
    {
      if(e.isPopupTrigger())
        popup.show(e.getComponent(),
                e.getX(), e.getY());
    }
  }
 
  
  class SAMRecordComparator implements Comparator<Object>
  {
    public int compare(Object o1, Object o2) 
    {
      SAMRecord pr1 = (SAMRecord) o1;
      SAMRecord pr2 = (SAMRecord) o2;
      
      int cmp = pr1.getReadName().compareTo(pr2.getReadName());
      
      if(cmp == 0)
      {
        if(pr1.getAlignmentStart() < pr2.getAlignmentStart())
          return -1;
        else
          return 1;
      }
      return cmp;
    }
  }


  public Dimension getPreferredScrollableViewportSize()
  {
    return getPreferredSize();
  }

  public int getScrollableBlockIncrement(Rectangle visibleRect,
      int orientation, int direction)
  {
    if (orientation == SwingConstants.HORIZONTAL)
      return visibleRect.width - maxUnitIncrement;
    else 
      return visibleRect.height - maxUnitIncrement;
  }

  public boolean getScrollableTracksViewportHeight()
  {
    return false;
  }

  public boolean getScrollableTracksViewportWidth()
  {
    return false;
  }

  public int getScrollableUnitIncrement(Rectangle visibleRect, int orientation,
                                        int direction)
  {
  //Get the current position.
    int currentPosition = 0;
    if (orientation == SwingConstants.HORIZONTAL) 
        currentPosition = visibleRect.x;
    else 
        currentPosition = visibleRect.y;

    //Return the number of pixels between currentPosition
    //and the nearest tick mark in the indicated direction.
    if (direction < 0)
    {
      int newPosition = currentPosition -
                        (currentPosition / maxUnitIncrement)
                         * maxUnitIncrement;
      return (newPosition == 0) ? maxUnitIncrement : newPosition;
    } 
    else 
    {
      return ((currentPosition / maxUnitIncrement) + 1)
              * maxUnitIncrement
              - currentPosition;
    }
  }
  
  /**
   * Artemis event notification
   */
  public void displayAdjustmentValueChanged(DisplayAdjustmentEvent event)
  {
    setDisplay(event.getStart(), event.getEnd(), event);
  }
  
  public static void main(String[] args)
  {
    String bam = args[0];
    int nbasesInView = 1000;
    String reference = null;
    
    for(int i=0;i<args.length; i++)
    {
      if(args[i].equals("-a"))
        bam = args[++i];
      else if(args[i].equals("-r"))
        reference = args[++i];
      else if(args[i].equals("-v"))
        nbasesInView = Integer.parseInt(args[++i]);
      else if(args[i].equals("-s"))
        System.setProperty("samtoolDir", args[++i]);
      else if(args[i].startsWith("-h"))
      { 
        System.out.println("-h\t show help");
        
        System.out.println("-a\t BAM/SAM file to display");
        System.out.println("-r\t reference file (optional)");
        System.out.println("-v\t number of bases to display in the view (optional)");
        System.out.println("-s\t samtool directory");

        System.exit(0);
      }
    }

    final JamView view = new JamView(bam, reference, nbasesInView);
    JFrame frame = new JFrame("JAM");
    
    frame.addWindowFocusListener(new WindowFocusListener()
    {
      public void windowGainedFocus(WindowEvent e)
      {
        view.requestFocus();
      }
      public void windowLostFocus(WindowEvent e){}
    });
    
    view.addJamToPanel((JPanel)frame.getContentPane(), false, null);
    frame.pack();

    view.jspView.getVerticalScrollBar().setValue(
        view.jspView.getVerticalScrollBar().getMaximum());
    frame.setVisible(true);
  }

  public void selectionChanged(SelectionChangeEvent event)
  {
    repaint();
  }
}
