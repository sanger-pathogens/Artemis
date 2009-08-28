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
import java.awt.Dimension;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
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
import java.awt.event.WindowEvent;
import java.awt.event.WindowFocusListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.StringReader;
import java.util.Collections;
import java.util.Comparator;
import java.util.Hashtable;
import java.util.List;
import java.util.Vector;

import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.Scrollable;
import javax.swing.SwingConstants;
import javax.swing.UIManager;

import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.SimpleEntryGroup;
import uk.ac.sanger.artemis.components.EntryFileDialog;
import uk.ac.sanger.artemis.components.MessageDialog;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.sequence.Bases;
import uk.ac.sanger.artemis.sequence.NoSequenceException;
import uk.ac.sanger.artemis.util.Document;
import uk.ac.sanger.artemis.util.DocumentFactory;
import uk.ac.sanger.artemis.util.OutOfRangeException;

public class JamView extends JPanel
                     implements Scrollable
{
  private static final long serialVersionUID = 1L;
  private List<Read> readsInView;
  private Hashtable<String, Integer> seqLengths = new Hashtable<String, Integer>();
  private Vector<String> seqNames = new Vector<String>();
  private String bam;

  private EntryGroup entryGroup;
  private JScrollPane jspView;
  private JComboBox combo;
  private JCheckBox checkBoxSingle;
  Ruler ruler = new Ruler();
  private int nbasesInView;
  private int laststart;
  private int lastend;
  private int maxUnitIncrement = 4;
  private int ALIGNMENT_PIX_PER_BASE;

 
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
      entryGroup = new SimpleEntryGroup();
      try
      {
        getEntry(reference,entryGroup);
      }
      catch (NoSequenceException e)
      {
        e.printStackTrace();
      }
    }
    
    readHeader();
    
    // set font size
    setFont(getFont().deriveFont(12.f));
    final javax.swing.plaf.FontUIResource font_ui_resource =
      new javax.swing.plaf.FontUIResource(getFont());
   //  Options.getOptions().getFontUIResource();

    java.util.Enumeration keys = UIManager.getDefaults().keys();
    while(keys.hasMoreElements()) 
    {
      Object key = keys.nextElement();
      Object value = UIManager.get(key);
      if(value instanceof javax.swing.plaf.FontUIResource) 
        UIManager.put(key, font_ui_resource);
    }
    FontMetrics fm  = getFontMetrics(getFont());
    ALIGNMENT_PIX_PER_BASE  = (int) (fm.stringWidth("A")*1.1);
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
      readsInView = new Vector<Read>();
    else
      readsInView.clear();
	RunSamTools samtools = new RunSamTools(cmd, null, null, readsInView);
		
	if(samtools.getProcessStderr() != null)
      System.out.println(samtools.getProcessStderr());
	    
	samtools.waitForStdout();
  }
  
  /**
   * Override
   */
  public void paintComponent(Graphics g)
  {
	super.paintComponent(g);
	Graphics2D g2 = (Graphics2D)g;

	String refName = (String) combo.getSelectedItem();
    int seqLength = seqLengths.get(refName);

	float pixPerBase = ((float)getWidth())/(float)(seqLength);
	
    double x = jspView.getViewport().getViewRect().getX();
    int start = (int) (seqLength * ( (float)x / (float)getWidth()));
    int end   = (int) (start + ((float)jspView.getViewport().getWidth() / 
                                (float)pixPerBase));

    if(laststart != start ||
       lastend   != end)
    {
      try
      {
        readFromBam(start, end);      
        if(pixPerBase < ALIGNMENT_PIX_PER_BASE)
          Collections.sort(readsInView, new ReadComparator());
      }
      catch(OutOfMemoryError ome)
      {
        JOptionPane.showMessageDialog(this, "Out of Memory");
        return;
      }
    }
    
    laststart = start;
    lastend   = end;
	if(pixPerBase >= ALIGNMENT_PIX_PER_BASE)
	  drawBaseAlignment(g2, seqLength, pixPerBase, start, end);
	else
      drawLineView(g2, seqLength, pixPerBase, start, end);
  }
  
  private void drawBaseAlignment(Graphics2D g2, int seqLength, 
                                 float pixPerBase, final int start, final int end)
  {
    FontMetrics fm =  getFontMetrics(getFont());
    int ypos = 0;
    
    ruler.start = start;
    ruler.end = end;
    ruler.repaint();
    
    //drawBaseScale(g2, start, end, ypos);
    boolean draw[] = new boolean[readsInView.size()];
    for(int i=0; i<readsInView.size(); i++)
      draw[i] = false;
    
    ypos+=6;
    
    for(int i=0; i<readsInView.size(); i++)
    {
      if (!draw[i])
      {
        Read thisRead = readsInView.get(i);
        ypos+=10;

        drawSequence(g2, thisRead, pixPerBase, ypos);
        draw[i] = true;
        
        int thisEnd = thisRead.pos+thisRead.seq.length();
        for(int j=i+1; j<readsInView.size(); j++)
        {
          if (!draw[j])
          {
            Read nextRead = readsInView.get(j);
            if(nextRead.pos > thisEnd+1)
            {
              drawSequence(g2, nextRead, pixPerBase, ypos);
              draw[j] = true;
              thisEnd = nextRead.pos+nextRead.seq.length();
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
  private void drawSequence(Graphics2D g2, Read read, float pixPerBase, int ypos)
  {
    if ((read.flag & 0x0001) != 0x0001 || 
        (read.flag & 0x0008) == 0x0008)
      g2.setColor(Color.black);
    else
      g2.setColor(Color.blue);
    
    int xpos;
    
    for(int i=0;i<read.seq.length(); i++)
    {
      xpos = ((read.pos-1) + i)*ALIGNMENT_PIX_PER_BASE;
      g2.drawString(read.seq.substring(i, i+1), xpos, ypos);
    }
    
  }
  
  private void drawLineView(Graphics2D g2, int seqLength, float pixPerBase, int start, int end)
  {   
    drawScale(g2, start, end, pixPerBase);
    
    Stroke originalStroke = g2.getStroke();
    Stroke stroke =
            new BasicStroke (1.2f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_ROUND);
    int scaleHeight = 15;
    
    for(int i=0; i<readsInView.size(); i++)
    {
      Read thisRead = readsInView.get(i);
      Read nextRead = null;      

      if( (thisRead.flag & 0x0001) != 0x0001 || // read is not paired in sequencing
          (thisRead.flag & 0x0008) == 0x0008 )  // mate is unmapped 
      {
        if(checkBoxSingle.isSelected())
        {
          System.out.println("HERE "+thisRead.qname+
              "\t\t pos "+thisRead.pos+"\tseq length "+thisRead.seq.length()+
              "\tisize "+thisRead.isize);
          int ypos = (getHeight() - scaleHeight) - thisRead.seq.length();
          g2.setColor(Color.orange);
          drawRead(g2, thisRead, pixPerBase, stroke, ypos);
        }
        continue;
      }
      
      int ypos = (getHeight() - scaleHeight) - ( Math.abs(thisRead.isize) );
      if(i < readsInView.size()-1)
      {
        nextRead = readsInView.get(++i);
        
        if(thisRead.qname.equals(nextRead.qname))
        {
          if( (thisRead.flag & 0x0010) == 0x0010 && // strand of the query (1 for reverse)
              (nextRead.flag & 0x0010) == 0x0010 )
            g2.setColor(Color.red);
          else
            g2.setColor(Color.blue);
 
          int thisEnd = drawRead(g2, thisRead, pixPerBase, stroke, ypos)[1];
          int nextStart = drawRead(g2, nextRead, pixPerBase, stroke, ypos)[0];
          
          if(thisEnd < nextStart && (nextStart-thisEnd)*pixPerBase > 2.f)
          {
            g2.setStroke(originalStroke);
            g2.setColor(Color.LIGHT_GRAY);
            g2.drawLine((int)(thisEnd*pixPerBase), ypos, (int)(nextStart*pixPerBase), ypos);
          }
        }
        else
        {
          drawLoneRead(g2, thisRead, ypos, pixPerBase, originalStroke, stroke);
          i--;
        }
      }
      else
      {
        drawLoneRead(g2, thisRead, ypos, pixPerBase, originalStroke, stroke);
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
  private void drawLoneRead(Graphics2D g2, Read thisRead, int ypos, 
      float pixPerBase, Stroke originalStroke, Stroke stroke)
  {
    boolean drawLine = true;
    g2.setColor(Color.blue); 
    if(ypos <= 0)
    {
      ypos = thisRead.seq.length();
      drawLine = false;
      g2.setColor(Color.orange); 
    }

    int thisEnd = drawRead(g2, thisRead, pixPerBase, stroke, ypos)[1];
    if(drawLine)
    {
      g2.setStroke(originalStroke);
      g2.setColor(Color.LIGHT_GRAY);
      int nextStart = (int) ((thisRead.mpos-1)*pixPerBase);
      g2.drawLine((int)(thisEnd*pixPerBase), ypos, nextStart, ypos);
    }
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
      drawTicks(g2, start, end, pixPerBase, 200);
    else
      drawTicks(g2, start, end, pixPerBase, 50);
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
  
  private int[] drawRead(Graphics2D g2, Read read,
		               float pixPerBase, Stroke stroke, int ypos)
  {
    int thisStart = read.pos-1;
    int thisEnd   = thisStart + read.seq.length();
    g2.setStroke(stroke);
    g2.drawLine((int)(thisStart*pixPerBase), ypos, (int)(thisEnd*pixPerBase), ypos);
    return new int[] { thisStart, thisEnd };
  }
  
  public void addToPanel(final JPanel panel)
  {
    JPanel topPanel = new JPanel(new GridBagLayout());
    GridBagConstraints gc = new GridBagConstraints();
    
    combo = new JComboBox(seqNames);
    combo.setEditable(false);
    combo.addItemListener(new ItemListener()
    {
      public void itemStateChanged(ItemEvent e)
      {
        laststart = -1;
        lastend   = -1;
        setZoomLevel(JamView.this.nbasesInView);
      }
    });
    gc.fill = GridBagConstraints.NONE;
    gc.anchor = GridBagConstraints.FIRST_LINE_START;
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
    
    panel.setPreferredSize(new Dimension(1000,500));
    setLength(nbasesInView);
    
    jspView = new JScrollPane(this, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
                                    JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        
    panel.setLayout(new BorderLayout());
    panel.add(topPanel, BorderLayout.NORTH);
    panel.add(jspView, BorderLayout.CENTER);

    jspView.getVerticalScrollBar().setValue(
        jspView.getVerticalScrollBar().getMaximum());
       
    addMouseListener(new MouseAdapter()
    {
      public void mouseClicked(MouseEvent e)
      {
        JamView.this.requestFocus();
      }
    });
    
    addKeyListener(new KeyAdapter()
    {
      public void keyPressed(final KeyEvent event)
      {
        switch(event.getKeyCode())
        {
          case KeyEvent.VK_UP:
            int startBase = getBaseAtStartOfView();
            setZoomLevel( (int) (JamView.this.nbasesInView*1.1) );
            setViewToBasePosition(startBase);
            repaint();
            break;
          case KeyEvent.VK_DOWN:
            startBase = getBaseAtStartOfView();
            setZoomLevel( (int) (JamView.this.nbasesInView*.9) );
            setViewToBasePosition(startBase);
            repaint();
            break;
          default:
            break;
        }
      }
    });
    
    setFocusable(true);
    requestFocusInWindow();
    addFocusListener(new FocusListener() 
    {
      public void focusGained(FocusEvent fe) {}
      public void focusLost(FocusEvent fe) {}
    });
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
    setLength(this.nbasesInView);
    revalidate();
    repaint();
  }
  
  /**
   * Set the ViewPort so it starts at the given base position.
   * @param base
   */
  private void setViewToBasePosition(int base)
  {
    Point p = jspView.getViewport().getViewPosition();
    
    String refName = (String) combo.getSelectedItem();
    int seqLength = seqLengths.get(refName);
    p.x = (int) ((getPreferredSize().width)*(((float)base)/(float)seqLength));
    jspView.getViewport().setViewPosition(p);
  }

  /**
   * Set the panel size based on the number of bases visible.
   * @param nbasesInView
   */
  private void setLength(int basesToShow)
  {
    String refName = (String) combo.getSelectedItem();
    int seqLength = seqLengths.get(refName);
    double pixPerBase = 1000.d/(double)(basesToShow);
    
    System.out.println(pixPerBase+"  "+ALIGNMENT_PIX_PER_BASE);
    if(pixPerBase > ALIGNMENT_PIX_PER_BASE)
    {
      pixPerBase = ALIGNMENT_PIX_PER_BASE;
      jspView.getVerticalScrollBar().setValue(0);
      jspView.setColumnHeaderView(ruler);
    }
    else if(jspView != null)
      jspView.setColumnHeaderView(null);
    Dimension d = new Dimension();
    d.setSize((seqLength*pixPerBase), 800.d);
    setPreferredSize(d);
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
  
  class Ruler extends JPanel
  {
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
  }
  
  class ReadComparator implements Comparator
  {
    public int compare(Object o1, Object o2) 
    {
      Read pr1 = (Read) o1;
      Read pr2 = (Read) o2;
      
      return pr1.qname.compareTo(pr2.qname);
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
    
    view.addToPanel((JPanel)frame.getContentPane());
    frame.pack();
    view.jspView.getVerticalScrollBar().setValue(
        view.jspView.getVerticalScrollBar().getMaximum());
    frame.setVisible(true);
  }
}
