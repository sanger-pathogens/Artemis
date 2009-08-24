
package uk.ac.sanger.artemis.components.alignment;

import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.Stroke;
import java.awt.event.FocusEvent;
import java.awt.event.FocusListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.StringReader;
import java.util.Collections;
import java.util.Comparator;
import java.util.Hashtable;
import java.util.List;
import java.util.Vector;

import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

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
{
  private static final long serialVersionUID = 1L;
  private List<Read> readsInView;
  private Hashtable<String, Integer> seqLengths = new Hashtable<String, Integer>();
  private Vector<String> seqNames = new Vector<String>();
  private String bam;
  private String sam;
  private EntryGroup entryGroup;
  private JScrollPane jspView;
  private JComboBox combo;
  private int nbasesInView;
 
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
  }
  
  private void readHeader()
  {
	String cmd[] = { "/Users/tjc/BAM/samtools-0.1.5c/samtools",  
			     "view", "-H", bam };
	
    RunSamTools samtools = new RunSamTools(cmd, null, null);
	
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
  
  private void readFromBam(int start, int end)
  {
    String refName = (String) combo.getSelectedItem();
	String cmd[] = { "/Users/tjc/BAM/samtools-0.1.5c/samtools",  
				     "view", bam, refName+":"+start+"-"+end };
		
	for(int i=0; i<cmd.length;i++)
	  System.out.print(cmd[i]+" ");
	System.out.println();
	RunSamTools samtools = new RunSamTools(cmd, null, null);
		
	if(samtools.getProcessStderr() != null)
      System.out.println(samtools.getProcessStderr());
	    
	sam = samtools.getProcessStdout();

	StringReader samReader = new StringReader(sam);
	BufferedReader buff = new BufferedReader(samReader);
	String line;
	try 
	{
	  if(readsInView == null)
	    readsInView = new Vector<Read>();
	  else
		readsInView.clear();
	  while((line = buff.readLine()) != null)
	  {
		String fields[] = line.split("\t");

		Read pread = new Read();
		pread.qname = fields[0];
		pread.flag  = Integer.parseInt(fields[1]);
		pread.rname = fields[2];
		pread.pos   = Integer.parseInt(fields[3]);
		pread.mapq  = Short.parseShort(fields[4]);
		pread.mpos  = Integer.parseInt(fields[7]);
		pread.isize = Integer.parseInt(fields[8]);
		pread.seq   = fields[9];
		readsInView.add(pread);
	  }
	   	  
	  Collections.sort(readsInView, new ReadComparator());
	} 
	catch (IOException e) 
	{
	  e.printStackTrace();
	}
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
	
	int start = (int) (seqLength * ( x / getWidth()));
	int end   = (int) (start + (jspView.getViewport().getWidth() / pixPerBase));
	
	int scaleHeight = 15;
	readFromBam(start, end);
	
	drawScale(g2, start, end, pixPerBase);
	
	Stroke originalStroke = g2.getStroke();
	Stroke stroke =
		    new BasicStroke (2.0f, BasicStroke.CAP_BUTT, BasicStroke.JOIN_ROUND);

	for(int i=0; i<readsInView.size(); i++)
	{
      Read thisRead = readsInView.get(i);
      Read nextRead = null;
      
      if( (thisRead.flag & 0x0001) != 0x0001 ||
    	  (thisRead.flag & 0x0008) == 0x0008 )  // not single read
    	continue;
      
      if(i < readsInView.size()-1)
      {
    	nextRead = readsInView.get(i+1);
    	i++;
      }
      
      if(nextRead != null &&
         (nextRead.flag & 0x0001) == 0x0001 &&
         (nextRead.flag & 0x0008) != 0x0008 )
      {         
    	if(thisRead.qname.equals(nextRead.qname))
    	{
          if( (thisRead.flag & 0x0010) == 0x0010 &&
              (nextRead.flag & 0x0010) == 0x0010 )
            g2.setColor(Color.red);
          else
            g2.setColor(Color.blue);
          
          int ypos = (getHeight() - scaleHeight) - ( Math.abs(thisRead.isize) );
          int thisStart = drawRead(g2, thisRead, pixPerBase, stroke, ypos);
    	  int nextStart = drawRead(g2, nextRead, pixPerBase, stroke, ypos);
    	  g2.setStroke(originalStroke);
    	  g2.setColor(Color.LIGHT_GRAY);
    	  g2.drawLine(thisStart, ypos, nextStart, ypos);
    	}
    	else
    	  i--;
      }
	}
  }
  
  private int getBaseAtStartOfView()
  {
    String refName = (String) combo.getSelectedItem();
    int seqLength = seqLengths.get(refName);
    double x = jspView.getViewport().getViewRect().getX();
    return (int) (seqLength * ( x / getWidth()));
  }
  
  private void drawScale(Graphics2D g2, int start, int end, float pixPerBase)
  {
    g2.setColor(Color.black);
    g2.drawLine( (int)(start*pixPerBase), getHeight()-14, (int)(end*pixPerBase), getHeight()-14);
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
  }
  
  private void drawTicks(Graphics2D g2, int start, int end, float pixPerBase, int division)
  {
    int markStart = (Math.round(start/division)*division)+division;
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
    }
  }
  
  private int drawRead(Graphics2D g2, Read read,
		                float pixPerBase, Stroke stroke, int ypos)
  {
    int thisStart = (int) (read.pos*pixPerBase);
    int thisEnd   = (int) (thisStart + (read.seq.length()*pixPerBase));
    g2.setStroke(stroke);
    g2.drawLine(thisStart, ypos, thisEnd, ypos);
    return thisStart;
  }



  
  public void showJFrame()
  {
    this.nbasesInView = nbasesInView;
    JFrame frame = new JFrame("JamTool");

    combo = new JComboBox(seqNames);
    combo.setEditable(false);
    combo.addItemListener(new ItemListener()
    {
      public void itemStateChanged(ItemEvent e)
      {
        setZoomLevel(JamView.this.nbasesInView);
      }
    });
    
    frame.setPreferredSize(new Dimension(1000,800));
    setLength(nbasesInView);
    
    jspView = new JScrollPane(this, JScrollPane.VERTICAL_SCROLLBAR_NEVER,
                                    JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
        
    JPanel panel = (JPanel) frame.getContentPane();

    panel.setLayout(new BorderLayout());
    panel.add(combo, BorderLayout.NORTH);
    panel.add(jspView, BorderLayout.CENTER);
    frame.pack();
    
    setFocusable(true);
    requestFocusInWindow();
    
    addKeyListener(new KeyAdapter()
    {
      public void keyPressed(final KeyEvent event)
      {
        switch(event.getKeyCode())
        {
          case KeyEvent.VK_UP:
            setZoomLevel( (int) (JamView.this.nbasesInView*1.1) );
            repaint();
            break;
          case KeyEvent.VK_DOWN:
            setZoomLevel( (int) (JamView.this.nbasesInView*.9) );
            repaint();
            break;
          default:
            break;
        }
      }
    });
    
    addFocusListener(new FocusListener() 
    {
      public void focusGained(FocusEvent fe) {}
      public void focusLost(FocusEvent fe) {}
    });

    frame.setVisible(true);
  }
  
  private void setZoomLevel(final int nbasesInView)
  {
    int startBase = getBaseAtStartOfView();
   
    this.nbasesInView = nbasesInView;
    setLength(this.nbasesInView);
    revalidate();
    repaint();
    System.out.println(startBase+"  "+getBaseAtStartOfView());
    
    String refName = (String) combo.getSelectedItem();
    int seqLength = seqLengths.get(refName);


    float pixPerBase = ((float)getWidth())/(float)(seqLength);
    Point p = jspView.getViewport().getViewPosition();
    p.x = (int)(startBase*pixPerBase);
    jspView.getViewport().setViewPosition(p);
    revalidate();
    repaint();
    
    System.out.println(startBase+"  "+getBaseAtStartOfView());
  }

  private void setLength(int basesToShow)
  {
    String refName = (String) combo.getSelectedItem();
    int seqLength = seqLengths.get(refName);
    float pixPerBase = 1000.f/(float)(basesToShow);
    setPreferredSize( new Dimension((int) (seqLength*pixPerBase), 
                                            getSize().height) );
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
  
  class Read
  {
    String qname; // query name
    int flag;     // bitwise flag
    String rname; // reference name
    int pos;      // leftmost coordinate
    short mapq;   // MAPping Quality
    String cigar; // extended CIGAR string
    String mrnm;  // Mate Reference sequence NaMe; Ò=Ó if the same as rname
    int mpos;     // leftmost Mate POSition
    int isize;    // inferred Insert SIZE
    String seq;   // query SEQuence; Ò=Ó for a match to the reference; n/N/. for ambiguity
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
  
  public static void main(String[] args)
  {
    String bam = args[0];
    int nbasesInView = Integer.parseInt(args[1]);
    String reference = null;
    if(args.length > 2)
      reference = args[2];
    
	JamView view = new JamView(bam, reference, nbasesInView);
	view.showJFrame();
  }
}
