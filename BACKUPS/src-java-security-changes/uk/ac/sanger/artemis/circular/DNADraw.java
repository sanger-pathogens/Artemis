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

import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Image;
import java.awt.Point;
import java.awt.RenderingHints;
import java.awt.event.*;
import java.awt.geom.AffineTransform;
import java.awt.image.BufferedImage;
import java.awt.image.RenderedImage;
import java.awt.print.PageFormat;
import java.awt.print.Printable;
import java.awt.print.PrinterException;
import java.awt.print.PrinterJob;

import java.util.Collections;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.List;
import java.util.Vector;
import java.awt.datatransfer.Transferable;
import java.awt.dnd.DnDConstants;
import java.awt.dnd.DragGestureEvent;
import java.awt.dnd.DragGestureListener;
import java.awt.dnd.DragSource;
import java.awt.dnd.DragSourceDragEvent;
import java.awt.dnd.DragSourceDropEvent;
import java.awt.dnd.DragSourceEvent;
import java.awt.dnd.DragSourceListener;
import java.awt.dnd.DropTarget;
import java.awt.dnd.DropTargetDragEvent;
import java.awt.dnd.DropTargetDropEvent;
import java.awt.dnd.DropTargetEvent;
import java.awt.dnd.DropTargetListener;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

import javax.swing.JButton;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.KeyStroke;
import javax.swing.SwingUtilities;

import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.editor.BrowserControl;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.sequence.Bases;


public class DNADraw extends ScrollPanel
                     implements Printable, DragGestureListener,
                     DragSourceListener, DropTargetListener
{
  private static final long serialVersionUID = 1L;
  public static JScrollPane jsp = new JScrollPane();
  private DNADraw current_dna;
  private JFrame mainFrame;

  private Point location    = new Point();
  private Dimension border  = new Dimension(150,100);
  private Dimension panelSize = new Dimension(690,690);
  private Dimension linearPanelSize = new Dimension(800,350);
  private Hashtable lineAttr;
  private Vector minorTicks;
  private Vector majorTicks;
  private Vector block;
  private Vector restrictionEnzyme;

  private int startTick = 0;
  private int minorTick = 100;
  private int majorTick = 500;

//
// store the tick positions -- there appears to be
// a bug in AffineTransform when it comes to using
// elements from the matrix when printing
//
  private int[] tickMajorXPositions;
  private int[] tickMajorYPositions;

  private int[] tickMinorXPositions;
  private int[] tickMinorYPositions;

  private boolean labelTicks = true;
  
  private int[] reXPositions;
  private int[] reYPositions;
  private boolean close = false;
  private EntryGroup artemisEntryGroup;
  private Bases bases;
  private Graph gcGraph;
  private Graph gcSkewGraph;
  private Vector<Graph> userGraphs = new Vector<Graph>();
  
  final JMenu userDraw = new JMenu("Draw");
  final JMenu userOptions = new JMenu("Options");
  
  protected static int Y_SHIFT = -35;
  protected static int THETA = -90;
  private TrackManager trackManager;
  private AffineTransform original;
  
  // linear plot variables
  private JMenuItem linearPlotOptions = new JMenuItem("Linear plot...");;
  private int numberOfLines;
  private int basesPerLine = 20000;
  private float lineHeight = 200.f;
  private float singleBaseWidth;
  private int borderWidth2;
  private int borderHeight2;
  
  public DNADraw()
  {
    super(new BorderLayout());
    location.setLocation(75.d, 75.d);
    
    current_dna = this;
    setBackground(Color.white);
    setPreferredSize(panelSize);
    setOpaque(false);
    setToolTipText("");
    
    if( System.getProperty("java.awt.headless") == null ||
       !System.getProperty("java.awt.headless").equals("true") )
    {
      DragSource dragSource = DragSource.getDefaultDragSource();
      dragSource.createDefaultDragGestureRecognizer(
         this,                             // component where drag originates
         DnDConstants.ACTION_COPY_OR_MOVE, // actions
         this);      
      setDropTarget(new DropTarget(this,this));
    }
    lineAttr = new Hashtable();
    lineAttr.put("start",new Integer(0));
    lineAttr.put("end",new Integer(4000));
    lineAttr.put("lsize",new Integer(5));
    lineAttr.put("circular",new Boolean(true));

    MouseListener mouseListener = new MouseAdapter()
    {
      public void mouseClicked(MouseEvent me)
      {
        if(me.getClickCount() == 2 &&
           !me.isPopupTrigger())
        {
          final Block b = getBlockAtLocation(me.getPoint());
          if(b != null)
          {
            final JFrame f = new JFrame("Properties");
            JButton butt = new JButton("Delete");
            butt.addActionListener(new ActionListener()
            {
              public void actionPerformed(ActionEvent e)
              {
                block.remove(b);
                current_dna.repaint();
                f.setVisible(false);
                f.dispose();
              }
            });
            b.showProperties(f, DNADraw.this, butt);
          }
        }
      }
    };
    this.addMouseListener(mouseListener);
  }
  
  public String getVersion()
  {
    final ClassLoader cl = this.getClass().getClassLoader();
    try
    {
      String line;
      InputStream in = cl.getResourceAsStream("etc/versions");
      BufferedReader reader = new BufferedReader(new InputStreamReader(in));
      while((line = reader.readLine()) != null)
      {
        if(line.startsWith("DNAPlotter"))
          return line.substring( "DNAPlotter".length() ).trim();
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
   * Add list of features to a track
   * @param features
   * @param track
   * @param drawLabel
   */
  public void addFeaturesToTrack(final List features, 
                                 final Track track,
                                 final boolean drawLabel)
  {
    for(int i=0; i<features.size(); i++)
    {
      addFeatureToTrack((uk.ac.sanger.artemis.circular.Feature)features.get(i),
          track, drawLabel);
    }
  }
  
  /**
   * Add a feature to a track
   * @param tmpf
   * @param track
   * @param drawLabel
   */
  public void addFeatureToTrack(final uk.ac.sanger.artemis.circular.Feature tmpf, 
                                final Track track,
                                final boolean drawLabel)
  {
    final uk.ac.sanger.artemis.Feature f = tmpf.getArtemisFeature();
    Vector ranges = f.getLocation().getRanges();
      
    for(int j=0; j<ranges.size(); j++)
    {
      Range range = (Range) ranges.get(j);

      Block drawBlock = new Block(f.getIDString(), 
            range.getStart(),
            range.getEnd(), 
            f.getColour(), 
            10.f, 
            track, this);
      drawBlock.setDrawLabel(drawLabel);
      drawBlock.setFeature(f);
      addBlock(drawBlock);
    }
  }
  
  public String getToolTipText(MouseEvent me) 
  {
    Block b = getBlockAtLocation(me.getPoint());
    if(b != null)
    {
      if(b.getFeature() != null)
        return b.getLabel()+" "+b.getFeature().getLocation().toStringShort();
      
      return b.getLabel()+" "+b.getBstart()+".."+b.getBend();
    }
    return null;
  }

  public DNADraw(Vector minorTicks, Vector majorTicks,
                 Vector block,
                 Vector restrictionEnzyme)
  {
    this();
    this.minorTicks = minorTicks;
    this.block = block;
    this.restrictionEnzyme = restrictionEnzyme;
  }


  public DNADraw(Vector block, Vector restrictionEnzyme,
                 Hashtable lineAttr, int startTick,
                 int minorTick, int majorTick)
  {
    this();
    this.block = block;
    this.restrictionEnzyme = restrictionEnzyme;
    this.lineAttr     = lineAttr;
    this.startTick    = startTick;
    this.minorTick    = minorTick;
    this.majorTick    = majorTick;

    if(!isCircular())
      setPreferredSize(linearPanelSize);
    
    basesPerLine = majorTick;
    calculateTickPosistions();
  }


  /**
  *
  * Get the width/diameter of the DNA map
  *
  */
  protected double getDiameter()
  {
    return ((double)getWidth())-border.getWidth();
  }


  protected Point getLocationPoint()
  {
    return location;
  }


  protected void zoomIn()
  {
    int wid = getWidth();
    wid     = wid+(int)(wid*0.1);
    int hgt = getHeight();
    if(isCircular())
      hgt = hgt+(int)(hgt*0.1);
    zoom(wid,hgt);
  }


  protected void zoomOut()
  {
    int wid = getWidth();
    wid     = wid-(int)(wid*0.1);
    int hgt = getHeight();  
    if(isCircular())
      hgt = hgt-(int)(hgt*0.1);
    zoom(wid,hgt);
  }

  
  private void zoom(int wid, int hgt)
  {
    if(isCircular())
    {
      panelSize = new Dimension(wid,hgt);
      setPreferredSize(panelSize);
      setSize(panelSize);
    }
    else
    {
      linearPanelSize = new Dimension(wid,hgt);
      setPreferredSize(linearPanelSize);
      setSize(linearPanelSize);
    }
    repaint();
  }

  
  protected void paintComponent(Graphics g)
  {
    super.paintComponent(g);
    Graphics2D g2 = (Graphics2D)g;  

    if(isCircular())
      drawCircularPanel(g2,true);
    else
      drawLinearPanel(g2);
  }


  protected boolean isCircular()
  {
    return ((Boolean)lineAttr.get("circular")).booleanValue();
  }


  protected void addBlock(Block b)
  {
    if(getGeneticMarker() == null)
      setGeneticMarker(new Vector());
    this.getGeneticMarker().add(b);
    validate();
  }


  protected void drawLinearPanel(Graphics2D g2)
  {
    RenderingHints qualityHints = new RenderingHints(RenderingHints.KEY_ANTIALIASING,
        RenderingHints.VALUE_ANTIALIAS_ON); 
    qualityHints.put(RenderingHints.KEY_RENDERING,               
    RenderingHints.VALUE_RENDER_QUALITY); 
    g2.setRenderingHints(qualityHints);
    
    FontMetrics fm = g2.getFontMetrics();
    double hgt = fm.getAscent();
    g2.setColor(Color.black);
    double widDash = 4;
    
    int lineSize = 5;
    try
    {
      lineSize = getLineSize();
    }
    catch(NullPointerException npe)
    {
      System.out.println("No line size specified using default!");
    }

    borderWidth2  = border.width/2;
    borderHeight2 = border.height/2;
    
    BasicStroke basicstroke = new BasicStroke(
        lineSize,
        BasicStroke.CAP_BUTT, 
        BasicStroke.JOIN_MITER);
    g2.setStroke(basicstroke);
    
    int lastBase  = getEnd();
    numberOfLines = Math.round( (lastBase/basesPerLine) + 0.5f);
    
    int height = (int) ((numberOfLines * lineHeight)+border.height);
    Dimension size = new Dimension(getWidth(),height);
    setPreferredSize(size);
    setSize(size);
    
    singleBaseWidth = (float)(getWidth()-border.width)/(float)basesPerLine;
    
    for(int i=0; i<numberOfLines; i++)
    {
      int ypos = (int) (borderHeight2+lineHeight+(lineHeight*i));
      
      if(i<numberOfLines-1)
        g2.drawLine(borderWidth2,ypos,getWidth()-borderWidth2,ypos);
      else
      {
        int lastLineWidth = (int) ((lastBase-(i*basesPerLine))*singleBaseWidth)+borderWidth2;
        g2.drawLine(borderWidth2,ypos,lastLineWidth,ypos);
      }
    }
    
    g2.setColor(Color.black);
    g2.setStroke(new BasicStroke(1.f));

    if(majorTicks == null || minorTicks == null)
      calculateTickPosistions();

    Enumeration enumTk = minorTicks.elements();
    while(enumTk.hasMoreElements())
    {
      int tick = ((Integer)enumTk.nextElement()).intValue();
      int lineNumber = Math.round((float)(tick-1)/((float)basesPerLine) + 0.5f)-1;
      
      int ypos = (int) (borderHeight2+lineHeight+(lineHeight*lineNumber));
      int xpos = (int) ((tick-(lineNumber*basesPerLine))*singleBaseWidth)+borderWidth2;
      int y = ypos+(int)((lineSize+widDash)/2);
      
      g2.drawLine(xpos,ypos,xpos,y);
    }

    enumTk = majorTicks.elements();
    while(enumTk.hasMoreElements())
    {
      int tick = ((Integer)enumTk.nextElement()).intValue();
      int lineNumber = Math.round((float)(tick-1)/((float)basesPerLine) + 0.5f)-1;
      
      if(lineNumber < 0)
        lineNumber = 0;
      
      int ypos = (int) (borderHeight2+lineHeight+(lineHeight*lineNumber));
      int xpos = (int) ((tick-(lineNumber*basesPerLine))*singleBaseWidth)+borderWidth2;
      int y = ypos+(int)((lineSize+widDash)/2);
      
      g2.drawLine(xpos,ypos,xpos,y);
      
      String label = Integer.toString(tick);
      xpos-=(fm.stringWidth(label)/2);
      y+=hgt;
      g2.drawString(label,xpos,y);
    }

    /*if(restrictionEnzyme != null)
    {
      enumTk = restrictionEnzyme.elements();
      while(enumTk.hasMoreElements())
      {
        Vector re = (Vector)enumTk.nextElement();
        String reLabel = (String)re.elementAt(0);
        int pos = ((Integer)re.elementAt(1)).intValue();
        g2.setColor((Color)re.elementAt(2));
        int x = ((diameter-location.x)*(pos-start)/(end-start))+location.x;
        int y = ymid-(lineSize/2)-(int)widDash;
        g2.drawLine(x,ymid,x,y);
        x-=(fm.stringWidth(reLabel)/2);
        y-=hgt;
        g2.drawString(reLabel,x,y);
      }
    }*/

    // draw features
    Vector markers = getGeneticMarker();
    for(int i=0; i<markers.size(); i++)
    {
      Block b = (Block)markers.get(i);
      b.drawLinear(g2);
    }
  }

  protected void drawCircularPanel(Graphics2D g2, boolean record)
  {
    RenderingHints qualityHints = new RenderingHints(RenderingHints.KEY_ANTIALIASING,
        RenderingHints.VALUE_ANTIALIAS_ON); 
    qualityHints.put(RenderingHints.KEY_RENDERING,               
    RenderingHints.VALUE_RENDER_QUALITY); 
    g2.setRenderingHints(qualityHints);
    g2.setColor(Color.black);
 
    FontMetrics fm = g2.getFontMetrics();
    double hgt = fm.getAscent();
    final double widthPanel  = getWidth(); 
    final double heightPanel = getHeight();

    double rad = 360.d;
    double pi  = Math.PI;
    double widDash = 4;

    double ddiameter  = widthPanel-border.getWidth();
    double ddiameter2 = ddiameter/2.d;
    int diameter = (int)ddiameter;

    int lineSize = 5;
    try
    {
      lineSize = getLineSize();
    }
    catch(NullPointerException npe)
    {
      System.out.println("No line size specified using default!");
    }

    original = g2.getTransform();
    AffineTransform origin = (AffineTransform) original.clone();
    
    origin.translate(0, Y_SHIFT);
    
    // rotate so that origin is at the top
    origin.translate(widthPanel/2, heightPanel/2);
    origin.rotate(Math.toRadians(THETA));
    origin.translate(-widthPanel/2, -heightPanel/2);

    g2.setTransform(origin);
    
    g2.setStroke(new BasicStroke((float)lineSize));
    g2.drawArc(location.x,location.y,
               diameter,diameter,0,360);

    /* draw track circle
     * int shift = (int)(diameter*(1.d-0.9d)/2);
    g2.drawArc(location.x+shift,location.y+shift,
        (int)(diameter*0.9),(int)(diameter*0.9),0,360);*/

    AffineTransform newOrig;

    if(restrictionEnzyme != null)
    {
      if(record)
      {
        int nsize = restrictionEnzyme.size();
        reXPositions = new int[nsize];
        reYPositions = new int[nsize];
      }
      Enumeration enumRes = restrictionEnzyme.elements();
      while(enumRes.hasMoreElements())
      {
        Vector re = (Vector)enumRes.nextElement();
        String reLabel = (String)re.elementAt(0);
        int pos = ((Integer)re.elementAt(1)).intValue();
        g2.setColor((Color)re.elementAt(2));
        double ang = getAngleFromPosition(pos,rad);

        newOrig = (AffineTransform)(origin.clone());
        newOrig.rotate(Math.toRadians(-ang),
                       widthPanel/2.d,heightPanel/2.d);

        int widLabel  = (lineSize+fm.stringWidth(reLabel))/2;
        int widREDash = (int)(widDash+widDash+lineSize)+widLabel;

        int x = 0;
        int y = 0;
        if(record)
        {
          x = (int)( ddiameter2 + (newOrig.getScaleX()*
                       (ddiameter2 + 10 + widLabel + widREDash) ) -
                       widLabel );
          y = (int)( ddiameter2 + (newOrig.getShearY()*
                       (ddiameter2 + 10 + widREDash + (hgt/2.d)) ) +
                       hgt/2.d );

          int index = restrictionEnzyme.indexOf(re);
          reXPositions[index] = x;
          reYPositions[index] = y;
        }
        else
        {
          int index = restrictionEnzyme.indexOf(re);
          x = reXPositions[index];
          y = reYPositions[index];
        }

        AffineTransform labelTransform = (AffineTransform)(origin.clone());
        labelTransform.translate(location.x+ddiameter2, location.y+ddiameter2);
        labelTransform.rotate(Math.toRadians(-DNADraw.THETA));
        labelTransform.translate(-(location.x+ddiameter2), -(location.y+ddiameter2));
        g2.setTransform(labelTransform);
        
        
        g2.drawString(reLabel,location.x+x,location.y+y);
        g2.setTransform(newOrig);
        g2.setStroke(new BasicStroke(1.f));
        int xLine = location.x+(int)(ddiameter);
        int yLine = location.y+(int)(ddiameter/2.d);
        g2.drawLine(xLine,yLine,(int)(xLine+widREDash),yLine);
        g2.setTransform(origin);
      }
    }

    if(majorTicks == null || minorTicks == null)
    {
      calculateTickPosistions();
      if(majorTick>0)
        basesPerLine = majorTick;
    }
    //major ticks
    drawCircularTicks(g2,ddiameter,ddiameter2,diameter,origin,
                      widthPanel,heightPanel,rad,pi,widDash,fm,
                      lineSize,record,majorTicks,false);

    //minor ticks
    drawCircularTicks(g2,ddiameter,ddiameter2,diameter,origin,
                      widthPanel,heightPanel,rad,pi,widDash/2,fm,
                      lineSize,record,minorTicks,true);

    // draw features
    Vector markers = getGeneticMarker();
    for(int i=0; i<markers.size(); i++)
    {
      Block b = (Block)markers.get(i);
      b.drawCircular(g2);
    }
    
  }


  private void drawCircularTicks(Graphics2D g2, double ddiameter,
            double ddiameter2, int diameter, AffineTransform origin,
            double widthPanel,double heightPanel, double rad, double pi,
            double widDash, FontMetrics fm, int lineSize,
            boolean record, Vector ticks, final boolean smallTicks)
  {

    double hgt = fm.getAscent();

    g2.setColor(Color.black);
    if(record)
    {
      int nsize = ticks.size();

      if(smallTicks)
      {
        tickMinorXPositions = new int[nsize];
        tickMinorYPositions = new int[nsize];
      }
      else
      {
        tickMajorXPositions = new int[nsize];
        tickMajorYPositions = new int[nsize];
      }
    }

    AffineTransform newOrig;
    Enumeration enumTk = ticks.elements();
    while(enumTk.hasMoreElements())
    {
      int tick = ((Integer)enumTk.nextElement()).intValue();
      double theta = Math.toRadians(-getAngleFromPosition(tick,rad));
      if(theta > pi)
        theta = theta - pi*2.d;

      newOrig = (AffineTransform)(origin.clone());

      // rotate and add tick mark
      newOrig.rotate(theta,widthPanel/2.d,heightPanel/2.d);
      String label = Integer.toString(tick);
      double wid = fm.stringWidth(label);

      int x = 0;
      int y = 0;
      if(record)
      {
        x = (int)( (ddiameter2) + (newOrig.getScaleX()*
                   (widDash+lineSize+3+(diameter+wid)/2.d)) - (wid/2.d));

        y = (int)( (ddiameter2) + (newOrig.getShearY()*
                   (widDash+lineSize+3+(diameter+hgt)/2.d)) + (hgt/2.d));

        int index = ticks.indexOf(new Integer(tick));

        if(smallTicks)
        {
          tickMinorXPositions[index] = x;
          tickMinorYPositions[index] = y;
        }
        else
        {
          tickMajorXPositions[index] = x;
          tickMajorYPositions[index] = y;
        }
      }
      else    // use stored positions for printing
      {
        int index = ticks.indexOf(new Integer(tick));
        if(smallTicks)
        {
          x = tickMinorXPositions[index];
          y = tickMinorYPositions[index];
        }
        else
        {
          x = tickMajorXPositions[index];
          y = tickMajorYPositions[index];
        }
      }

      if(labelTicks && !smallTicks)        // add tick label
      {
        AffineTransform labelTransform = (AffineTransform)(origin.clone());
        labelTransform.translate(location.x+ddiameter2, location.y+ddiameter2);
        labelTransform.rotate(Math.toRadians(-THETA));
        labelTransform.translate(-(location.x+ddiameter2), -(location.y+ddiameter2));
        g2.setTransform(labelTransform);
        g2.drawString(label,
                    location.x+x,
                    location.y+y);
        g2.setTransform(origin);
      }
      
      g2.setTransform(newOrig);

      g2.setStroke(new BasicStroke(1.f));
      int xLine = location.x+(int)(ddiameter);
      int yLine = location.y+(int)(ddiameter/2.d);
      g2.drawLine(xLine,yLine,(int)(xLine+lineSize+widDash),yLine);

/*
      System.out.println("THETA "+Math.toDegrees(theta));
      System.out.println("m00 "+newOrig.getScaleX()+
                         " m01 "+newOrig.getShearX()+
                         " m02 "+newOrig.getTranslateX());
      System.out.println("m10 "+newOrig.getScaleY()+
                         " m12 "+newOrig.getTranslateY());
*/
      g2.setTransform(origin);
    }

    return;
  }


  /**
  *
  * Calculate the tick marks to be drawn
  *
  */
  protected void calculateTickPosistions() 
  {
    minorTicks = new Vector();
    majorTicks = new Vector();
    int start = getStart();
    int end   = getEnd();

    if(majorTick == 0)
      return;

    for(int i=startTick; i<end; i+=majorTick)
      if(i >= start)
        majorTicks.add(new Integer(i));

    if(minorTick == 0)
      return;

    for(int i=startTick; i<end; i+=minorTick)
    {
      Integer tk = new Integer(i);
      if(i >= start && !majorTicks.contains(tk))
        minorTicks.add(tk);
    }
  }


  /**
  *
  * Return the position tick marks start at
  *
  */
  protected int getStartTick()
  { 
    return startTick;
  }

  
  /**
  *
  * Set the position tick marks start at
  *
  */
  protected boolean setStartTick(int startTick)
  {
    this.startTick = startTick;
    if((startTick >= getStart()) && (startTick < getEnd()))
      return true;
    
    return false;
  }

  
  /**
  *
  * Return the interval for the tick marks
  *
  */
  protected int getTickInterval()
  {
    return majorTick;
  }


  /**
  *
  * Set the interval for the tick marks
  *
  */
  public boolean setTickInterval(int majorTick)
  {
    if(majorTick < (getEnd()-getStart()))
    {
      this.majorTick = majorTick;
      return true;
    }
    return false;
  }


  /**
  *
  * Return the interval for the tick marks
  *
  */
  protected int getMinorTickInterval()
  {
    return minorTick;
  }


  /**
  *
  * Set the interval for the tick marks
  *
  */
  public boolean setMinorTickInterval(int minorTick)
  {
    if(minorTick < (getEnd()-getStart()))
    {
      this.minorTick = minorTick;
      return true;
    }
    return false;
  }


  /**
  *
  * Return an angle in degrees
  *
  */
  protected double getAngleFromPosition(int pos,double rad)
  {
    int start = getStart();
    int end   = getEnd();
    return - ((pos-start)*rad)/(end-start);
  }


  /**
  *
  * The method @print@ must be implemented for @Printable@ interface.
  * Parameters are supplied by system.
  *
  */
  public int print(Graphics g, PageFormat pf, int pageIndex)
                                       throws PrinterException
  {
    Graphics2D g2 = (Graphics2D)g;
    g2.setColor(Color.black);    //set default foreground color to black

    //RepaintManager.currentManager(this).setDoubleBufferingEnabled(false);
    Dimension d = this.getSize();    //get size of document
    double panelWidth  = d.width;    //width in pixels
    double panelHeight = d.height;   //height in pixels
    double pageHeight = pf.getImageableHeight();   //height of printer page
    double pageWidth  = pf.getImageableWidth();    //width of printer page
    double scale = pageWidth/panelWidth;
    int totalNumPages = (int)Math.ceil(scale * panelHeight / pageHeight);

    // Make sure not print empty pages
    if(pageIndex >= totalNumPages)
    {
      System.out.println("NO SUCH PAGE "+pageIndex); 
      return Printable.NO_SUCH_PAGE;
    }
    
    // Shift Graphic to line up with beginning of print-imageable region
    g2.translate(pf.getImageableX(), pf.getImageableY());
    // Shift Graphic to line up with beginning of next page to print
    g2.translate(0f, -pageIndex*pageHeight);
    // Scale the page so the width fits...
    g2.scale(scale, scale);
    
    try
    {
      drawAll(g2,false);
    }
    catch(Exception e){ e.printStackTrace(); }
    return Printable.PAGE_EXISTS;
  }

  public void drawAll(Graphics2D g2, boolean l)
  {
    if(((Boolean)lineAttr.get("circular")).booleanValue())
      drawCircularPanel(g2,l);   //repaint the page for printing
    else
      drawLinearPanel(g2);

    if(gcGraph != null && containsGraph(gcGraph))
    {
      if(isCircular())
        gcGraph.draw(g2);
      else
        gcGraph.drawLinear(g2);
    }
    if(gcSkewGraph != null && containsGraph(gcSkewGraph))
    {
      if(isCircular())
        gcSkewGraph.draw(g2);
      else
        gcSkewGraph.drawLinear(g2);
    }
    if(userGraphs.size() > 0)
    {
      for(int i=0; i<userGraphs.size(); i++)
      {
        Graph userGraph = userGraphs.get(i);
        if(containsGraph(userGraph))
        {
          if(isCircular())
            userGraph.draw(g2);
          else
            userGraph.drawLinear(g2);
        }
      }
    }
  }

  /**
   * Check if a Graph object is contained by the panel.
   * @param g
   * @return
   */
  protected boolean containsGraph(Graph g)
  {
    int ncomponents = getComponentCount();
    for(int i=0; i<ncomponents; i++)
      if(getComponent(i).equals(g))
        return true;
    
    return false;
  }

  public void doPrintActions()
  {
    final PrinterJob pj=PrinterJob.getPrinterJob();
    pj.setPrintable(this);
    pj.printDialog();
    try
    {
      setCursor(new Cursor(Cursor.WAIT_CURSOR));
      pj.print();
      setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
    }
    catch (Exception pe) { pe.printStackTrace(); }
  }


  public void setRestrictionEnzyme(Vector restrictionEnzyme)
  {
    this.restrictionEnzyme = restrictionEnzyme;
  }

  public void setGeneticMarker(Vector block)
  {
    this.block = block;
  }

  protected Hashtable getLineAttributes()
  {
    return lineAttr;
  }

  public Hashtable getFeaturePoints()
  {
    Collections.sort(block, new BlockComparator());
    final Hashtable h = new Hashtable(block.size());
    int len = 0;
    double lastAngle = -10.d;
    
    for(int i=0; i<block.size(); i++)
    {
      Block b = (Block)block.get(i);
      double thisAngle = b.getMidAngle();
      
      if( Math.abs(thisAngle - lastAngle ) < 2 )
        len = len+10;
      else
        len = 0;
      
      h.put(b.getFeature().getIDString()+";"+b.getBstart()+".."+b.getBend(),
            b.getLinePoints( len ));
      
      lastAngle = thisAngle;
    }
    return h;
  }

  public void setLineAttributes(Hashtable lineAttr)
  {
    this.lineAttr = lineAttr;
    
    if(isCircular())
    {
      setSize(panelSize);
      setPreferredSize(panelSize);
    }
    
    linearPlotOptions.setEnabled(!isCircular());
    revalidate();
    repaint();
  }


  protected void setLineSize(int lineSize)
  {
    lineAttr.put("lsize",new Integer(lineSize));
  }


  protected int getLineSize()
  {
    return ((Integer)lineAttr.get("lsize")).intValue();
  }

  public void setPlasmidLocation(int x,int y)
  {
    location.setLocation(x,y);
  } 

  
  public JMenuBar createMenuBar()
  {
    JMenuBar menuBar = new JMenuBar();

// file menu
    JMenu fileMenu = new JMenu("File");
    fileMenu.setMnemonic(KeyEvent.VK_F);
    menuBar.add(fileMenu);

    final JMenuItem openMenu = new JMenuItem("Open");
    openMenu.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        Wizard.getDNADrawFromFile(DNADraw.this);
        majorTicks = null;

        if(gcGraph != null)
          gcGraph = new GCGraph(DNADraw.this);
        if(gcSkewGraph != null)
          gcSkewGraph = new GCSkewGraph(DNADraw.this);
        
        repaint();
        
        /*EmbossCirdnaReader dnaRead = new EmbossCirdnaReader();
        block = dnaRead.getBlock();
        restrictionEnzyme = dnaRead.getRestrictionEnzyme();

        lineAttr.put("start",new Integer(dnaRead.getStart()));
        lineAttr.put("end",new Integer(dnaRead.getEnd()));
   
        current_dna = new DNADraw(block,restrictionEnzyme,
                                  lineAttr,0,100,100);
        jsp.setViewportView(current_dna);*/
      }
    });
    fileMenu.add(openMenu);
    fileMenu.add(new JSeparator());
 
    final JMenuItem readInEntry = new JMenuItem("Read In Entry on Separate Track...");
    readInEntry.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        Wizard.readEntry(DNADraw.this, getBases());
        trackManager.refresh();
        repaint();
      }
    });
    fileMenu.add(readInEntry);
    
// print
    final JMenuItem printMenu = new JMenuItem("Print...");
    fileMenu.add(printMenu);

    printMenu.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        doPrintActions();
      }
    });

    JMenuItem printImage = new JMenuItem("Save As jpeg/png Image...");
    printImage.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        PrintDNAImage pdnai = new PrintDNAImage(current_dna);
        pdnai.printAsSinglePage();
      }
    });
    fileMenu.add(printImage);

// print preview
    JMenuItem printPreview = new JMenuItem("Print Preview...");
    printPreview.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        PrintDNAImage pdnai = new PrintDNAImage(current_dna);
        pdnai.printPreview();
      }
    });
    fileMenu.add(printPreview);
    fileMenu.add(new JSeparator());

    fileMenu.add(TrackManager.getExportTrackTemplateMenuItem(null, this));
    
    if(trackManager == null)
      trackManager = new TrackManager(DNADraw.this);
    fileMenu.add(trackManager.getImportTrackTemplateMenuItem(null));
    fileMenu.add(new JSeparator());
    
    final JMenuItem fileMenuExit;
    if(!close)
      fileMenuExit = new JMenuItem("Exit");
    else
      fileMenuExit = new JMenuItem("Close");
    fileMenuExit.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        if(!close)
          System.exit(0);
        else
        {
          mainFrame.setVisible(false);
          mainFrame.dispose();
        }
      }
    });
    fileMenu.add(fileMenuExit);

// view menu
    JMenu viewMenu = new JMenu("View");
    menuBar.add(viewMenu);
  
    JMenuItem zoomIn = new JMenuItem("Zoom In");
    zoomIn.setAccelerator(KeyStroke.getKeyStroke(
                    KeyEvent.VK_I, ActionEvent.CTRL_MASK));
    zoomIn.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        zoomIn();
      }
    });
    viewMenu.add(zoomIn);

    JMenuItem zoomOut = new JMenuItem("Zoom Out");
    zoomOut.setAccelerator(KeyStroke.getKeyStroke(
                    KeyEvent.VK_O, ActionEvent.CTRL_MASK));
    zoomOut.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        zoomOut();
      }
    });
    viewMenu.add(zoomOut);

    
    JMenu graph_menu = new JMenu("Graph");
    graph_menu.setMnemonic(KeyEvent.VK_G);
    menuBar.add(graph_menu);
    
    JMenu gc = new JMenu("GC plot");
    graph_menu.add(gc);
    
    final JCheckBoxMenuItem gcPlot = new JCheckBoxMenuItem("Draw");
    
    if(gcGraph == null)
      gcPlot.setState(false);
    else
      gcPlot.setState(true);
    gcPlot.addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent event) 
      {
        if(gcPlot.isSelected())
        {
          if(gcGraph == null)
          {
            setCursor(new Cursor(Cursor.WAIT_CURSOR));
            gcGraph = new GCGraph(DNADraw.this);
            gcGraph.setTrack(0.4);
            setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
          }
          add(gcGraph);
          revalidate();
        }
        else
          remove(gcGraph);
        repaint();
      }
    });
    gc.add(gcPlot);
    
    JMenuItem gcOptions = new JMenuItem("Options...");
    gc.add(gcOptions);
    gcOptions.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        gcGraph.showOptions();
      } 
    });
    
    JMenu gcSkew = new JMenu("GC Skew");
    graph_menu.add(gcSkew);
    
    final JCheckBoxMenuItem gcSkewPlot = new JCheckBoxMenuItem("Draw");
    if(gcSkewGraph == null)
      gcSkewPlot.setState(false);
    else
      gcSkewPlot.setState(true);

    gcSkewPlot.addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent event) 
      {
        if(gcSkewPlot.isSelected())
        {
          if(gcSkewGraph == null)
          {
            setCursor(new Cursor(Cursor.WAIT_CURSOR));
            gcSkewGraph = new GCSkewGraph(DNADraw.this);
            gcSkewGraph.setTrack(0.2);
            setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
          }
          add(gcSkewGraph);
          revalidate();
        }
        else
          remove(gcSkewGraph);
        repaint();
      }
    });
    gcSkew.add(gcSkewPlot);

    JMenuItem gcSkewOptions = new JMenuItem("Options...");
    gcSkew.add(gcSkewOptions);
    gcSkewOptions.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        gcSkewGraph.showOptions();
      } 
    });

    
    //
    // User graph
    
    JMenu user = new JMenu("User");
    graph_menu.add(user);
    
    final JMenuItem userPlot = new JMenuItem("Add");
    userPlot.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent ae) 
      {
        final uk.ac.sanger.artemis.components.StickyFileChooser dialog =
            new uk.ac.sanger.artemis.components.StickyFileChooser ();

        dialog.setDialogTitle ("Select a data file name ...");
        dialog.setDialogType (JFileChooser.OPEN_DIALOG);

        final int status = dialog.showOpenDialog (null);

        if(status != JFileChooser.APPROVE_OPTION ||
           dialog.getSelectedFile () == null) 
        {
          return;
        }

        final java.io.File file =
            new java.io.File (dialog.getCurrentDirectory (),
                      dialog.getSelectedFile ().getName ());

        if(file.length() == 0)
          return;
          
        //final uk.ac.sanger.artemis.util.Document document =
        //    new uk.ac.sanger.artemis.util.FileDocument (file);
          
        UserGraph userGraph;
        try
        {
          userGraph = 
            new UserGraph(DNADraw.this, file.getAbsolutePath());
        }
        catch(IOException e)
        {
          e.printStackTrace();
          return;
        }
          
        if(userGraph == null)
        {
          setCursor(new Cursor(Cursor.WAIT_CURSOR));
          try
          {
            userGraph = new UserGraph(DNADraw.this, file.getAbsolutePath());
          }
          catch(IOException e)
          {
            e.printStackTrace();
          }
          userGraph.setTrack(0.2);
          setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
        }

        setUserGraph(userGraph);
        add(userGraph);
        
        revalidate();
        repaint();
      }
    });
    user.add(userPlot);
    user.add(userDraw);
    user.add(userOptions);
    
    
// options menu
    JMenu optionMenu = new JMenu("Options");
    menuBar.add(optionMenu);


    JMenuItem wizard = new JMenuItem("DNA Wizard...");
    wizard.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {   
        new Wizard(current_dna);
      }
    });
    optionMenu.add(wizard);
    optionMenu.add(new JSeparator());

    JMenuItem tracksMenu = new JMenuItem("Track Manager...");

    if(getArtemisEntryGroup() != null)
    {
      if(trackManager == null)
        trackManager = new TrackManager(DNADraw.this);
      tracksMenu.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          trackManager.setVisible(true);
        }
      });
    }
    else
      tracksMenu.setEnabled(false);
    
    optionMenu.add(tracksMenu);
    
    
    /*JMenuItem line = new JMenuItem("DNA Properties...");
    line.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        JFrame f = new JFrame("DNA Properties");
        LineAttribute la = new LineAttribute(current_dna);
        JScrollPane laScroll = new JScrollPane(la);
        JPanel laPane = (JPanel)f.getContentPane();
        laPane.add(laScroll,BorderLayout.CENTER);
        f.setJMenuBar(la.createMenuBar(f));
        f.pack();
        f.setVisible(true);
      }
    });
    optionMenu.add(line);*/

    
    JMenuItem tickMarks = new JMenuItem("Tick marks...");
    tickMarks.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        JFrame f = new JFrame("Tick Marks");
        Ticks tk = new Ticks(current_dna,true);
        JScrollPane tkScroll = new JScrollPane(tk);
        JPanel tkPane = (JPanel)f.getContentPane();
        tkPane.add(tkScroll,BorderLayout.CENTER);
        f.setJMenuBar(tk.createMenuBar(f));
        f.pack();
        f.setVisible(true);
      }
    });
    optionMenu.add(tickMarks);


    JMenuItem gmarker = new JMenuItem("Features...");
    gmarker.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        JFrame f = new JFrame("Features");
        GeneticMarker gm = new GeneticMarker(current_dna,
                                             block);
        JScrollPane gmScroll = new JScrollPane(gm);
        JPanel gmPane = (JPanel)f.getContentPane();
        gmPane.add(gmScroll,BorderLayout.CENTER);
        f.setJMenuBar(gm.createMenuBar(f));
        f.pack();
        f.setVisible(true);
      }
    });
    optionMenu.add(gmarker);

    JMenuItem reSites = new JMenuItem("Restriction Enzyme...");
    reSites.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        JFrame f = new JFrame("Restriction Enzyme");
        RestrictionEnzyme re = new RestrictionEnzyme(current_dna,
                                              restrictionEnzyme);
        JScrollPane reScroll = new JScrollPane(re);
        JPanel rePane = (JPanel)f.getContentPane();
        rePane.add(reScroll,BorderLayout.CENTER);
        f.setJMenuBar(re.createMenuBar(f));
        f.pack();
        f.setVisible(true);
      }
    });
    optionMenu.add(reSites);
    
    linearPlotOptions.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        showLinearPlotOptions();
      }
    });
    linearPlotOptions.setEnabled(!isCircular());
    optionMenu.add(linearPlotOptions);  
    optionMenu.addSeparator();
    
    final JCheckBoxMenuItem labelTick = new JCheckBoxMenuItem("Show tick numbers", labelTicks);
    labelTick.addItemListener(new ItemListener()
    {

      public void itemStateChanged(ItemEvent e)
      {
        labelTicks = labelTick.isSelected();
        repaint();
      }
    });
    optionMenu.add(labelTick);  

    final JCheckBoxMenuItem asCircular = new JCheckBoxMenuItem("Display as Circular Plot", isCircular());
    asCircular.addItemListener(new ItemListener()
    {

      public void itemStateChanged(ItemEvent e)
      {
        lineAttr.put("circular",new Boolean(asCircular.isSelected()));
        setLineAttributes(lineAttr);
        linearPlotOptions.setEnabled(!isCircular());
      }
    });
    optionMenu.add(asCircular);

// help manu
    JMenu helpMenu = new JMenu("Help");
    menuBar.add(helpMenu);
   
    JMenuItem aboutMenu = new JMenuItem("About");
    helpMenu.add(aboutMenu);
    aboutMenu.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        /*ClassLoader cl = DNADraw.this.getClass().getClassLoader();
        URL url = cl.getResource("etc/readmeDNADraw.html");*/
        BrowserControl.displayURL("http://www.sanger.ac.uk/Software/Artemis/circular/");
      }
    });

    return menuBar;
  }

  /**
   * Provide options for the linear plot
   */
  private void showLinearPlotOptions()
  {
    GridBagLayout grid = new GridBagLayout();
    GridBagConstraints c = new GridBagConstraints();
    c.ipady = 3;
    c.ipadx = 5;

    JPanel optionBox = new JPanel(grid);
    c.gridx = 1; 
    c.gridy = 0;
    c.anchor = GridBagConstraints.WEST;
    TextFieldFloat lineHeightField = new TextFieldFloat();
    lineHeightField.setValue(getLineHeight());
    optionBox.add(lineHeightField, c);
    c.gridx = 0;
    c.anchor = GridBagConstraints.EAST;
    optionBox.add(new JLabel("Line Height"), c);
    
    c.gridx = 1; 
    c.gridy = 1;
    c.anchor = GridBagConstraints.WEST;
    TextFieldInt basesPerLineField = new TextFieldInt();
    basesPerLineField.setValue(getBasesPerLine());
    optionBox.add(basesPerLineField, c);
    c.gridx = 0;
    c.anchor = GridBagConstraints.EAST;
    optionBox.add(new JLabel("Bases Per Line"), c);
    
    JOptionPane.showMessageDialog(null, optionBox,
        "Linear Plot Options", JOptionPane.QUESTION_MESSAGE);
    
    setLineHeight( (float) lineHeightField.getValue() );
    setBasesPerLine(basesPerLineField.getValue());
    revalidate();
    repaint();
  }
  
  public void setCloseAndDispose(boolean close, JFrame mainFrame)
  {
    this.mainFrame = mainFrame;
    this.close = close;
  }

  protected Vector getGeneticMarker()
  {
    return block;
  }


  protected Vector getRestrictionEnzyme()
  {
    return restrictionEnzyme;
  }


  protected int getStart()
  {
    return ((Integer)lineAttr.get("start")).intValue();
  }

 
  protected int getEnd()
  {
    return ((Integer)lineAttr.get("end")).intValue();
  }

  
  protected void setStart(int start)
  {
    lineAttr.put("start",new Integer(start));
    calculateTickPosistions();
  }


  protected void setEnd(int end)
  {
    lineAttr.put("end",new Integer(end));
    calculateTickPosistions();
  }

  public Block getBlockAtLocation(Point loc)
  {
    for(int i=0; i<block.size(); i++)
    {
      Block b = (Block) block.get(i);
      if(b.isOverMe(loc.x, loc.y))
        return b;
    }
    return null;
  }
  
  public Block getBlockAtBasePosition(int bend)
  {
    for(int i=0; i<block.size(); i++)
    {
      Block b = (Block) block.get(i);
      if(b.getBend() == bend)
        return b;
    }
    return null;
  }
  
////////////////////
// DRAG AND DROP
////////////////////
  public RenderedImage createImage(Block b)
  {
    // Create a buffered image in which to draw
    BufferedImage bufferedImage = new BufferedImage(getWidth(), getHeight(),
                                      BufferedImage.TYPE_INT_ARGB);
  
    // Create a graphics contents on the buffered image
    Graphics2D g2d = bufferedImage.createGraphics();
    RenderingHints qualityHints = new RenderingHints(RenderingHints.KEY_ANTIALIASING,
        RenderingHints.VALUE_ANTIALIAS_ON); 
    qualityHints.put(RenderingHints.KEY_RENDERING,               
    RenderingHints.VALUE_RENDER_QUALITY); 
    g2d.setRenderingHints(qualityHints);
    g2d.setTransform( ((Graphics2D)getGraphics()).getTransform() );
    // Draw graphics
    //b.drawCircular(g2d);
    
    g2d.setStroke(new BasicStroke(2.f));
    g2d.setColor(b.getColour());
    g2d.drawArc(0, 0, 10, 10, 0, 360);
    g2d.dispose();
  
    return bufferedImage;
  }
  
// drag source
  public void dragGestureRecognized(DragGestureEvent e)
  {
    // ignore if mouse popup trigger
    InputEvent ie = e.getTriggerEvent();
    if(ie instanceof MouseEvent)
      if(((MouseEvent)ie).isPopupTrigger())
        return;

    Point loc = e.getDragOrigin();
    Block b = getBlockAtLocation(loc);
    
    if(b != null)
    {
      Image image = (Image)createImage(b);
      e.startDrag(DragSource.DefaultCopyDrop,  // cursor
          image, new Point(-1, -1),
          (Transferable)b,             // transferable data
          this);
    }
  }

  public void dragDropEnd(DragSourceDropEvent e) {}
  public void dragEnter(DragSourceDragEvent e) {}
  public void dragExit(DragSourceEvent e) {}
  public void dragOver(DragSourceDragEvent e) {}
  public void dropActionChanged(DragSourceDragEvent e) {}

// drop sink
  public void dragEnter(DropTargetDragEvent e)
  {
    if(e.isDataFlavorSupported(Block.BLOCK))
      e.acceptDrag(DnDConstants.ACTION_COPY_OR_MOVE);
  }

  public void drop(DropTargetDropEvent e)
  {
    Transferable t = e.getTransferable();
    if(t.isDataFlavorSupported(Block.BLOCK))
    {
      try
      {
        Point loc = e.getLocation();
        Block b   = (Block)t.getTransferData(Block.BLOCK);
        b.setBlockLocation(loc.x,loc.y,trackManager);
        DNADraw.this.repaint();
      }
      catch(Exception ufe){} 
    }
  }

  public void dragOver(DropTargetDragEvent e)
  {
  }

  public void dropActionChanged(DropTargetDragEvent e) {}
  public void dragExit(DropTargetEvent e){}
  
  public EntryGroup getArtemisEntryGroup()
  {
    return artemisEntryGroup;
  }

  public void setArtemisEntryGroup(EntryGroup artemisEntryGroup)
  {
    this.artemisEntryGroup = artemisEntryGroup;
  }
  
  public Bases getBases()
  {
    if(getArtemisEntryGroup() != null)
      return getArtemisEntryGroup().getSequenceEntry().getBases();
    return this.bases;
  }
  
  public void setBases(Bases bases)
  {
    this.bases = bases;
  }
  
  protected Vector getBlock()
  {
    return block;
  }


  public int getBasesPerLine()
  {
    return basesPerLine;
  }


  public void setBasesPerLine(int basesPerLine)
  {
    this.basesPerLine = basesPerLine;
  }


  public int getNumberOfLines()
  {
    return numberOfLines;
  }


  public void setNumberOfLines(int numberOfLines)
  {
    this.numberOfLines = numberOfLines;
  }


  public float getLineHeight()
  {
    return lineHeight;
  }


  public void setLineHeight(float lineHeight)
  {
    this.lineHeight = lineHeight;
  }


  public float getSingleBaseWidth()
  {
    return singleBaseWidth;
  }


  public void setSingleBaseWidth(float singleBaseWidth)
  {
    this.singleBaseWidth = singleBaseWidth;
  }


  public int getBorderWidth2()
  {
    return borderWidth2;
  }


  public void setBorderWidth2(int borderWidth2)
  {
    this.borderWidth2 = borderWidth2;
  }


  public int getBorderHeight2()
  {
    return borderHeight2;
  }


  public void setBorderHeight2(int borderHeight2)
  {
    this.borderHeight2 = borderHeight2;
  }

  public TrackManager getTrackManager()
  {
    return trackManager;
  }

  public void setTrackManager(TrackManager trackManager)
  {
    this.trackManager = trackManager;
  }

  public Graph getGcGraph()
  {
    return gcGraph;
  }
  
  public Vector<Graph> getUserGraphs()
  {
    return userGraphs;
  }

  public void setGcGraph(Graph gcGraph)
  {
    this.gcGraph = gcGraph;
  }

  public Graph getGcSkewGraph()
  {
    return gcSkewGraph;
  }

  public void setGcSkewGraph(Graph gcSkewGraph)
  {
    this.gcSkewGraph = gcSkewGraph;
  }


  protected void setUserGraph(final UserGraph thisUserGraph)
  {
    userGraphs.add(thisUserGraph);

    JMenuItem thisGraphOptions = new JMenuItem(thisUserGraph.getFileName());
    thisGraphOptions.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        thisUserGraph.showOptions();
      } 
    });
    userOptions.add(thisGraphOptions);
    
    final JCheckBoxMenuItem thisGraphDraw = 
      new JCheckBoxMenuItem(thisUserGraph.getFileName());
    thisGraphDraw.setSelected(true);
    thisGraphDraw.addItemListener(new ItemListener()
    {
      public void itemStateChanged(ItemEvent e)
      {
        if(thisGraphDraw.isSelected())
          thisUserGraph.setVisible(true);
        else
          thisUserGraph.setVisible(false);
      }
    });
    userDraw.add(thisGraphDraw);
  }
  
  public static void main(final String arg[])
  {
    SwingUtilities.invokeLater(new Runnable() 
    {
      public void run() 
      {
        Wizard wiz = null;
        
        if(arg.length > 0 && arg[0].equals("-t"))
        {
          try 
        	  {
        	    wiz = new Wizard(arg[1]);
        	  } 
        	  catch (Exception e) 
        	  {
            System.err.println("\nUnable to load template file: " + arg[1]);
        	    System.exit(1);
          }
        }
        else
        {
          wiz = new Wizard((DNADraw)null);
        }
        
        final DNADraw dna = wiz.getDNADraw();
        
        //
        final String version = dna.getVersion();
        final JFrame f = new JFrame();
        if(version == null)
          f.setTitle("DNAPlotter");
        else
          f.setTitle("DNAPlotter :: "+version);
          
        Dimension d = f.getToolkit().getScreenSize();

        jsp.setViewportView(dna);
        jsp.getViewport().setBackground(Color.white);
        f.getContentPane().add(jsp);
        f.setJMenuBar(dna.createMenuBar());
        
        //dna.add(new Graph(dna));
        f.pack();
        f.setLocation(((int)d.getWidth()-f.getWidth())/4,
                      ((int)d.getHeight()-f.getHeight())/2);

        
        //dna.add(dna.new BlockPanel(dna));
        f.setVisible(true);
        
        if(wiz.getWorkerGraph() != null)
          wiz.getWorkerGraph().start();
      }
    });
  }

}

