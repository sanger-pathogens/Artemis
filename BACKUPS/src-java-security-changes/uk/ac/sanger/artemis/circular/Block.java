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

import uk.ac.sanger.artemis.Feature;

import javax.swing.*;

import java.awt.geom.AffineTransform;
import java.awt.geom.Line2D;
import java.awt.*;
import java.awt.datatransfer.*;
import java.awt.event.*;
import java.io.IOException;
import java.util.Vector;

import javax.swing.event.*;
import java.awt.geom.Arc2D;

public class Block implements Transferable
{
  private static final long serialVersionUID = 1L;
  private DNADraw current_dna;
  private double angStart;
  private double angEnd;

  private Vector rect;
  private boolean drawLabel = false;
  final public static DataFlavor BLOCK =
         new DataFlavor(Block.class, "Block");
  static DataFlavor blockFlavors[] = { BLOCK };
  
  private String label;
  private int bstart;
  private int bend;
  private Feature feature;
  private Color colour;
  private float strokeSize;
  private boolean arrowHead;
  private boolean arrowTail;
  
  /** track to place the block on - as a fraction of the radii*/
  private Track track;
  
  
  public Block(final String label, 
               int bstart,
               int bend, 
               Color colour, 
               float strokeSize, 
               Track track)
  {
    super();
    this.label = label;
    this.bstart = bstart;
    this.bend = bend;
    this.colour = colour;
    this.strokeSize = strokeSize;
    this.track = track;
  }


  public Block(final String label, 
               int bstart,
               int bend, 
               Color colour, 
               float strokeSize, 
               Track fracRadii,
               DNADraw current_dna)
  {
    this(label, bstart, bend, colour, strokeSize, fracRadii);
    this.current_dna = current_dna;
  }


  protected void draw(Graphics2D g2)
  {
    if(current_dna.isCircular())
      drawCircular(g2);
    else
      drawLinear(g2);   
  }


  protected void showProperties(final JFrame f, final DNADraw draw,
                                JButton delete)
  {
//  final JFrame f = new JFrame("Properties");
//set up menu bar
    JMenuBar menuBar = new JMenuBar();
    JMenu fileMenu = new JMenu("File");
    fileMenu.setMnemonic(KeyEvent.VK_F);
    menuBar.add(fileMenu);
                                                                             
    JMenuItem closeMenu = new JMenuItem("Close");
    closeMenu.setAccelerator(KeyStroke.getKeyStroke(
              KeyEvent.VK_E, ActionEvent.CTRL_MASK));
                                                                             
    closeMenu.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        f.dispose();
      }
    });
    fileMenu.add(closeMenu);

//set up window
    f.setJMenuBar(menuBar);
    JPanel pane = (JPanel)f.getContentPane();
    String markerLabel = getLabel();
    int bstart = getBstart();
    int bend   = getBend();
    Color colour = getColour();
    float strokeSize  = getStrokeSize();
    boolean arrowHead = isArrowHead();
    boolean arrowTail = isArrowTail();
         
    Box bdown = Box.createVerticalBox();
    pane.add(bdown);
     
    bdown.add(Box.createVerticalStrut(4));
    Dimension d = new Dimension(200,30); 
    Box bacross = Box.createHorizontalBox(); 
    final JTextField annotation = new JTextField();
    annotation.setPreferredSize(d);
    annotation.setMaximumSize(d);
    annotation.setText(markerLabel);
    annotation.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        setLabel(annotation.getText());
        if(draw != null)
          draw.repaint();
      }
    });
    bacross.add(annotation);
    bacross.add(new JLabel(" Label"));
    bacross.add(Box.createHorizontalGlue());
    bdown.add(bacross);

    d = new Dimension(65,30);
    bacross = Box.createHorizontalBox();  
    final TextFieldInt start = new TextFieldInt();
    start.setPreferredSize(d);
    start.setMaximumSize(d);
    start.setValue(bstart);
    start.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        setBstart(start.getValue()); 
        if(draw != null)
          draw.repaint();
      }
    });
    bacross.add(start); 
    bacross.add(new JLabel(" Start"));
    bacross.add(Box.createHorizontalGlue());
    bdown.add(bacross);

    bacross = Box.createHorizontalBox();
    final TextFieldInt end = new TextFieldInt();
    end.setPreferredSize(d);
    end.setMaximumSize(d);
    end.setValue(bend);
    end.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        setBend(end.getValue());
        if(draw != null)
          draw.repaint();
      }
    });
    bacross.add(end);
    bacross.add(new JLabel(" End"));
    bacross.add(Box.createHorizontalGlue());
    bdown.add(bacross);

//Set up the dialog that the button brings up.
    final JButton colourButton = GeneticMarker.setUpColorButton(colour);
    final JButton applyButton = new JButton("Apply");

    applyButton.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        setColour(colourButton.getBackground());
        setLabel(annotation.getText());
        if(draw != null)
          draw.repaint();
      }
    });

    bacross = Box.createHorizontalBox();
    bacross.add(new JLabel("Pick a Colour: "));
    bacross.add(colourButton);
    bacross.add(applyButton);
    
    bacross.add(Box.createHorizontalGlue());
    bdown.add(bacross);

    bacross = Box.createHorizontalBox();    
    final JSlider slider;
    
    if(strokeSize <= 25)
      slider = new JSlider(1,25,(int)strokeSize);
    else
      slider = new JSlider(1,25+(int)strokeSize,(int)strokeSize);
    
    bacross.add(slider);
    bacross.add(new JLabel(" Line width"));
    bacross.add(Box.createHorizontalGlue());
    slider.addChangeListener(new ChangeListener()
    {
      public void stateChanged(ChangeEvent e)
      {
        setStrokeSize((float)slider.getValue());
        if(draw != null)
          draw.repaint();
      }
    });
    bacross.add(Box.createHorizontalGlue());
    bdown.add(bacross);

    bacross = Box.createHorizontalBox();
    bacross.add(new JLabel("Arrow :"));
    bacross.add(Box.createHorizontalGlue());
    bdown.add(bacross);

    bacross = Box.createHorizontalBox();
    final JRadioButton head = new JRadioButton("Head");
    final JRadioButton tail = new JRadioButton("Tail");
    final JRadioButton none = new JRadioButton("None");
    head.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        if(head.isSelected())
        {
          setArrowHead(true);   
          setArrowTail(false);
        }
        else
        {
          setArrowHead(false);
          setArrowTail(true);
        }
        if(draw != null)
          draw.repaint();
      }
    });
    bacross.add(head);

    tail.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        if(tail.isSelected())
        {
          setArrowTail(true);
          setArrowHead(false);
        }
        else
        {
          setArrowTail(false);
          setArrowHead(true);
        }
        if(draw != null)
          draw.repaint();
      }
    });
    bacross.add(tail);

    none.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        if(none.isSelected())
        {
          setArrowTail(false);
          setArrowHead(false);
        }
        if(draw != null)
          draw.repaint();
      }
    });
    bacross.add(none);
    bacross.add(Box.createHorizontalGlue());
    bdown.add(bacross);

    bacross = Box.createHorizontalBox();
    final JCheckBox label = new JCheckBox("Show Label", drawLabel);
    label.addItemListener(new ItemListener()
    {
      public void itemStateChanged(ItemEvent e)
      {
        drawLabel = label.isSelected();
        setLabel(annotation.getText());
        draw.repaint();
      }  
    });
    bacross.add(label);
    bacross.add(Box.createHorizontalGlue());
    bdown.add(bacross);
    
    
    bacross = Box.createHorizontalBox();
    bacross.add(delete);
    bacross.add(Box.createHorizontalGlue());
    bdown.add(bacross);

    ButtonGroup group = new ButtonGroup();
    group.add(head);
    group.add(tail);
    group.add(none);
    if(arrowHead)
      head.setSelected(true);
    else if(arrowTail)
      tail.setSelected(true);

    f.pack();
    f.setVisible(true);
  }

  protected void drawLinear(final Graphics2D g2)
  {
    RenderingHints qualityHints = new RenderingHints(RenderingHints.KEY_ANTIALIASING,
        RenderingHints.VALUE_ANTIALIAS_ON); 
    qualityHints.put(RenderingHints.KEY_RENDERING,               
    RenderingHints.VALUE_RENDER_QUALITY); 
    g2.setRenderingHints(qualityHints);

    String markerLabel = getLabel();
    int bstart = getBstart();
    int bend   = getBend();
    
   
    FontMetrics fm = g2.getFontMetrics();
    g2.setColor(getColour());

    //int shift = (int)getTrack().getPosition();

    int basesPerLine    = current_dna.getBasesPerLine();
    int lineNumberStart = Math.round((float)(bstart-1)/((float)basesPerLine) + 0.5f)-1;
    int lineNumberEnd   = Math.round((float)(bend-1)/((float)basesPerLine) + 0.5f)-1;
    float singleBaseWidth = current_dna.getSingleBaseWidth();
     
    int borderWidth2 = current_dna.getBorderWidth2();
    int borderHeight2 = current_dna.getBorderHeight2();
    
    int ypos = (int) ( (lineNumberStart*current_dna.getLineHeight()) - (strokeSize/2.f) -
                        ((1-track.getPosition())*current_dna.getLineHeight())  +
                        borderHeight2+current_dna.getLineHeight());

    int xstart = (int) ((bstart-(lineNumberStart*basesPerLine))*singleBaseWidth)+borderWidth2;
    int xend   = (int) ((bend-(lineNumberEnd*basesPerLine))*singleBaseWidth)+borderWidth2;
    
    if(rect == null)
      rect = new Vector();
    
    BasicStroke basicstroke;
    if(bstart == bend)
    {
      // SNP / variation feature
      basicstroke = new BasicStroke(
          1.f,
          BasicStroke.CAP_BUTT, 
          BasicStroke.JOIN_MITER);
      g2.setStroke(basicstroke);
      g2.drawLine(xstart,ypos,xstart,(int) (ypos+strokeSize));
      Rectangle r = new Rectangle();
      r.setLocation(xstart,ypos);
      r.setSize(current_dna.getWidth()-borderWidth2-xstart,2);
      rect.add(r);
      return;
    }
    
    basicstroke = new BasicStroke(
        strokeSize,
        BasicStroke.CAP_BUTT, 
        BasicStroke.JOIN_MITER);
    g2.setStroke(basicstroke);

    if(arrowTail)
    {
      int[] xPoints = {xstart+(int)strokeSize,xstart+(int)strokeSize,xstart};
      int[] yPoints = {ypos+(int)strokeSize,ypos-(int)strokeSize,ypos};
      xstart+=strokeSize/2;
      g2.fillPolygon(xPoints,yPoints,3);
    }


    Rectangle r;
    
    if(lineNumberEnd > lineNumberStart)
    {
      g2.drawLine(xstart,ypos,current_dna.getWidth()-borderWidth2,ypos);
      r = new Rectangle();
      r.setLocation(xstart,ypos);
      r.setSize(current_dna.getWidth()-borderWidth2-xstart,(int)strokeSize);
      rect.add(r);
      
      for(int i=lineNumberStart+1; i<lineNumberEnd; i++)
      {
        ypos = (int) ( (i*current_dna.getLineHeight()) - (strokeSize/2.f) -
            ((1-track.getPosition())*current_dna.getLineHeight())  +
            borderHeight2+current_dna.getLineHeight());
        g2.drawLine(borderWidth2,ypos,current_dna.getWidth()-borderWidth2,ypos);
        
        r = new Rectangle();
        r.setLocation(borderWidth2,ypos);
        r.setSize(current_dna.getWidth()-borderWidth2-borderWidth2,(int)strokeSize);
        rect.add(r);
      }
      
      ypos = (int) ( (lineNumberEnd*current_dna.getLineHeight()) - (strokeSize/2.f) -
          ((1-track.getPosition())*current_dna.getLineHeight())  +
          borderHeight2+current_dna.getLineHeight());
      xstart = borderWidth2;
    }
    
    if(arrowHead)
    {
      int[] xPoints = {xend-(int)strokeSize,xend-(int)strokeSize,xend};
      int[] yPoints = {ypos+(int)strokeSize,ypos-(int)strokeSize,ypos};
      xend-=strokeSize/2;
      g2.fillPolygon(xPoints,yPoints,3);
    }
    
    g2.drawLine(xstart,ypos,xend,ypos);
    
    r = new Rectangle();
    r.setLocation(xstart,ypos);
    r.setSize(xend-xstart,(int)strokeSize);
    rect.add(r);
    
    
    if(drawLabel)
    {
      int xmid = xstart+(int)(((xend-xstart)/2.)-(fm.stringWidth(markerLabel)/2.));
      g2.drawString(markerLabel,xmid,ypos-strokeSize);
    } 
  }


  protected void drawCircular(Graphics2D g2)
  {
    RenderingHints qualityHints = new RenderingHints(RenderingHints.KEY_ANTIALIASING,
                                                  RenderingHints.VALUE_ANTIALIAS_ON); 
    qualityHints.put(RenderingHints.KEY_RENDERING,               
                     RenderingHints.VALUE_RENDER_QUALITY); 
    g2.setRenderingHints(qualityHints);

    
    int bstart = getBstart();
    int bend   = getBend();
    Color colour = getColour();

    FontMetrics fm = g2.getFontMetrics();
    double hgt = fm.getAscent();
    g2.setColor(colour);
  
    final double ddiameter   = current_dna.getDiameter();
    final double widthPanel  = current_dna.getWidth();
    final double heightPanel = current_dna.getHeight();
    final double dradii      = ddiameter/2.d;
    final double rad = 360.d;
    Point location = current_dna.getLocationPoint();

    AffineTransform origin  = g2.getTransform();
    AffineTransform newOrig = (AffineTransform)(origin.clone());

    angStart = current_dna.getAngleFromPosition(bstart,rad) ;
    angEnd   = current_dna.getAngleFromPosition(bend,rad)  - angStart;

    double shift = dradii*(1.d-getTrack().getPosition());
    double bdiameter = ddiameter*getTrack().getPosition();
 
    // 
    float stroke  = getStrokeSize();
    float stroke2 = getStrokeSize()/2.f;
    int strokeInt = (int)stroke;
    
    double d2 = (bdiameter/2.d);
    double d  = d2*2.d;
    
    double locx = location.getX()+shift;
    double locy = location.getY()+shift;
    
    if(arrowHead)
    {
      newOrig.rotate(Math.toRadians(-angStart-angEnd),
          locx+d2,locy+d2);
      
      int xmid = (int)(locx+d);
      int ymid = location.y+(int)(dradii);
      int[] xPoints = {xmid-strokeInt,xmid+strokeInt,xmid};
      int[] yPoints = {ymid-strokeInt,ymid-strokeInt,ymid};
      g2.setTransform(newOrig);
      g2.fillPolygon(xPoints,yPoints,3);
      angEnd +=  Math.toDegrees(Math.atan(((float)(strokeInt)*0.8f)/(float)d2));
    }
    else if(arrowTail)
    {
      newOrig.rotate(Math.toRadians(-angStart),
          locx+d2,locy+d2);

      int xmid = (int)(locx+d);
      int ymid = location.y+(int)(dradii);
      int[] xPoints = {xmid-strokeInt,xmid+strokeInt,xmid};
      int[] yPoints = {ymid+strokeInt,ymid+strokeInt,ymid};
      g2.setTransform(newOrig);
      g2.fillPolygon(xPoints,yPoints,3);
      
      double adj = Math.toDegrees(Math.atan(((float)strokeInt*0.8f)/(float)d2));
      angStart -=  adj;
      angEnd += adj;
    }
 

    g2.setTransform(origin);
    // double arcLength = Math.toRadians(angEnd)*(bdiameter/2.);
    
    // if too small draw arc as a line
    if(Math.abs(angEnd) < 0.4d)
    {
      newOrig = (AffineTransform) (origin.clone());
      
      double xLine = d2-stroke2 ; 
      newOrig.translate(locx+d2,locy+d2);
      newOrig.rotate(Math.toRadians(-angStart));

      BasicStroke basicstroke = new BasicStroke(
          1,
          BasicStroke.CAP_BUTT, 
          BasicStroke.JOIN_MITER);

      g2.setStroke(basicstroke);
      g2.setTransform(newOrig);
      //g2.drawLine((int)xLine, 0, (int)(xLine+strokeInt), 0);
      g2.draw(new Line2D.Double(xLine, 0.d, xLine+stroke, 0.d));
    }
    else
    {
      BasicStroke basicstroke = new BasicStroke(
          stroke,
          BasicStroke.CAP_BUTT, 
          BasicStroke.JOIN_MITER);

      g2.setStroke(basicstroke);
      /*g2.drawArc(locx,locy, d, d,
               Math.round(Math.round(angStart)),
               Math.round(Math.round(angEnd)));*/
      
      g2.draw(new Arc2D.Double(locx,locy, d, d, angStart, angEnd,
          Arc2D.OPEN));
    }
    
    g2.setTransform(origin);

    if(drawLabel)
    {
      newOrig = (AffineTransform)(origin.clone());
      newOrig.rotate(Math.toRadians(-angStart-(angEnd/2.d)),
                   widthPanel/2.d,heightPanel/2.d);
   
      int xblock = 0;
      int yblock = 0;

      int widMarker = fm.stringWidth(getLabel())/2;
      xblock = (int)( dradii - (shift/2.d) + (newOrig.getScaleX()*
                     (dradii - 3 - widMarker - ((shift+getStrokeSize())/2.d))) -
                      widMarker );

      xblock = (int)( dradii + (newOrig.getScaleX()*
                     (dradii - shift - 3 - widMarker - (getStrokeSize()/2.d))) -
                      widMarker );

      yblock = (int)( dradii + (newOrig.getShearY()*
                     (dradii - shift - 3 - ((getStrokeSize() + hgt)/2.d))) +
                      hgt/2.d );
    
      newOrig = (AffineTransform)(origin.clone());
      
      
      AffineTransform labelTransform = (AffineTransform)(origin.clone());
      labelTransform.translate(location.x+dradii, location.y+dradii);
      labelTransform.rotate(Math.toRadians(-DNADraw.THETA));
      labelTransform.translate(-(location.x+dradii), -(location.y+dradii));
      g2.setTransform(labelTransform);
      
      g2.drawString(getLabel(),
                  location.x+xblock,
                  location.y+yblock);
      g2.setTransform(origin);
    }
  }


  /**
  * Routine to set the position of the block
  * @param x    x position 
  * @param y    y position
  */
  public void setBlockLocation(int x, int y, final TrackManager viewer)
  {
    if(current_dna.isCircular())
    {
      double dradii = current_dna.getDiameter()/2.d;
      Point location = current_dna.getLocationPoint();

      double x_origin = location.x+dradii;
      double y_origin = location.y+dradii+DNADraw.Y_SHIFT;
      double len = Math.sqrt(Math.pow((y_origin-y),2)+
                             Math.pow((x_origin-x),2));   
      getTrack().setPosition(len/dradii);
    }
    else
      getTrack().setPosition(y);
    viewer.refresh();
  }

  protected double getMidAngle()
  {
    return Math.abs(getAngStart() + (getAngEnd()/2.d));
  }

  protected Point[] getLinePoints(final int extraLength)
  {
    if(current_dna.isCircular())
    {
      double dradii = current_dna.getDiameter()/2.d;
      Point location = current_dna.getLocationPoint();

      double x_origin = location.x+dradii;
      double y_origin = location.y+dradii+DNADraw.Y_SHIFT;

      double midAngle = getMidAngle();
      double trackPos = (track.getPosition()*dradii)+(track.getSize()/2.f);
      
      double x  = 0.d;
      double y  = 0.d;
      double x2 = 0.d;
      double y2 = 0.d;
      double x3 = 0.d;
      int lineLength = 100+extraLength;
      
      if(midAngle <= 90.d)
      {
        x = (Math.sin(Math.toRadians(midAngle)) * trackPos)+x_origin;
        y = y_origin-(Math.cos(Math.toRadians(midAngle)) * trackPos);
        
        x2 = (Math.sin(Math.toRadians(midAngle)) * (trackPos+lineLength))+x_origin;
        y2 = y_origin-(Math.cos(Math.toRadians(midAngle)) * (trackPos+lineLength));
        
        x3 = x2+10;
      }
      else if(midAngle <= 180.d)
      {
        midAngle-=90;
        x = (Math.cos(Math.toRadians(midAngle)) * trackPos)+x_origin;
        y = (Math.sin(Math.toRadians(midAngle)) * trackPos)+y_origin;
        
        x2 = (Math.cos(Math.toRadians(midAngle)) * (trackPos+lineLength))+x_origin;
        y2 = (Math.sin(Math.toRadians(midAngle)) * (trackPos+lineLength))+y_origin;
        
        x3 = x2+10;
      }
      else if(midAngle <= 270.d)
      {
        midAngle-=180;
        x = x_origin-(Math.sin(Math.toRadians(midAngle)) * trackPos);
        y = (Math.cos(Math.toRadians(midAngle)) * trackPos)+y_origin;
        
        x2 = x_origin-(Math.sin(Math.toRadians(midAngle)) * (trackPos+lineLength));
        y2 = (Math.cos(Math.toRadians(midAngle)) * (trackPos+lineLength))+y_origin;
        
        x3 = x2-10;
      }
      else if(midAngle <= 360.d)
      {
        midAngle-=270;
        x = x_origin-(Math.cos(Math.toRadians(midAngle)) * trackPos);
        y = y_origin-(Math.sin(Math.toRadians(midAngle)) * trackPos);
        
        x2 = x_origin-(Math.cos(Math.toRadians(midAngle)) * (trackPos+lineLength));
        y2 = y_origin-(Math.sin(Math.toRadians(midAngle)) * (trackPos+lineLength));
        
        x3 = x2-10;
      }
      
      final Point[] line = new Point[3];
      line[0] = new Point((int)x, (int)y);
      line[1] = new Point((int)x2, (int)y2);
      line[2] = new Point((int)x3, (int)y2);
      
      return line;
    }
    else
      return null;
  }
  
  /**
  *
  * Routine to check whether a given location is over
  * the drawn block.
  * @param x	x position 
  * @param y	y position
  *
  */
  public boolean isOverMe(int x, int y)
  {
    if(current_dna.isCircular())
    {
      double dradii = current_dna.getDiameter()/2.d;
      Point location = current_dna.getLocationPoint();

      double x_origin = location.x+dradii;
      double y_origin = location.y+dradii+DNADraw.Y_SHIFT;

      double ttheta = Math.abs(y_origin-y)/Math.abs(x_origin-x);
      double angDegrees = Math.toDegrees(Math.atan(ttheta));
      
      if(x < x_origin && y > y_origin)
        angDegrees = (180 - angDegrees);
      else if(x < x_origin && y < y_origin)
        angDegrees = (angDegrees + 180);
      else if(x > x_origin && y < y_origin)
        angDegrees = (360 - angDegrees);

      angDegrees = angDegrees - DNADraw.THETA;
      if(angDegrees > 360)
        angDegrees-=360;
      
      if( angDegrees >= -angStart-.05 && angDegrees <= -(angStart+angEnd-0.05) )
      {
        double len = Math.sqrt(Math.pow((y_origin - y), 2)
            + Math.pow((x_origin - x), 2));
        //double rat = (len / dradii);

        /*System.out.println(getLabel() + "   " + getFeature().getIDString());
        System.out.println("ANGLE " + angDegrees );
        System.out.println("START " + angStart + "  " + angStart);
        System.out.println("END   " + (angStart + angEnd) + "  " + angEnd);*/
        double trackPos = track.getPosition()*dradii;
        
        if(len <= trackPos+(track.getSize()/2.f) && 
           len >= trackPos-(track.getSize()/2.f))
          return true;
      }
    }
    else
    {
      for(int i=0; i<rect.size(); i++)
      {
        Rectangle r = (Rectangle)rect.get(i);
        if(r.contains(x,y))
          return true;
      }
    }
    
    return false;
  }

  // Transferable
  public DataFlavor[] getTransferDataFlavors()
  {
    return blockFlavors;
  }

  public boolean isDataFlavorSupported(DataFlavor f)
  {
    if(f.equals(BLOCK))
      return true;
    return false;
  }

  public Object getTransferData(DataFlavor d)
      throws UnsupportedFlavorException, IOException
  {
    if(d.equals(BLOCK))
      return this;
    else throw new UnsupportedFlavorException(d);
  }


  public Track getTrack()
  {
    return track;
  }
  
  
  public void setTrack(Track track)
  {
    this.track = track;
  }


  public double getAngEnd()
  {
    return angEnd;
  }


  public void setAngEnd(double angEnd)
  {
    this.angEnd = angEnd;
  }


  public double getAngStart()
  {
    return angStart;
  }


  public void setAngStart(double angStart)
  {
    this.angStart = angStart;
  }


  public boolean isArrowHead()
  {
    return arrowHead;
  }


  public void setArrowHead(boolean arrowHead)
  {
    this.arrowHead = arrowHead;
  }


  public boolean isArrowTail()
  {
    return arrowTail;
  }


  public void setArrowTail(boolean arrowTail)
  {
    this.arrowTail = arrowTail;
  }


  public int getBend()
  {
    return bend;
  }


  public void setBend(int bend)
  {
    this.bend = bend;
  }


  public int getBstart()
  {
    return bstart;
  }


  public void setBstart(int bstart)
  {
    this.bstart = bstart;
  }


  public Color getColour()
  {
    if(track.getColour() != null)
      return track.getColour();
    return colour;
  }


  public void setColour(Color colour)
  {
    this.colour = colour;
  }


  public String getLabel()
  {
    return label;
  }

  public void setLabel(String label)
  {
    this.label = label;
  }


  public float getStrokeSize()
  {
    return strokeSize;
  }


  public void setStrokeSize(float strokeSize)
  {
    this.strokeSize = strokeSize;
  }


  public Feature getFeature()
  {
    return feature;
  }


  public void setFeature(Feature feature)
  {
    this.feature = feature;
  }


  public boolean isDrawLabel()
  {
    return drawLabel;
  }


  public void setDrawLabel(boolean drawLabel)
  {
    this.drawLabel = drawLabel;
  }

}

