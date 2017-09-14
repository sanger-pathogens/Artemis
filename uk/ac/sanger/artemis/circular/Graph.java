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
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.RenderingHints;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.geom.AffineTransform;

import javax.swing.*;

import uk.ac.sanger.artemis.sequence.Bases;

public abstract class Graph extends JPanel
{
  private static final long serialVersionUID = 1L;
  
  private DNADraw currentDna;
  private int windowSize = 10000;
  private int baseStepSize = 200;
  private float[] value_array = null;
  private float gcAverage = 0;
  private float maxValue = Float.MIN_VALUE;
  private float minValue = Float.MAX_VALUE;
  private float graphHeight = 0.2f;
  private int strokeSize = 1;
  /** track to place the graph on - as a fraction of the radii*/
  private double track = .4d;
  private Color minusColour = new Color(0.6f, 0.f, 0.6f);
  private Color plusColour  = new Color(0.7f, 0.7f, 0.1f);
  
  public Graph(DNADraw currentDna)
  {
    this.currentDna = currentDna;
    setOpaque(false);
    setPreferredSize(currentDna.getPreferredSize());
  }
  
  
  protected void paintComponent(Graphics g)
  {
    super.paintComponent(g);
    Graphics2D g2  = (Graphics2D)g;
    if(getCurrentDna().isCircular())
      draw(g2);
    else
      drawLinear(g2);
  }
  
  /**
   *  Recalculate the values in value_array_array, step_size, min_value and
   *  max_value.
   **/
  protected abstract float calculateValue(int start, int end);
  
  public void draw(Graphics2D g2)
  {
    RenderingHints qualityHints = new RenderingHints(RenderingHints.KEY_ANTIALIASING,
        RenderingHints.VALUE_ANTIALIAS_ON); 
    qualityHints.put(RenderingHints.KEY_RENDERING,               
    RenderingHints.VALUE_RENDER_QUALITY); 
    g2.setRenderingHints(qualityHints);

    AffineTransform origin  = g2.getTransform();
    AffineTransform newOrig = (AffineTransform)(origin.clone());
    newOrig = (AffineTransform) (origin.clone());

    double ddiameter   = getCurrentDna().getDiameter();
    double widthPanel  = getCurrentDna().getWidth();
    double heightPanel = getCurrentDna().getHeight();
    
    Bases bases = getBases();
    /*if(gcAverage < 1.f)
      gcAverage = calculateValue(1,bases.getLength()); */ 
    
    int nvalues = bases.getLength()/getBaseStepSize();
    
    if(value_array == null)
      calcGraphValues();
    
    int minPos = (int)((ddiameter/2.d)*getTrack());
    int maxPos = minPos + (int)((ddiameter/2.d)*getGraphHeight());
    
    BasicStroke basicstroke = new BasicStroke(
        getStrokeSize(),
        BasicStroke.CAP_BUTT, 
        BasicStroke.JOIN_MITER);

    g2.setStroke(basicstroke);

    newOrig.translate(widthPanel/2.d,heightPanel/2.d);
    g2.setTransform(newOrig);
    
    basicstroke = new BasicStroke(
      getStrokeSize(),
        BasicStroke.CAP_BUTT, 
        BasicStroke.JOIN_MITER);

    g2.setStroke(basicstroke);
    
    float theta = ((float)getBaseStepSize()*360.f)/((float)bases.getLength());
    
    int avVal = minPos + (int)( ((gcAverage-minValue)/(maxValue-minValue))*(maxPos-minPos) );
    for(int i=0;i<nvalues-1; i++)
    {
      //int midBase = (i*baseStepSize)+(getWindowSize()/2);
      newOrig.rotate(Math.toRadians(theta));
      g2.setTransform(newOrig);
      
      float xA = value_array[i];
      if(xA == 0.f)
        continue;
      if(value_array[i] >= gcAverage)
        g2.setColor(getPlusColour());
      else
        g2.setColor(getMinusColour());

      int val = minPos + (int)( ((xA-minValue)/(maxValue-minValue))*(maxPos-minPos) );
      g2.drawLine(avVal, 0, val, 0);
    }
    g2.setTransform(origin);
  }
  
  
  public void drawLinear(Graphics2D g2)
  {
    RenderingHints qualityHints = new RenderingHints(RenderingHints.KEY_ANTIALIASING,
        RenderingHints.VALUE_ANTIALIAS_ON); 
    qualityHints.put(RenderingHints.KEY_RENDERING,               
    RenderingHints.VALUE_RENDER_QUALITY); 
    g2.setRenderingHints(qualityHints);

    int basesPerLine = getCurrentDna().getBasesPerLine();
    float lineHeight = getCurrentDna().getLineHeight();
    float singleBaseWidth = getCurrentDna().getSingleBaseWidth();
    int borderWidth2 = getCurrentDna().getBorderWidth2();
    int borderHeight2 = getCurrentDna().getBorderHeight2();
    
    Bases bases = getBases();
    int nvalues = bases.getLength()/getBaseStepSize();
    
    if(value_array == null)
      calcGraphValues();
    
    int minPos = (int)(lineHeight*(1-getTrack()));
    int maxPos = minPos + (int)(lineHeight*getGraphHeight());
    
    BasicStroke basicstroke = new BasicStroke(
        getStrokeSize(),
        BasicStroke.CAP_BUTT, 
        BasicStroke.JOIN_MITER);
    g2.setStroke(basicstroke);
    
    int avVal = minPos - (int)( ((gcAverage-minValue)/(maxValue-minValue))*(maxPos-minPos) );
    for(int i=0;i<nvalues-1; i++)
    {
      float xA = value_array[i];
       
      if(value_array[i] >= gcAverage)
        g2.setColor(getPlusColour());
      else
        g2.setColor(getMinusColour());

      int basePos = (i*getBaseStepSize())+1;
      int lineNumber = Math.round((float)(basePos-1)/((float)basesPerLine) + 0.5f)-1;
      
      int ypos = (int)(borderHeight2 + (lineNumber*lineHeight) 
          + (lineHeight*getGraphHeight()));
      
      int val = minPos - (int)( ((xA-minValue)/(maxValue-minValue))*(maxPos-minPos) );
      int xpos = (int) ((basePos-(lineNumber*basesPerLine))*singleBaseWidth)+borderWidth2;
      
      g2.drawLine(xpos, avVal+ypos, xpos, val+ypos);
    }
  }

  protected void calcGraphValues()
  {
    Bases bases = getBases();
    int nvalues = bases.getLength()/getBaseStepSize();
    value_array = new float [nvalues];
    gcAverage   = 0;
    for(int i=0; i<nvalues; i++)
    {
      int start = (i*getBaseStepSize())+1;
      int end = start+getWindowSize();
      
      if(end > bases.getLength())
        end = bases.getLength();
      
      value_array[i] = calculateValue(start, end);
      if(value_array[i] > maxValue)
        maxValue = value_array[i];
      if(value_array[i] < minValue)
        minValue = value_array[i];
      gcAverage += value_array[i];
    }
    gcAverage = gcAverage/nvalues;
  }
  
  protected int getWindowSize()
  {
    return windowSize;
  }
  
  protected void  setWindowSize(int windowSize)
  {
    value_array = null;
    this.windowSize = windowSize;
  }
  
  protected Bases getBases()
  {
  	return currentDna.getBases();
  }

  public int getStrokeSize() 
  {
	return strokeSize;
  }

  public void setStrokeSize(int strokeSize) 
  {
  	this.strokeSize = strokeSize;
  }

  public double getTrack()
  {
    return track;
  }

  public void setTrack(double track)
  {
    this.track = track;
  }

  public int getBaseStepSize()
  {
    return baseStepSize;
  }

  public void setBaseStepSize(int baseStepSize)
  {
    value_array = null;
    this.baseStepSize = baseStepSize;
  }

  public DNADraw getCurrentDna()
  {
    return currentDna;
  }


  public float getMaxValue()
  {
    return maxValue;
  }


  public void setMaxValue(float maxValue)
  {
    this.maxValue = maxValue;
  }


  public float getMinValue()
  {
    return minValue;
  }


  public void setMinValue(float minValue)
  {
    this.minValue = minValue;
  }


  public float getGraphHeight()
  {
    return graphHeight;
  }


  public void setGraphHeight(float graphHeight)
  {
    this.graphHeight = graphHeight;
  }
  

  public Color getMinusColour()
  {
    return minusColour;
  }


  public void setMinusColour(Color minusColour)
  {
    this.minusColour = minusColour;
  }


  public Color getPlusColour()
  {
    return plusColour;
  }


  public void setPlusColour(Color plusColour)
  {
    this.plusColour = plusColour;
  }
  
  
  /**
   * Used to change the graph default options
   */
  protected void showOptions()
  {
    GridBagLayout grid = new GridBagLayout();
    GridBagConstraints c = new GridBagConstraints();
    c.ipady = 3;
    c.ipadx = 5;

    JPanel optionBox = new JPanel(grid);
    
    // GRAPH HEIGHT
    c.gridx = 0;
    c.gridy = 0;
    c.anchor = GridBagConstraints.EAST;
    
    optionBox.add(new JLabel("Graph Height"), c);
    
    c.gridx = 1;
    c.anchor = GridBagConstraints.WEST;
    //final JSlider slider = new JSlider(1,60,(int)(getGraphHeight()*100));
    final TextFieldFloat graphHeightField = new TextFieldFloat();
    graphHeightField.setValue(getGraphHeight());
    graphHeightField.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        setGraphHeight((float) graphHeightField.getValue());
        repaint();
      }
    });
    optionBox.add(graphHeightField, c);
    
    
    // WINDOW SIZE
    c.gridx = 0;
    c.gridy = 1;
    c.anchor = GridBagConstraints.EAST;
    optionBox.add(new JLabel("Window Size"), c);
    
    c.gridx = 1;
    c.anchor = GridBagConstraints.WEST;
    final TextFieldInt winSize = new TextFieldInt();
    winSize.setValue(getWindowSize());
    
    winSize.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        setWindowSize(winSize.getValue());
        repaint();
      }
    });
    optionBox.add(winSize, c);


    // STEP SIZE
    c.gridx = 0;
    c.gridy = 2;
    c.anchor = GridBagConstraints.EAST;
    optionBox.add(new JLabel("Step Size"), c);

    c.gridx = 1;
    c.anchor = GridBagConstraints.WEST;
    final TextFieldInt stepSize = new TextFieldInt();
    stepSize.setValue(getBaseStepSize());
    
    stepSize.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        setBaseStepSize(stepSize.getValue());
        repaint();
      }
    });
    optionBox.add(stepSize, c);   
    
    // TRACK POSITION
    c.gridx = 0;
    c.gridy = 4;
    c.anchor = GridBagConstraints.EAST;
    optionBox.add(new JLabel("Track"), c);   
    
    
    c.gridx = 1;
    c.anchor = GridBagConstraints.WEST;
    final TextFieldFloat trackField = new TextFieldFloat();
    trackField.setValue(getTrack());
    
    trackField.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        setTrack(trackField.getValue());
        repaint();
      }
    });
    optionBox.add(trackField, c);
    
    // COLOUR
    final JButton button = new JButton();
    final JColorChooser colorChooser = new JColorChooser();
    
    ActionListener okListener = new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        button.setBackground(colorChooser.getColor());
        setMinusColour(colorChooser.getColor());
        repaint();
      }
    };
    setColorButton(button, getMinusColour(), okListener, colorChooser);
    c.gridx = 0;
    c.gridy = 5;
    c.anchor = GridBagConstraints.EAST;
    optionBox.add(new JLabel("Below Average"), c); 
    c.gridx = 1;
    c.anchor = GridBagConstraints.WEST;
    optionBox.add(button, c); 
    
    
    // COLOUR
    final JButton button2 = new JButton();
    final JColorChooser colorChooser2 = new JColorChooser();
    
    ActionListener okListener2 = new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        button2.setBackground(colorChooser2.getColor());
        setPlusColour(colorChooser2.getColor());
        repaint();
      }
    };
    setColorButton(button2, getPlusColour(), okListener2, colorChooser2);
    c.gridx = 0;
    c.gridy = 6;
    c.anchor = GridBagConstraints.EAST;
    optionBox.add(new JLabel("Above Average"), c);
    c.gridx = 1;
    c.anchor = GridBagConstraints.WEST;
    optionBox.add(button2, c); 
    
    JOptionPane.showMessageDialog(null, 
        optionBox, "Graph Options", 
        JOptionPane.PLAIN_MESSAGE);
    
    try
    {
      setGraphHeight((float) graphHeightField.getValue());
      setWindowSize(winSize.getValue());
      setBaseStepSize(stepSize.getValue());
      setTrack(trackField.getValue());
      //setMinusColour(colorChooser.getColor());
      //setPlusColour(colorChooser2.getColor());
    }
    catch(Exception e){ e.printStackTrace(); }
    repaint();
  }
  
  /**
   * Used to write out options to template file
   * @return
   */
  protected String getOptionsStr()
  {
    return "height="+getGraphHeight()+" window_size="+getWindowSize()+
           " base_step_size="+getBaseStepSize()+" track="+getTrack()+
           " minus_colour="+getMinusColour().getRed()+":"+
                            getMinusColour().getGreen()+":"+
                            getMinusColour().getBlue()+
           " plus_colour="+getPlusColour().getRed()+":"+
                           getPlusColour().getGreen()+":"+
                           getPlusColour().getBlue();
  }
  
  /**
   * Used when reading in a template file
   * @param options
   */
  protected void setOptionsStr(final String options[])
  {
    for(int i=0; i<options.length; i++)
    {
      if(options[i].startsWith("height"))
        setGraphHeight(Float.parseFloat(options[i+1]));
      else if(options[i].startsWith("window_size"))
        setWindowSize(Integer.parseInt(options[i+1]));
      else if(options[i].startsWith("base_step_size"))
        setBaseStepSize(Integer.parseInt(options[i+1]));
      else if(options[i].startsWith("track"))
        setTrack(Double.parseDouble(options[i+1]));
      else if(options[i].startsWith("minus_colour"))
      {
        String col[] = options[i+1].split(":");
        setMinusColour(new Color(Integer.parseInt(col[0]),
                                 Integer.parseInt(col[1]),
                                 Integer.parseInt(col[2])));
      }
      else if(options[i].startsWith("plus_colour"))
      {
        String col[] = options[i+1].split(":");
        setPlusColour(new Color(Integer.parseInt(col[0]),
                                 Integer.parseInt(col[1]),
                                 Integer.parseInt(col[2])));
      }
    }
  }
  
  private void setColorButton(final JButton button,
                              final Color col,
                              final ActionListener okListener,
                              final JColorChooser colorChooser)
  { 
    button.setOpaque(true);
    button.setBackground(col);
    button.setBorderPainted(false);
    button.setMargin(new Insets(0,0,0,0));
    Dimension d = new Dimension(25,25);
    button.setPreferredSize(d);

   //Set up the dialog that the button brings up.
    final JDialog dialog = JColorChooser.createDialog(button,
                                        "Pick a Color",
                                        true,
                                        colorChooser,
                                        okListener,
                                        null); 

    button.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        dialog.setVisible(true);
      }
    });
  }

}