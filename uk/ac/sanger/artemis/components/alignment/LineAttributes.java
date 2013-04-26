/* LineAttributes.java
 *
 * created: 2010
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2010  Genome Research Limited
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
 **/


package uk.ac.sanger.artemis.components.alignment;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.List;
import java.util.Random;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JColorChooser;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
import javax.swing.ListCellRenderer;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

/**
 * Define line attributes for graphs lines
 */
public class LineAttributes
{
  /** defines the colour */
  private Color lineColour = Color.black;
  /** defines the line Stroke - size and style */
  private BasicStroke stroke;
  
  private static float dotDash[] = {10.f, 5.f, 3.f, 5.f};
  private static float dash[] = {10.f};
  
  private static BasicStroke style1 = 
      new BasicStroke(1.f);
  private static BasicStroke style2 = 
      new BasicStroke(1.f, BasicStroke.CAP_BUTT,
                      BasicStroke.JOIN_MITER, 
                      1.f, dash, 0.f);
  private static BasicStroke style3 = 
      new BasicStroke(1.f, BasicStroke.CAP_BUTT, 
                      BasicStroke.JOIN_MITER, 
                      1.f, dotDash, 0.f);
  
  /** available stroke types */
  private static BasicStroke[] STROKES = 
       new BasicStroke[]{ style1, style2, style3 };
  
  /** fill in underneath wiggle plots */
  public static String PLOT_TYPES[]    = 
          { "Open", "Filled" };
  private String plotType = PLOT_TYPES[0];

  /**
   * Contruct a LineAttributes instance
   * @param lineColour
   */
  public LineAttributes(Color lineColour)
  {
    this.lineColour = lineColour;
    this.stroke = new BasicStroke(1.f);
  }
  
  public Color getLineColour()
  {
    return lineColour;
  }

  public void setLineColour(Color lineColour)
  {
    this.lineColour = lineColour;
  }

  public BasicStroke getStroke()
  {
    return stroke;
  }

  public void setStroke(BasicStroke stroke)
  {
    this.stroke = stroke;
  }
  
  protected static BasicStroke[] getStrokes()
  {
    return STROKES;
  }
  
  /**
   * Given a stroke return the index for it based on the style.
   * @param stroke
   * @return
   */
  private static int getStyleIndex(BasicStroke stroke)
  {
    float myDash[] = stroke.getDashArray();
    if(myDash != null && 
       myDash.length == dotDash.length)
      return 2;
    else if(myDash != null && 
            myDash.length == dash.length &&
            myDash[0] == dash[0])
      return 1;
    else
      return 0;
  }
  
  public String getPlotType()
  {
    return plotType;
  }

  public static LineAttributes[]  init(int numPlots)
  {
    final Color frameColour[] = { 
        Color.red,
        Color.blue,
        Color.black,
        new Color(0,200,0), 
        Color.magenta,
        new Color(100, 100, 0),
        new Color(0, 100, 100),
        new Color(100, 0, 100),
        new Color(255, 100, 0),
        new Color(0, 255, 255),
        new Color(100, 0, 255),
        new Color(150, 0, 75),
        new Color(255, 200, 0),
        new Color(255, 0, 150)};

    final LineAttributes lines[] = new LineAttributes[numPlots];
    
    if(numPlots == 1)
    {
      lines[0] = new LineAttributes(Color.black);
      return lines;
    }
    
    final Random rand = new Random();
    for(int i=0; i<numPlots; i++)
    {
      if(i < frameColour.length)
        lines[i] = new LineAttributes(frameColour[i]);
      else
        lines[i] = new LineAttributes(
            new Color(rand.nextFloat(), rand.nextFloat(), rand.nextFloat()));
    }
    return lines;
  }

  
  /**
   * Utility used by uk.ac.sanger.artemis.components.Plot
   * @param numPlots
   * @param lines
   * @param plot
   * @return
   */
  public static LineAttributes[] configurePlots(
                                    final List<String> bamList, 
                                    final LineAttributes[] lines,
                                    final JPanel covergePanel)
  {
    JPanel panel = new JPanel(new GridBagLayout());
    GridBagConstraints c = new GridBagConstraints();

    final int numPlots = bamList.size();
    final LineAttributes[] thislines;
    if(lines.length < numPlots)
      thislines = init(numPlots);
    else
      thislines = lines;
    
    int gridx = 0;
    int gridy = 0;
    c.gridy = gridy++;
    c.gridx = gridx;
    // open / filled    
    panel.add(new JLabel("Plot Options:"), c);
    final JComboBox openPlot = new JComboBox(PLOT_TYPES);
    openPlot.setSelectedItem(thislines[0].plotType);
    openPlot.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        for(int i=0; i<numPlots; i++)
          thislines[i].plotType =    
            (String) openPlot.getSelectedItem();
        covergePanel.repaint();
      }
    });
    c.gridy = gridy++;
    panel.add(openPlot, c);
    
    c.gridy = gridy++;
    panel.add(new JLabel("Line Options:"), c);
    c.gridx = gridx++;
    c.gridy = gridy++;
    gridx++;
    
    final JCheckBox applyToAll = new JCheckBox("Apply colour to all");

    panel.add(new JLabel("Colour"), c);
    c.gridx = gridx++;
    panel.add(new JLabel("Line style"), c);
   
    c.gridx = gridx;
    panel.add(new JLabel("Line size"), c);
    
    for(int i=0; i<numPlots; i++)
    {
      c.gridy = gridy++;
      final int colourNumber = i;
      gridx = 0;

      c.gridx = 0;
      c.gridwidth = GridBagConstraints.REMAINDER;
      c.anchor = GridBagConstraints.NORTHWEST;
      final JLabel bamFileName = new JLabel(bamList.get(i));
      panel.add(bamFileName,c);
      
      c.gridwidth = 1;
      c.gridy = gridy++;
      final JLabel colourLabel = new JLabel("   ");
      colourLabel.setBackground(thislines[i].getLineColour());
      colourLabel.setOpaque(true);
      c.gridx = gridx++;
      panel.add(colourLabel,c);

      JButton butt = new JButton("Select");
      butt.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent _)
        {
          Color newColour = JColorChooser.showDialog(null, "Colour Chooser",
              thislines[colourNumber].getLineColour());
          
          if(applyToAll.isSelected())
          {
            for(int iplot=0; iplot<numPlots; iplot++)
              thislines[iplot].setLineColour(newColour);
          }
          else
            thislines[colourNumber].setLineColour(newColour); 
          colourLabel.setBackground(thislines[colourNumber].getLineColour());
          covergePanel.repaint();
        }
      });
      c.gridx = gridx++;
      panel.add(butt, c);
      
      // line style
      final JSlider slider = new JSlider(1, 10, 
          (int)lines[colourNumber].getStroke().getLineWidth());
      Integer index[] = new Integer[STROKES.length];
      for(int j=0; j<index.length; j++)
        index[j] = j;
      final JComboBox lineStyle = new JComboBox(index);
      lineStyle.setRenderer(new LineStyleListRenderer());
      lineStyle.setSelectedIndex(
          LineAttributes.getStyleIndex(thislines[colourNumber].getStroke()));
      lineStyle.addItemListener(new ItemListener()
      {
        public void itemStateChanged(ItemEvent e)
        {
          thislines[colourNumber].setStroke(
              STROKES[lineStyle.getSelectedIndex()]);
          setLineSize(covergePanel, slider, thislines, colourNumber);
        }
      });
      c.gridx = gridx++;
      panel.add(lineStyle, c);
      
      
      // line size
      slider.addChangeListener(new ChangeListener()
      {
        public void stateChanged(ChangeEvent e)
        {
          setLineSize(covergePanel, slider, thislines, colourNumber);
        }
      });
      c.gridx = gridx;
      panel.add(slider, c);
      
    }
    
    c.gridx = 0;
    c.gridy = gridy++;
    c.gridwidth = GridBagConstraints.REMAINDER;
    c.anchor = GridBagConstraints.WEST;
    panel.add(applyToAll, c);

    JScrollPane jsp = new JScrollPane(panel, 
        JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
        JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
    Dimension d = Toolkit.getDefaultToolkit().getScreenSize();
    
    if(panel.getPreferredSize().height > d.height/2)
    {
      d.setSize(jsp.getPreferredSize().width, d.height/2);
      jsp.setPreferredSize(d);
    }
    
    String config_options[] = { "OK" };
    JOptionPane.showOptionDialog(null,
                                 jsp, "Configure Lines",
                                 JOptionPane.DEFAULT_OPTION,
                                 JOptionPane.QUESTION_MESSAGE,
                                 null, config_options, config_options[0]);
    return thislines;
  }
  
  
  private static void setLineSize(JPanel covergePanel, JSlider slider,
                           LineAttributes[] thislines, int number)
  {
    BasicStroke oldStroke = thislines[number].getStroke();
    BasicStroke newStroke = new BasicStroke(slider.getValue(),
        BasicStroke.CAP_SQUARE, BasicStroke.JOIN_MITER, 
        1.f, oldStroke.getDashArray(), 0.f);
    thislines[number].setStroke(newStroke);
    covergePanel.repaint();
  }
}


/**
 * Renderer for the JComboBox to define different line styles.
 */
class LineStyleListRenderer extends JComponent implements ListCellRenderer
{
  private static final long serialVersionUID = 1L;
  private int selectedIndex = 0;
  public LineStyleListRenderer()
  {
    super();
    setMinimumSize(new Dimension(100, 20));
    setPreferredSize(new Dimension(100, 20));
  }

  /**
   * This method finds the selected index.
   */
  public Component getListCellRendererComponent(JList list, Object value,
      int index, boolean isSelected, boolean cellHasFocus)
  {
    selectedIndex = ((Integer) value).intValue();
    return this;
  }
  
  public void paintComponent(Graphics g)
  {
    super.paintComponent(g);
    Graphics2D g2 = (Graphics2D)g;
    g2.setStroke(LineAttributes.getStrokes()[selectedIndex]);
    g2.drawLine(10, 10, 90, 10);
  }
}