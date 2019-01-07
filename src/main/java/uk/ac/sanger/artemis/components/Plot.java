/* Plot.java
 *
 * created: Thu Dec 17 1998
 *
 * This file is part of Artemis
 *
 * Copyright (C) 1998,1999,2000  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/Plot.java,v 1.26 2009-08-18 09:01:44 tjc Exp $
 **/

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.circular.TextFieldFloat;
import uk.ac.sanger.artemis.plot.*;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.GridLayout;
import java.awt.Image;
import java.awt.Scrollbar;
import java.awt.event.*;
import java.util.Vector;

import javax.swing.JMenu;
import javax.swing.JPanel;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JSplitPane;
import javax.swing.JTextField;
import javax.swing.JOptionPane;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JMenuItem;
import javax.swing.JScrollBar;
import javax.swing.JPopupMenu;

import org.apache.batik.svggen.SVGGraphics2D;

/**
 *  This class implements a simple plot component.
 *  @author Kim Rutherford
 **/

public abstract class Plot extends JPanel 
{
  private static final long serialVersionUID = 1L;

  /** scroll bar for changing the window size. */
  private JScrollBar window_changer = null;

  /** height of the font used in this component. */
  private int font_height;

  /** off screen image used for double buffering when drawing */
  private Image offscreen;

  /**
   *  The object that will generate the value we plot in this component.
   **/
  private Algorithm algorithm;

  /**
   *  If true then a scale line will be drawn at the bottom of the graph when
   *  drawScaleLine() is called.
   **/
  private boolean draw_scale;
  
  private final int SCROLL_NOB_SIZE = 10;

  /**
   *  Set to true if drawMultiValueGraph() should call recalculateValues().
   *  It is reset to false by recalculateValues().
   **/
  protected boolean recalculate_flag = true;
  
  private boolean showAverage = true;

  /**
   *  The x position of the last click or -1 if the user hasn't clicked
   *  anywhere yet or the user clicked outside the graph.
   **/
  private int cross_hair_position = -1;

  /**
   *  The x position of the start of the last mouse drag or -1 if the user
   *  hasn't clicked anywhere yet or the user clicked outside the graph.
   **/
  private int drag_start_position = -1;

  /**
   *  A vector of those objects listening for PlotMouse events.
   **/
  final private Vector<PlotMouseListener> listener_list = new Vector<PlotMouseListener>();

  /**
   *  Recalculate the values all the state that is used for drawing the plot
   **/
  protected abstract void recalculateValues();

  /**
   *  Get the position in the Feature or Sequence of the given x canvas
   *  position.  This is the label used when the user clicks the mouse in on
   *  the canvas (see drawCrossHair()).
   **/
  protected abstract int getPointPosition(final int canvas_x_position);
  
  protected abstract void calculateFeatures(boolean fromPeak);  

  protected abstract void showAveragesForRange();
  
  /** number of graph lines to be drawn */
  private int numPlots;

  protected LineAttributes lines[];
  
  private int lastPaintHeight = getHeight();
  
  /** the minimum distance in pixels between the labels */
  private final static int MINIMUM_LABEL_SPACING = 50;
 
  /**
   *  Create a new plot component.
   *  @param algorithm The object that will generate the values we plot in
   *    this component.
   *  @param draw_scale If true then a scale line will be drawn at the bottom
   *    of the graph.
   **/
  public Plot(Algorithm algorithm, boolean draw_scale) 
  {
    super(new BorderLayout());
    this.algorithm = algorithm;
    this.draw_scale = draw_scale;

    setFont(Options.getOptions().getFont());
    FontMetrics fm = getFontMetrics(getFont());
    font_height = fm.getHeight();

    final int MAX_WINDOW;
    if(getAlgorithm().getDefaultMaxWindowSize() != null) 
      MAX_WINDOW = getAlgorithm().getDefaultMaxWindowSize().intValue();
    else 
      MAX_WINDOW = 500;

    final int MIN_WINDOW;
    if(getAlgorithm().getDefaultMinWindowSize() != null) 
      MIN_WINDOW = getAlgorithm().getDefaultMinWindowSize().intValue();
    else 
      MIN_WINDOW = 5;

    final int START_WINDOW;
    if(getAlgorithm().getDefaultWindowSize() == null) 
      START_WINDOW = 10;
    else 
      START_WINDOW = getAlgorithm().getDefaultWindowSize().intValue();

    window_changer = new JScrollBar(Scrollbar.VERTICAL);
    window_changer.setValues(START_WINDOW, SCROLL_NOB_SIZE,
                             MIN_WINDOW, MAX_WINDOW + SCROLL_NOB_SIZE);
    if(MAX_WINDOW >= 50) 
      window_changer.setBlockIncrement(MAX_WINDOW/50);
    else 
      window_changer.setBlockIncrement(1);

    window_changer.addAdjustmentListener(new AdjustmentListener()
    {
      public void adjustmentValueChanged(AdjustmentEvent e) 
      {
        recalculate_flag = true;
        repaint();
      }
    });

    addComponentListener(new ComponentAdapter() 
    {
      public void componentShown(ComponentEvent e) 
      {
        recalculate_flag = true;
        repaint();
      }
    });

    add(window_changer, "East");
    addMouseListener(mouse_listener);
    addMouseMotionListener(mouse_motion_listener);
  }

  /**
   *  Return the algorithm that was passed to the constructor.
   **/
  public Algorithm getAlgorithm() 
  {
    return algorithm;
  }

  /**
   *  Return the current value of the window size, as set by the
   *  window_changer scrollbar.
   **/
  public int getWindowSize()
  {
    return window_changer.getValue();
  }

  final MouseListener mouse_listener = new MouseAdapter()
  {
    /**
     *  Listen for mouse press popup menu and crosshair events.
     **/
    public void mousePressed(MouseEvent event) 
    {
      if(event.isPopupTrigger() || event.isMetaDown()) 
      {
        final JComponent parent = (JComponent)event.getSource();
        final JPopupMenu popup  = new JPopupMenu("Plot Options");

        // configure colours for multiple graph plots
        final JMenuItem config = new JMenuItem("Configure...");
        config.addActionListener(new ActionListener()
        {
          public void actionPerformed(ActionEvent actionEvent)
          {
            lines = LineAttributes.configurePlots(numPlots, lines, Plot.this);
          }
        });
        popup.add(config);


        final JMenuItem setScale = new JMenuItem("Set the Window Size...");
        setScale.addActionListener(new ActionListener()
        {
          public void actionPerformed(ActionEvent actionEvent)
          {
            final JTextField newWinSize = new JTextField(Integer.toString(getWindowSize()));
            String window_options[] = { "Set Window Size", "Cancel" };
            int select = JOptionPane.showOptionDialog(null,
                                        newWinSize,
                                        "Set Window Size",
                                         JOptionPane.DEFAULT_OPTION,
                                         JOptionPane.QUESTION_MESSAGE,
                                         null, window_options, window_options[0]);
            final int value;
            try
            {
              value = Integer.parseInt(newWinSize.getText().trim());
            }
            catch(NumberFormatException nfe)
            {
              return;
            }
            if(value > window_changer.getMaximum() ||
               value < window_changer.getMinimum())
            {
              window_options[0] = "Continue";
              select = JOptionPane.showOptionDialog(null,
                                        "Value selected: " + value +
                                        " is outside the range\n"+
                                        " Min: "+window_changer.getMinimum() +
                                        " Max: "+window_changer.getMaximum(),
                                        "Set Window Size",
                                         JOptionPane.DEFAULT_OPTION,
                                         JOptionPane.WARNING_MESSAGE,
                                         null, window_options, window_options[1]);
              if(select == 1)
                return;

              if(value > window_changer.getMaximum())
                window_changer.setMaximum(value+10);
              else
                window_changer.setMinimum(value);
            }

            if(select == 0)
            {
              recalculate_flag = true;
              window_changer.setValue(value);
              repaint();
            }
          }
        });
        popup.add(setScale);


        final JCheckBoxMenuItem scaling_toggle =
          new JCheckBoxMenuItem("Scaling",getAlgorithm().scalingFlag());
        scaling_toggle.addItemListener(new ItemListener() 
        {
          public void itemStateChanged(ItemEvent actionEvent) 
          {
            getAlgorithm().setScalingFlag(scaling_toggle.getState());
            recalculate_flag = true;
            repaint();
          }
        });
        popup.add(scaling_toggle);

        final JCheckBoxMenuItem showAverageLn = new JCheckBoxMenuItem("Show average", showAverage);
        showAverageLn.addItemListener(new ItemListener() 
        {
          public void itemStateChanged(ItemEvent itemEvent) 
          {
            showAverage = showAverageLn.isSelected();
            repaint();
          }
        });
        popup.add(showAverageLn);

        if(Plot.this instanceof BasePlot)
        {
          final JMenuItem showMinMaxValues =
            new JMenuItem("Set Min/Max Values...");
          popup.add(showMinMaxValues);
          showMinMaxValues.addActionListener(new ActionListener() 
          {
            public void actionPerformed(ActionEvent e) 
            {
              JPanel gridPane = new JPanel(new GridLayout(2,2));
              TextFieldFloat minValue = new TextFieldFloat();
              minValue.setValue( ((BasePlot)Plot.this).getMin_value() );
              gridPane.add(new JLabel("Min:"));
              gridPane.add(minValue);
              
              TextFieldFloat maxValue = new TextFieldFloat();
              maxValue.setValue( ((BasePlot)Plot.this).getMax_value() );
              gridPane.add(new JLabel("Max:"));
              gridPane.add(maxValue);

              String window_options[] = { "Set", "Cancel" };
              int select = JOptionPane.showOptionDialog(null, gridPane,
                  "Set Min/Max Plot Values", JOptionPane.DEFAULT_OPTION,
                  JOptionPane.QUESTION_MESSAGE, null, window_options,
                  window_options[0]);
              if(select == 1)
                return;
              
              getAlgorithm().setUserMaxMin(true);
              getAlgorithm().setUserMin((float) minValue.getValue());
              getAlgorithm().setUserMax((float) maxValue.getValue());
              ((BasePlot)Plot.this).setMin_value((float) minValue.getValue());
              ((BasePlot)Plot.this).setMax_value((float) maxValue.getValue());
              getAlgorithm().setScalingFlag(false);
              repaint();
            }
          });
        }

        popup.addSeparator();

        final JMenu max_window_size =
              new JMenu("Maximum Window Size");
        popup.add(max_window_size);

        final int[] window_sizes = 
        {
          100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000,
          200000, 500000, 1000000
        };

        JMenuItem window_size_item;
        for(int i = 0 ; i < window_sizes.length ; ++i) 
        {
          final int size = i;
          window_size_item = new JMenuItem(" " + window_sizes[i]);

          window_size_item.addActionListener(new ActionListener() 
          {
            public void actionPerformed(ActionEvent actionEvent) 
            {
              final int new_maximum = window_sizes[size];
              if(new_maximum > window_changer.getMinimum()) 
              {
                window_changer.setMaximum(new_maximum + SCROLL_NOB_SIZE);
                recalculate_flag = true;
                repaint();
              }
            }
          });

          max_window_size.add(window_size_item);
        }
        
        if(numPlots == 1 && getAlgorithm() instanceof BaseAlgorithm)
        {
          popup.addSeparator();
          
          final JMenu createFeaturesFrom = new JMenu("Create features from graph");
          popup.add(createFeaturesFrom);
          
          final JMenuItem createFeaturesFromPeak =
            new JMenuItem("peaks...");
          createFeaturesFrom.add(createFeaturesFromPeak);
          
          createFeaturesFromPeak.addActionListener(new ActionListener() 
          {
            public void actionPerformed(ActionEvent e) 
            {
              calculateFeatures(true);
            }
          });
          
          final JMenuItem createFeaturesFromDip =
            new JMenuItem("trough...");
          createFeaturesFrom.add(createFeaturesFromDip);
          
          createFeaturesFromDip.addActionListener(new ActionListener() 
          {
            public void actionPerformed(ActionEvent e) 
            {
              calculateFeatures(false);
            }
          });
        }

        ///       
        if(Plot.this instanceof BasePlot)
        { 
          final JMenuItem showAverages =
            new JMenuItem("Values and average(s) for selected range...");
          popup.add(showAverages);
        
          showAverages.addActionListener(new ActionListener() 
          {
            public void actionPerformed(ActionEvent e) 
            {
              showAveragesForRange();
            }
          });
        }

        final JSplitPane splitPane = getJSplitPane();
        if(splitPane == null)
        {
          popup.addSeparator();
          final JMenu graphHeight = new JMenu("Graph Height");
          popup.add(graphHeight);
          final JMenuItem smaller = new JMenuItem("smaller");
          final JMenuItem larger = new JMenuItem("larger");
          final JMenuItem setHeight = new JMenuItem("set...");

          graphHeight.add(smaller);
          graphHeight.add(larger);
          graphHeight.add(setHeight);

          smaller.addActionListener(new ActionListener()
          {
            public void actionPerformed(ActionEvent e)
            {
              rescale((int) (getSize().height * 0.9f));
            }
          });

          larger.addActionListener(new ActionListener()
          {
            public void actionPerformed(ActionEvent e)
            {
              rescale((int) (getSize().height * 1.1f));
            }
          });

          setHeight.addActionListener(new ActionListener()
          {
            public void actionPerformed(ActionEvent e)
            {
              final JTextField newGraphHgt = new JTextField(Integer
                  .toString(getSize().height));
              String window_options[] = { "Set Window Size", "Cancel" };
              int select = JOptionPane.showOptionDialog(null, newGraphHgt,
                  "Set Window Size", JOptionPane.DEFAULT_OPTION,
                  JOptionPane.QUESTION_MESSAGE, null, window_options,
                  window_options[0]);

              if(select == 1)
                return;
              try
              {
                rescale(Integer.parseInt(
                    newGraphHgt.getText().trim()));
              }
              catch(NumberFormatException nfe)
              {
              }
            }
          });
        }

        parent.add(popup);
        popup.show(parent, event.getX(), event.getY());
      } 
      else 
      {
        final int point_x = event.getPoint().x;
        final int point_y = event.getPoint().y;

        if(point_y > getLabelHeight()) 
        {
          cross_hair_position = point_x;
          drag_start_position = point_x;
        } 
        else
          cancelCrossHairs();

        if(event.getClickCount() == 2) 
          fireDoubleClickEvent();
        else
          fireClickEvent();

        repaint();
      }
    }
  };

  final MouseMotionListener mouse_motion_listener =
      new MouseMotionAdapter() 
  {
    public void mouseDragged(MouseEvent event) 
    {
      if(isMenuTrigger(event))
        return;

      final int point_x = event.getPoint().x;
      final int point_y = event.getPoint().y;

      if(point_y > getLabelHeight()) 
        cross_hair_position = point_x;
      else 
        cancelCrossHairs();

      fireDragEvent();
      repaint();
    }
  };

  /**
   * Find JSplitPane parent container or null if it does not
   * belong to one.
   * @return
   */
  private JSplitPane getJSplitPane()
  {
    JComponent child = Plot.this;
    JSplitPane splitPane = null;
    int count = 0;
    try
    {
      while(splitPane == null && count < 10)
      {
        JComponent plotParent = (JComponent) child.getParent();
        if(plotParent instanceof JSplitPane)
          splitPane = (JSplitPane) plotParent;
        child = plotParent;
        count++;
      }
    }
    catch(Exception e){}
    return splitPane;
  }
  
  /**
   * Reset graph height
   * @param hgt
   */
  private void rescale(int hgt)
  {
    setSize(getSize().width, hgt);
    
    if(Plot.this instanceof BasePlot)
      BasePlot.HEIGHT = getSize().height;
    else if(Plot.this instanceof FeaturePlot)
      FeaturePlot.HEIGHT = getSize().height;
    
    offscreen = null;
    revalidate();
  }

  /**
   *  Return true if and only if the given MouseEvent (a mouse press) should
   *  pop up a JPopupMenu.
   **/
  private boolean isMenuTrigger(final MouseEvent event) 
  {
    if(event.isPopupTrigger() ||
      (event.getModifiers() & InputEvent.BUTTON3_MASK) != 0) 
      return true;
    else 
      return false;
  }

  /**
   *  Call mouseClick() on each of the PlotMouseListener objects in the
   *  listener list.
   **/
  private void fireClickEvent() 
  {
    PlotMouseListener listener;
    for(int i = 0; i < listener_list.size(); ++i)
    {
      listener = listener_list.elementAt(i);
      listener.mouseClick(getPointPosition(cross_hair_position));
    }
  }

  /**
   *  Call mouseDragged() on each of the PlotMouseListener objects in the
   *  listener list.
   **/
  private void fireDragEvent()
  {
    PlotMouseListener listener;
    for(int i = 0; i < listener_list.size(); ++i) 
    {
      listener = listener_list.elementAt(i);
      listener.mouseDrag(getPointPosition(drag_start_position),
                         getPointPosition(cross_hair_position));
    }
  }

  /**
   *  Call mouseDoubleClick() on each of the PlotMouseListener objects in the
   *  listener list.
   **/
  private void fireDoubleClickEvent() 
  {
    PlotMouseListener listener;
    for(int i = 0; i < listener_list.size(); ++i) 
    {
      listener = listener_list.elementAt(i);
      listener.mouseDoubleClick(getPointPosition(cross_hair_position));
    }
  }

  /**
   *  Adds the given listener to receive mouse events from this object.
   *  @param l the listener.
   **/
  public void addPlotMouseListener(final PlotMouseListener listener) 
  {
    listener_list.addElement(listener);
  }

  /**
   *  Removes the given listener from the list of those objects that receive
   *  mouse events from this object.
   *  @param l the listener.
   **/
  public void removePlotMouseListener(final PlotMouseListener listener) 
  {
    listener_list.removeElement(listener);
  }

  /**
   *  The main paint function for the canvas.  An off screen image used for
   *  double buffering when drawing the canvas.
   *  @param g The Graphics object of the canvas.
   **/
  protected void paintComponent(final Graphics g) 
  {
    super.paintComponent(g);
    if(!isVisible()) 
      return;

    final int width  = getWidth() - window_changer.getWidth();
    final int height = getHeight();

    // there is no point painting a zero width canvas
    if(height <= 0 || width <= 0) 
      return;

    if(offscreen == null || lastPaintHeight != height)
      offscreen = createImage(width, height);

    final Graphics og;
    if(g instanceof SVGGraphics2D)
      og = g;
    else
    {
      og = offscreen.getGraphics();
      og.setClip(0, 0, width, height);
      og.setColor(Color.WHITE);
      og.fillRect(0, 0, width, height);
    }
    
    // Redraw the graph on the canvas using the algorithm from the
    // constructor.

    if(lines == null && getAlgorithm() instanceof BaseAlgorithm)
    {
      final int get_values_return_count =
       ((BaseAlgorithm)getAlgorithm()).getValueCount();
      
      if(getAlgorithm() instanceof UserDataAlgorithm)
        lines = ((UserDataAlgorithm)getAlgorithm()).getLineAttributes();
      if(lines == null)
        lines = LineAttributes.init(get_values_return_count);
    }
    
    numPlots = drawMultiValueGraph(og,lines);
    drawLabels(og,numPlots);
    
    if( !(g instanceof SVGGraphics2D) )
    {
      g.drawImage(offscreen, 0, 0, null);
      og.dispose();
    }
    lastPaintHeight = height;
  }

  protected void resetOffscreenImage()
  {
    offscreen = null;
  }

  /**
   *  Return the canvas x position of the last click or -1 if the user hasn't
   *  clicked anywhere yet.
   **/
  protected int getCrossHairPosition() 
  {
    if(cross_hair_position >= getSize().width) 
      return -1;
    else 
      return cross_hair_position;
  }

  /**
   *  Force this component to stop drawing crosshairs.
   **/
  protected void cancelCrossHairs() 
  {
    cross_hair_position = -1;
    drag_start_position = -1;
  }

  /**
   *  Draw the scale line at the bottom of the graph (used by FeaturePlot).
   *  @param start The base on the left
   *  @param end The base on the right
   **/
  protected void drawScaleLine(final Graphics g,
                               final int start, final int end) 
  {
    final int hgt = getHeight();
    final int scale_number_y_pos =  hgt - 1;
    final float bases_per_pixel = 1.0F;

    final int possible_index_of_first_label = start / MINIMUM_LABEL_SPACING;
    final int index_of_first_label;

    if(possible_index_of_first_label == 0) 
      index_of_first_label = 1;
    else 
      index_of_first_label = possible_index_of_first_label;

    final int index_of_last_label = end / MINIMUM_LABEL_SPACING;
    for(int i = index_of_first_label; i <= index_of_last_label; i++)
    {
      final int scale_number_x_pos =
        (int)((i * MINIMUM_LABEL_SPACING - start) / bases_per_pixel);

      g.drawString(String.valueOf((int)(i * MINIMUM_LABEL_SPACING)),
                   scale_number_x_pos + 2,
                   scale_number_y_pos);

      g.drawLine(scale_number_x_pos, hgt - getScaleHeight() / 2,
                 scale_number_x_pos, hgt - getScaleHeight());
    }
  }

  /**
   *  Plot the given points onto a Graphics object.
   *  @param min_value The minimum of the plot_values.
   *  @param max_value The maximum of the plot_values.
   *  @param step_size The current step size for this algorithm.  This is
   *    never greater than window_size.
   *  @param window_size The window size used in calculating plot_values.
   *  @param total_unit_count The maximum number of residues/bases we can
   *    show.  This is used to draw the scale line and to calculate the
   *    distance (in pixels) between plot points.
   *  @param start_position The distance from the edge of the canvas (measured
   *    in residues/bases) to start drawing the plot.
   *  @param plot_values The values to plot.
   **/
  protected void drawPoints(final Graphics g,
                            final float min_value, final float max_value,
                            final int step_size, final int window_size,
                            final int total_unit_count,
                            final int start_position,
                            final float [] plot_values,
                            final int value_index,
                            final int numberPlots,
                            final boolean isWiggle,
                            final boolean isBlast) 
  {
    final float residues_per_pixel =
      (float) total_unit_count / getSize().width;

    // this is the height of the graph (slightly smaller than the canvas for
    // ease of viewing).
    final int graph_height = getSize().height -
      getLabelHeight() -       // leave room for the algorithm name
      getScaleHeight() -       // leave room for the scale
      2;

    // too small to draw
    if(graph_height < 5) 
      return;

    Color definedColours[] = null;
    String plotType = null;
    if(lines != null)
    {
      plotType = lines[value_index].getPlotType();
    
      int NUMBER_OF_SHADES = 254;
      if(plotType.equals(LineAttributes.PLOT_TYPES[2]))
      {
        definedColours = makeColours(lines[value_index].getLineColour(),
                                     NUMBER_OF_SHADES);
      }
    }

    final int number_of_values = plot_values.length;
    int start_residue;
    int end_residue;
    int start_x;
    int end_x;

    for(int i = 0; i<number_of_values - 1; ++i) 
    {
      if( (isBlast || isWiggle) && plot_values[i] == 0)
      {
        if( !(isBlast && plotType.equals(LineAttributes.PLOT_TYPES[0])) )
          continue;
      }
      start_residue = window_size / 2 + i * step_size + start_position;
      start_x = (int)(start_residue / residues_per_pixel);
      
      if(isWiggle)
      {
        int span = ((UserDataAlgorithm)getAlgorithm()).getWiggleSpan(value_index);
        end_residue   = start_residue + span;
      }
      else
        end_residue   = start_residue + step_size;
      end_x = (int)(end_residue / residues_per_pixel);

      // this is a number between 0.0 and 1.0
      final float scaled_start_value =
        (plot_values[i] - min_value) / (max_value - min_value);
      final int start_y =
        graph_height - (int)(scaled_start_value * graph_height) +
        getLabelHeight() + 1;

      final float scaled_end_value =
        (plot_values[i+1] - min_value) / (max_value - min_value);
      final int end_y =
        graph_height - (int)(scaled_end_value * graph_height) +
        getLabelHeight() + 1;
      
      if(plotType == null ||
         plotType.equals(LineAttributes.PLOT_TYPES[0]))
      {
        if(isWiggle)
        {
          g.drawLine(start_x, graph_height+getLabelHeight() + 1, start_x, start_y);
          g.drawLine(start_x, start_y, end_x, start_y);
          g.drawLine(end_x, graph_height+getLabelHeight() + 1, end_x, start_y);
        }
        else
          g.drawLine(start_x, start_y, end_x, end_y);
      }
      else if(plotType.equals(LineAttributes.PLOT_TYPES[1]))
      {
        if(isWiggle)
          g.fillRect(start_x, start_y, end_x-start_x, graph_height+getLabelHeight() + 1);
        {
          int xPoints[] = { start_x, end_x, end_x, start_x };
          int yPoints[] = { start_y, end_y,  
              graph_height+getLabelHeight() + 1,
              graph_height+getLabelHeight() + 1};
          g.fillPolygon(xPoints, yPoints, 4);
        }
      }
      else
      {
        int ytop = getLabelHeight() + 1 + 
                   (value_index*(graph_height/numberPlots));
        int ybtm = (graph_height/numberPlots);

        // set color based on value
        int colourIndex = 
          (int)(definedColours.length * 0.999 * scaled_start_value);

        if(colourIndex > definedColours.length - 1)
          colourIndex = definedColours.length - 1;
        else if (colourIndex < 0)
          colourIndex = 0;
        
        g.setColor(definedColours[ colourIndex ]);
        g.fillRect(start_x, ytop, end_x-start_x, ybtm);
      }
    }
  }
  
  /**
   * Generate the colours for heat map plots.
   * @param col
   * @param NUMBER_OF_SHADES
   * @return
   */
  public static Color[] makeColours(Color col, int NUMBER_OF_SHADES)
  {
    Color definedColour[] = new Color[NUMBER_OF_SHADES];
    for(int i = 0; i < NUMBER_OF_SHADES; ++i)
    {
      int R = col.getRed();
      int G = col.getGreen();
      int B = col.getBlue();
      int scale = NUMBER_OF_SHADES-i;
      
      if((R+scale) <= 255)
        R += scale;
      else
        R = 254;
      if((G+scale) <= 255)
        G += scale;
      else
        G = 254;
      if((B+scale) <= 255)
        B += scale;
      else
        B = 254;

      definedColour[i] = new Color(R,G,B);
    }
    return definedColour;
  }

  /**
   *  Redraw the graph on the canvas using the algorithm.
   *  @param g The object to draw into.
   **/
  protected abstract int drawMultiValueGraph(final Graphics g, LineAttributes[] lines);

  /**
   *  Draw a line representing the average of the algorithm over the feature.
   *  @param g The object to draw into.
   *  @param min_value The minimum value of the function for the range we are
   *    viewing
   *  @param max_value The maximum value of the function for the range we are
   *    viewing
   **/
  protected void drawGlobalAverage(final Graphics g,
                                    final float min_value,
                                    final float max_value) 
  {
    // if a heatmap do not show the average
    if(!showAverage || (
        lines != null && lines[0].getPlotType().equals(LineAttributes.PLOT_TYPES[2])))
        return;
    
    final Float average = getAlgorithm().getAverage();

    if(average != null) 
    {
      g.setColor(Color.gray);

      // this is the height of the graph (slightly smaller than the canvas for
      // ease of viewing).
      final int graph_height =
        getSize().height - getFontHeight();

      // this is a number between 0.0 and 1.0
      final float scaled_average =
        (average.floatValue() - min_value) / (max_value - min_value);

      final int position =
        graph_height -
        (int)(scaled_average * graph_height) +
        getFontHeight() + 1;

      g.drawLine(0, position,
                 getSize().width, position);

      final FontMetrics fm = g.getFontMetrics();

      final int width = getSize().width;

      final String average_string =
        String.valueOf(Math.round(average.floatValue() * 100.0) / 100.0);

      g.drawString(average_string,
                   width - fm.stringWidth(average_string) - 
                   window_changer.getWidth() - 1,
                   position);
    }
  }

  /**
   *  Put the algorithm name in the top left corner of the canvas and the
   *  window size in the bottom left.
   *  @param g The object to draw into.
   **/
  private void drawLabels(final Graphics g, final int numPlots) 
  {
    g.setColor(Color.black);

    String desc = getAlgorithm().getAlgorithmName() + "  Window size: " +
                  String.valueOf(window_changer.getValue());

    g.drawString(desc, 2, font_height);

    if(numPlots < 2 || numPlots > 10)
      return;

    final FontMetrics fm = g.getFontMetrics();
    int font_width = fm.stringWidth("2");

    int width = 0;
    for(LineAttributes ln : lines)
      width += ln.getLabelWidth(fm);
    width = getWidth() - window_changer.getWidth() - width;

    g.translate(width,0);
    ((BaseAlgorithm)getAlgorithm()).drawLegend(g,font_height,
                                               font_width,lines, numPlots);
    g.translate(-width,0);
  }

  /**
   *  The method converts the min_value and max_value to String objects and
   *  then draws them onto the canvas.  The min_value is drawn at the bottom
   *  right, max_value at the top right.
   **/
  protected void drawMinMax(final Graphics g,
                            final float min_value, final float max_value)
  {
    g.setColor(Color.black);

    final int width  = getWidth() - window_changer.getWidth();
    final int height = getHeight();

    g.drawLine(0, height - getScaleHeight(),
               width, height - getScaleHeight());

    g.drawLine(0, getLabelHeight(),
               width, getLabelHeight());

    final FontMetrics fm = g.getFontMetrics();

    final String min_string =
      String.valueOf(((int)(min_value * 100)) / 100.0);

    g.drawString(min_string,
                 width - fm.stringWidth(min_string) - 1,
                 height - 1 - getScaleHeight());

    final String max_string =
      String.valueOf(((int)(max_value * 100)) / 100.0);

    g.drawString(max_string,
                 width - fm.stringWidth(max_string) - 1,
                 1 + getFontHeight() * 2);
  }

  /**
   *  Draw a vertical line at the given position.
   *  @param label The label to use on the crosshair
   *  @param label_pos The position on the line at which the label should be
   *    drawn (0 is nearest the top).
   **/
  protected void drawCrossHair(final Graphics g, final int x_position,
                               final String label, final int label_pos) 
  {
    if(x_position >= 0) 
    {
      g.drawLine(x_position, getLabelHeight(),
                  x_position, getSize().height);

      g.drawString(label, x_position + 2,
                   getFontHeight() * (2 + label_pos) + 2);
    }
  }

  /**
   *  Return the amount of vertical space (in pixels) to use for the scale.
   **/
  protected int getScaleHeight() 
  {
    if(draw_scale) 
      return getFontHeight() + 2;
    else 
      return 0;
  }

  /**
   *  Return the height in algorithm name and label line (returns the font
   *  height plus a small amount).
   **/
  protected int getLabelHeight() 
  {
    return getFontHeight() + 2;
  }

  /**
   *  Return the height in pixels of the current font.
   **/
  private int getFontHeight() 
  {
    return font_height;
  }
  
  /**
   *  Used to get the Y coordinate for the tooltip text.
   *  @param step_size The current step size for this algorithm.  This is
   *    never greater than window_size.
   *  @param window_size The window size used in calculating plot_values.
   *  @param start_position The distance from the edge of the canvas (measured
   *    in residues/bases) to start drawing the plot.
   *  @param plot_values The values to plot.
   *  @param base_pos    The base (from getXCoordinate) position.
   **/
  protected float getYCoordinate(
      final int step_size, final int window_size,
      final int start_position,
      final float plot_values[], int base_pos)
  {
    int ypos = (int)((base_pos - start_position - (window_size/2))/step_size);
    if(ypos < 0)
      ypos = 0;
    else if(ypos > plot_values.length-1)
      ypos = plot_values.length-1;
    
    return plot_values[ypos];
  }
}
