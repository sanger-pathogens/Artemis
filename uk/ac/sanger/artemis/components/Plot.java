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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/Plot.java,v 1.3 2004-10-05 09:24:37 tjc Exp $
 **/

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.plot.*;

import java.awt.*;
import java.awt.event.*;

import javax.swing.JPanel;
import javax.swing.JComponent;
import javax.swing.JScrollPane;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JMenuItem;
import javax.swing.JScrollBar;
import javax.swing.JPopupMenu;

/**
 *  This class implements a simple plot component.
 *
 *  @author Kim Rutherford
 *  @version $Id: Plot.java,v 1.3 2004-10-05 09:24:37 tjc Exp $
 **/

public abstract class Plot extends JPanel 
{

  /**
   *  The drawing area for this component.
   **/
  private JComponent canvas = null;

  /**
   *  A scroll bar for changing the window size.
   **/
  private JScrollBar window_changer = null;

  /**
   *  The height of the font used in this component.
   **/
  private int font_height;

  /**
   *  Off screen image used for double buffering when drawing the canvas.
   **/
  private Image offscreen;

  /**
   *  The object that will generate the value we plot in this component.
   **/
  private Algorithm algorithm;

  /**
   *  If true then a scale line will be drawn at the bottom of the graph when
   *  drawScaleLine () is called.
   **/
  private boolean draw_scale;

  /**
   *  Set to true if drawMultiValueGraph () should call recalculateValues ().
   *  It is reset to false by recalculateValues ().
   **/
  protected boolean recalculate_flag = true;

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
  final private java.util.Vector listener_list = new java.util.Vector ();

  /**
   *  Create a new plot component.
   *  @param algorithm The object that will generate the values we plot in
   *    this component.
   *  @param draw_scale If true then a scale line will be drawn at the bottom
   *    of the graph.
   **/
  public Plot (Algorithm algorithm, boolean draw_scale) 
  {
    this.algorithm = algorithm;
    this.draw_scale = draw_scale;

    final Font font = Options.getOptions ().getFont ();

    setFont (font);
    FontMetrics fm = getFontMetrics (font);
    font_height = fm.getHeight();

    setLayout(new BorderLayout());

    final int MAX_WINDOW;

    if(getAlgorithm().getDefaultMaxWindowSize () != null) 
      MAX_WINDOW = getAlgorithm().getDefaultMaxWindowSize().intValue();
    else 
      MAX_WINDOW = 500;

    final int MIN_WINDOW;

    if(getAlgorithm ().getDefaultMinWindowSize () != null) 
      MIN_WINDOW = getAlgorithm ().getDefaultMinWindowSize ().intValue ();
    else 
      MIN_WINDOW = 5;

    final int START_WINDOW;

    if(getAlgorithm().getDefaultWindowSize() == null) 
      START_WINDOW = 10;
    else 
      START_WINDOW = getAlgorithm().getDefaultWindowSize().intValue();

    window_changer = new JScrollBar (Scrollbar.VERTICAL);
    window_changer.setValues (START_WINDOW, SCROLL_NOB_SIZE,
                              MIN_WINDOW, MAX_WINDOW + SCROLL_NOB_SIZE);
    if (MAX_WINDOW >= 50) 
      window_changer.setBlockIncrement (MAX_WINDOW/50);
    else 
      window_changer.setBlockIncrement (1);

    window_changer.addAdjustmentListener(new AdjustmentListener()
    {
      public void adjustmentValueChanged(AdjustmentEvent e) 
      {
        recalculate_flag = true;
        repaintCanvas ();
      }
    });

    addComponentListener(new ComponentAdapter() 
    {
      public void componentShown(ComponentEvent e) 
      {
        recalculate_flag = true;
        repaintCanvas ();
      }
    });

    add(window_changer, "East");

    createCanvas();

    getCanvas().setBackground (Color.white);
  }

  final int SCROLL_NOB_SIZE = 10;

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
    return window_changer.getValue ();
  }

  /**
   *  Create the canvas object for this BasePlot and add it to the component.
   **/
  private void createCanvas()
  {
    canvas = new JComponent() 
    {
      /**
       *  Set the offscreen buffer to null as part of invalidation.
       **/
      public void invalidate() 
      {
        super.invalidate();
        offscreen = null;
      }

      /**
       *  Override update to *not* erase the background before painting
       */
      public void update(Graphics g) 
      {
        paint(g);
      }

      /**
       *  Paint the canvas.
       */
      public void paint(Graphics g) 
      {
        paintCanvas (g);
      }
    };

    final MouseListener mouse_listener = new MouseAdapter()
    {
        /**
         *  Listen for mouse press events so that we can do a popup menu and a
         *  crosshair.
         **/
        public void mousePressed(MouseEvent event) 
        {
          if(event.isPopupTrigger() || event.isMetaDown()) 
          {
            final JComponent parent = (JComponent)event.getSource();
            final JPopupMenu popup  = new JPopupMenu("Plot Options");

            final JCheckBoxMenuItem scaling_toggle =
              new JCheckBoxMenuItem ("Scaling");

            scaling_toggle.setState(getAlgorithm ().scalingFlag ());
            scaling_toggle.addItemListener(new ItemListener () 
            {
              public void itemStateChanged (ItemEvent _) 
              {
                getAlgorithm ().setScalingFlag (scaling_toggle.getState ());
                recalculate_flag = true;
                repaintCanvas();
              }
            });

            popup.add(scaling_toggle);
            popup.addSeparator();

            final JMenuItem max_window_size =
              new JMenuItem ("Maximum Window Size:");

            popup.add (max_window_size);

            final int [] window_sizes = 
            {
              100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000,
              200000, 500000, 1000000
            };

            for(int i = 0 ; i < window_sizes.length ; ++i) 
            {
              final int size = i;

              JMenuItem window_size_item = new JMenuItem(" " + window_sizes[i]);

              window_size_item.addActionListener(new ActionListener() 
              {
                public void actionPerformed (ActionEvent _) 
                {
                  final int new_maximum = window_sizes[size];
                  if(new_maximum > window_changer.getMinimum ()) 
                  {
                    window_changer.setMaximum(new_maximum + SCROLL_NOB_SIZE);
                    recalculate_flag = true;
                    repaintCanvas ();
                  }
                }
              });

              popup.add (window_size_item);
            }

            parent.add (popup);
            popup.show (parent, event.getX (), event.getY ());
          } 
          else 
          {
            final int point_x = event.getPoint ().x;
            final int point_y = event.getPoint ().y;

            if(point_y > getLabelHeight ()) 
            {
              cross_hair_position = point_x;
              drag_start_position = point_x;
            } 
            else
              cancelCrossHairs ();

            if(event.getClickCount() == 2) 
              fireDoubleClickEvent();
            else
              fireClickEvent();

            repaintCanvas();
          }
        }
      };

    canvas.addMouseListener (mouse_listener);

    final MouseMotionListener mouse_motion_listener =
      new MouseMotionAdapter() 
    {
      public void mouseDragged (MouseEvent event) 
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
        repaintCanvas();
      }
    };

    canvas.addMouseMotionListener (mouse_motion_listener);

    add (canvas, "Center");
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
   *  Call mouseClick () on each of the PlotMouseListener objects in the
   *  listener list.
   **/
  private void fireClickEvent() 
  {
    for (int i = 0 ; i < listener_list.size() ; ++i)
    {
      final PlotMouseListener listener =
        (PlotMouseListener)listener_list.elementAt(i);
      listener.mouseClick(getPointPosition(cross_hair_position));
    }
  }

  /**
   *  Call mouseDragged () on each of the PlotMouseListener objects in the
   *  listener list.
   **/
  private void fireDragEvent()
  {
    for(int i = 0; i < listener_list.size(); ++i) 
    {
      final PlotMouseListener listener =
        (PlotMouseListener)listener_list.elementAt(i);
      listener.mouseDrag(getPointPosition(drag_start_position),
                         getPointPosition(cross_hair_position));
    }
  }

  /**
   *  Call mouseDoubleClick () on each of the PlotMouseListener objects in the
   *  listener list.
   **/
  private void fireDoubleClickEvent() 
  {
    for(int i = 0; i < listener_list.size(); ++i) 
    {
      final PlotMouseListener listener =
        (PlotMouseListener)listener_list.elementAt (i);
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
  protected void paintCanvas(final Graphics g) 
  {
    if(!isVisible()) 
      return;

    final int canvas_width  = canvas.getSize().width;
    final int canvas_height = canvas.getSize().height;

    if(canvas_height <= 0 || canvas_width <= 0) 
    {
      // there is no point painting a zero width canvas
      return;
    }

    if(offscreen == null) 
      offscreen = createImage(canvas_width,
                              canvas_height);

    Graphics og = offscreen.getGraphics ();
    og.setClip(0, 0, canvas_width, canvas_height);

    og.setColor(new Color (240, 240, 240));

    og.fillRect(0, 0, canvas.getSize().width, canvas.getSize().height);

    // Redraw the graph on the canvas using the algorithm from the
    // constructor.
    drawMultiValueGraph(og);

    drawLabels(og);

    g.drawImage(offscreen, 0, 0, null);
    og.dispose();
  }

  /**
   *  Return the canvas x position of the last click or -1 if the user hasn't
   *  clicked anywhere yet.
   **/
  protected int getCrossHairPosition() 
  {
    if(cross_hair_position >= getCanvas().getSize().width) 
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

  // the minimum distance in pixels between the labels
  private final static int MINIMUM_LABEL_SPACING = 50;

  /**
   *  Draw the scale line at the bottom of the graph.
   *  @param start The base on the left
   *  @param end The base on the right
   **/
  protected void drawScaleLine(final Graphics g,
                               final int start, final int end) 
  {
    final int canvas_width = canvas.getSize ().width;
    final int canvas_height = canvas.getSize ().height;

    final int scale_number_y_pos =  canvas_height - 1;

    final float bases_per_pixel = 1.0F;

    // set the spacing so that the labels are at multiples of 10
    final int base_label_spacing = MINIMUM_LABEL_SPACING;

    final int label_spacing = (int)(base_label_spacing / bases_per_pixel);

    final int possible_index_of_first_label = start / base_label_spacing;

    final int index_of_first_label;

    if (possible_index_of_first_label == 0) 
      index_of_first_label = 1;
    else 
      index_of_first_label = possible_index_of_first_label;

    final int index_of_last_label = end / base_label_spacing;

    for (int i = index_of_first_label ;
         i <= index_of_last_label ;
         ++i) 
    {
      final String label_string =
        String.valueOf ((int)(i * base_label_spacing));

      final int scale_number_x_pos =
        (int) ((i * base_label_spacing - start) / bases_per_pixel);

      g.drawString (label_string,
                    scale_number_x_pos + 2,
                    scale_number_y_pos);

      g.drawLine (scale_number_x_pos, canvas_height - getScaleHeight () / 2,
                  scale_number_x_pos, canvas_height - getScaleHeight ());
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
  protected void drawPoints (final Graphics g,
                             final float min_value, final float max_value,
                             final int step_size, final int window_size,
                             final int total_unit_count,
                             final int start_position,
                             final float [] plot_values) 
  {
    final float residues_per_pixel =
      (float) total_unit_count / canvas.getSize ().width;

    // this is the height of the graph (slightly smaller than the canvas for
    // ease of viewing).
    final int graph_height = canvas.getSize ().height -
      getLabelHeight () -       // leave room for the algorithm name
      getScaleHeight () -       // leave room for the scale
      2;

    if (graph_height < 5) {
      // too small to draw
      return;
    }

    final int number_of_values = plot_values.length;

    for(int i = 0; i<number_of_values - 1; ++i) 
    {
      final int start_residue =
        window_size / 2 + i * step_size + start_position;
      final int end_residue =
        start_residue + step_size;

      final int start_x = (int) (start_residue / residues_per_pixel);
      final int end_x = (int) (end_residue / residues_per_pixel);

      // this is a number between 0.0 and 1.0
      final float scaled_start_value =
        (plot_values[i] - min_value) / (max_value - min_value);
      final int start_y =
        graph_height - (int) (scaled_start_value * graph_height) +
        getLabelHeight () + 1;

      final float scaled_end_value =
        (plot_values[i+1] - min_value) / (max_value - min_value);
      final int end_y =
        graph_height - (int) (scaled_end_value * graph_height) +
        getLabelHeight () + 1;

      g.drawLine (start_x, start_y, end_x, end_y);
    }
  }

  /**
   *  Redraw the graph on the canvas using the algorithm.
   *  @param g The object to draw into.
   **/
  protected abstract void drawMultiValueGraph (final Graphics g);

  /**
   *  Draw a line representing the average of the algorithm over the feature.
   *  @param g The object to draw into.
   *  @param min_value The minimum value of the function for the range we are
   *    viewing
   *  @param max_value The maximum value of the function for the range we are
   *    viewing
   **/
  protected void drawGlobalAverage (final Graphics g,
                                    final float min_value,
                                    final float max_value) 
  {
    final Float average = getAlgorithm ().getAverage ();

    if (average != null) {
      g.setColor (Color.gray);

      // this is the height of the graph (slightly smaller than the canvas for
      // ease of viewing).
      final int graph_height =
        canvas.getSize ().height - getFontHeight ();

      // this is a number between 0.0 and 1.0
      final float scaled_average =
        (average.floatValue () - min_value) / (max_value - min_value);

      final int position =
        graph_height -
        (int) (scaled_average * graph_height) +
        getFontHeight () + 1;

      g.drawLine (0, position,
                  canvas.getSize ().width, position);

      final FontMetrics fm = g.getFontMetrics ();

      final int canvas_width = canvas.getSize ().width;

      final String average_string =
        String.valueOf (Math.round (average.floatValue () * 100.0) / 100.0);

      g.drawString (average_string,
                    canvas_width - fm.stringWidth (average_string) - 1,
                    position);
    }
  }

  /**
   *  Put the algorithm name in the top left corner of the canvas and the
   *  window size in the bottom left.
   *  @param g The object to draw into.
   **/
  private void drawLabels (final Graphics g) 
  {
    g.setColor (Color.black);
    g.drawString (getAlgorithm ().getAlgorithmName () +
                  "  Window size: " +
                  String.valueOf (window_changer.getValue ()),
                  2,
                  font_height);
  }

  /**
   *  The method converts the min_value and max_value to String objects and
   *  then draws them onto the canvas.  The min_value is drawn at the bottom
   *  right, max_value at the top right.
   **/
  protected void drawMinMax(final Graphics g,
                            final float min_value, final float max_value)
  {
    g.setColor (Color.black);

    final int canvas_width = canvas.getSize ().width;
    final int canvas_height = canvas.getSize ().height;

    g.drawLine (0, canvas_height - getScaleHeight (),
                canvas_width, canvas_height - getScaleHeight ());

    g.drawLine (0, getLabelHeight (),
                canvas_width, getLabelHeight ());

    final FontMetrics fm = g.getFontMetrics ();

    final String min_string =
      String.valueOf (((int) (min_value * 100)) / 100.0);

    g.drawString (min_string,
                  canvas_width - fm.stringWidth (min_string) - 1,
                  canvas_height - 1 - getScaleHeight ());

    final String max_string =
      String.valueOf (((int) (max_value * 100)) / 100.0);

    g.drawString (max_string,
                  canvas_width - fm.stringWidth (max_string) - 1,
                  1 + getFontHeight () * 2);
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
    if (x_position >= 0) 
    {
      g.drawLine (x_position, getLabelHeight (),
                  x_position, canvas.getSize ().height);

      g.drawString (label, x_position + 2,
                    getFontHeight () * (2 + label_pos) + 2);
    }
  }

  /**
   *  Return the JComponent of this Plot.
   **/
  protected JComponent getCanvas() 
  {
    return canvas;
  }

  /**
   *  Call repaint () on the canvas object.
   **/
  protected void repaintCanvas() 
  {
    canvas.repaint ();
  }

  /**
   *  Recalculate the values all the state that is used for drawing the plot
   **/
  protected abstract void recalculateValues ();

  /**
   *  Get the position in the Feature or Sequence of the given x canvas
   *  position.  This is the label used when the user clicks the mouse in on
   *  the canvas (see drawCrossHair ()).
   **/
  protected abstract int getPointPosition (final int canvas_x_position);

  /**
   *  Return the amount of vertical space (in pixels) to use for the scale.
   **/
  private int getScaleHeight() 
  {
    if (draw_scale) 
      return getFontHeight() + 2;
    else 
      return 0;
  }

  /**
   *  Return the height in algorithm name and label line (returns the font
   *  height plus a small amount).
   **/
  private int getLabelHeight() 
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

}
