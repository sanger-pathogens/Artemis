/* BasePlot.java
 *
 * created: Tue Dec 15 1998
 *
 * This file is part of Artemis
 *
 * Copyright (C) 1998,1999,2000,2001  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/BasePlot.java,v 1.3 2004-11-25 15:24:24 tjc Exp $
 **/

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.sequence.*;
import uk.ac.sanger.artemis.plot.*;
import uk.ac.sanger.artemis.util.OutOfRangeException;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

/**
 *  A component for plotting functions over the base sequence.  Scrolling and
 *  scale is tied to a FeatureDisplay component.
 *
 *  @author Kim Rutherford
 *  @version $Id: BasePlot.java,v 1.3 2004-11-25 15:24:24 tjc Exp $
 **/

public class BasePlot extends Plot
    implements DisplayAdjustmentListener, SelectionChangeListener 
{
  /**
   *  Create a new FeatureDisplay object.
   *  @param algorithm The object that will generate the value we plot in
   *    this component.
   *  @param selection Used to set and display the current selection in the
   *    BasePlot.
   *  @param goto_event_source The object the we will call gotoBase() on.
   *    This allows the user to double click on a base in a BasePlot and have
   *    the FeatureDisplay follow.
   **/
  public BasePlot(final BaseAlgorithm algorithm,
                  final Selection selection,
                  final GotoEventSource goto_event_source) 
  {
    super(algorithm, false);   // false means don't draw the scale line

    this.selection = selection;
    this.goto_event_source = goto_event_source;
    this.bases = getBaseAlgorithm().getBases();

    getSelection().addSelectionChangeListener(this);

    addPlotMouseListener(new PlotMouseListener() 
    {
      /**
       *  Set the selection to be a range from start_base to end_base.
       **/
      private void setSelectionRange(final int start_base,
                                     final int end_base) 
      {
        final Strand strand = bases.getForwardStrand();

        try 
        {
          final MarkerRange marker_range =
            strand.makeMarkerRangeFromPositions(start_base, end_base);
          getSelection().setMarkerRange(marker_range);
        } 
        catch(uk.ac.sanger.artemis.util.OutOfRangeException e) 
        {
          getSelection().clear();
        }
      }

      /**
       *  Called when the user clicks somewhere on the plot canvas.
       *  @param position the base/amino acid position of the click.  This is
       *    -1 if and only if the click was outside the graph (eg. in the
       *    label at the top)
       **/
      public void mouseClick(final int position) 
      {
      }

      /**
       *  Called when the user drags the mouse over the plot.
       *  @param drag_start_position The base/amnino acid position where the
       *    drag started or -1 if the drag was started outside the graph.
       *  @param current_position the base/amino acid position of the click.
       *    This is -1 if and only if the user has dragged the mouse out of
       *    the graph (eg. in the label at the top)
       **/
      public void mouseDrag(int drag_start_position,
                            int current_position) 
      {
        if(rev_comp_display) 
        {
          drag_start_position =
            bases.getComplementPosition(drag_start_position);
          current_position =
            bases.getComplementPosition(current_position);
        }
        setSelectionRange(drag_start_position,
                          current_position);
      }

      /**
       *  Called when the user double-clicks somewhere on the plot.
       *  @param position the base/amino acid position of the click.  This is
       *    -1 if and only if the click was outside the graph (eg. in the
       *    label at the top)
       **/
      public void mouseDoubleClick(int position)
      {
        if(rev_comp_display) 
          position = bases.getComplementPosition(position);
        
        setSelectionRange(position, position);
        getGotoEventSource().gotoBase(position);
      }
    });
  }

  /**
   *  Used by getPreferredSize() and getMinimumSize();
   **/
  private final static int HEIGHT;

  static 
  {
    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();

    final Integer base_plot_height = 
      Options.getOptions().getIntegerProperty("base_plot_height");

    if(base_plot_height == null) 
    {
      if(screen.height <= 600) 
        HEIGHT = 100;
      else 
        HEIGHT = 150;
    }
    else
      HEIGHT = base_plot_height.intValue();
  }

  /**
   *  Overridden to set the component height to 150.
   **/
  public Dimension getPreferredSize() 
  {
    return(new Dimension(getSize().width, HEIGHT));
  }

  /**
   *  Overridden to set the component height to 150.
   **/
  public Dimension getMinimumSize() 
  {
    return (new Dimension(getSize().width, HEIGHT));
  }

  /**
   *  Implementation of the DisplayAdjustmentListener interface.  Invoked when
   *  a component or changes the scale.
   **/
  public void displayAdjustmentValueChanged(DisplayAdjustmentEvent event) 
  {
    start_base = event.getStart();
    end_base = event.getEnd();
    width_in_bases = event.getWidthInBases();
    rev_comp_display = event.isRevCompDisplay();
    recalculate_flag = true;

    if(event.getType() == DisplayAdjustmentEvent.ALL_CHANGE_ADJUST_EVENT) 
    {
      selection_start_marker = null;
      selection_end_marker = null;
    }

    repaint();
  }

  /**
   *  We listen for SelectionChange events so that we can update the
   *  crosshairs.
   **/
  public void selectionChanged(SelectionChangeEvent event) 
  {
    selection_start_marker = null;
    selection_end_marker = null;
    repaint();
  }

  /**
   *  Return the algorithm that was passed to the constructor.
   **/
  public BaseAlgorithm getBaseAlgorithm()
  {
    return (BaseAlgorithm)super.getAlgorithm();
  }

  /**
   *  Return the new start base to display, from the last event.
   **/
  private int getStart()
  {
    return start_base;
  }

  /**
   *  Return the new end base to display, from the last event.
   **/
  private int getEnd() 
  {
    return end_base;
  }

  /**
   *  Return the width in bases of the display, from the last event.
   **/
  private int getWidthInBases() 
  {
    return width_in_bases;
  }

  /**
   *  This array is used by drawMultiValueGraph().  It is reallocated when
   *  the scale changes.
   **/
  private float[][] value_array_array = null;

  /**
   *  The number of bases to step before each evaluation of the algorithm.
   *  (Set by recalculateValues()).
   **/
  private int step_size = 0;

  /**
   *  The maximum of the values in value_array_array.
   **/
  private float min_value = Float.MAX_VALUE;

  /**
   *  The minimum of the values in value_array_array.
   **/
  private float max_value = Float.MIN_VALUE;

  /**
   *  Recalculate the values in value_array_array, step_size, min_value and
   *  max_value.
   **/
  protected void recalculateValues()
  {
    final Float algorithm_minimum = getAlgorithm().getMinimum();
    final Float algorithm_maximum = getAlgorithm().getMaximum();

    // use the Algorithm specified maximum if there is one - otherwise
    // calculate it
    if(algorithm_maximum == null) 
      max_value = Float.MIN_VALUE;
    else
      max_value = algorithm_maximum.floatValue();

    // use the Algorithm specified minimum if there is one - otherwise
    // calculate it
    if(algorithm_minimum == null) 
      min_value = Float.MAX_VALUE;
    else
      min_value = algorithm_minimum.floatValue();

    final int window_size = getWindowSize();

    final Integer default_step_size =
      getAlgorithm().getDefaultStepSize(window_size);

    if(default_step_size == null) 
      step_size = 1;
    else
    {
      if(default_step_size.intValue() < window_size) 
        step_size = default_step_size.intValue();
      else
        step_size = window_size;
    }

    int real_start = getStart();

    if(real_start < 1)
      real_start = 1;

    final int unit_count = getEnd() - real_start;

    // the number of plot points in the graph
    final int number_of_values =
      (unit_count - (window_size - step_size)) / step_size;

    if(number_of_values < 2) 
    {
      // there is nothing to plot
      value_array_array = null;
      return;
    }

    getBaseAlgorithm().setRevCompDisplay(rev_comp_display);

    // the number of values that getValues() will return
    final int get_values_return_count =
      getBaseAlgorithm().getValueCount();

    if(value_array_array == null) 
      value_array_array = new float [get_values_return_count][];

    if(value_array_array[0] == null ||
       value_array_array[0].length != number_of_values) 
    {
      for(int i = 0 ; i < value_array_array.length ; ++i) 
        value_array_array[i] = new float [number_of_values];
    }
    else 
    {
      // reuse the previous arrays
    }

    float [] temp_values = new float [get_values_return_count];

    for(int i = 0 ; i < number_of_values ; ++i) 
    {
      getBaseAlgorithm().getValues(real_start + i * step_size,
                                   real_start + i * step_size +
                                   window_size - 1,
                                   temp_values);

      for(int value_index = 0 ;
          value_index < get_values_return_count ;
          ++value_index) 
      {
        final float current_value = temp_values[value_index];

        value_array_array[value_index][i] = current_value;

        // use the Algorithm specified maximum if there is one - otherwise
        // calculate it
        if(algorithm_maximum == null)
        {
          if (current_value > max_value) 
            max_value = current_value;
        }

        // use the Algorithm specified minimum if there is one - otherwise
        // calculate it
        if(algorithm_minimum == null) 
        {
          if(current_value < min_value) 
            min_value = current_value;
        }
      }
    }

    recalculate_flag = false;
  }

  /**
   *  Redraw the graph on the canvas using the algorithm, start_base and
   *  end_base.  This method plots BaseWindowAlgorithm objects only.
   *  @param g The object to draw into.
   **/
  public void drawMultiValueGraph(Graphics g) 
  {
    if(recalculate_flag)
      recalculateValues();

    if(value_array_array == null) 
    {
      // there is nothing to draw - probably because the sequence is too short
      drawMinMax(g, 0, 1);
      return;
    }

    final int window_size = getWindowSize();

    // the number of values to plot at each x position
    final int get_values_return_count =
      getBaseAlgorithm().getValueCount();

    final int number_of_values = value_array_array[0].length;

    if(number_of_values > 1) 
      drawGlobalAverage(g, min_value, max_value);

    for(int value_index = 0; value_index < get_values_return_count;
        ++value_index)
    {
      if(get_values_return_count == 1) 
        g.setColor(Color.black);
      else 
      {
        switch(value_index) 
        {
          case 0:
            g.setColor (new Color (255, 0, 0));
            break;
          case 1:
            g.setColor (new Color (0, 200, 0));
            break;
          case 2:
            g.setColor (new Color (0, 0, 255));
            break;
          default:
            g.setColor (Color.black);
        }
        if(rev_comp_display)
          System.out.println("TRUE "+value_index+" 0 = R, 1 = G, 2 = B");
        else
          System.out.println("FALSE "+value_index+" 0 = R, 1 = G, 2 = B");
      }

      final int offset;

      if(getStart() < 1) 
        offset = 1 - getStart();
      else
        offset = 0;

      drawPoints(g, min_value, max_value, step_size, window_size,
                 getWidthInBases(),
                 offset,
                 value_array_array[value_index]);
    }

    drawMinMax(g, min_value, max_value);

    if(getCrossHairPosition() >= 0)
    {
      final int cross_hair_position = getCrossHairPosition();
      final int selection_base = getPointPosition(cross_hair_position);

      if(selection_base >= 1)
      {
        if(selection_base > end_base) 
          cancelCrossHairs();
        else 
        {
          final String label_string;

          if(rev_comp_display) 
          {
            label_string =
              String.valueOf(bases.getLength() - selection_base + 1);
          } 
          else
            label_string = String.valueOf(selection_base);

          drawCrossHair(g, cross_hair_position, label_string, 0);
        }
      }
    }

    if(getCrossHairPosition() >= 0 && getSelectionStartMarker() != null) 
    {
      final int selection_first_base =
        getSelectionStartMarker().getRawPosition();

      final String label_string;

      if(rev_comp_display) 
      {
        label_string =
          String.valueOf(bases.getLength() - selection_first_base + 1);
      } 
      else
        label_string = String.valueOf(selection_first_base);

      if(Math.abs(selection_first_base -
                  getPointPosition(getCrossHairPosition())) > 3) 
      {
        drawCrossHair(g, getCanvasPosition(selection_first_base),
                      label_string, 1);
      } 
      else 
      {
        // don't draw - too close to main cross hair
      }
    }

    if(getCrossHairPosition() >= 0 && getSelectionEndMarker() != null)
    {
      if(getSelectionStartMarker() != null &&
         Math.abs((getSelectionEndMarker().getRawPosition() -
                   getSelectionStartMarker().getRawPosition())) >= 3)
      {
        final int selection_last_base =
          getSelectionEndMarker().getRawPosition();

        final String label_string;

        if(rev_comp_display)
        {
          label_string =
            String.valueOf(bases.getLength() - selection_last_base + 1);
        } 
        else 
          label_string = String.valueOf(selection_last_base);

        if(Math.abs(selection_last_base -
                    getPointPosition(getCrossHairPosition())) > 3)
        {
          drawCrossHair(g, getCanvasPosition (selection_last_base),
                        label_string, 2);
        }
        else
        {
          // don't draw - too close to main cross hair
        }
      }
    }
  }

  /**
   *  Get the position in the Feature of the given canvas x position.  This
   *  base position is the label used when the user clicks the mouse in on the
   *  canvas (see drawCrossHair()).
   **/
  protected int getPointPosition(final int canvas_x_position) 
  {
    return (int)((1.0 * canvas_x_position / getSize().width) *
                 getWidthInBases()) + getStart();
  }

  /**
   *  Return the canvas position of the given base,
   **/
  private int getCanvasPosition(final int base)
  {
    return (int)((1.0 * base - getStart()) / getWidthInBases() *
                 getSize().width);
  }

  /**
   *  Return the Marker of the start base of the Selection or null if there is
   *  nothing selected.
   **/
  private Marker getSelectionStartMarker()
  {
    if(selection_start_marker == null) 
    {
      selection_start_marker = getSelection().getLowestBaseOfSelection();

      if(selection_start_marker != null &&
         rev_comp_display) 
      {
        final Strand strand = bases.getReverseStrand();
        final int orig_position = selection_start_marker.getRawPosition();
        final int rev_comp_position =
          bases.getComplementPosition(orig_position);
        try 
        {
          selection_start_marker = strand.makeMarker(orig_position);
        }
        catch(OutOfRangeException e) 
        {
          throw new Error("internal error - unexpected exception: " + e);
        }
      }
    }

    return selection_start_marker;
  }

  /**
   *  Return the Marker of the end base of the Selection or null if there is
   *  nothing selected.
   **/
  private Marker getSelectionEndMarker() 
  {
    if(selection_end_marker == null) 
    {
      selection_end_marker = getSelection().getHighestBaseOfSelection();

      if(selection_end_marker != null &&
         rev_comp_display) 
      {
        final Strand strand = bases.getReverseStrand();
        final int orig_position = selection_end_marker.getRawPosition();
        final int rev_comp_position =
          bases.getComplementPosition(orig_position);
        try 
        {
          selection_end_marker = strand.makeMarker(orig_position);
        }
        catch(OutOfRangeException e) 
        {
          throw new Error("internal error - unexpected exception: " + e);
        }
      }
    }

    return selection_end_marker;
  }

  /**
   *  Return the Selection object that was passed to the constructor.
   **/
  private Selection getSelection() 
  {
    return selection;
  }

  /**
   *  Return the GotoEventSource object that was passed to the constructor.
   **/
  private GotoEventSource getGotoEventSource()
  {
    return goto_event_source;
  }

  /**
   *  The start base to plot, as obtained from the DisplayAdjustmentEvent.
   **/
  private int start_base;

  /**
   *  The end base to plot, as obtained from the DisplayAdjustmentEvent.
   **/
  private int end_base;

  /**
   *  The width in bases of the display, as obtained from the
   *  DisplayAdjustmentEvent.
   **/
  private int width_in_bases;

  /**
   *  True if and only if the FeatureDisplay is drawing in reverse complement
   *  mode.
   **/
  private boolean rev_comp_display;

  /**
   *  The Bases that this BasePlot is graphing.
   **/
  private Bases bases;

  /**
   *  The Marker of the start of the selection.  This is a cache used by
   *  getSelectionStartMarker().
   **/
  private Marker selection_start_marker = null;

  /**
   *  The Marker of the end of the selection.  This is a cache used by
   *  getSelectionEndMarker()
   **/
  private Marker selection_end_marker = null;

  /**
   *  The Selection that was passed to the constructor.
   **/
  private Selection selection;

  /**
   *  The GotoEventSource that was passed to the constructor.
   **/
  private GotoEventSource goto_event_source;
}
