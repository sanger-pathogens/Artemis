/* FeaturePlot.java
 *
 * created: Wed Dec 16 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/FeaturePlot.java,v 1.11 2009-07-20 15:11:17 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;

import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureChangeEvent;
import uk.ac.sanger.artemis.FeatureChangeListener;
import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.plot.FeatureAlgorithm;
import uk.ac.sanger.artemis.plot.LineAttributes;
import uk.ac.sanger.artemis.sequence.AminoAcidSequence;


/**
 *  The components of this class display a plot of a FeatureAlgorithm for a
 *  particular feature.
 *  @author Kim Rutherford
 **/

public class FeaturePlot extends Plot
    implements DisplayAdjustmentListener, FeatureChangeListener {

  private static final long serialVersionUID = 1L;

  /**
   *  Create a new FeatureDisplay object.
   *  @param algorithm The object that will generate the values we plot in
   *    this component.
   **/
  public FeaturePlot (FeatureAlgorithm algorithm) {
    super (algorithm, true);    // true means draw a scale at the bottom of
                                // the graph

    getFeature ().addFeatureChangeListener (this);

    recalculate_flag = true;
  }

  /**
   *  Used by getPreferredSize () and getMinimumSize ();
   **/
  protected static int HEIGHT;

  static {
    final Integer feature_plot_height =
      Options.getOptions ().getIntegerProperty ("feature_plot_height");

    if (feature_plot_height == null) {
      HEIGHT = 160;
    } else {
      HEIGHT = feature_plot_height.intValue ();
    }
  }

  /**
   *  Overridden to set the component height to 150.
   **/
  public Dimension getPreferredSize() {
    return (new Dimension(getSize ().width, HEIGHT));
  }
  
  /**
   *  Overridden to set the component height to 150.
   **/
  public Dimension getMinimumSize() {
    return (new Dimension(getSize ().width, HEIGHT));
  }
  
  /**
   *  Remove this object as a feature change listener.  (Called by
   *  FeaturePlotGroup)
   **/
  void stopListening () {
    getFeature ().removeFeatureChangeListener (this);
  }

  /**
   *  Implementation of the DisplayAdjustmentListener interface.  Invoked when
   *  a component scrolls or changes the scale.
   **/
  public void displayAdjustmentValueChanged (DisplayAdjustmentEvent event) {
    start_base = event.getStart ();
    end_base = event.getEnd ();
    //width_in_bases = event.getWidthInBases ();

    recalculate_flag = true;
    repaint();
  }

  /**
   *  Implementation of the FeatureChangeListener interface.
   *  @param event The change event.
   **/
  public void featureChanged (FeatureChangeEvent event) {
    recalculate_flag = true;
  }

  /**
   *  Return the algorithm that was passed to the constructor.
   **/
  public FeatureAlgorithm getFeatureAlgorithm () {
    return (FeatureAlgorithm) super.getAlgorithm ();
  }

  /**
   *  Return the new start base to display, from the last event.
   **/
  private int getStart () {
    return start_base;
  }

  /**
   *  Return the new end base to display, from the last event.
   **/
  private int getEnd () {
    return end_base;
  }

  /**
   *  This array is used by drawGraph ().  It is reallocated when the scale
   *  changes.
   **/
  private float [][] value_array_array = null;

  /**
   *  The number of bases to step before each evaluation of the algorithm.
   *  (Set by recalculateValues ()).
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
  protected void recalculateValues () {
    final Float algorithm_minimum = getAlgorithm ().getMinimum ();
    final Float algorithm_maximum = getAlgorithm ().getMaximum ();

    // use the Algorithm specified maximum if there is one - otherwise
    // calculate it
    if (algorithm_maximum == null) {
      max_value = Float.MIN_VALUE;
    } else {
      max_value = algorithm_maximum.floatValue ();
    }

    // use the Algorithm specified minimum if there is one - otherwise
    // calculate it
    if (algorithm_minimum == null) {
      min_value = Float.MAX_VALUE;
    } else {
      min_value = algorithm_minimum.floatValue ();
    }

    final int window_size = getWindowSize ();

    final Integer default_step_size =
      getAlgorithm ().getDefaultStepSize (window_size);

    if (default_step_size == null) {
      step_size = window_size;
    } else {
      if (default_step_size.intValue () < window_size) {
        step_size = default_step_size.intValue ();
      } else {
        step_size = window_size;
      }
    }

    final int unit_count = getEnd () - getStart ();

    // the number of plot points in the graph
    final int number_of_values =
      (unit_count - (window_size - step_size)) / step_size;

    if (number_of_values < 2) {
      // there is nothing to plot
      value_array_array = null;
      return;
    }

    // the number of values that getValues () will return
    final int get_values_return_count =
      getFeatureAlgorithm ().getValueCount ();

    if (value_array_array == null) {
      value_array_array = new float [get_values_return_count][];
    }

    if (value_array_array[0] == null ||
        value_array_array[0].length != number_of_values) {
      for (int i = 0 ; i < value_array_array.length ; ++i) {
        value_array_array[i] = new float [number_of_values];
      }
    } else {
      // reuse the previous arrays
    }

    if (!isVisible ()) {
      return;
    }

    float [] temp_values = new float [get_values_return_count];

    for (int i = 0 ; i < number_of_values ; ++i) {
      getFeatureAlgorithm ().getValues (getStart () + i * step_size,
                                        getStart () + i * step_size +
                                        window_size - 1,
                                        temp_values);

      for (int value_index = 0 ;
           value_index < get_values_return_count ;
           ++value_index) {
        final float current_value = temp_values[value_index];

        value_array_array[value_index][i] = current_value;

        // use the Algorithm specified maximum if there is one - otherwise
        // calculate it
        if (algorithm_maximum == null) {
          if (current_value > max_value) {
            max_value = current_value;
          }
        }

        // use the Algorithm specified minimum if there is one - otherwise
        // calculate it
        if (algorithm_minimum == null) {
          if (current_value < min_value) {
            min_value = current_value;
          }
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
  protected int drawMultiValueGraph (Graphics g, LineAttributes[] lines) {
    if (recalculate_flag) {
      recalculateValues ();
    }
    
    if (value_array_array == null) {
      // there is nothing to draw - probably because the sequence is too short
      drawMinMax (g, 0, 1);
      
      return 0;
    }
    
    final int window_size = getWindowSize ();

    // the number of values to plot at each x position
    final int get_values_return_count =
      getFeatureAlgorithm ().getValueCount ();

    // the number of x positions to plot points at
    final int number_of_values = value_array_array[0].length;

    if (number_of_values > 1) {
      drawGlobalAverage (g, min_value, max_value);
    }

    for (int value_index = 0 ;
         value_index < get_values_return_count ;
         ++value_index) {
      if (get_values_return_count == 1) {
        g.setColor (Color.black);
      } else {
        switch (value_index) {
        case 0:
          g.setColor (new Color (255, 0, 0));
          break;
        case 1:
          g.setColor (new Color (100, 255, 100));
          break;
        case 2:
          g.setColor (new Color (0, 0, 255));
          break;
        default:
          g.setColor (Color.black);
        }
      }

      drawPoints (g, min_value, max_value, step_size, window_size,
                  getSize ().width,
                  0, // no offset.
                  value_array_array[value_index], value_index, 
                  get_values_return_count, false, false);
    }

    drawMinMax (g, min_value, max_value);

    drawScaleLine (g, getStart (), getEnd ());

    final int cross_hair_position = getCrossHairPosition ();

    if (cross_hair_position >= 0) {
      if (cross_hair_position >= getTranslation ().length ()) {
        cancelCrossHairs ();
      } else {
        drawCrossHair (g, cross_hair_position,
                       String.valueOf (getPointPosition (cross_hair_position +
                                                         getStart () - 1)),
                       0);
      }
    }

    return get_values_return_count;
  }

  /**
   *  Get the position in the Feature of the given canvas x position.  This
   *  amino acid position is the label used when the user clicks the mouse in
   *  on the canvas (see drawCrossHair ()).
   **/
  protected int getPointPosition (final int canvas_x_position) {
    return canvas_x_position + 1;
  }

  /**
   *  Return the feature that this component is plotting.
   **/
  private Feature getFeature () {
    return getFeatureAlgorithm ().getFeature ();
  }

  /**
   *  Return the translation of the bases of the feature we are plotting.
   **/
  private AminoAcidSequence getTranslation () {
    return getFeature ().getTranslation ();
  }

  /**
   *  The start base to plot, as obtained from the DisplayAdjustmentEvent.
   **/
  private int start_base;

  /**
   *  The end base to plot, as obtained from the DisplayAdjustmentEvent.
   **/
  private int end_base;

  protected void calculateFeatures(boolean fromPeak)
  {
    // TODO Auto-generated method stub
  }

  protected void showAveragesForRange() {}
}
