/* GCSDWindowAlgorithm.java
 *
 * created: Mon Oct 25 1999
 *
 * This file is part of Artemis
 *
 * Copyright (C) 1999  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/plot/GCSDWindowAlgorithm.java,v 1.4 2009-03-17 17:47:42 tjc Exp $
 */

package uk.ac.sanger.artemis.plot;

import uk.ac.sanger.artemis.util.*;
import uk.ac.sanger.artemis.io.Range;

import uk.ac.sanger.artemis.sequence.*;

/**
 *  Objects of this class have one useful method - getValues (), which takes a
 *  range of bases and returns a single floating point number, which is the
 *  percent GC content of the range if the GC content of the range is more
 *  than 2.5 standard deviations from the average GC content of the sequence.
 *  It returns the average GC content of the sequence otherwise.  The Strand
 *  to use is set in the constructor.
 *
 *  @author Kim Rutherford
 *  @version $Id: GCSDWindowAlgorithm.java,v 1.4 2009-03-17 17:47:42 tjc Exp $
 **/

public class GCSDWindowAlgorithm extends BaseAlgorithm {
  /**
   *  Create a new GCSDWindowAlgorithm object.
   *  @param strand The strand to do the calculation on.
   **/
  public GCSDWindowAlgorithm (final Strand strand) {
    super (strand, "GC Content (%) With A 2.5 SD Cutoff", "sd_gc_content");

    setScalingFlag (false);
  }

  /**
   *  Return the percent GC between a pair of bases if the GC content of the
   *  range is more than 2.5 standard deviations from the average GC content
   *  of the sequence.  It returns the average GC content of the sequence
   *  otherwise.
   *  @param start The start base (included in the range).
   *  @param end The end base (included in the range).
   *  @param values The one return value for this algorithm is returned in
   *    this array.
   **/
  public void getValues (int start, int end, float [] values) {
    final int window_size = end - start + 1;

    final float standard_deviation;

    if (window_size > getDefaultMaxWindowSize ().intValue ()) {
      standard_deviation = calculateSD (window_size);
    } else {
      if (standard_deviations[window_size - 1] < 0) {
        // set the cached value
        standard_deviation = calculateSD (window_size);
        standard_deviations[window_size - 1] = standard_deviation;

        //       System.err.println ("SD: " + standard_deviation);
      } else {
        standard_deviation = standard_deviations[window_size - 1];
      }
    }

    final String sequence;

    try {
      sequence = getStrand ().getSubSequence (new Range (start, end));
    } catch (OutOfRangeException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    }

    float gc_count = 0;

    for (int i = 0 ; i < sequence.length () ; ++i) {
      final char this_char = sequence.charAt (i);
//      System.out.println (this_char);

      if (this_char == 'g' || this_char == 'c') {
        ++gc_count;
      }
    }

    final float gc_content = gc_count/sequence.length () * 100;

    final float gc_average =
      getStrand ().getBases ().getAverageGCPercent ();

    if (Math.abs (gc_content - gc_average) < standard_deviation * 2.5) {
      values[0] = gc_average;
    } else {
      values[0] = gc_content;
    }
  }

  /**
   *  Calculate and return the standard deviation of the GC content of the
   *  Bases object of the Strand that was passed to the constructor.
   **/
  private float calculateSD (final int window_size) {
//      System.err.println ("calculateSD:");

    final int sequence_length = getStrand ().getBases ().getLength ();

    final String bases = getStrand ().getBases ().toString ();

    // the number of windows to search
    final int window_count = sequence_length - window_size;

//      System.err.println ("window_size: " + window_size);
//      System.err.println ("window_count: " + window_count);

    // this is the sum over all the windows of:
    // (GC content) * (percent GC) / window_size * window_size
    double sum_so_far = 0;

    // this is the sum over all the windows of (GC content)/window_size
    double gc_so_far = 0;

    double current_gc_count = 0;

    for (int i = - window_size ; i < window_count ; ++i) {
      if (i > 0) {
        final char previous_char = bases.charAt (i - 1);

        if (previous_char == 'g' || previous_char == 'c') {
          --current_gc_count;
        }
      }

      final char new_char = bases.charAt (i + window_size);

      if (new_char == 'g' || new_char == 'c') {
        ++current_gc_count;
      }

      if (i >= 0) {
        final double this_value = current_gc_count / window_size;

        sum_so_far += this_value * this_value;

        gc_so_far += this_value;
      }
    }

//      System.err.println ("sum_so_far: " + sum_so_far);
//      System.err.println ("gc_so_far: " + gc_so_far);

    final double gc_average = gc_so_far / window_count;

//      System.err.println ("gc_average: " + gc_average);
//      System.err.println ("sum_so_far/window_count: " +
//                          (sum_so_far/window_count));
//      System.err.println ("gc_average*gc_average: " + (gc_average * gc_average));

//      System.err.println ("sd*sd: " + (sum_so_far / window_count -
//                                       gc_average * gc_average));

    return (float) Math.sqrt (sum_so_far / window_count -
                              gc_average * gc_average) * 100;
  }

  /**
   *  Return the number of values a call to getValues () will return - one
   *  in this case.
   **/
  public int getValueCount () {
    return 1;
  }

  /**
   *  Return the default or optimal window size.
   *  @return null is returned if this algorithm doesn't have optimal window
   *    size.
   **/
  public Integer getDefaultWindowSize () {
    final Integer super_window_size = super.getDefaultWindowSize ();
    if (super_window_size != null) {
      // the superclass version of getDefaultWindowSize () returns non-null
      // iff the user has set the window size in the options file
      return super_window_size;
    }
    return new Integer (1000);
  }

  /**
   *  Return the default maximum window size for this algorithm.
   *  @return null is returned if this algorithm doesn't have maximum window
   *    size.
   **/
  public Integer getDefaultMaxWindowSize () {
    final Integer super_max_window_size = super.getDefaultMaxWindowSize ();
    if (super_max_window_size != null) {
      // the superclass version of getDefaultMaxWindowSize () returns non-null
      // iff the user has set the max window size in the options file
      return super_max_window_size;
    }
    return new Integer (5000);
  }

  /**
   *  Return the default minimum window size for this algorithm.
   *  @return null is returned if this algorithm doesn't have minimum window
   *    size.
   **/
  public Integer getDefaultMinWindowSize () {
    final Integer super_min_window_size = super.getDefaultMinWindowSize ();
    if (super_min_window_size != null) {
      // the superclass version of getDefaultMinWindowSize () returns non-null
      // iff the user has set the minimum window size in the options file
      return super_min_window_size;
    }
    return new Integer (100);
  }

  /**
   *  Return the default or optimal step size.
   *  @return null is returned if this algorithm doesn't have optimal step
   *    size.
   **/
  public Integer getDefaultStepSize (int window_size) {
    return new Integer (1);
  }

  /**
   *  Return the average value of function over the whole strand.
   *  @return null is returned if this algorithm doesn't have an average or if
   *    the average can't be calculated.
   **/
  public Float getAverage () {
    return new Float (getStrand ().getBases ().getAverageGCPercent ());
  }

  /**
   *  Return the maximum value of this algorithm.
   *  @return The maximum is 100.
   **/
  protected Float getMaximumInternal () {
    return new Float (100);
  }

  /**
   *  Return the minimum value of this algorithm.
   *  @return The minimum is 0.
   **/
  protected Float getMinimumInternal () {
    return new Float (0);
  }

  /**
   *  This array contains a cache of the standard deviations of the GC content
   *  of the Bases object of the Strand that was passed to the constructor.
   *  (Set by the constructor).  The array is indexed by window size.
   **/
  float [] standard_deviations;

  {
    final int max_window_size = getDefaultMaxWindowSize ().intValue ();
    standard_deviations = new float [max_window_size];

    for (int i = 0 ; i < max_window_size ; ++i) {
      standard_deviations[i] = -1;
    }
  }
}
