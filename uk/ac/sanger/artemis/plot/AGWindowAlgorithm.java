/* AGWindowAlgorithm.java
 *
 * created: Wed Apr 14 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/plot/AGWindowAlgorithm.java,v 1.4 2009-03-17 17:47:42 tjc Exp $
 */

package uk.ac.sanger.artemis.plot;

import uk.ac.sanger.artemis.util.*;
import uk.ac.sanger.artemis.io.Range;

import uk.ac.sanger.artemis.sequence.*;

/**
 *  Objects of this class have one useful method - getValues (), which takes a
 *  range of bases and returns a single floating point number, which is the
 *  percent AG content of the range.  The Strand to use is set in the
 *  constructor.
 *
 *  @author Kim Rutherford
 *  @version $Id: AGWindowAlgorithm.java,v 1.4 2009-03-17 17:47:42 tjc Exp $
 **/

public class AGWindowAlgorithm extends BaseAlgorithm {
  /**
   *  Create a new GCWindowAlgorithm object.
   *  @param strand The strand to do the calculation on.
   **/
  public AGWindowAlgorithm (final Strand strand) {
    super (strand, "AG Content (%)", "ag_content");

    setScalingFlag (true);
  }

  /**
   *  Return the percent GC between a pair of bases.
   *  @param start The start base (included in the range).
   *  @param end The end base (included in the range).
   *  @param values The one return value for this algorithm is returned in
   *    this array.
   **/
  public void getValues (int start, int end, final float [] values) {
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

      if (this_char == 'g' || this_char == 'a') {
        ++gc_count;
      }
    }

//    System.out.println ("start: " + start + " end: " + end + " returning: " + gc_count/sequence.length ()); 

    values[0] = gc_count/sequence.length () * 100;
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
    return new Integer (120);
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
    return new Integer (500);
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
    return new Integer (24);
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
   *  Return the default or optimal step size.
   *  @return null is returned if this algorithm doesn't have optimal step
   *    size.
   **/
  public Integer getDefaultStepSize (int window_size) {
    if (window_size > 10) {
      return new Integer (window_size / 10);
    } else {
      return null;
    }
  } 

  /**
   *  Return the average value of function over the whole strand.
   *  @return null is returned if this algorithm doesn't have an average or if
   *    the average can't be calculated.
   **/
  public Float getAverage () {
    return new Float (getStrand ().getBases ().getAverageAGPercent ());
  }
}
