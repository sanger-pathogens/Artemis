/* CumulativeATSkewAlgorithm.java
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 2004  Genome Research Limited
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

package uk.ac.sanger.artemis.plot;

import uk.ac.sanger.artemis.util.*;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.sequence.*;

/**
 *
 * Grigoriev A (1999)  Strand-specific compositional asymmetries in 
 * double-stranded DNA viruses.
 * Virus Research 60, 1-19.
 *
 *  @author Derek Gatherer
 **/

public class CumulativeATSkewAlgorithm extends BaseAlgorithm 
{
  /**
   *  Create a new CumulativeATSkewAlgorithm object
   *  @param strand The strand to do the calculation on.
   **/
  public CumulativeATSkewAlgorithm(final Strand strand)
  {
    super(strand, "Cumulative AT Skew, (A-T)/(A+T)", "at_skew");
    setScalingFlag(true);
  }

  /**
   *  @param start The start base (included in the range).
   *  @param end The end base (included in the range).
   *  @param values The one return value for this algorithm is returned in
   *    this array.
   **/
  public void getValues(int start, int end, final float [] values) 
  {
    final String sequence;
    String subseq;

// sliding window code here

    int leap = end-start;    
    values[0] = 0;
 
    float a_count = 0;
    float t_count = 0;
   
    for(int window = 0 ; window < end ; window += leap) 
    {
      try
      {
        subseq = getStrand().getSubSequence(new Range(window, window+leap));
      }
      catch(OutOfRangeException e) 
      {
        throw new Error("internal error - unexpected exception: " + e);
      }

      a_count = 0;
      t_count = 0;

      for(int i = 0 ; i < subseq.length() ; ++i) 
      {
        final char this_char = subseq.charAt(i);

        if(this_char == 'a')
          ++a_count;

        if(this_char == 't') 
          ++t_count;
      }

      if(a_count + t_count > 0) 
        values[0] += (a_count - t_count) / (a_count + t_count); 
    }
  }

  /**
   *  Return the number of values a call to getValues() will return - one
   *  in this case.
   **/
  public int getValueCount() 
  {
    return 1;
  }

  /**
   *  Return the default or optimal window size.
   *  @return null is returned if this algorithm doesn't have optimal window
   *    size.
   **/
  public Integer getDefaultWindowSize() 
  {
    final Integer super_window_size = super.getDefaultWindowSize();
    if(super_window_size != null) 
    {
      // the superclass version of getDefaultWindowSize() returns non-null
      // iff the user has set the window size in the options file
      return super_window_size;
    }
    return new Integer(60);
  }

  /**
   *  Return the default maximum window size for this algorithm.
   *  @return null is returned if this algorithm doesn't have maximum window
   *    size.
   **/
  public Integer getDefaultMaxWindowSize() 
  {
    final Integer super_max_window_size = super.getDefaultMaxWindowSize();
    if(super_max_window_size != null) 
    {
      // the superclass version of getDefaultMaxWindowSize() returns non-null
      // iff the user has set the max window size in the options file
      return super_max_window_size;
    }
    return new Integer(5000);
  }

  /**
   *  Return the default minimum window size for this algorithm.
   *  @return null is returned if this algorithm doesn't have minimum window
   *    size.
   **/
  public Integer getDefaultMinWindowSize() 
  {
    final Integer super_min_window_size = super.getDefaultMinWindowSize();
    if(super_min_window_size != null) 
    {
      // the superclass version of getDefaultMinWindowSize() returns non-null
      // iff the user has set the minimum window size in the options file
      return super_min_window_size;
    }
    return new Integer(10);
  }

  /**
   *  Return the default or optimal step size.
   *  @return null is returned if this algorithm doesn't have optimal step
   *    size.
   **/
  public Integer getDefaultStepSize(int window_size) 
  {
    if(window_size > 10) 
      return new Integer(window_size / 10);
    else
      return null;
  }

  /**
   *  Return the maximum value of this algorithm.
   *  @return The maximum is 100.
   **/
  protected Float getMaximumInternal() 
  {
    return new Float(1);
  }

  /**
   *  Return the minimum value of this algorithm.
   *  @return The minimum is 0.
   **/
  protected Float getMinimumInternal() 
  {
    return new Float(-1);
  }

  /**
   *  Return the average value of function over the whole strand.
   *  @return null is returned if this algorithm doesn't have an average or if
   *    the average can't be calculated.
   **/
  public Float getAverage() 
  {
    final float a_count = getStrand().getACount();
    final float t_count = getStrand().getTCount();

    final float a_minus_t = a_count - t_count;
    final float a_plus_t = a_count + t_count;

    if(a_plus_t > 0) 
      return new Float(a_minus_t/a_plus_t);
    else 
      return new Float(0.0);
  }
}
