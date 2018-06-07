/* GCFrameAlgorithm.java
 *
 * created: Tue Dec 15 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/plot/GCFrameAlgorithm.java,v 1.8 2009-06-26 15:52:48 tjc Exp $
 */

package uk.ac.sanger.artemis.plot;

import uk.ac.sanger.artemis.sequence.*;
import uk.ac.sanger.artemis.util.*;
import uk.ac.sanger.artemis.io.Range;

import java.awt.*;

/**
 *  Objects of this class have one useful method - getValues (), which takes a
 *  range of bases and returns three floating point numbers, which are the
 *  percent GC content in each frame.  The Strand to use is set in the
 *  constructor.
 *
 *  @author Kim Rutherford
 *  @version $Id: GCFrameAlgorithm.java,v 1.8 2009-06-26 15:52:48 tjc Exp $
 **/
public class GCFrameAlgorithm extends BaseAlgorithm 
{
  /**
   *  Create a new GCFrameAlgorithm object.
   *  @param strand The Strand to do the calculation on.
   **/
  public GCFrameAlgorithm(final Strand strand) 
  {
    super(strand, makeName(strand), "gc_frame");
    setScalingFlag(true);
  }

  /**
   *  This is used as temporary storage by getValues().
   **/
  private int gc_counts [] = new int [getValueCount()];

  /**
   *  Return the percent gc between a pair of bases in each of the three
   *  frames.
   *  @param start The start base (included in the range).
   *  @param end The end base (included in the range).  If the start/end pair
   *    doesn't give a multiple of three bases end is moved down so that it is
   *    a multiple of three.
   *  @param values The three results are returned in this array, hence it
   *    should be three long.  The first value is the gc content of positions
   *    start, start+3, start+6, etc, the second value is start+1, start+4,
   *    etc.
   **/
  public void getValues(int start, int end, final float [] values) 
  {
    if(isRevCompDisplay())
    {
      final int new_end =
        getStrand().getBases().getComplementPosition(start);
      final int new_start =
        getStrand().getBases().getComplementPosition(end);

      end = new_end;
      start = new_start;
    }

    // add 1 or 2 if necessary to make the range a multiple of 3
    if(getStrand().isForwardStrand())
      end -= (end - start + 1) % 3;
    else
      start += (end - start + 1) % 3;

    for(int i = 0; i < getValueCount(); ++i)
      gc_counts[i] = 0;
    
    final char[] sub_sequence;

    try 
    {
      sub_sequence = getStrand().getRawSubSequenceC(new Range(start, end));
    } 
    catch(OutOfRangeException e) 
    {
      throw new Error("internal error - unexpected exception: " + e);
    }

    final int sub_sequence_length = sub_sequence.length;

    if(getStrand().isForwardStrand()) 
    {
      for(int i = 0 ; i < sub_sequence_length ; i += 3) 
      {
        for(int frame = 0 ; frame < 3 ; ++frame) 
        {
          final char this_char = sub_sequence[i + frame];

          if(this_char == 'g' || this_char == 'c') 
            ++gc_counts[(frame + start) % 3];
        }
      }
    } 
    else
    {
      final int whole_sequence_length = getStrand().getSequenceLength();
      final int whole_sequence_length_mod3 = whole_sequence_length % 3;

      for(int i = 0; i < sub_sequence_length; i += 3)
      {
        for(int frame = 0; frame < 3; ++frame) 
        {
          final char this_char = sub_sequence[i + frame];

          if(this_char == 'g' || this_char == 'c') 
            ++gc_counts[(frame + start + 3 - whole_sequence_length_mod3) % 3];
        }
      }
    }

    // multiply by 3 because we are taking every third base
    for(int frame = 0 ; frame < 3 ; ++frame) 
    {
      values[frame] = 1.0F * gc_counts[frame]/sub_sequence_length * 3 * 100;
    }
  }

  /**
  *  Override drawLegend() 
  */
  public void drawLegend(Graphics g, int font_height,
                         int font_width, LineAttributes[] lines, int numPlots)
  {
    return;
  }


  /**
   *  Return the number of values a call to getValues() will return - three
   *  in this case.
   **/
  public int getValueCount()
  {
    return 3;
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
    return new Integer(120);
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
    return new Integer(500);
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
    return new Integer(24);
  }

  /**
   *  Return the default or optimal step size.
   *  @return null is returned if this algorithm doesn't have optimal step
   *    size.
   **/
  public Integer getDefaultStepSize(int window_size) 
  {
    if(window_size > 8)
      return new Integer(window_size / 8);
    else
      return null;
  }

  /**
   *  Return the maximum value of this algorithm.
   *  @return The maximum is 100.
   **/
  protected Float getMaximumInternal()
  {
    return new Float(100);
  }

  /**
   *  Return the minimum value of this algorithm.
   *  @return The minimum is 0.
   **/
  protected Float getMinimumInternal() 
  {
    return new Float(0);
  }

  /**
   *  Return the average value of function over the whole strand.
   *  @return null is returned if this algorithm doesn't have an average or if
   *    the average can't be calculated.
   **/
  public Float getAverage() 
  {
    return new Float(getStrand().getBases().getAverageGCPercent());
  }

  /**
   *  Returns "GC Frame Plot" if the given strand is a forward strand
   *  otherwise returns "Reverse GC Frame Plot".
   **/
  private static String makeName(final Strand strand) 
  {
    if(strand.isForwardStrand())
      return "GC Frame Plot";
    else
      return "Reverse GC Frame Plot";
  }
}
