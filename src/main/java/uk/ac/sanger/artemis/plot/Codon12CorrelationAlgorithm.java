/* Codon12CorrelationAlgorithm.java
 *
 * created: Mon Feb  1 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/plot/Codon12CorrelationAlgorithm.java,v 1.4 2006-06-23 10:40:14 tjc Exp $
 */

package uk.ac.sanger.artemis.plot;

import uk.ac.sanger.artemis.util.*;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.sequence.*;

import java.awt.*;

/**
 *  Objects of this class have one useful method - getValues(), which takes a
 *  range of bases and returns three floating point numbers, which are the
 *  codon position 1 and 2 correlation scores in each of the three frames for
 *  the given strand.  The Strand to use is set in the constructor.
 *
 *  @author Kim Rutherford
 *  @version $Id: Codon12CorrelationAlgorithm.java,v 1.4 2006-06-23 10:40:14 tjc Exp $
 **/
public class Codon12CorrelationAlgorithm extends BaseAlgorithm 
{
  /**
   *  Create a new Codon12CorrelationAlgorithm object.
   *  @param strand The strand to do the calculation on.
   **/
  public Codon12CorrelationAlgorithm(final Strand strand)
  {
    super(strand, makeName(strand), "correlation_score");
    setScalingFlag(true);
  }

  /**
   *  Magic numbers used by getValues() and Feature.get12CorrelationScore()
   *  to calculate the correlation score for codon positions 1 and 2.  These
   *  are the values for t, c, a, g in the first position.
   **/
  public static double[] correlation_score_factors_1 = 
  {
    17.7, 21.1, 27.7, 33.6
  };

  /**
   *  Magic numbers used by getValues() and Feature.get12CorrelationScore()
   *  to calculate the correlation score for codon positions 1 and 2.  These
   *  are the values for t, c, a, g in the second position.
   **/
  public static double[] correlation_score_factors_2 = 
  {
    27.1, 23.8, 31.0, 18.2
  };

  /**
   *  Return the codon position 1 and 2 correlation scores between a pair of
   *  bases in each of the three frames.
   *  @param start The start base (included in the range).
   *  @param end The end base (included in the range).  If the start/end pair
   *    doesn't give a multiple of three bases end is moved down so that it is
   *    a multiple of three.
   *  @param values The three results are returned in this array, hence it
   *    should be three long.  There is one value for each frame.
   **/
  public void getValues(int start, int end, final float[] values)
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

    final char[] sub_sequence_raw;

    try 
    {
      sub_sequence_raw = getStrand().getRawSubSequenceC(new Range(start, end));
    } 
    catch(OutOfRangeException e) 
    {
      throw new Error("internal error - unexpected exception: " + e);
    }

    final float gc_counts[] = new float[3];

    // the first index is the position the second is the base (t,c,a,g)
    final int[][] positional_base_counts = new int[4][3];
    final int sub_sequence_length = sub_sequence_raw.length;

    if(getStrand().isForwardStrand())
    {
      for(int i = 0 ; i < sub_sequence_length ; ++i)
      {
        final int base_index = Bases.getIndexOfBase(sub_sequence_raw[i]);
        if(base_index < 4) 
          ++positional_base_counts[base_index][i % 3];
      }
    } 
    else 
    {
      final char[] complement_sub_sequence_raw = Bases.complement(sub_sequence_raw);
      
      for(int i = 0 ; i < sub_sequence_length ; ++i) 
      {
        final int base_index =
          Bases.getIndexOfBase(complement_sub_sequence_raw[i]);
        if(base_index < 4) 
        {
          final int position_index = i % 3;
          ++positional_base_counts[base_index][position_index];
        }
      }
    }
    
    final int whole_sequence_length = getStrand().getSequenceLength();
    final int whole_sequence_length_mod3 = whole_sequence_length % 3;
    
    for(int frame = 0 ; frame < 3 ; ++frame) 
    {
      final double cor1_2_score =
        3.0 * (1.0 * positional_base_counts[0][frame]/sub_sequence_length *
               correlation_score_factors_1[0] +
               
               1.0 * positional_base_counts[1][frame]/sub_sequence_length * 
               correlation_score_factors_1[1] +
               
               1.0 * positional_base_counts[2][frame]/sub_sequence_length * 
               correlation_score_factors_1[2] +
               
               1.0 * positional_base_counts[3][frame]/sub_sequence_length * 
               correlation_score_factors_1[3] +
               
               1.0 * positional_base_counts[0][(frame + 1)%3]/
               sub_sequence_length * 
               correlation_score_factors_2[0] +
               
               1.0 * positional_base_counts[1][(frame + 1)%3]/
               sub_sequence_length * 
               correlation_score_factors_2[1] +
               
               1.0 * positional_base_counts[2][(frame + 1)%3]/
               sub_sequence_length * 
               correlation_score_factors_2[2] +
               
               1.0 * positional_base_counts[3][(frame + 1)%3]/
               sub_sequence_length * 
               correlation_score_factors_2[3]) +
        0.5;         // add 0.5 because that is what the old uk.ac.sanger.artemis did
      
      // add 2 brings the frame colouring in-line with codon usage
      if(getStrand().isForwardStrand())
        values [(start + frame + 2) % 3] = (float)cor1_2_score;
      else
        values [(start + frame + whole_sequence_length_mod3 + 2) % 3] = (float)cor1_2_score;
    }
    
    /*
    for(int frame = 0 ; frame < 3 ; ++frame) 
    {
      if(getStrand().isForwardStrand())
        System.out.println("FWD "+frame+"  "+((start + frame + 2) % 3));
      else
        System.out.println("BWD "+frame+"  "+((start + frame + whole_sequence_length_mod3 -1) % 3)+
                                        "  "+whole_sequence_length_mod3);
    }
    */
  }

   
  /**
   *  Return the number of values a call to getValues () will return - three
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
    return new Integer(240);
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
    return new Integer(600);
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
    return new Integer(48);
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
    return new Float(52.7);
  }

  /**
   *  Returns "Codon 1 and 2 Scores" if the given strand is a forward strand
   *  otherwise returns "Reverse Codon 1 and 2 Scores".
   **/
  private static String makeName(final Strand strand) 
  {
    if(strand.isForwardStrand())
      return "Correlation Scores";
    else 
      return "Reverse Correlation Scores";
  }
}
