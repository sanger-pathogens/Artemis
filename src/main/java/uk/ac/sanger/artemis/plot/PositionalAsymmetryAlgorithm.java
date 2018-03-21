/* PositionalAsymmetryAlgorithm.java
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
 **/



package uk.ac.sanger.artemis.plot;

import uk.ac.sanger.artemis.util.*;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.sequence.*;
import java.lang.Math.*;

/**
 *
 *  Positional asymmetry (Shulman et al. 1981)
 *  Shulman MJ, Steinberg CM, Westmoreland N (1981) The coding function of 
 *  nucleotide sequences can be discerned by statistical analysis. J Theor Biol 
 *  88:409-20
 *
 *  @author Derek Gatherer
 **/

public class PositionalAsymmetryAlgorithm extends BaseAlgorithm 
{
  /**
   *  Create a new GCWindowAlgorithm object.
   *  @param strand The strand to do the calculation on.
   **/
  public PositionalAsymmetryAlgorithm(final Strand strand) 
  {
    super(strand, makeName(strand), "positional_asymmetry"); 
    setScalingFlag(true);
  }

  /**
   *  Positional asymmetry
   *  @param start The start base (included in the range).
   *  @param end The end base (included in the range).
   *  @param values The one return value for this algorithm is returned in
   *    this array.
   **/
  public void getValues(int start, int end, final float [] values) 
  {
    if(getStrand().isForwardStrand())  // rather than isRevCompDisplay()
    {
//    System.out.println("isRevCompDisplay does not activate here");
    }
    else
    {
      final int new_end =
        getStrand().getBases().getComplementPosition(start);
      final int new_start =
        getStrand().getBases().getComplementPosition(end);

      end = new_end;
      start = new_start;
//      System.out.println("Revcomp, so new start:"+start+"new end:"+end);
    }
    
    final String sequence;

    try 
    {
      sequence = getStrand().getSubSequence(new Range(start, end));
    }
    catch(OutOfRangeException e) 
    {
      throw new Error("internal error - unexpected exception: " + e);
    }

    float g_count = 0;
    float c_count = 0;
    float a_count = 0;
    float t_count = 0;
    float[] a_pos = {0,0,0};
    float[] c_pos = {0,0,0};
    float[] g_pos = {0,0,0};
    float[] t_pos = {0,0,0};
    
    float a_chi = 0; 
    float c_chi = 0; 
    float g_chi = 0; 
    float t_chi= 0; 
    float chi_total = 0;

    for(int i = 0; i < sequence.length(); ++i) 
    {
      char this_char = sequence.charAt(i);

      if(this_char == 'g') 
        ++g_count;
      else if(this_char == 'c') 
        ++c_count;
      else if(this_char == 't') 
        ++a_count;
      else if(this_char == 'a') 
        ++t_count;
    }

    float exp_a = a_count/3;
    float exp_c = c_count/3;
    float exp_g = g_count/3;
    float exp_t = t_count/3;
    
    for(int frame = 0; frame < 3; ++frame) 
    { 
      for(int i = frame; i < sequence.length(); i=i+3) 
      {
        char this_char = sequence.charAt(i);
	if(this_char == 'g') 
          ++g_pos[frame];
        else if(this_char == 'c') 
          ++c_pos[frame];
        else if(this_char == 't')
          ++t_pos[frame];
        else if(this_char == 'a')
          ++a_pos[frame];
      }
    }
    
    for(int frame = 0; frame < 3; ++frame) 
    {
      if(exp_a >= 1) 
      {
        if(exp_a <= 5) 
          a_chi += (Math.pow((Math.abs(exp_a-a_pos[frame])-0.5),2))/exp_a;
        else
          a_chi += (Math.pow((exp_a-a_pos[frame]),2))/exp_a;
      }
      if(exp_c >= 1) 
      {
        if(exp_c <= 5)
          c_chi += (Math.pow((Math.abs(exp_c-c_pos[frame])-0.5),2))/exp_c;
        else
          c_chi += (Math.pow((exp_c-c_pos[frame]),2))/exp_c;
      }
      if(exp_g >= 1) 
      {
        if(exp_g <= 5) 
          g_chi += (Math.pow((Math.abs(exp_g-g_pos[frame])-0.5),2))/exp_g;
        else
          g_chi += (Math.pow((exp_g-g_pos[frame]),2))/exp_g;
      }
      if(exp_t >= 1) 
      {
        if(exp_t <= 5) 
          t_chi += (Math.pow((Math.abs(exp_t-t_pos[frame])-0.5),2))/exp_t;
        else 
          t_chi += (Math.pow((exp_t-t_pos[frame]),2))/exp_t;
      }
    }
    chi_total = a_chi + c_chi + g_chi + t_chi;
         
    
//    System.out.println ("start: " + start + " end: " + end + " returning: " + gc_count/sequence.length()); 

    values[0] = chi_total;
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
    return new Integer(500);
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
    return new Float(10000);
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
  
  private static String makeName(final Strand strand) 
  {
    if(strand.isForwardStrand()) 
      return "Positional Asymmetry";
    else
      return "Reverse Positional Asymmetry";
  }
}
