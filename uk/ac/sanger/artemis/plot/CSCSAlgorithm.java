/* CSCSAlgorithm.java
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
 */

package uk.ac.sanger.artemis.plot;

import uk.ac.sanger.artemis.util.*;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.sequence.*;

/**
 *
 * Corrected Scaled Chi Square (Mathur and Tuli 1991)
 * Mathur M, Tuli R (1991) Analysis of codon usage in genes for nitrogen 
 * fixation from phylogenetically diverse organisms. J Mol Evol 32:364-373
 *
 *  @author Derek Gatherer
 *  original version 10-09-03
 *  revised 01-12-04
 **/

public class CSCSAlgorithm extends BaseAlgorithm 
{
  /**
   *  Create a new CSCSAlgorithm object.
   *  @param strand The Strand to do the calculation on.
   **/
  public CSCSAlgorithm(final Strand strand) 
  {
    super(strand, makeName(strand), "Corrected Scaled Chi Square");
    setScalingFlag(true);
  }

  /**
   *  @param start The start base (included in the range).
   *  @param end The end base (included in the range).  If the start/end pair
   *    doesn't give a multiple of three bases end is moved down so that it is
   *    a multiple of three.
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

    final String sub_sequence;
    
    // add 1 or 2 if necessary to make the range a multiple of 3
    if(getStrand().isForwardStrand())
      end -= (end - start + 1) % 3;
    else
      start += (end - start + 1) % 3;

    try 
    {
      sub_sequence = getStrand().getSubSequence(new Range(start, end));
    }
    catch(OutOfRangeException e) 
    {
      throw new Error("internal error - unexpected exception: " + e);
    }

    final float[][][][] exp_value = new float[4][4][4][4];  // 3D for bases, 1 for frame

    char [] sequence_raw;
    sequence_raw = sub_sequence.toCharArray();

    int [][][][] obs_value = new int [4][4][4][4];
    float [] chi_square = new float [3];
    
    for(int x = 0 ; x < 4 ; ++x) 
    {
      for(int y = 0 ; y < 4 ; ++y) 
      {
        for(int z = 0 ; z < 4 ; ++z) 
        {
          for(int a = 0 ; a < 4 ; ++a) 
          {
            obs_value[x][y][z][a] = 0;
            exp_value[x][y][z][a] = 0;
          }
        }
      }
    }

    char this_f_base = 0;
    char next_f_base = 0;
    char last_f_base = 0;
    int this_f_base_index = 0;
    int next_f_base_index = 0;
    int last_f_base_index = 0;

    int GC = 0;  // GC count
    float GCfreq = 0; // GC/length

// get GC content
    for(int c = 0 ; c < sequence_raw.length; c++) 
    {
      char this_base = sequence_raw[c]; 
      if(this_base == 'g' || this_base == 'c') 
        GC++; 
    }
    GCfreq = (float)GC/sequence_raw.length;

    chi_square = new float [3];
    for(int i = 0 ; i < sequence_raw.length - 5 ; i+=3) 
    {
      for(int frame = 0; frame < 3; frame++) 
      {
        this_f_base = sequence_raw[i + frame];
        next_f_base = sequence_raw[i + 1 + frame];
        last_f_base = sequence_raw[i + 2 + frame];

        this_f_base_index = Bases.getIndexOfBase(this_f_base);
        next_f_base_index = Bases.getIndexOfBase(next_f_base);
        last_f_base_index = Bases.getIndexOfBase(last_f_base);
        
        // ignore Ns
        if(this_f_base_index < 4 && next_f_base_index < 4 && last_f_base_index < 4)
          ++obs_value[this_f_base_index][next_f_base_index][last_f_base_index][(frame+start)%3];
      } 
    }

 // having collected the observed values in all 3 frames, formulate the expected
 
    for(int frame_index = 0 ; frame_index < 3 ; ++frame_index) 
    {
      int frame = frame_index; //(frame_index+start)%3;   // so it corresponds to the above
 // Phe 001 000
 // 2-mer AT         exp = total aa * (1-GC)
      exp_value[0][0][0][frame] = (obs_value[0][0][0][frame]+obs_value[0][0][1][frame])*(1-GCfreq);
 // 2-mer GC         exp = total aa * GC
      exp_value[0][0][1][frame] = (obs_value[0][0][0][frame]+obs_value[0][0][1][frame])*GCfreq;
 // Leu = 002 003 100 101 102 103
 // arg leu SNW      exp = total aa * GC * (1-GC) /1.5     100 003 102
 // arg leu SNS      exp = total aa * GC * GC /1.5         101 103
 // arg leu WNW      exp = total aa * (1-GC) * (1-GC)/1.5  002
      exp_value[0][0][2][frame] = (obs_value[0][0][2][frame]+obs_value[0][0][3][frame]+obs_value[1][0][0][frame]
         +obs_value[1][0][1][frame]+obs_value[1][0][2][frame]+obs_value[1][0][3][frame])* ((1-GCfreq) * (1-GCfreq))/(float)1.5;
      exp_value[0][0][3][frame] = (obs_value[0][0][2][frame]+obs_value[0][0][3][frame]+obs_value[1][0][0][frame]
         +obs_value[1][0][1][frame]+obs_value[1][0][2][frame]+obs_value[1][0][3][frame])* (GCfreq * (1-GCfreq)) /(float)1.5;
      exp_value[1][0][0][frame] = (obs_value[0][0][2][frame]+obs_value[0][0][3][frame]+obs_value[1][0][0][frame]
         +obs_value[1][0][1][frame]+obs_value[1][0][2][frame]+obs_value[1][0][3][frame])* (GCfreq * (1-GCfreq)) /(float)1.5;
      exp_value[1][0][1][frame] = (obs_value[0][0][2][frame]+obs_value[0][0][3][frame]+obs_value[1][0][0][frame]
         +obs_value[1][0][1][frame]+obs_value[1][0][2][frame]+obs_value[1][0][3][frame])* (GCfreq * GCfreq) /(float)1.5;
      exp_value[1][0][2][frame] = (obs_value[0][0][2][frame]+obs_value[0][0][3][frame]+obs_value[1][0][0][frame]
         +obs_value[1][0][1][frame]+obs_value[1][0][2][frame]+obs_value[1][0][3][frame])* ((1-GCfreq) * (1-GCfreq))/(float)1.5;
      exp_value[1][0][3][frame] = (obs_value[0][0][2][frame]+obs_value[0][0][3][frame]+obs_value[1][0][0][frame]
         +obs_value[1][0][1][frame]+obs_value[1][0][2][frame]+obs_value[1][0][3][frame])* (GCfreq * GCfreq) /(float)1.5;
// Ile = 200 201 202
// ile A/T          exp = total aa * (1-GC)/1.5
      exp_value[2][0][0][frame] = (obs_value[2][0][0][frame]+obs_value[2][0][1][frame]+obs_value[2][0][2][frame])* (1-GCfreq)/(float)1.5;
      exp_value[2][0][1][frame] = (obs_value[2][0][0][frame]+obs_value[2][0][1][frame]+obs_value[2][0][2][frame])* GCfreq/(float)1.5;
      exp_value[2][0][2][frame] = (obs_value[2][0][0][frame]+obs_value[2][0][1][frame]+obs_value[2][0][2][frame])* (1-GCfreq)/(float)1.5;
// Val = 300 301 302 303
// 4-mer AT         exp = total aa * (1-GC)/2
      exp_value[3][0][0][frame] = (obs_value[3][0][0][frame]+obs_value[3][0][1][frame]+obs_value[3][0][2][frame]+obs_value[3][0][3][frame])*(1-GCfreq)/2;
// 4-mer GC         exp = total aa * GC/2
      exp_value[3][0][1][frame] = (obs_value[3][0][0][frame]+obs_value[3][0][1][frame]+obs_value[3][0][2][frame]+obs_value[3][0][3][frame])*GCfreq/2;
      exp_value[3][0][2][frame] = (obs_value[3][0][0][frame]+obs_value[3][0][1][frame]+obs_value[3][0][2][frame]+obs_value[3][0][3][frame])*(1-GCfreq)/2;
      exp_value[3][0][3][frame] = (obs_value[3][0][0][frame]+obs_value[3][0][1][frame]+obs_value[3][0][2][frame]+obs_value[3][0][3][frame])*GCfreq/2;
// Ser = 010 011 012 013 230 231
// ser AT           exp = total aa * (1-GC)/3
      exp_value[0][1][0][frame] = (obs_value[0][1][0][frame]+obs_value[0][1][1][frame]+obs_value[0][1][2][frame]
         +obs_value[0][1][3][frame]+obs_value[2][3][0][frame]+obs_value[2][3][1][frame])*(1-GCfreq)/3;
// ser GC           exp = total aa * GC/3
      exp_value[0][1][1][frame] = (obs_value[0][1][0][frame]+obs_value[0][1][1][frame]+obs_value[0][1][2][frame]
         +obs_value[0][1][3][frame]+obs_value[2][3][0][frame]+obs_value[2][3][1][frame])*GCfreq/3;
      exp_value[0][1][2][frame] = (obs_value[0][1][0][frame]+obs_value[0][1][1][frame]+obs_value[0][1][2][frame]
         +obs_value[0][1][3][frame]+obs_value[2][3][0][frame]+obs_value[2][3][1][frame])*(1-GCfreq)/3;
      exp_value[0][1][3][frame] = (obs_value[0][1][0][frame]+obs_value[0][1][1][frame]+obs_value[0][1][2][frame]
        +obs_value[0][1][3][frame]+obs_value[2][3][0][frame]+obs_value[2][3][1][frame])*GCfreq/3;
      exp_value[2][3][0][frame] =  (obs_value[0][1][0][frame]+obs_value[0][1][1][frame]+obs_value[0][1][2][frame]
        +obs_value[0][1][3][frame]+obs_value[2][3][0][frame]+obs_value[2][3][1][frame])*(1-GCfreq)/3;
      exp_value[2][3][1][frame] = (obs_value[0][1][0][frame]+obs_value[0][1][1][frame]+obs_value[0][1][2][frame]
        +obs_value[0][1][3][frame]+obs_value[2][3][0][frame]+obs_value[2][3][1][frame])*GCfreq/3;
// Pro = 110 111 112 113    
      exp_value[1][1][0][frame] = (obs_value[1][1][0][frame]+obs_value[1][1][1][frame]+obs_value[1][1][2][frame]+obs_value[1][1][3][frame])*(1-GCfreq)/2;
      exp_value[1][1][1][frame] = (obs_value[1][1][0][frame]+obs_value[1][1][1][frame]+obs_value[1][1][2][frame]+obs_value[1][1][3][frame])*GCfreq/2;
      exp_value[1][1][2][frame] = (obs_value[1][1][0][frame]+obs_value[1][1][1][frame]+obs_value[1][1][2][frame]+obs_value[1][1][3][frame])*(1-GCfreq)/2;
      exp_value[1][1][3][frame] = (obs_value[1][1][0][frame]+obs_value[1][1][1][frame]+obs_value[1][1][2][frame]+obs_value[1][1][3][frame])*GCfreq/2;
// Thr = 210 211 212 213    
      exp_value[2][1][0][frame] = (obs_value[2][1][0][frame]+obs_value[2][1][1][frame]+obs_value[2][1][2][frame]+obs_value[2][1][3][frame])*(1-GCfreq)/2;
      exp_value[2][1][1][frame] = (obs_value[2][1][0][frame]+obs_value[2][1][1][frame]+obs_value[2][1][2][frame]+obs_value[2][1][3][frame])*GCfreq/2;
      exp_value[2][1][2][frame] = (obs_value[2][1][0][frame]+obs_value[2][1][1][frame]+obs_value[2][1][2][frame]+obs_value[2][1][3][frame])*(1-GCfreq)/2;
      exp_value[2][1][3][frame] = (obs_value[2][1][0][frame]+obs_value[2][1][1][frame]+obs_value[2][1][2][frame]+obs_value[2][1][3][frame])*GCfreq/2;
// Ala = 310 311 312 313   
      exp_value[3][1][0][frame] = (obs_value[3][1][0][frame]+obs_value[3][1][1][frame]+obs_value[3][1][2][frame]+obs_value[3][1][3][frame])*(1-GCfreq)/2;
      exp_value[3][1][1][frame] = (obs_value[3][1][0][frame]+obs_value[3][1][1][frame]+obs_value[3][1][2][frame]+obs_value[3][1][3][frame])*GCfreq/2;
      exp_value[3][1][2][frame] = (obs_value[3][1][0][frame]+obs_value[3][1][1][frame]+obs_value[3][1][2][frame]+obs_value[3][1][3][frame])*(1-GCfreq)/2;
      exp_value[3][1][3][frame] = (obs_value[3][1][0][frame]+obs_value[3][1][1][frame]+obs_value[3][1][2][frame]+obs_value[3][1][3][frame])*GCfreq/2;
// Tyr = 020 021             
      exp_value[0][2][0][frame] = (obs_value[0][2][0][frame]+obs_value[0][2][1][frame])*(1-GCfreq);
      exp_value[0][2][1][frame] = (obs_value[0][2][0][frame]+obs_value[0][2][1][frame])*GCfreq;
// His = 120 121
      exp_value[1][2][0][frame] = (obs_value[1][2][0][frame]+obs_value[1][2][1][frame])*(1-GCfreq);
      exp_value[1][2][1][frame] = (obs_value[1][2][0][frame]+obs_value[1][2][1][frame])*GCfreq;
// Gln = 122 123    
      exp_value[1][2][2][frame] = (obs_value[1][2][2][frame]+obs_value[1][2][3][frame])*(1-GCfreq);
      exp_value[1][2][3][frame] = (obs_value[1][2][2][frame]+obs_value[1][2][3][frame])*GCfreq;
// Asn = 220 221    
      exp_value[2][2][0][frame] = (obs_value[2][2][0][frame]+obs_value[2][2][1][frame])*(1-GCfreq);
      exp_value[2][2][1][frame] = (obs_value[2][2][0][frame]+obs_value[2][2][1][frame])*GCfreq;
// Lys = 222 223     
      exp_value[2][2][2][frame] = (obs_value[2][2][2][frame]+obs_value[2][2][3][frame])*(1-GCfreq);
      exp_value[2][2][3][frame] = (obs_value[2][2][2][frame]+obs_value[2][2][3][frame])*GCfreq;
// Asp = 320 321     
      exp_value[3][2][0][frame] = (obs_value[3][2][0][frame]+obs_value[3][2][1][frame])*(1-GCfreq);
      exp_value[3][2][1][frame] = (obs_value[3][2][0][frame]+obs_value[3][2][1][frame])*GCfreq;
// Glu = 322 323   
      exp_value[3][2][2][frame] = (obs_value[3][2][2][frame]+obs_value[3][2][3][frame])*(1-GCfreq);
      exp_value[3][2][3][frame] = (obs_value[3][2][2][frame]+obs_value[3][2][3][frame])*GCfreq;
//Cys = 030 031      
      exp_value[0][3][0][frame] = (obs_value[0][3][0][frame]+obs_value[0][3][1][frame])*(1-GCfreq);
      exp_value[0][3][1][frame] = (obs_value[0][3][0][frame]+obs_value[0][3][1][frame])*GCfreq;
// Arg = 130 131 132 133 232 233
// arg leu SNW      exp = total aa * GC * (1-GC) /1.5     130 233 132
// arg leu SNS      exp = total aa * GC * GC /1.5         131 133
// arg leu WNW      exp = total aa * (1-GC) * (1-GC)/1.5  232
      exp_value[1][3][0][frame] = (obs_value[1][3][0][frame]+obs_value[1][3][1][frame]+obs_value[1][3][2][frame]
         +obs_value[1][3][3][frame]+obs_value[2][3][2][frame]+obs_value[2][3][3][frame])* (GCfreq * (1-GCfreq)) /(float)1.5;
      exp_value[1][3][1][frame] = (obs_value[1][3][0][frame]+obs_value[1][3][1][frame]+obs_value[1][3][2][frame]
         +obs_value[1][3][3][frame]+obs_value[2][3][2][frame]+obs_value[2][3][3][frame])* (GCfreq * GCfreq) /(float)1.5;
      exp_value[1][3][2][frame] = (obs_value[1][3][0][frame]+obs_value[1][3][1][frame]+obs_value[1][3][2][frame]
         +obs_value[1][3][3][frame]+obs_value[2][3][2][frame]+obs_value[2][3][3][frame])* (GCfreq * (1-GCfreq)) /(float)1.5;
      exp_value[1][3][3][frame] = (obs_value[1][3][0][frame]+obs_value[1][3][1][frame]+obs_value[1][3][2][frame]
         +obs_value[1][3][3][frame]+obs_value[2][3][2][frame]+obs_value[2][3][3][frame])* (GCfreq * GCfreq) /(float)1.5;
      exp_value[2][3][2][frame] = (obs_value[1][3][0][frame]+obs_value[1][3][1][frame]+obs_value[1][3][2][frame]
         +obs_value[1][3][3][frame]+obs_value[2][3][2][frame]+obs_value[2][3][3][frame])* ((1-GCfreq) * (1-GCfreq))/(float)1.5;
      exp_value[2][3][3][frame] = (obs_value[1][3][0][frame]+obs_value[1][3][1][frame]+obs_value[1][3][2][frame]
         +obs_value[1][3][3][frame]+obs_value[2][3][2][frame]+obs_value[2][3][3][frame])* (GCfreq * (1-GCfreq)) /(float)1.5;
// Gly = 330 331 332 333   
      exp_value[3][3][0][frame] = (obs_value[3][3][0][frame]+obs_value[3][3][1][frame]+obs_value[3][3][2][frame]+obs_value[3][3][3][frame])*(1-GCfreq)/2;
      exp_value[3][3][1][frame] = (obs_value[3][3][0][frame]+obs_value[3][3][1][frame]+obs_value[3][3][2][frame]+obs_value[3][3][3][frame])*GCfreq/2;
      exp_value[3][3][2][frame] = (obs_value[3][3][0][frame]+obs_value[3][3][1][frame]+obs_value[3][3][2][frame]+obs_value[3][3][3][frame])*(1-GCfreq)/2;
      exp_value[3][3][3][frame] = (obs_value[3][3][0][frame]+obs_value[3][3][1][frame]+obs_value[3][3][2][frame]+obs_value[3][3][3][frame])*GCfreq/2;

      chi_square[frame] = 0;
      exp_value[2][0][3][frame] = 0;  // don't count Met
      exp_value[0][3][3][frame] = 0;  // don't count Trp
      exp_value[0][2][2][frame] = 0;  // don't count STOP
      exp_value[0][2][3][frame] = 0;  // don't count STOP
      exp_value[0][3][2][frame] = 0;  // don't count STOP

// having calculated expected, now do chi_square
      chi_square[0]=0; chi_square[1]=0; chi_square[2]=0;
      for(int first_base = 0; first_base < 4; ++first_base)
      {
        for(int second_base = 0; second_base < 4; ++second_base) 
        {
          for(int third_base = 0; third_base < 4; ++third_base) 
          {
            final float obs = obs_value[first_base][second_base][third_base][frame];
            final float exp = exp_value[first_base][second_base][third_base][frame];
            if(exp >= 1) 
            {
              if(exp <= 5)   
              {
                if(Math.abs(exp-obs)>0.25) 
                    chi_square[frame] += (Math.pow((Math.abs(exp-obs)-0.5),2))/exp;
                else 
                  chi_square[frame] += (Math.pow((exp-obs),2))/exp;
              }
              else 
                chi_square[frame] += (Math.pow((exp-obs),2))/exp;
            }
          }
        } 
      }
      values[frame] = chi_square[frame];
    }  
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
    return new Integer(24);
  }

  /**
   *  Return the default or optimal step size.
   *  @return null is returned if this algorithm doesn't have optimal step
   *    size.
   **/
  public Integer getDefaultStepSize(int window_size) 
  {
    if(window_size > 10) 
    {
      Integer step = new Integer(24);
//      Integer step = new Integer(window_size/10);
//      int step_int = step.intValue();
//      step_int+=(step_int % 3);  // step is multiple of 3
//      Integer step_out = new Integer(step_int);
//      return step_out;
      return step; 

    }
    else 
      return null;
  }

  /**
   *  Return the maximum value of this algorithm.
   *  @return The maximum is 2.
   **/
  protected Float getMaximumInternal() 
  {
    return new Float(1000);
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
    return null;
  }  
  
  private static String makeName(final Strand strand) 
  {
    if(strand.isForwardStrand()) 
      return "Corrected Scaled Chi Square";
    else 
      return "Reverse Corrected Scaled Chi Square";
  }
}  // end of class
