/* EntropyAlgorithm.java
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

import java.lang.Math.*;
import java.util.*;
import java.text.*;


import uk.ac.sanger.artemis.util.*;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.sequence.*;

/**
 *
 *  Informational Entropy (Konopka 1984)
 *  Konopka A (1984) Is the information content of DNA evolutionarily 
 *  significant? J Theor Biol 107:697-704
 *
 *  @author Derek Gatherer 
 **/

public class EntropyAlgorithm extends BaseAlgorithm 
{
  /**
   *  Create a new EntropyAlgorithm object.
   *  @param strand The strand to do the calculation on.
   **/
  public EntropyAlgorithm(final Strand strand) 
  {
    super(strand, makeName(strand), "entropy"); 
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
    
    Collator co = Collator.getInstance();
    TreeMap wordSet = new TreeMap(co);
    int total = 0;

    try 
    {
      sequence = getStrand().getSubSequence(new Range(start, end));
    }
    catch(OutOfRangeException e)
    {
      throw new Error("internal error - unexpected exception: " + e);
    }

//    System.out.println("start:"+start+"end:"+end);
//    System.out.println(sequence);

// calculates the overlapping triplet entropy

    for(int i = 0; i < sequence.length()-3; ++i)   // increment by 1, therefore overlapping
    {
      String codon = sequence.substring(i,i+3);   // codon String is 3 in length
      Integer number = (Integer)wordSet.get(codon);
      if(number == null)
	number = new Integer(0);
      wordSet.put(codon, new Integer(number.intValue()+1)); // count them up
      total++;	      
    }

    Set mappings = wordSet.entrySet();
    double ent = 0;

    for(Iterator i = mappings.iterator(); i.hasNext();) 
    {
      Map.Entry e = (Map.Entry)i.next();
      String in_hash = e.getValue().toString();
      float as_num = Float.parseFloat(in_hash);
      float freq = as_num/total;
      ent -= freq*Math.log(freq)/Math.log(2);  //  H = sigma(p* log(2) p)
    }
    values[0] = (float)ent;
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
    return new Integer(25);
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
  
  private static String makeName(final Strand strand) 
  {
    if(strand.isForwardStrand()) 
      return "Informational Entropy";
    else
      return "Reverse Informational Entropy";
  }
}
