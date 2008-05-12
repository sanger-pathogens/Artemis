/*
 * Copyright (C) 2008  Genome Research Limited
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
 *  @author: Tim Carver
 */

package uk.ac.sanger.artemis.circular;

import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.util.OutOfRangeException;

public class GCGraph extends Graph
{
  private static final long serialVersionUID = 1L;
  
  public GCGraph(DNADraw currentDna)
  {
    super(currentDna);
  }
  
  /**
   *  Recalculate the values in value_array_array, step_size, min_value and
   *  max_value.
   **/
  protected float calculateValue(int start, int end)
  {
    String sequence;
    try 
    {
      if(end<=getBases().getLength())
        sequence = 
          getBases().getForwardStrand().getSubSequence (new Range (start, end));
      else 
      {
        sequence = getBases().getForwardStrand().getSubSequence (
                   new Range (start, getBases().getLength()));
      sequence = sequence +
                   getBases().getForwardStrand().getSubSequence (
               new Range (1, getWindowSize()-(getBases().getLength()-start)));
      }
    } 
    catch (OutOfRangeException e) 
    {
      throw new Error ("internal error - unexpected exception: " + e);
    }

    
    float gc_count = 0;

    for (int i = 0 ; i < sequence.length() ; ++i) 
    {
      final char this_char = sequence.charAt (i);
      if (this_char == 'g' || this_char == 'c')
        ++gc_count;
    }
    return gc_count/sequence.length() * 100;
  }
}