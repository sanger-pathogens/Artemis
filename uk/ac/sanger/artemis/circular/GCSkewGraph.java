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

import uk.ac.sanger.artemis.io.Sequence;

public class GCSkewGraph extends Graph
{
  private static final long serialVersionUID = 1L;

  public GCSkewGraph(DNADraw currentDna)
  {
    super(currentDna);
  }
 
 /**
  * Return the value of (G content - C content)/(G content + C content)
  *  between the given pair of bases.
  */
  protected float calculateValue(int start, int end)
  {
    char[] sequence;
    if(end<=getBases().getLength())
      sequence = 
        getBases().getSequence().getCharSubSequence(start, end);
    else 
    {
      final Sequence s = getBases().getSequence();
      char[] seq1 = s.getCharSubSequence(start, getBases().getLength());
      char[] seq2 = s.getCharSubSequence(1, getWindowSize()-(getBases().getLength()-start));
      sequence = new char[seq1.length+seq2.length];
      System.arraycopy(seq1, 0, sequence, 0, seq1.length);
      System.arraycopy(seq2, 0, sequence, seq1.length-1, seq2.length);
    }

    
    float g_count = 0;
    float c_count = 0;

    for (int i = 0 ; i < sequence.length ; ++i) 
    {
      final char this_char = sequence[i];

      if (this_char == 'g') 
        ++g_count;

      if (this_char == 'c')
        ++c_count;
    }

    return (g_count - c_count) / (g_count + c_count);
  }
}