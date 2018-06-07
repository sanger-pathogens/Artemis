/* CodonUsageAlgorithm.java
 *
 * created: Tue Apr 13 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/plot/CodonUsageAlgorithm.java,v 1.2 2006-06-23 10:40:14 tjc Exp $
 */

package uk.ac.sanger.artemis.plot;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;

import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.sequence.*;

import uk.ac.sanger.artemis.util.*;
import uk.ac.sanger.artemis.io.Range;

/**
 *  Objects of this class have one useful method - getValues (), which takes a
 *  range of bases and returns three floating point numbers, which are the
 *  codon usage scores in each frame.  The Strand to use is set in the
 *  constructor.
 *  See Gribskov et al. (Nucl. Acids Res. 12(1); 539-549 (1984)).
 *
 *  @author Kim Rutherford
 *  @version $Id: CodonUsageAlgorithm.java,v 1.2 2006-06-23 10:40:14 tjc Exp $
 **/

public class CodonUsageAlgorithm extends BaseAlgorithm {
  /**
   *  Create a new CodonUsageAlgorithm object.  See getValues () for detail on
   *  how the weightings object is used.
   *  @param strand The Strand to do the calculation on.
   *  @param usage_data This object is the source of the codon usage data that
   *    will be used by getValues ().
   **/
  public CodonUsageAlgorithm (final Strand strand,
                              final CodonUsageWeight usage_data) {
    super (strand, makeName (strand, usage_data), "codon_usage");

    this.usage_data = usage_data;

    setScalingFlag (true);
  }

  /**
   *  Return the average of the user weightings for the codons in the given
   *  window.
   *  @param start The start base (included in the range).
   *  @param end The end base (included in the range).
   *  @param values The three return values for this algorithm are returned in
   *    this array.  There is one value for each frame and the value is the
   *    the average of the weightings for the codons in the range.
   **/
  public void getValues(int start, int end, final float[] values) 
  {
    if (isRevCompDisplay ()) 
    {
      final int new_end =
        getStrand ().getBases ().getComplementPosition (start);
      final int new_start =
        getStrand ().getBases ().getComplementPosition (end);

      end = new_end;
      start = new_start;
    }

    // add 1 or 2 if necessary to make the range a multiple of 3
    if(getStrand ().isForwardStrand ()) 
      end -= (end - start + 1) % 3;
    else 
      start += (end - start + 1) % 3;

    final char[] sequence;

    try 
    {
      sequence = getStrand ().getRawSubSequenceC(new Range (start, end));
    }
    catch (OutOfRangeException e)
    {
      throw new Error ("internal error - unexpected exception: " + e);
    }

    float[] totals = { 0, 0, 0 };

    // a count of the number of codons we have seen
    int codon_count = 0;

    final int sub_sequence_length = sequence.length;

    if(getStrand ().isForwardStrand ()) 
    {
      for (int frame = 0 ; frame < 3 ; ++frame) 
      {
        final int real_frame = (frame + start + 2) % 3;
        for (int i = frame ; i < sub_sequence_length - 3 ; i += 3) 
        {
          final char base1 = sequence[i];
          final char base2 = sequence[i + 1];
          final char base3 = sequence[i + 2];

          final float this_weight =
            usage_data.getCodonValue(base1, base2, base3);

          ++codon_count;
          totals[real_frame] += Math.log(this_weight);
        }
      }
    }
    else
    {
      for(int frame = 2; frame >= 0 ; --frame) 
      {
        final int real_frame = (frame + start + 2) % 3;
        for (int i = frame ; i < sub_sequence_length - 3 ; i += 3) 
        {
          final char base1 = Bases.complement(sequence[i + 2]);
          final char base2 = Bases.complement(sequence[i + 1]);
          final char base3 = Bases.complement(sequence[i]);

          final float this_weight =
            usage_data.getCodonValue(base1, base2, base3);

          ++codon_count;
          totals[real_frame] += Math.log(this_weight);
        }
      }
    }

    for (int frame = 0 ; frame < 3 ; ++frame) 
    {
      if (codon_count == 0) 
        values[frame] = 0;
      else 
        values[frame] = (float)Math.exp(totals[frame] / codon_count);
    }
  }

  
  /**
   *  Return the number of values a call to getValues () will return - three
   *  in this case.
   **/
  public int getValueCount () {
    return 3;
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
   *  Return the default or optimal step size.
   *  @return null is returned if this algorithm doesn't have optimal step
   *    size.
   **/
  public Integer getDefaultStepSize (int window_size) {
    if (window_size > 8) {
      return new Integer (window_size / 8);
    } else {
      return null;
    }
  }

  /**
   *  Return the maximum value of this algorithm.
   *  @return The maximum is 100.
   **/
  protected Float getMaximumInternal () {
    return new Float (2);
  }

  /**
   *  Return the minimum value of this algorithm.
   *  @return The minimum is 0.
   **/
  protected Float getMinimumInternal () {
    return new Float (0);
  }

  /**
   *  Return the average value of function over the whole strand.
   *  @return null is returned if this algorithm doesn't have an average or if
   *    the average can't be calculated.
   **/
  public Float getAverage () {
    return null;
  }

  /**
   *  Calculate the codon usage score for the given Feature.
   **/
  public float getFeatureScore (final Feature feature) {
    final String sequence = feature.getTranslationBases ();

    float total = 0F;

    for (int i = 0 ; i < sequence.length () ; i += 3) {

      final char base1 = sequence.charAt (i);
      final char base2 = sequence.charAt (i + 1);
      final char base3 = sequence.charAt (i + 2);

      final float this_weight =
        usage_data.getCodonValue (base1, base2, base3);

      total += Math.log (this_weight);
    }

    final int codon_count = sequence.length () / 3;

    return (float) Math.exp (total / codon_count);
  }
  
  /**
   *  Returns "Codon Usage Scores" if the given strand is a forward strand
   *  otherwise returns "Reverse Codon Usage Scores".
   *  @param strand The Strand to do the calculation on.
   *  @param codon_weight This is used to get the name of the file that the
   *    usage information came from.
   **/
  private static String makeName (final Strand strand,
                                  final CodonWeight codon_weight) {
    if (strand.isForwardStrand ()) {
      return "Codon Usage Scores from " + codon_weight.getName ();
    } else {
      return "Reverse Codon Usage Scores from " + codon_weight.getName ();
    }
  }

  /**
   *  The CodonWeight reference that was passed to the constructor.
   **/
  final private CodonWeight usage_data;
}
