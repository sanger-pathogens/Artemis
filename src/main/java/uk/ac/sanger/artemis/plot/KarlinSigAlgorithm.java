/* KarlinSigAlgorithm.java
 *
 * created: Thu Aug  9 2001
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/plot/KarlinSigAlgorithm.java,v 1.4 2008-06-27 10:01:30 tjc Exp $
 */

package uk.ac.sanger.artemis.plot;

import uk.ac.sanger.artemis.util.*;
import uk.ac.sanger.artemis.io.Range;

import uk.ac.sanger.artemis.sequence.*;

/**
 *  This class contains the code for calculating the absolute "genomic
 *  signature" difference.
 *  For details see "Global dinucleotide signatures and analysis of genomic
 *  heterogeneity" Samuel Karlin - Current Opinion in Microbiology 1998,
 *  1:598-610.
 *  Objects of this class have one useful method - getValues (), which takes a
 *  range of bases and returns three floating point numbers, which are the
 *  codon usage scores in each frame.  The Strand to use is set in the
 *  constructor.
 *
 *  @author Kim Rutherford
 *  @version $Id: KarlinSigAlgorithm.java,v 1.4 2008-06-27 10:01:30 tjc Exp $
 **/

public class KarlinSigAlgorithm extends BaseAlgorithm {
  /**
   *  Create a new KarlinSigAlgorithm object.
   *  @param strand The Strand to do the calculation on.
   **/
  public KarlinSigAlgorithm (final Strand strand) {
    super (strand, "Karlin Signature Difference", "karlin_sig");

    setScalingFlag (true);
  }

  /**
   *  Return the Karlin genomic signature for the given range of bases.
   *  @param start The start base (included in the range).
   *  @param end The end base (included in the range).  If the start/end pair
   *    doesn't give a multiple of three bases end is moved down so that it is
   *    a multiple of three.
   *  @param values The one return value for this algorithm is returned in
   *    this array.
   **/
  public void getValues (int start, int end, final float [] values) {
    // add 1 or 2 if necessary to make the range a multiple of 3
    end -= (end - start + 1) % 3;

    final char[] sub_sequence;

    try {
      // reset
      if(end-start > 1000)
        ((uk.ac.sanger.artemis.io.StreamSequence)(getStrand().getBases().getSequence())).forceReset();
      sub_sequence = getStrand().getRawSubSequenceC(new Range (start, end));
    } catch (OutOfRangeException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    }

    final float [][] global_relative_abundance_values =
      getGlobalRelativeAbundance ();
    final float [][] subseq_relative_abundance_values =
      getRelativeAbundance (sub_sequence);

    float signature_difference = 0;

    for (int first_base = 0 ; first_base < 4 ; ++first_base) {
      for (int second_base = 0 ; second_base < 4 ; ++second_base) {
        final float global_value =
          global_relative_abundance_values[first_base][second_base];
        final float subseq_value =
          subseq_relative_abundance_values[first_base][second_base];
        signature_difference += Math.abs (global_value - subseq_value);
      }
    }

    values [0] = (float) signature_difference / 16f ;
  }

  /**
   *  Return the number of values a call to getValues () will return - three
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
    return new Integer (240);
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
    return new Integer (5000);
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
    if (window_size > 10) {
      return new Integer (window_size / 10);
    } else {
      return null;
    }
  }

  /**
   *  Return the maximum value of this algorithm.
   *  @return The maximum is 2.
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
   *  Return a 4x4 array containing the relative abundance values for each
   *  dinucleotide pair.  Indexed by base (t,c,a,g).  The value for the
   *  dinucleotide "TT" is stored in global_signature[0][0], "TC" is stored in
   *  [0][1], etc.
   **/
  private float [][] getRelativeAbundance (final char [] sequence_forward_raw) {
    final float [][] return_value = new float [4][4];

    final char [] sequence_reverse_raw =
      Bases.reverseComplement (sequence_forward_raw);

    final int [] base_counts = new int [4];
    final int [][] dinucleotide_base_counts = new int [4][4];

    int this_f_base_index = Bases.getIndexOfBase (sequence_forward_raw[0]);
    int next_f_base_index = 0;
    int this_r_base_index = Bases.getIndexOfBase (sequence_reverse_raw[0]);
    int next_r_base_index = 0;

    for (int i = 0 ; i < sequence_forward_raw.length - 1 ; ++i) 
    {
      next_f_base_index = Bases.getIndexOfBase (sequence_forward_raw[i + 1]);

      if (this_f_base_index < 4 && next_f_base_index < 4) {
        ++base_counts[this_f_base_index];
        ++dinucleotide_base_counts[this_f_base_index][next_f_base_index];
      } else {
        // ignore Ns
      }

      next_r_base_index = Bases.getIndexOfBase (sequence_reverse_raw[i + 1]);

      if (this_r_base_index < 4 && next_r_base_index < 4) {
        ++base_counts[this_r_base_index];
        ++dinucleotide_base_counts[this_r_base_index][next_r_base_index];
      } else {
        // ignore Ns
      }
      
      this_f_base_index = next_f_base_index;
      this_r_base_index = next_r_base_index;
    }

    // remember to add the last base
    if (next_f_base_index < 4) 
      ++base_counts[next_f_base_index];

    if (next_r_base_index < 4) 
      ++base_counts[next_r_base_index];


    for (int first_base_index = 0 ;
         first_base_index < 4 ;
         ++first_base_index) 
    {
      for (int second_base_index = 0 ;
           second_base_index < 4 ;
           ++second_base_index) 
      {
        final float dinucleotide_frequency =
          1f * dinucleotide_base_counts[first_base_index][second_base_index] /
          (sequence_reverse_raw.length - 1) / 2;
        final float first_base_frequency =
          1f * base_counts[first_base_index] /
          sequence_reverse_raw.length / 2;
        final float second_base_frequency =
          1f * base_counts[second_base_index] /
          sequence_reverse_raw.length / 2;

        return_value[first_base_index][second_base_index] =
          dinucleotide_frequency /
          (first_base_frequency * second_base_frequency);
      }
    }

    return return_value;
  }

  /**
   *  Return the relative abundance values for the complete sequence.  Indexed
   *  by base (t,c,a,g).  The value for the dinucleotide "TT" is stored in
   *  global_signature[0][0], "TC" is stored in [0][1], etc.
   **/
  private float [][] getGlobalRelativeAbundance () {
    if (global_relative_abundance_values == null) {
      try {
        final Range whole_range =
          new Range (1, getStrand ().getSequenceLength ());
        final char[] sequence =
          getStrand().getRawSubSequenceC(whole_range);
        global_relative_abundance_values = getRelativeAbundance (sequence);
      } catch (OutOfRangeException e) {
        throw new Error ("internal error - unexpected exception: " + e);
      }
    }

    return global_relative_abundance_values;
  }

  /**
   *  Storage for relative abundance values for the complete sequence.
   *  Indexed by base (t,c,a,g).  The value for the dinucleotide "TT" is
   *  stored in global_signature[0][0], "TC" is stored in [0][1], etc.
   **/
  private float [][] global_relative_abundance_values = null;
}
