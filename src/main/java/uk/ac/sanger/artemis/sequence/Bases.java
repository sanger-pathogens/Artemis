/* Bases.java
 *
 * created: Sun Oct 11 1998
 *
 * This file is part of Artemis
 *
 * Copyright (C) 1998-2005  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/sequence/Bases.java,v 1.26 2009-03-27 14:00:51 tjc Exp $
 */

package uk.ac.sanger.artemis.sequence;

import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.util.*;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.io.EmblStreamSequence;
import uk.ac.sanger.artemis.io.Sequence;
import uk.ac.sanger.artemis.io.StreamSequence;

import org.biojava.bio.symbol.IllegalSymbolException;

import java.util.WeakHashMap;
import java.util.Iterator;

/**
 *  This class is a wrapper for the uk.ac.sanger.artemis.io.Sequence class
 *  that allows us to control what is done to the sequence and to send events
 *  to interested objects when changes happen.  Note: a '@' character is used
 *  as a marker when we don't have a base letter, for example complementing a
 *  non-base letter returns '@'.
 *
 *  @author Kim Rutherford
 *  @version $Id: Bases.java,v 1.26 2009-03-27 14:00:51 tjc Exp $ */

public class Bases 
{
  /**
   *  Indicates the bases should be read in the forward direction for a
   *  particular operation.
   **/
  static public final int FORWARD = 1;

  /**
   *  Indicates the bases should be read in the reverse direction for a
   *  particular operation.
   **/
  static public final int REVERSE = 2;

  /**
   *  The lowest possible value for use with addSequenceChangeListener ().
   **/
  static public final int MIN_PRIORITY = -5;

  /**
   *  An arbitrary value for use with addSequenceChangeListener ().
   **/
  static public final int MEDIUM_PRIORITY = 0;

  /**
   *  The highest possible value for use with addSequenceChangeListener ().
   **/
  static public final int MAX_PRIORITY = 5;

  /**
   *  A cache of the forward & reverse stop codon positions.
   *  0 means not set/cached yet, 1 not a stop codon, 2 and 3 are a 
   *  stop codon on fwd or reverse strand respectively.
   **/
  private byte [] stop_codon_cache = null;

  /**
   *  A cache of the forward & reverse start codon positions.
   *  0 means not set/cached yet, 1 not a start codon, 2 and 3 are a 
   *  start codon on fwd or reverse strand repectively.
   **/
  private byte [] start_codon_cache = null;
  
  /**
   *  Create a new Bases object.
   *  @param sequence The raw sequence that the new object will use.
   **/
  public Bases(final Sequence sequence) 
  {
    this.embl_sequence = sequence;

    stop_codon_cache = null;

    forward_strand = new Strand(this, FORWARD);
    reverse_strand = new Strand(this, REVERSE);

    for(int i = 0 ; i < listener_hash_map_array.length ; ++i) 
      listener_hash_map_array [i] = new WeakHashMap();
  }

  /**
   *  Return the object representing the forward sequence of bases for this
   *  object.
   **/
  public Strand getForwardStrand() 
  {
    return forward_strand;
  }

  /**
   *  Return the object representing the reverse complemented sequence of
   *  bases for this Bases objects.
   **/
  public Strand getReverseStrand() 
  {
    return reverse_strand;
  }

  /**
   *  Returns the length of the sequence in bases.
   **/
  public int getLength() 
  {
    return embl_sequence.length();
  }

  /**
   *  Return a String representation of the sequence.
   **/
  public String toString() 
  {
    return embl_sequence.getSubSequence(1,getLength());
  }

  /**
   *  Reverse and complement both of the Strand objects (by swapping them and
   *  reverse complementing the sequence).
   *  @exception ReadOnlyException If the Bases cannot be changed.
   **/
  public void reverseComplement()
      throws ReadOnlyException 
  {
    stop_codon_cache = null;

    final Strand temp = forward_strand;
    forward_strand = reverse_strand;
    reverse_strand = temp;

//  final String new_sequence =
//    reverseComplement(getSequence().getSubSequence(1, getLength()));

    final char[] new_sequence =
      reverseComplement(getSequence().getCharSubSequence(1, getLength()));

    try 
    {
//    getSequence().setFromChar(new_sequence.toCharArray());
      getSequence().setFromChar(new_sequence);
    } 
    catch (IllegalSymbolException e) 
    {
      throw new Error ("internal error - unexpected exception: " + e);
    }

    final SequenceChangeEvent event =
      new SequenceChangeEvent (this, SequenceChangeEvent.REVERSE_COMPLEMENT);

    fireSequenceChangeEvent (event);
  }

  /**
   *  This array is used to convert between bases and indices.  See
   *  getIndexOfBase()
   **/
  public final static char[] letter_index = 
  {
    't', 'c', 'a', 'g', 'n'
  };

  /**
   *  Given a base letter return its index where t = 0, c = 1, a = 2, g = 3, 4
   *  otherwise.
   *  See letter_index.
   **/
  public final static int getIndexOfBase(final char base) 
  {
    switch(base) 
    {
      case 'c':
        return 1;
      case 'a':
        return 2;
      case 'g':
        return 3;
      case 't':
      case 'u':
        return 0;
    }
  
    return 4;
  }

  /**
   *  Return the complement of the given Range.  eg. if the sequence length is
   *  100 and the Range is 1..10 then the return value will be 90..100.
   **/
  private Range complementRange (final Range range) {
    final int real_start = getComplementPosition (range.getEnd ());
    final int real_end   = getComplementPosition (range.getStart ());

    try {
      final Range real_range = new Range (real_start, real_end);

      return real_range;
    } catch (OutOfRangeException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    }
  }

  /**
   *  Return the complement of the given position on the sequence.  eg. if the
   *  sequence length is 100 and the position is 10 then the return value will
   *  be 90.
   **/
  public int getComplementPosition (final int position) {
    return getLength () - position + 1;
  }

  /**
   *  Return the raw of a base position on this object.  The raw position of a
   *  base on the forward strand is the same as the position itself.  The raw
   *  position of a base on the reverse strand is position of the
   *  corresponding bases on the forward strand.
   *  @param position The position of the base.
   *  @param direction The direction (strand) that the position refers to.
   **/
  public int getRawPosition (final int position, final int direction) {
    if (direction == FORWARD) {
      return position;
    } else {
      return getComplementPosition (position);
    }
  }

  /**
   *  Translate a sequence of bases into the corresponding single letter amino
   *  acid codes.
   *  @param range The range of the bases to translated.  If the range.start
   *    - range.end + 1 is not a multiple of three the last codon is
   *    incomplete and will not be translated.  If the range is out of range
   *    ie. it has a start or end less than one or greater than the length of
   *    the sequence, then the out of range codons will be translated as '.'.
   *  @param direction The direction of the translation.  If FORWARD the
   *    translation will happen as expected, if REVERSE the translation will
   *    be done on the reverse complement.
   *  @param unknown_is_x If this parameter is true codons that contain
   *    ambiguous bases will be translated as 'x', if false they will be
   *    translated as '.'
   *  @return The translated sequence in one letter abbreviated form.
   **/
  public AminoAcidSequence getTranslation(final Range range,
                                          final int direction,
                                          final boolean unknown_is_x) 
  {
    // getSubSequenceC() will return a sequence going in the right direction
    // so we don't have to worry.
    final char[] sub_sequence = getSubSequenceC(range, direction);
    return AminoAcidSequence.getTranslation(sub_sequence, unknown_is_x);
  }


  public AminoAcidSequence getSpacedTranslation(final Range range,
                                           final int direction,
                                           final boolean unknown_is_x) 
  {
    // getSubSequenceC() will return a sequence going in the right direction
    // so we don't have to worry.
    final char[] sub_sequence = getSubSequenceC(range, direction);
    return AminoAcidSequence.getSpacedTranslation(sub_sequence, unknown_is_x);
  }

  /**
   *  Return an array containing the positions of the codons that match the
   *  strings given by the query_codons argument.  Only those codons that are
   *  in the same frame as the first base of the range are returned.
   *  @param range The inclusive range of bases to get the codons from.
   *  @param direction The direction of the translation.  REVERSE means
   *    translate the reverse complement bases (the positions in the range
   *    argument are complemented first.)
   *  @param query_codons The codons to search for.  Each element of this
   *    vector should be a string that is 3 characters long.
   *  @return An array containing the positions of the first base of the
   *    codons.  This array is padded with zeros at the end.
   **/
  public int [] getMatchingCodons (final Range range, final int direction,
                                   final StringVector query_codons) {
    final Range real_range;

    if(direction == FORWARD)
      real_range = range;
    else
      real_range = complementRange(range);

    // guess the number of codons in getCount () bases - there are
    // query_codons.size() search codons in every 64 codons if G+C is 50%
    // and we have getCount()/3 codons to look at.

    float at_content = (100 - getAverageGCPercent()) / 100;

    int array_start_size =
      (int) (range.getCount () *
             at_content * at_content * (2-at_content) *
             query_codons.size () / 64);

    if(array_start_size < 20)
      array_start_size = 20;

    // this array will be resized as necessary
    int[] return_positions = new int[array_start_size];

    int current_return_array_index = 0;

    final String sequence_string =
      getSequence ().getSubSequence (1, getLength ());

    final int range_start_index = real_range.getStart () - 1;
    final int range_end_index = real_range.getEnd () - 1;

    if(direction == FORWARD) 
    {
      for (int i = range_start_index ; i < range_end_index - 2 ; i += 3) {
        if (i < 0 || i >= sequence_string.length () - 2) {
          continue;
        }

        boolean is_matching_codon =
          isMatchingCodon (sequence_string, i, direction, query_codons);

        if (is_matching_codon) {
          if (current_return_array_index == return_positions.length) {
            // first reallocate the array
            final int [] new_array =
              new int [return_positions.length * 3 / 2 + 1];

            System.arraycopy (return_positions, 0,
                              new_array, 0,
                              return_positions.length);
            return_positions = new_array;
          }

          return_positions[current_return_array_index] = i + 1;

          ++current_return_array_index;
        }
      }
    } else {

      for (int i = range_end_index ; i > range_start_index + 2 ; i -= 3) {
        if (i < 2 || i >= sequence_string.length ()) {
          continue;
        }

        boolean is_matching_codon =
          isMatchingCodon (sequence_string, i, direction, query_codons);

        if (is_matching_codon) {
          if (current_return_array_index == return_positions.length) {
            // first reallocate the array
            final int [] new_array =
              new int [return_positions.length * 3 / 2 + 1];

            System.arraycopy (return_positions, 0,
                              new_array, 0,
                              return_positions.length);
            return_positions = new_array;
          }

          // return the complemented base position
          return_positions[current_return_array_index] =
            sequence_string.length () - i;

          ++current_return_array_index;
        }
      }
    }

    return return_positions;

  }

  /**
   *  Check a three character substring and return true if and only if the
   *  three bases match an element of the query_codons argument.  If the
   *  direction is REVERSE then the three bases to check are at start_index,
   *  start_index - 1 and start_index - 2.  In that case true is returned if
   *  and only the complement of those three bases matches.
   **/
  private boolean isMatchingCodon (final String sequence_string,
                                   final int start_index,
                                   final int direction,
                                   final StringVector query_codons) {
    for (int query_codon_index = 0 ;
         query_codon_index < query_codons.size () ;
         ++query_codon_index) {
      if (isMatchingCodon (sequence_string, start_index, direction,
                           (String)query_codons.elementAt (query_codon_index))) {
        return true;
      }
    }

    return false;
  }

  /**
   *  Check a three character substring and return true if and only if the
   *  three bases match the query_codon argument.  If the direction is
   *  REVERSE then the three bases to check are at start_index, start_index -
   *  1 and start_index - 2.  In that case true is returned if and only the
   *  complement of those three bases matches.
   **/
  private boolean isMatchingCodon (final String sequence_string,
                                   final int start_index,
                                   final int direction,
                                   final String query_codon) {
    if (direction == FORWARD) {
      if (query_codon.charAt (0) == sequence_string.charAt (start_index) &&
          query_codon.charAt (1) == sequence_string.charAt (start_index + 1) &&
          query_codon.charAt (2) == sequence_string.charAt (start_index + 2)) {
        return true;
      }
    } else {
      final char first_letter =
        complement (sequence_string.charAt (start_index));
      final char second_letter =
        complement (sequence_string.charAt (start_index - 1));
      final char third_letter =
        complement (sequence_string.charAt (start_index - 2));

      if (query_codon.charAt (0) == first_letter &&
          query_codon.charAt (1) == second_letter &&
          query_codon.charAt (2) == third_letter) {
        return true;
      }
    }

    return false;
  }

  /**
   *  Returns stop_codon_cache after allocating it (if it is null).
   **/
  private byte[] getStopCodonCache() 
  {
    if (stop_codon_cache == null) 
    { 
      final int nbytes = getLength() >> 1 >> 1;
      stop_codon_cache = new byte[nbytes+1];
    }

    return stop_codon_cache;
  }

  
  /**
   *  Returns start_codon_cache after allocating it (if it is null).
   **/
  private byte[] getStartCodonCache() 
  {
    if (start_codon_cache == null) 
    { 
      final int nbytes = getLength() >> 1 >> 1;
      start_codon_cache = new byte[nbytes+1];
    }

    return start_codon_cache;
  }

  /**
   *  Clear stop codon cache (forward and reverse).
   **/
  public void clearCodonCache()
  {
    stop_codon_cache = null;
    start_codon_cache = null;
  }
 

  /**
   *  Return an array containing the positions of the stop codons.  Only those
   *  codons that are in the same frame as the first base of the range are
   *  returned.
   *  @param range The inclusive range of bases to get the stop codons from.
   *  @param direction The direction of the translation.  REVERSE means
   *    translate the reverse complement bases (the positions in the range
   *    argument are complemented first.)
   *  @return An array containing the positions of the first base of the stop
   *    codons.  This array is padded with zeros at the end.
   **/
  protected int[] getStopCodons(final Range range, final int direction)
  {
    final Range real_range;

    if(direction == FORWARD)
      real_range = range;
    else
      real_range = complementRange (range);

    // guess the number of stop codons in getCount() bases - there are 3
    // stop codons in every 64 codons if G+C is 50% and we have getCount()/3
    // codons to look at.

    float at_content = (100 - getAverageGCPercent()) / 100;

    int array_start_size =
      (int)(range.getCount() *
             at_content * at_content * (2-at_content) * 3 / 64);

    if(array_start_size < 20)
      array_start_size = 20;

    // this array will be resized as necessary
    int[] return_positions = new int[array_start_size];

    int current_return_array_index = 0;
    int range_start_index = real_range.getStart();
    int range_end_index   = real_range.getEnd();

    final int sequence_length = getLength();

    if(range_start_index < 1)
    {
      if(direction == FORWARD)
        range_start_index = 3 + (range_start_index % 3);
      else
        range_start_index =  1;
    }
    if(range_end_index > sequence_length)
      range_end_index = sequence_length;

    final char sequence_string[] =
      getSequence().getCharSubSequence(range_start_index, range_end_index);

    range_start_index--;
    range_end_index--;

    // whether a codon is a stop codon or not is cached in
    // 2 bit chunks (i.e. 4 per byte)
    int ncurrent_byte;
    int bit_position;
    byte bitty;

    final byte[] this_stop_codon_flags = getStopCodonCache();
    if(direction == FORWARD)
    {
    
      for(int i = range_start_index; i < range_end_index + 1; i += 3)
      {
        if(i < 0 || i >= sequence_length-1)
          continue;

        ncurrent_byte = i >> 1 >> 1;
        bit_position  = i % 4;

        // determine if codon type is cached or not
        bitty = (byte) ((this_stop_codon_flags[ncurrent_byte]
                              >> (2*bit_position) ) & 0x0003);
        if(bitty == 0)
        {
          // not cached yet
          setCache(range_start_index, range_end_index, sequence_string, i,
              null, this_stop_codon_flags, ncurrent_byte, bit_position);
        }

        bitty = (byte) ((this_stop_codon_flags[ncurrent_byte]
                           >> (2*bit_position) ) & 0x0003);
        if( bitty == 1 || bitty == 3 ) 
          continue;

        // if we reach here this is a stop codon
        if(current_return_array_index == return_positions.length)
        {
            // first reallocate the array
          final int[] new_array =
            new int[return_positions.length * 3 / 2 + 1];

          System.arraycopy(return_positions, 0,
                           new_array, 0,
                           return_positions.length);
          return_positions = new_array;
        }

        return_positions[current_return_array_index] = i + 1;
        ++current_return_array_index;
      }
    }
    else
    {
      for (int i = range_end_index ; i > range_start_index + 2 ; i -= 3)
      {
        if(i < 2 || i >= sequence_length)
          continue;

        ncurrent_byte = i >> 1 >> 1;
        bit_position = i % 4;
        bitty = (byte) ((this_stop_codon_flags[ncurrent_byte]
                              >> (2*bit_position) ) & 0x0003);

        if(bitty == 0)
        {
          // not cached yet
          setCache(range_start_index, range_end_index, sequence_string, i,
              null, this_stop_codon_flags, ncurrent_byte, bit_position);
        }

        bitty = (byte) ((this_stop_codon_flags[ncurrent_byte]
                           >> (2*bit_position) ) & 0x0003);

        if( bitty == 1 || bitty != 3 )  
          continue;

        // if we reach here this is a stop codon
        if(current_return_array_index == return_positions.length)
        {
          // first reallocate the array
          final int[] new_array =
            new int[return_positions.length * 3 / 2 + 1];

          System.arraycopy(return_positions, 0,
                           new_array, 0,
                           return_positions.length);
          return_positions = new_array;
        }

        return_positions[current_return_array_index] =
            sequence_length - i;
        ++current_return_array_index;
      }
    }

    return return_positions;
  }

  /**
   * Return an 2D array containing the stop or start codons in a range for
   *  all 3 frames of the strand. 
   *  @param range The inclusive range of bases to get the codons from.
   *  @param direction The direction of the translation.  REVERSE means
   *    translate the reverse complement bases (the positions in the range
   *    argument are complemented first.)
   *  @param query_codons if this is NULL then this assumes we are looking
   *    for stop codons, otherwise this is used to look for start codons.
   *  @return An array containing the positions of the first base of the stop
   *    codons.  This array is padded with zeros at the end.
   **/
  protected int[][] getStopOrStartCodons(final Range range, 
                                         final int direction,                            
                                         final StringVector query_codons) 
  {
    final Range real_range;

    if(direction == FORWARD)
      real_range = range;
    else
      real_range = complementRange(range);

    // guess the number of stop codons in getCount() bases - there are 3
    // stop codons in every 64 codons if G+C is 50% and we have getCount()/3
    // codons to look at.

    float at_content = (100 - getAverageGCPercent()) / 100;

    int array_start_size =
      (int)(range.getCount() *
            at_content * at_content * (2-at_content) * 3 / 64);

    if(array_start_size < 20)
      array_start_size = 20;
    // this array will be resized as necessary
    int[][] return_positions = new int[3][array_start_size];

    int[] current_return_array_index = new int[3];
    current_return_array_index[0] = 0;
    current_return_array_index[1] = 0;
    current_return_array_index[2] = 0;

    int range_start_index = real_range.getStart();
    int range_end_index   = real_range.getEnd();

    final int sequence_length = getLength();

    if(range_start_index < 1)
    {
      if(direction == FORWARD)
        range_start_index = 3 + (range_start_index % 3);
      else
        range_start_index =  1;
    }
    
    if(range_end_index > sequence_length)
      range_end_index = sequence_length;
 
    range_start_index--;
    range_end_index--;
    char[] sequence_string = null;

    // whether a codon is a stp codon or not is cached in
    // 2 bit chunks (i.e. 4 per byte)
    int ncurrent_byte;
    int bit_position;
    int nframe = 0;
    byte bitty;

    final byte[] this_forward_codon_flags;
    // if this is null then searching for stop codons
    if(query_codons == null)
      this_forward_codon_flags = getStopCodonCache();
    else
      this_forward_codon_flags = getStartCodonCache();
      
    for(int i = range_start_index; i < range_end_index+1; i += 1)
    {
      if(i < 0 || i >= sequence_length)
        continue;

      ncurrent_byte = i >> 1 >> 1;
      bit_position  = i % 4;

      // determine if codon type is cached or not
      bitty = (byte) ((this_forward_codon_flags[ncurrent_byte]
                                >> (2*bit_position) ) & 0x0003);
        
      if(bitty == 0)  // not cached yet
      {
        if(sequence_string == null)
          sequence_string = getSequence().getCharSubSequence(range_start_index+1,
                                                               range_end_index+1);

        setCache(range_start_index, range_end_index, sequence_string, i,
                 query_codons, this_forward_codon_flags, ncurrent_byte, 
                 bit_position);
        bitty = (byte) ((this_forward_codon_flags[ncurrent_byte]
                                 >> (2*bit_position) ) & 0x0003);
      }

      if(  bitty == 1 ||                         // not a stop/start codon
          (direction == FORWARD && bitty == 3) ||
          (direction != FORWARD && bitty != 3 ))  
        continue;

      if(direction == FORWARD) 
        nframe = (i-range_start_index) % 3;
      else
        nframe = (range_end_index-i) % 3;
      
      // if we reach here this is a stop/start codon
      if(current_return_array_index[nframe] == return_positions[nframe].length) 
      {
        // first reallocate the array
        final int[][] new_array =
            new int[3][return_positions[nframe].length * 3 / 2 + 1];

        for(int j=0; j<3; j++)
          System.arraycopy(return_positions[j], 0,
                           new_array[j], 0,
                           return_positions[j].length);
        return_positions = new_array;
      }

      if(direction == FORWARD)
      {
        if(i==0)
          return_positions[nframe][current_return_array_index[nframe]] = i + 1;
        else
          return_positions[nframe][current_return_array_index[nframe]] = i;
      }
      else
        return_positions[nframe][current_return_array_index[nframe]] =
              sequence_length - i;
      ++current_return_array_index[nframe];
    }
    
    return return_positions;
  }

  /**
   * Set the codon cache for forward and reverse strand.
   * @param range_start_index
   * @param range_end_index
   * @param sequence_string
   * @param i
   * @param query_codons
   * @param this_codon_flags
   * @param ncurrent_byte
   * @param bit_position
   */
  private void setCache(int range_start_index, 
                        int range_end_index, 
                        char[] sequence_string, 
                        int i,
                        final StringVector query_codons,
                        final byte[] this_codon_flags,
                        int ncurrent_byte,
                        int bit_position)
  {
    // test if stop (or start) codon
    boolean ismatch = false;
    
    // forward codon
    if(i < range_end_index-1)
      if(query_codons == null)
        ismatch = isStopCodon(sequence_string[i-range_start_index],
                              sequence_string[i-range_start_index+1],
                              sequence_string[i-range_start_index+2]);
      else
        ismatch = isCodon(sequence_string[i-range_start_index],
                          sequence_string[i-range_start_index+1],
                          sequence_string[i-range_start_index+2],
                          query_codons);

    if(ismatch)
    {
      this_codon_flags[ncurrent_byte] =                // forward strand stop/start = 2
             (byte)(this_codon_flags[ncurrent_byte] 
                     | (0x0002 << 2*bit_position));
    }
    else
    {
      this_codon_flags[ncurrent_byte] =                // cached no stop/start = 1
        (byte)(this_codon_flags[ncurrent_byte] 
               | (0x0001 << 2*bit_position));
    }

    // reverse codon
    ismatch = false;
    if(i-range_start_index > 1 && i-range_start_index < sequence_string.length)
      if(query_codons == null)
        ismatch = isStopCodon(complement(sequence_string[i-range_start_index]),
                              complement(sequence_string[i-range_start_index-1]),
                              complement(sequence_string[i-range_start_index-2]));
      else
        ismatch = isCodon(complement(sequence_string[i-range_start_index]),
                          complement(sequence_string[i-range_start_index-1]),
                          complement(sequence_string[i-range_start_index-2]),
                          query_codons);
    if(ismatch)
      this_codon_flags[ncurrent_byte] =                // reverse strand stop/start = 3
             (byte)(this_codon_flags[ncurrent_byte] 
                     | (0x0003 << 2*bit_position));
  }

  /**
   *  Return the base at the given position.
   **/
  public char getBaseAt (final int position)
      throws OutOfRangeException 
  {
    if(position > getLength()) 
      throw new OutOfRangeException(position + " > " + getLength());
    
    if(position < 1) 
      throw new OutOfRangeException(position + " < " + 1);
    
    return getSequence().charAt(position);
  }

  /**
   *  Return a sub sequence of the bases from this object.
   *  @param range The range of the bases to be extracted.
   *  @param direction The direction of the returned sequence.  If FORWARD the
   *    sub sequence will be as expected, if REVERSE it will be reverse
   *    complemented.
   *  @return The extracted sequence, which will include the end bases of the
   *    range.
   **/
  public String getSubSequence (final Range range, final int direction) {
    final Range real_range;

    if(direction == FORWARD) 
      real_range = range;
    else 
      real_range = complementRange (range);
    
    // we need to make sure that we pass in-range coordinates to
    // Sequence.getSubSequence()
    final int sub_seq_start_index;
    final int sub_seq_end_index;

    if(real_range.getStart () < 1)
      sub_seq_start_index = 1;
    else
      sub_seq_start_index = real_range.getStart ();

    if(real_range.getEnd () > getLength ()) 
      sub_seq_end_index = getLength ();
    else
      sub_seq_end_index = real_range.getEnd ();

    String sub_sequence = 
      getSequence().getSubSequence(sub_seq_start_index, sub_seq_end_index);  

    // sanity checks - if the user asks for more bases than we
    // have, we return the symbol "@" for the out-of-range bases.
    if (real_range.getStart () < 1) {
      final int dummy_base_count = 1 - real_range.getStart ();
      final char [] dummy_bases = new char [dummy_base_count];

      for (int i = 0 ; i < dummy_base_count ; ++i) {
        dummy_bases[i] = '@';
      }

      sub_sequence = new String (dummy_bases) + sub_sequence;
    }

    if (real_range.getEnd () > getLength ()) {
      final int dummy_base_count = real_range.getEnd () - getLength ();
      final char [] dummy_bases = new char [dummy_base_count];

      for (int i = 0 ; i < dummy_base_count ; ++i) {
        dummy_bases[i] = '@';
      }

      sub_sequence = sub_sequence + new String (dummy_bases);
    }

    if (FORWARD == direction) {
      return sub_sequence;
    } else {
      return reverseComplement (sub_sequence);
    }
  }

  public char[] getSubSequenceC(final Range range, final int direction)
  {
    final Range real_range;

    if(direction == FORWARD)  
      real_range = range;
    else   
      real_range = complementRange (range);
    
    // we need to make sure that we pass in-range coordinates to
    // Sequence.getSubSequence()
    final int sub_seq_start_index;
    final int sub_seq_end_index;

    if(real_range.getStart () < 1) 
      sub_seq_start_index = 1;
    else 
      sub_seq_start_index = real_range.getStart ();

    if(real_range.getEnd () > getLength ()) 
      sub_seq_end_index = getLength ();
    else
      sub_seq_end_index = real_range.getEnd ();

    char[] sub_sequence = 
      getSequence().getCharSubSequence(sub_seq_start_index, sub_seq_end_index);

    if(real_range.getStart() < 1) 
    {
      final int dummy_base_count = 1 - real_range.getStart();
      final char[] dummy_bases = new char[dummy_base_count+sub_sequence.length];

      for(int i = 0; i < dummy_base_count; ++i)
        dummy_bases[i] = '@';

      System.arraycopy(sub_sequence, 0, dummy_bases, dummy_base_count, sub_sequence.length);
      sub_sequence = dummy_bases;
    }

    if(real_range.getEnd() > getLength()) 
    {
      final int dummy_base_count = real_range.getEnd() - getLength();
      final char[] dummy_bases = new char[dummy_base_count+sub_sequence.length];

      for(int i = sub_sequence.length; i < dummy_bases.length; ++i)
        dummy_bases[i] = '@';

      System.arraycopy(sub_sequence, 0, dummy_bases, 0, sub_sequence.length);
      sub_sequence = dummy_bases;
    }


    if(FORWARD == direction)
      return sub_sequence;
    else
      return reverseComplement(sub_sequence);
  }

  /**
   *  This method truncates the sequence use the start and end of the argument.
   *  @param constraint This contains the start and end base of the new
   *    sequence.
   *  @return the Bases truncated into the new coordinate system.
   **/
  public Bases truncate (final Range constraint) {
    final String bases_string = getSubSequence (constraint, FORWARD);

    final Sequence new_sequence = new EmblStreamSequence (bases_string);

    return new Bases (new_sequence);
  }


  /**
  *
  * Reverse complement a range of the sequence.
  *
  */
  public void reverseComplement(final Feature feature)
              throws ReadOnlyException 
  {
    stop_codon_cache = null;
   
    final Range range = feature.getMaxRawRange();
    final int range_start_index = range.getStart();
    final int range_end_index   = range.getEnd();

    // ensure we just get subsequence of interest
    ((StreamSequence)getSequence()).forceReset();
    // sequence to reverse complement
    final char[] sub_sequence = reverseComplement(getSequence().getCharSubSequence(
                                              range_start_index, range_end_index));  
    final char[] new_sequence = new char[getLength()];
    final char[] old_sequence = ((StreamSequence)getSequence()).getCharSequence();

//  System.out.println("range_start_index  "+range_start_index);
//  System.out.println("range_end_index    "+range_end_index);
//  System.out.println("getLength          "+getLength());
//  System.out.println("sub_sequence.length "+sub_sequence.length);
//  System.out.println(feature.getEntry().getEMBLEntry().toString());
//  System.out.println(new String(sub_sequence));

    // if not first contig
    if(range_start_index != 1)
      System.arraycopy(old_sequence, 0, new_sequence, 0, range_start_index-1);

    // copy in new sequence fragment that has been reverse complemented
    System.arraycopy(sub_sequence, 0, new_sequence, range_start_index-1, 
                                                    sub_sequence.length);

    // if not last contig
    if(range_end_index != getLength())
      System.arraycopy(old_sequence, range.getEnd(), new_sequence, range_end_index,
                                                    getLength()-range_end_index);
                                                     
    try 
    {
      embl_sequence.setFromChar(new_sequence);
    } 
    catch (IllegalSymbolException e) 
    {
      throw new Error ("internal error - unexpected exception: " + e);
    }

    final SequenceChangeEvent event =
      new SequenceChangeEvent(this, SequenceChangeEvent.CONTIG_REVERSE_COMPLEMENT,
                              range, sub_sequence.length);

    fireSequenceChangeEvent(event);
  }

  
  public void contigRearrange(final Feature feature, final int new_base_pos)
              throws ReadOnlyException
  {
    stop_codon_cache = null;
  
    final Range range = feature.getMaxRawRange();
    final int range_start_index = range.getStart();
    final int range_end_index   = range.getEnd();

    if(new_base_pos == range_start_index)
      return;

    final char[] new_sequence = new char[getLength()];
    final char[] old_sequence = ((StreamSequence)getSequence()).getCharSequence();

    int contig_length = 0;
    if(new_base_pos < range_start_index)
    {
      // if not first contig
      if(new_base_pos != 1)
        System.arraycopy(old_sequence, 0, new_sequence, 0, new_base_pos-1);

      contig_length = range_end_index - range_start_index + 1;
      // copy in new sequence fragment that has been reverse complemented
      System.arraycopy(old_sequence, range_start_index-1, 
                       new_sequence, new_base_pos-1, contig_length);

      System.arraycopy(old_sequence, new_base_pos-1,
                       new_sequence, new_base_pos+contig_length-1, 
                       range_start_index-new_base_pos);

      // if not last contig
      if(new_base_pos < getLength()+1)
        System.arraycopy(old_sequence, range_end_index,
                         new_sequence, range_end_index,
                         getLength()-range_end_index);
    }
    else 
    {
      System.arraycopy(old_sequence, 0, new_sequence, 0, range_start_index-1);

      System.arraycopy(old_sequence, range_end_index,
                       new_sequence, range_start_index-1, 
                       new_base_pos-range_end_index-1);

      System.arraycopy(old_sequence, range_start_index-1,
                       new_sequence, (range_start_index-1)+(new_base_pos-range_end_index)-1,
                       range_end_index - range_start_index + 1);

      // if not last contig
      if(new_base_pos < getLength()+1)
        System.arraycopy(old_sequence, new_base_pos-1,
                         new_sequence, new_base_pos-1,
                         getLength()-new_base_pos+1);
    }
    

    try 
    { 
      embl_sequence.setFromChar(new_sequence);
    }
    catch (IllegalSymbolException e)  
    {
      throw new Error ("internal error - unexpected exception: " + e);
    }

    final SequenceChangeEvent event =
      new SequenceChangeEvent(SequenceChangeEvent.CONTIG_REORDER,
                              new_base_pos, range);

    fireSequenceChangeEvent(event);
  }


  /**
   *  Delete the bases in the given range and send out a SequenceChange event
   *  to all the listeners.
   *  @param range The inclusive range of bases to delete.
   *  @return A String containing the deleted bases.
   *  @exception ReadOnlyException If this Bases object cannot be changed.
   **/
  public String deleteRange (final Range range)
      throws ReadOnlyException {
    stop_codon_cache = null;

    final String removed_bases =
      getSequence ().getSubSequence (range.getStart (), range.getEnd ());

    final String new_sequence =
      getSequence ().getSubSequence (1, range.getStart () - 1) +
      getSequence ().getSubSequence (range.getEnd () + 1,
                                     embl_sequence.length ());

    try {
      embl_sequence.setFromChar(new_sequence.toCharArray());
    } catch (IllegalSymbolException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    } 

    final SequenceChangeEvent event =
      new SequenceChangeEvent (this,
                               SequenceChangeEvent.DELETION,
                               range.getStart (),
                               removed_bases);

    fireSequenceChangeEvent (event);

    return removed_bases;
  }

  /**
   *  Insert the given bases at the given base position and send out a
   *  SequenceChange event to all the listeners.
   *  @param position The bases are inserted just before this base position if
   *    direction is FORWARD or just after if direction is REVERSE.
   *  @param direction If this is FORWARD, then the bases is the bases String
   *    will be inserted just before the base given by position.  If this is
   *    REVERSE the bases will be reversed, complemented and inserted just
   *    after the position.
   *  @param bases The bases to add (or the reverse complement of the bases to
   *    add if direction is REVERSE).
   *  @exception ReadOnlyException If this Bases object cannot be changed.
   **/
  public void addBases (final int position, final int direction,
                        final String bases)
      throws ReadOnlyException, IllegalSymbolException {
    stop_codon_cache = null;

    final String new_sequence;
    final int real_position;
    final String real_bases;

    if (direction == FORWARD) {
      real_position = position;
      real_bases = bases.toLowerCase ();
    } else {
      real_position = position + 1;
      real_bases = reverseComplement (bases.toLowerCase ());
    }

    new_sequence =
      getSequence ().getSubSequence (1, real_position - 1) +
      real_bases +
      getSequence ().getSubSequence (real_position, getLength ());

    getSequence ().setFromChar(new_sequence.toCharArray());

    final SequenceChangeEvent event =
      new SequenceChangeEvent (this,
                               SequenceChangeEvent.INSERTION,
                               real_position,
                               real_bases);

    fireSequenceChangeEvent (event);

    return;
  }

  /**
   *  There is one element in this array for each possible
   *  SequenceChangeListener priority.  This array is changed by
   *  addSequenceChangeListener() and removeSequenceChangeListener().
   **/
  final private WeakHashMap listener_hash_map_array [] =
    new WeakHashMap [MAX_PRIORITY - MIN_PRIORITY + 1];

  /**
   *  Adds the specified event listener to the list of object that receive
   *  sequence change events from this object.
   *  @param l the event change listener.
   *  @param priority The listeners are stored in a priority queue using this
   *    value.  Larger priority means that the listener will receive the event
   *    sooner (than lower priority listeners).  Values less than MIN_PRIORITY
   *    will be treated like MIN_PRIORITY values higher than MAX_PRIORITY will
   *    be treated like MAX_PRIORITY.
   **/
  public void addSequenceChangeListener (final SequenceChangeListener l,
                                         int priority) {
    if (priority < MIN_PRIORITY) {
      priority = MIN_PRIORITY;
    }

    if (priority > MAX_PRIORITY) {
      priority = MAX_PRIORITY;
    }

    listener_hash_map_array [priority - MIN_PRIORITY].put (l, null);
  }

  /**
   *  Removes the specified event listener so that it no longer receives
   *  sequence change events from this object.
   *  @param l the event change listener.
   **/
  public void removeSequenceChangeListener (final SequenceChangeListener l) {
    for (int i = 0 ; i < listener_hash_map_array.length ; ++i) {
      final WeakHashMap this_hash_map = listener_hash_map_array [i];

      if (this_hash_map.containsKey (l)) {
        this_hash_map.remove (l);
        return;
      }
    }
  }

  /**
   *  Send a SequenceChangeEvent to each object that is listening for it.
   **/
  private void fireSequenceChangeEvent (final SequenceChangeEvent event) {
    for (int i = listener_hash_map_array.length - 1 ; i >= 0 ; --i) {
      final WeakHashMap this_hash_map = listener_hash_map_array [i];

      if (this_hash_map != null) {
        final Iterator iter = this_hash_map.keySet ().iterator ();

        while (iter.hasNext()) 
        {
          final SequenceChangeListener this_listener =
            (SequenceChangeListener) iter.next();
          this_listener.sequenceChanged (event);
        }
      }
    }
  }

  /**
   *  Return the average gc percent for the sequence.
   **/
  public float getAverageGCPercent () {
    return ((float)(getSequence ().getCCount () +
                    getSequence ().getGCount ())) /
      getSequence ().length () * 100;
  }

  /**
   *  Return the average AG percent for the sequence as a percentage.
   **/
  public float getAverageAGPercent () {
    return ((float)(getSequence ().getACount () +
                    getSequence ().getGCount ())) /
      getSequence ().length () * 100;
  }

  /**
   *  Return the number of 'A's in this Bases object.
   **/
  public int getACount () {
    return getSequence ().getACount ();
  }

  /**
   *  Return the number of 'T's in this Bases object.
   **/
  public int getTCount () {
    return getSequence ().getTCount ();
  }

  /**
   *  Return the number of 'G's in this Bases object.
   **/
  public int getGCount () {
    return getSequence ().getGCount ();
  }

  /**
   *  Return the number of 'C's in this Bases object.
   **/
  public int getCCount () {
    return getSequence ().getCCount ();
  }

  /**
   *  Return a String containing the reverse complement of the argument
   *  String.  For example an argument of "aatc" will result in "gatt".
   **/
  public static String reverseComplement (final String sequence_string) {
    StringBuffer return_buffer = new StringBuffer (sequence_string.length ());

    for (int i = sequence_string.length () - 1 ; i >= 0 ; --i) {
      return_buffer.append (complement (sequence_string.charAt (i)));
    }

    return return_buffer.toString ();
  }

  /**
   *  Return a char[] containing the reverse complement of the argument
   *  String.  For example an argument of "aatc" will result in "gatt".
   **/
  public static char[] reverseComplement (final char[] sequence_char) 
  {
    final int length = sequence_char.length;
    final char[] return_sequence = new char[length];
    int j = 0;

    for(int i = length - 1 ; i >= 0 ; --i) 
    {
      return_sequence[j] = complement(sequence_char[i]);
      j++;
    }

    return return_sequence;
  }


  /**
   *  Return a String containing the complement of the argument String.  For
   *  example an argument of "aatc" will result in "ttag".
   **/
  public static String complement (final String sequence_string) {
    StringBuffer return_buffer = new StringBuffer (sequence_string.length ());

    for (int i = 0 ; i < sequence_string.length () ; ++i) {
      return_buffer.append (complement (sequence_string.charAt (i)));
    }

    return return_buffer.toString ();
  }

  /**
   *  Return a char array containing the complement of the argument char[].  For
   *  example an argument of "aatc" will result in "ttag".
   **/
  public static char[] complement (final char sequence[])
  {
    final char[] seq_comp = new char[sequence.length];

    for (int i = 0 ; i < sequence.length; ++i)
      seq_comp[i] = complement(sequence[i]);

    return seq_comp;
  }
  
  /**
   *  Returns the complement base of it's argument - c for g, a for t etc.
   *  The argument may be upper or lower case, but the result is always lower
   *  case.  This also works for IUB base codes: the complement of 'y' is 'r'
   *  because 'y' is 'c' or 't' and 'r' is 'a' or 'g', the complement of 'n'
   *  or 'x' (any base) is 'n'.
   **/
  public final static char complement (final char base) {

    switch (base) {
    case 'a': case 'A': return 't';
    case 't': case 'T': case 'u': case 'U': return 'a';
    case 'g': case 'G': return 'c';
    case 'c': case 'C': return 'g';
    case 'r': case 'R': return 'y';
    case 'y': case 'Y': return 'r';
    case 'k': case 'K': return 'm';
    case 'm': case 'M': return 'k';
    case 's': case 'S': return 's';
    case 'w': case 'W': return 'w';
    case 'b': case 'B': return 'v';
    case 'd': case 'D': return 'h';
    case 'h': case 'H': return 'd';
    case 'v': case 'V': return 'b';
    case 'n': case 'N': return 'n';
    case 'x': case 'X': return 'x';
    default:
      return '@';
//      throw new Error ("in Bases.complement - tried to complement a letter " +
//                       "that isn't a base");
    }
  }

  /**
   *  Return the Sequence object that was passed to the constructor.
   **/
  public Sequence getSequence () 
  {
    return embl_sequence;
  }

  
  /**
   *  Check a three character substring and return true if and only if the
   *  three bases translate to a stop codon.  If the direction is REVERSE
   *  then the three bases to check are at start_index, start_index - 1 and
   *  start_index - 2.  In that case true is returned if and only the
   *  complement of those three bases is a stop codon.
   *  Codons that contain an X are considered to be stop codons.
   **/
  private static boolean isCodon(char first_letter, char second_letter, char third_letter,
                                 final StringVector query_codons)
  {
    char[] tran = {first_letter, second_letter, third_letter };
    
    if(query_codons.contains( new String(tran) ))
      return true;
    
    return false;
  }
  
  /**
   *  Check a three character substring and return true if and only if the
   *  three bases translate to a stop codon.  If the direction is REVERSE
   *  then the three bases to check are at start_index, start_index - 1 and
   *  start_index - 2.  In that case true is returned if and only the
   *  complement of those three bases is a stop codon.
   *  Codons that contain an X are considered to be stop codons.
   **/
  private static boolean isStopCodon(char first_letter, char second_letter, char third_letter)
  {
    // codons that contain an X are considered to be stop codons.
    if(first_letter == 'x' || second_letter == 'x' || third_letter == 'x')
      return true;

    final char translation = AminoAcidSequence.getCodonTranslation(first_letter,
                                                                  second_letter,
                                                                  third_letter);

    if(translation == '+' || translation == '*' || translation == '#') 
      return true;
    else 
      return false;
  }


  /**
   *  Check a three character substring and return true if and only if the
   *  three bases are legal (see isLegalBase ()).
   **/
  /*private static boolean isLegalCodon (final String sequence_string,
                                       final int start_index,
                                       final int direction) {
    if (direction == FORWARD) {
      if (isLegalBase (sequence_string.charAt (start_index)) &&
          isLegalBase (sequence_string.charAt (start_index + 1)) &&
          isLegalBase (sequence_string.charAt (start_index + 2))) {
        return true;
      }
    } else {
      if (isLegalBase (sequence_string.charAt (start_index)) &&
          isLegalBase (sequence_string.charAt (start_index - 1)) &&
          isLegalBase (sequence_string.charAt (start_index - 2))) {
        return true;
      }
    }

    // this isn't a stop codon
    return false;

  }*/

  /**
   *  Return true if and only if the given base character is one of 'a', 't',
   *  'c', 'g' or 'u'.
   **/
  public final static boolean isLegalBase (final char base_char) {
    switch (base_char) {
    case 'a': case 'A': return true;
    case 't': case 'T': return true;
    case 'u': case 'U': return true;
    case 'g': case 'G': return true;
    case 'c': case 'C': return true;
    default:
      return false;
    }
  }

  /**
   *  The underlying sequence object that holds the data for this object.
   *  This is the same object that was passed to the constructor.
   **/
  private Sequence embl_sequence;

  /**
   *  The object representing the forward sequence of bases.
   **/
  private Strand forward_strand;

  /**
   *  The object representing the reverse (reverse complemented)
   *  sequence of bases.
   **/
  private Strand reverse_strand;
}
