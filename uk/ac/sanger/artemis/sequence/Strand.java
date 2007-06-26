/* Strand.java
 *
 * created: Sun Oct 11 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/sequence/Strand.java,v 1.7 2007-06-26 13:59:04 tjc Exp $
 */

package uk.ac.sanger.artemis.sequence;

import uk.ac.sanger.artemis.util.*;
import uk.ac.sanger.artemis.io.Range;

import org.biojava.bio.symbol.IllegalSymbolException;

/**
 *  This represents one strand of DNA.  The bases of the strand must be passed
 *  in to the constructor as a Bases object.  All operations where the
 *  direction of reading matters should be handled by this class.  After the
 *  constructor has been called with the correct direction (by the Bases
 *  constructor) other classes should be able to use this class without caring
 *  what direction the Strand represents.
 *
 *  @author Kim Rutherford
 *  @version $Id: Strand.java,v 1.7 2007-06-26 13:59:04 tjc Exp $
 **/

public class Strand {
  /**
   *  When passed to the constructor, FORWARD indicates the bases should be
   *  read in the forward direction.
   **/
  static public final int FORWARD = Bases.FORWARD;

  /**
   *  When passed to the constructor, REVERSE indicates the bases should be
   *  read in the forward direction.
   **/
  static public final int REVERSE = Bases.REVERSE;

  /**
   *  Create a new Strand object with an associated direction.  Typically this
   *  constructor is called twice with the sames Bases object, but with a
   *  different value for direction.
   *  @param bases The bases of this strand.
   *  @param direction The direction the bases should read.  If REVERSE then
   *    the bases will be automatically read in the reverse direction and
   *    complemented before use.
   **/
  Strand (final Bases bases, final int direction) {
    this.bases = bases;
  }

  /**
   *  Return the Bases reference that was passed to the constructor.
   **/
  public Bases getBases () {
    return bases;
  }

  /**
   *  Return a String containing the bases of this strand.
   **/
  public String getStrandBases () {
    if (getDirection () == FORWARD) {
      return getBases ().toString ();
    } else {
      return Bases.reverseComplement (getBases ().toString ());
    }
  }

  /**
   *  Return the direction of this strand.
   **/
  public int getDirection () {
    if (getBases ().getForwardStrand () == this) {
      return FORWARD;
    } else {
      return REVERSE;
    }
  }

  /**
   *  Return true if and only if this Strand is a forward strand.
   **/
  public boolean isForwardStrand () {
    if (getDirection () == FORWARD) {
      return true;
    } else {
      return false;
    }
  }

  /**
   *  Return the length of this Strand in bases.
   **/
  public int getSequenceLength () {
    return getBases ().getLength ();
  }

  /**
   *  Return an 2D array containing the stop or start codons in a range for
   *  all 3 frames of the strand.  
   *  @param range The inclusive range of bases to get the stop codons from.
   *  @param query_codons if this is NULL then this assumes we are looking
   *    for stop codons, otherwise this is used to look for start codons.
   **/
  public int [][] getStopOrStartCodons(final Range range, final StringVector query_codons) {
    return bases.getStopOrStartCodons(range, getDirection (), query_codons);
  }

  private int [] getStopCodons (final Range range) {
    return bases.getStopCodons (range, getDirection ());
  }


  /**
   *  Return an array containing the positions of the codons that match the
   *  strings given by the query_codons argument.  Only those codons that are
   *  in the same frame as the first base of the range are returned.
   *  @param range The inclusive range of bases to get the codons from.
   *  @param query_codons The codons to search for.  Each element of this
   *    vector should be a string that is 3 characters long.
   *  @return An array containing the positions of the first base of the
   *    codons.  This array is padded with zeros at the end.
   **/
  public int [] getMatchingCodons (final Range range,
                                   final StringVector query_codons) {
    return bases.getMatchingCodons (range, getDirection (), query_codons);
  }

  /**
   *  Return a String containing the bases (base letters) of the codon that
   *  starts at the given Marker.
   **/
  public static String getCodonAtMarker (final Marker marker) {
    try {
      final Range codon_range =
        new Range (marker.getPosition (), marker.getPosition () + 2);

      return marker.getStrand ().getSubSequence (codon_range);
    } catch (OutOfRangeException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    }
  }

  /**
   *  Return an array of MarkerRange objects - one object for each open
   *  reading frame bigger than the given size.
   *  @param search_range All the returned MarkerRanges must overlap this
   *    range.
   *  @param minimum_size All the returned MarkerRanges must be at least
   *    this many residues long.
   **/
  public static MarkerRange [] getOpenReadingFrameRanges (final MarkerRange search_range,
                                                          final int minimum_size,
                                                          final int sequence_end,
                                                          final int sequence_start) {
    final Strand search_strand = search_range.getStrand ();

    final MarkerRange [] [] frame_ranges = new MarkerRange [3] [];

    // get an array of MarkerRange objects for each frame

    for (int frame_offset = 0 ; frame_offset < 3 ; ++frame_offset) {
      final Range sequence_range;

      try {
        sequence_range =
          new Range (sequence_start + frame_offset, search_strand.getSequenceLength ());
      } catch (OutOfRangeException e) {
        throw new Error ("internal error - unexpected exception: " + e);
      }

      final int [] frame_stop_codons =
        search_strand.getStopCodons (sequence_range);

      frame_ranges [frame_offset] =
        search_strand.getORFsFromStopCodons (frame_stop_codons,
                                             minimum_size,
                                             frame_offset,
                                             search_range,
                                             sequence_end, 
                                             sequence_start);
    }

    // now copy the MarkerRange objects into a single array to return

    final int max_range_count =
      frame_ranges[0].length + frame_ranges[1].length + frame_ranges[2].length;

    final MarkerRange [] temp_range_array = new MarkerRange [max_range_count];

    int temp_range_array_index = 0;

    for (int i = 0 ; i < 3 ; ++i) {
      final MarkerRange [] frame_range_array = frame_ranges[i];

      for (int range_index = 0 ;
           range_index < frame_range_array.length &&
             frame_ranges[i][range_index] != null ;
           ++range_index) {

        final MarkerRange this_range = frame_ranges[i][range_index];

        if (search_range.overlaps (this_range)) {
          temp_range_array[temp_range_array_index] = this_range;

          ++temp_range_array_index;
        } else {
          // ignore this range
          continue;
        }
      }
    }

    final MarkerRange [] return_array =
      new MarkerRange [temp_range_array_index];

    System.arraycopy (temp_range_array, 0,
                      return_array, 0,
                      temp_range_array_index);

    return return_array;
  }

  /**
   *  Find an ORF around the given Marker object.  The Strand and frame to
   *  search are read from the Marker.
   *  @param search_marker Search from this marker.  The frame to search is
   *    the frame in which this marker point to the first base of a codon.
   *  @param ignore_illegal_codons If true this method will read through
   *    illegal bases (like N and the other ambiguity codes).  If false codons
   *    containing illegal bases will be treated like stop codons.
   *  @return A MarkerRange for the ORF containing search_marker or null is
   *    the marker points to a stop codon.
   **/
  public static MarkerRange getORFAroundMarker (final Marker search_marker,
                                                final boolean
                                                  ignore_illegal_codons) {
    // first search backwards

    final Range search_codon_range;

    try {
      search_codon_range = new Range (search_marker.getPosition (),
                                      search_marker.getPosition () + 2);
    } catch (OutOfRangeException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    }

    final Strand search_strand = search_marker.getStrand ();

    final String search_codon_sequence =
      search_strand.getSubSequence (search_codon_range);

    final char search_codon_char =
      AminoAcidSequence.getCodonTranslation (search_codon_sequence);

    if (AminoAcidSequence.isStopCodon (search_codon_char)) {
      return null;
    }

    final Marker start_marker = getStartOfORF (search_marker,
                                               ignore_illegal_codons);

    final Marker end_marker = getEndOfORF (search_marker,
                                           ignore_illegal_codons);

    try {
      return new MarkerRange (search_strand,
                              start_marker.getPosition (),
                              end_marker.getPosition ());
    } catch (OutOfRangeException e) {
      throw new Error ("internal error - unexpected OutOfRangeException");
    }
  }

  /**
   *  Return the Marker of the start of the ORF that contains the given
   *  Marker.  This method works by stepping backwards three bases at a time
   *  until the next codon is stop codon or until the start of sequence is
   *  reached.
   *  @param ignore_illegal_codons If true this method will read through
   *    illegal bases (like N and the other ambiguity codes).  If false codons
   *    containing illegal bases will be treated like stop codons.
   *  @return A Marker that points to a position just after a stop codon or
   *    to the first base of the Strand (in this frame).
   **/
  private static Marker getStartOfORF (final Marker search_marker,
                                       final boolean ignore_illegal_codons) {
    final Strand search_strand = search_marker.getStrand ();

    Marker current_marker = search_marker;

    Marker previous_marker = search_marker;

    while (true) {
      try {
        // this should be optimised to avoid using moveBy() which is slow

        current_marker = current_marker.moveBy (-3);

        final Range search_codon_range =
          new Range (current_marker.getPosition (),
                     current_marker.getPosition () + 2);

        final String codon_sequence =
          search_strand.getSubSequence (search_codon_range);

        if (codon_sequence.charAt (0) == 'x' ||
            codon_sequence.charAt (1) == 'x' ||
            codon_sequence.charAt (2) == 'x') {
          return previous_marker;
        }

        final char this_codon_char =
          AminoAcidSequence.getCodonTranslation (codon_sequence);

        if (AminoAcidSequence.isStopCodon (this_codon_char) ||
            !ignore_illegal_codons &&
            !AminoAcidSequence.isLegalCodon (this_codon_char)) {
          return previous_marker;
        }

        previous_marker = current_marker;
      } catch (OutOfRangeException e) {
        return previous_marker;
      }
    }
  }

  /**
   *  Return the Marker of the end of the ORF that contains the given
   *  Marker.  This method works by stepping forwards three bases at a time
   *  until the next codon is stop codon or until the end of sequence is
   *  reached.
   *  @param ignore_illegal_codons If true this method will read through
   *    illegal bases (like N and the other ambiguity codes).  If false codons
   *    containing illegal bases will be treated like stop codons.
   *  @return A Marker that points to a position just before a stop codon or
   *    to the last base of the Strand (in this frame).
   **/
  private static Marker getEndOfORF (final Marker search_marker,
                                     final boolean ignore_illegal_codons) {
    final Strand search_strand = search_marker.getStrand ();

    Marker current_marker = null;

    try {
      // use the end of the codon as the start position
      current_marker = search_marker.moveBy (2);
    } catch (OutOfRangeException e) {
      // end of sequence
      return search_marker;
    }

    Marker previous_marker = current_marker;

    while (true) {
      try {
        // this should be optimised to avoid using moveBy() which is slow

        current_marker = current_marker.moveBy (3);

        final Range search_codon_range =
          new Range (current_marker.getPosition () - 2,
                     current_marker.getPosition ());

        final String codon_sequence =
          search_strand.getSubSequence (search_codon_range);

        if (codon_sequence.charAt (0) == 'x' ||
            codon_sequence.charAt (1) == 'x' ||
            codon_sequence.charAt (2) == 'x') {
          return previous_marker;
        }

        final char this_codon_char =
          AminoAcidSequence.getCodonTranslation (codon_sequence);

        if (AminoAcidSequence.isStopCodon (this_codon_char) ||
            !ignore_illegal_codons &&
            !AminoAcidSequence.isLegalCodon (this_codon_char)) {
          return previous_marker;
        }

        previous_marker = current_marker;
      } catch (OutOfRangeException e) {
        return previous_marker;
      }
    }
  }

  /**
   *  Return an array of MarkerRange objects - one object for each open
   *  reading frame bigger than the given size.
   *  @param stop_codons An array containing the positions of the stop codons
   *    for one of the three translation frames.
   *  @param minimum_size All the returned ORFs will be at least this many
   *    amino acids long.
   *  @param frame_offset The frame to that the stop codons are in (0, 1 or 2)
   *  @param test_range All the returned MarkerRanges must overlap this
   *    range.
   *  @return An array of MarkerRange objects.  The array may contain some
   *    null references at the end.
   **/
  private MarkerRange [] getORFsFromStopCodons (final int [] stop_codons,
                                                final int minimum_size,
                                                final int frame_offset,
                                                final MarkerRange test_range,
                                                final int sequence_end,
                                                final int sequence_start) {
    // this array is the maximum possible size - the last elements of the
    // array will almost certainly be empty when we return
    final MarkerRange [] return_array =
      new MarkerRange [stop_codons.length + 1];

    int return_array_index = 0;

    // special case for the first ORF, the one between the first base of the
    // sequence and the first stop codon

    for (int i = -1 ;
         i < stop_codons.length && (i == -1 || stop_codons[i] != 0) ;
         ++i) {
      final int first_base_of_range;

      // an index of -1 indicates that the start of this ORF is the first base
      // of the sequence (offset to be in the correct frame).
      if (i == -1) {
        first_base_of_range = sequence_start + frame_offset;
      } else {
        // use the next base after the stop codon as the start of the range
        first_base_of_range = stop_codons[i] + 3;
      }

      if (first_base_of_range >= sequence_end ||
          first_base_of_range < sequence_start) {
        continue;
      }

      int last_base_of_range;

      // the last base in the range is the last base of the next stop codon or
      // the end of sequence (if there are no more stop codons).
      if (i + 1 == stop_codons.length || stop_codons[i + 1] == 0) {
        last_base_of_range = sequence_end;
      } else {
        // use the last base of the next stop codon as the end of the range
        last_base_of_range = stop_codons[i + 1] + 2;
      }

      if (last_base_of_range >= sequence_end) {
        last_base_of_range = sequence_end;
      }
      
      final int aa_count = (last_base_of_range - first_base_of_range) / 3;

      if (aa_count >= minimum_size &&
          last_base_of_range >= test_range.getStart ().getPosition () &&
          first_base_of_range <= test_range.getEnd ().getPosition ()) {
        try {

          return_array[return_array_index] =
            makeMarkerRangeFromPositions (first_base_of_range,
                                          last_base_of_range);

          ++return_array_index;
        } catch (OutOfRangeException e) {
          throw new Error ("internal error - unexpected OutOfRangeException");
        }
      }
    }

    return return_array;
  }

  /**
   *  Create and return a Marker on this Strand at the given position.  The
   *  position should refer to this strand not the underlying Bases object.
   **/
  public Marker makeMarker (int position)
      throws OutOfRangeException {
    final Marker new_marker = new Marker (this, position);

    return new_marker;
  }

  /**
   *  Make a Marker on this Strand from a raw position.  "Raw" means that if
   *  this is the REVERSE Strand then the position refers to a base in the the
   *  underlying sequence.
   **/
  public Marker makeMarkerFromRawPosition (int position)
      throws OutOfRangeException {
    if (getDirection () == FORWARD) {
      return makeMarker (position);
    } else {
      return makeMarker (getBases ().getComplementPosition (position));
    }
  }

  /**
   *  Make a MarkerRange on this Strand from two positions.
   **/
  public MarkerRange makeMarkerRangeFromPositions (int start_position,
                                                   int end_position)
      throws OutOfRangeException {
    return new MarkerRange (this, start_position, end_position);
  }

  /**
   *  Make a MarkerRange on this Strand from two raw positions.  The
   *  positions should be positions on the underlying Bases object rather than
   *  positions on this Strand.
   **/
  public MarkerRange makeMarkerRangeFromRawPositions (int raw_start_position,
                                                      int raw_end_position)
      throws OutOfRangeException {
    if (getDirection () == FORWARD) {
      return new MarkerRange (this, raw_start_position, raw_end_position);
    } else {
      final int real_start_position =
        getBases ().getComplementPosition (raw_start_position);
      final int real_end_position =
        getBases ().getComplementPosition (raw_end_position);
      return new MarkerRange (this, real_start_position, real_end_position);
    }
  }

  /**
   *  Return the position on the Bases object of the argument base.
   **/
  public int getRawPosition (int strand_position) {
    return getBases ().getRawPosition (strand_position, getDirection ());
  }

  /**
   *  Delete the bases in the given MarkerRange and send out a SequenceChange
   *  event to all the listeners.
   *  @exception ReadOnlyException If this Strand object cannot be changed.
   **/
  public static void deleteRange (final MarkerRange range)
      throws ReadOnlyException {
    range.getStrand ().getBases ().deleteRange (range.getRawRange ());
  }

  /**
   *  Add the given bases just before the given Marker position and send out a
   *  SequenceChange event to all the listeners.
   *  @param position The new bases are inserted just before this Marker
   *    position.
   *  @param bases_string The bases to insert.
   *  @exception ReadOnlyException If this Strand object cannot be changed.
   **/
  public static void addBases (final Marker position,
                               final String bases_string)
      throws ReadOnlyException, IllegalSymbolException {
    final Bases bases = position.getStrand ().getBases ();

    if (position.getStrand ().isForwardStrand ()) {
      bases.addBases (position.getRawPosition (), FORWARD, bases_string);
    } else {
      bases.addBases (position.getRawPosition (), REVERSE, bases_string);
    }
  }

  /**
   *  Translate a sequence of bases into the corresponding single letter amino
   *  acid codes.
   *  @param range The range of the bases to translated.  If the range.start
   *    - range.end + 1 is not a multiple of three the last codon is
   *    incomplete and will not be translated.
   *  @param unknown_is_x If this parameter is true codons that contain
   *    ambiguous bases will be translated as 'x', if false they will be
   *    translated as '.'
   *  @return The translated sequence in one letter abbreviated form.
   **/
  public AminoAcidSequence getTranslation(final Range range,
                                          final boolean unknown_is_x) {
    return getBases ().getTranslation (range, getDirection (), unknown_is_x);
  }

  public AminoAcidSequence getSpacedTranslation(final Range range,
                                          final boolean unknown_is_x) {
    return getBases ().getSpacedTranslation (range, getDirection (), unknown_is_x);
  }


  /**
   *  Return a sub-sequence of bases from this strand.
   *  @param range The inclusive range of bases to return.
   **/
  public String getSubSequence (Range range) {
    return getBases ().getSubSequence (range, getDirection ());
  }

  /**
   *  Return a sub-sequence of bases from this Bases object that underlies
   *  this Strand object.  
   *  @param range The inclusive range of bases to return.
   **/
  public char[] getRawSubSequenceC (Range range) {
    return getBases ().getSubSequenceC (range, FORWARD);
  }
  
  /**
   *  Return a sub-sequence of bases from this Bases object that underlies
   *  this Strand object.  This returns the same as getSubSequence () for
   *  FORWARD strands.
   *  @param range The inclusive range of bases to return.
   **/
  public String getRawSubSequence (Range range) {
    return getBases ().getSubSequence (range, FORWARD);
  }

  /**
   *  Return the bases referenced by the given MarkerRange.
   **/
  public static String markerRangeBases (final MarkerRange marker_range) {
    return marker_range.getStrand ().getSubSequence (marker_range.getRange ());
  }

  /**
   *  Return the number of 'A's on this Strand.
   **/
  public int getACount () {
    if (isForwardStrand ()) {
      return getBases ().getACount ();
    } else {
      return getBases ().getTCount ();
    }
  }

  /**
   *  Return the number of 'T's on this Strand.
   **/
  public int getTCount () {
    if (isForwardStrand ()) {
      return getBases ().getTCount ();
    } else {
      return getBases ().getACount ();
    }
  }

  /**
   *  Return the number of 'G's on this Strand.
   **/
  public int getGCount () {
    if (isForwardStrand ()) {
      return getBases ().getGCount ();
    } else {
      return getBases ().getCCount ();
    }
  }

  /**
   *  Return the number of 'C's on this Strand.
   **/
  public int getCCount () {
    if (isForwardStrand ()) {
      return getBases ().getCCount ();
    } else {
      return getBases ().getGCount ();
    }
  }

  /**
   *  The reference the was passed to the constructor.  This object does most
   *  of the real work for us, such as reverse complementing and translating.
   **/
  private Bases bases;
}


