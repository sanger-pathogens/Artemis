/* MarkerRange.java
 *
 * created: Tue Dec 29 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/sequence/MarkerRange.java,v 1.1 2004-06-09 09:52:20 tjc Exp $
 */

package uk.ac.sanger.artemis.sequence;

import uk.ac.sanger.artemis.util.*;

import uk.ac.sanger.artemis.io.Range;

/**
 *  This class implements a range that has end points specified by a pair of
 *  positions in a Strand object.
 *
 *  @author Kim Rutherford
 *  @version $Id: MarkerRange.java,v 1.1 2004-06-09 09:52:20 tjc Exp $
 **/

public class MarkerRange {
  /**
   *  Create a new MarkerRange object.  If the start position is greater than
   *  the end position the two positions will be silently swapped.
   *  @param strand The Strand that contains this range.
   *  @param start The start position.
   *  @param end The end position.
   **/
  public MarkerRange (final Strand strand, final int start, final int end)
      throws OutOfRangeException {
    if (start <= end) {
      this.start = strand.makeMarker (start);
      this.end = strand.makeMarker (end);
    } else {
      this.start = strand.makeMarker (end);
      this.end = strand.makeMarker (start);
    }
  }

  /**
   *  Create a new MarkerRange object that covers one base (given by a Marker
   *  object).
   *  @param marker The position of the new MarkerRange.
   **/
  public MarkerRange (final Marker marker) {
    this.start = marker;
    this.end = marker;
  }

  /**
   *  Return the Strand object that was passed to the constructor.
   **/
  public Strand getStrand () {
    return getStart ().getStrand ();
  }

  /**
   *  Return true if and only if this MarkerRange is a MarkerRange on the
   *  forward strand.
   **/
  public boolean isForwardMarker () {
    return getStrand ().isForwardStrand ();
  }

  /**
   *  Return the start Marker of this range.
   **/
  public Marker getStart () {
    return start;
  }

  /**
   *  Return the end Marker of this range.
   **/
  public Marker getEnd () {
    return end;
  }

  /**
   *  Return the number of bases in the range, inclusive of the end points.
   **/
  public int getCount () {
    return
      getRawEnd ().getRawPosition () - getRawStart ().getRawPosition ()  + 1;
  }
  
  /**
   *  Return the end marker if it marks a lower position on the underlying
   *  Bases of the Strand than the start Marker, otherwise return the start
   *  Marker.
   **/
  public Marker getRawStart () {
    if (getStart ().getRawPosition () <= getEnd ().getRawPosition ()) {
      return getStart ();
    } else {
      return getEnd ();
    }
  }
  
  /**
   *  Return the start marker if it marks a lower position on the underlying
   *  Bases of the Strand than the end Marker, otherwise return the end
   *  Marker.
   **/
  public Marker getRawEnd () {
    if (getEnd ().getRawPosition () <= getStart ().getRawPosition ()) {
      return getStart ();
    } else {
      return getEnd ();
    }
  }

  /**
   *  Return an embl.Range object that corresponds a range of bases on the
   *  underlying Bases object of this MarkerRange.
   **/
  public Range getRawRange () {
    try {
      if (getStrand ().isForwardStrand ()) {
        return new Range (getStart ().getRawPosition (),
                          getEnd ().getRawPosition ());
      } else {
        return new Range (getEnd ().getRawPosition (),
                          getStart ().getRawPosition ());
      }
    } catch (OutOfRangeException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    }
  }

  /**
   *  Return an embl.Range object for this MarkerRange.  The start and end
   *  positions of the Range will be the same as the start and end that was
   *  passed to the constructor.
   **/
  public Range getRange () {
    try {
      return new Range (getStart ().getPosition (), getEnd ().getPosition ());
    } catch (OutOfRangeException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    }
  }

  /**
   *  Return true if and only if this MarkerRange and the given MarkerRange
   *  overlap.
   **/
  public boolean overlaps (final MarkerRange arg_range) {
    if (getStart ().getPosition () < arg_range.getStart ().getPosition () &&
        getEnd ().getPosition () < arg_range.getStart ().getPosition ()) {
      return false;
    }
    if (getStart ().getPosition () > arg_range.getEnd ().getPosition () &&
        getEnd ().getPosition () > arg_range.getEnd ().getPosition ()) {
      return false;
    }
    return true;
  }

  /**
   *  Return a new MarkerRange object that is the same as this object but has
   *  been extended to contain the given range.
   *  @param arg_range The returned MarkerRange will cover this range too.
   **/
  public MarkerRange extendRange (final MarkerRange arg_range) {
    return combineRanges (arg_range, true);
  }

  /**
   *  Return a new MarkerRange object that is the same as this object but is
   *  guaranteed to contain the given range.
   *  @param arg_range The returned MarkerRange will cover this range too.
   *  @param truncate_if_contained If the arg_range is completely contained in
   *    this range then the returned range will be truncated if and only if
   *    this is true.  If true then one extreme of the returned range will be
   *    moved so that it coincides with the closest extreme of the arg_range.
   *    The returned range will be only as small as necessary to achieve this.
   *    eg. if this is 1..100, the arg is 20..30 and truncate_if_contained is
   *    true then the returned range will be 20..100
   **/
  public MarkerRange combineRanges (final MarkerRange arg_range,
                                    final boolean truncate_if_contained) {
    if (getStrand () != arg_range.getStrand ()) {
      throw new Error ("internal error - strands do not match");
    }

    final MarkerRange return_range;

    try {
       return_range =
        new MarkerRange (getStrand (),
                         getStart ().getPosition (),
                         getEnd ().getPosition ());
    } catch (OutOfRangeException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    }

    if (return_range.contains (arg_range)) {
      if (truncate_if_contained) {
        // truncate the return_range
        if ((arg_range.getStart ().getPosition () -
             return_range.getStart ().getPosition ()) >
            (return_range.getEnd ().getPosition () -
             arg_range.getEnd ().getPosition ())) {
          return_range.end = arg_range.getEnd ();
        } else {
          return_range.start = arg_range.getStart ();
        }
      } else {
        // do nothing
      }
    } else {
      if (return_range.getStart ().getPosition () >
          arg_range.getStart ().getPosition ()) {
        return_range.start = arg_range.getStart ();
      }

      if (return_range.getEnd ().getPosition () <
          arg_range.getEnd ().getPosition ()) {
        return_range.end = arg_range.getEnd ();
      }
    }

    return return_range;
  }

  /**
   *  Return true if and only if the given MarkerRange is contained within
   *  this range.
   **/
  public boolean contains (final MarkerRange arg_range) {
    if (getStart ().getPosition () <= arg_range.getStart ().getPosition () &&
        getEnd ().getPosition () >= arg_range.getEnd ().getPosition ()) {
      return true;
    } else {
      return false;
    }
  }

  /**
   *  Create and return a Location object from the embl package that describes
   *  this MarkerRange object.  For example if the sequence is 1000 bases long
   *  and the MarkerRange object covers the reverse strand bases 100..200 then
   *  the Location object will be "complement(800..900)".
   **/
  public uk.ac.sanger.artemis.io.Location createLocation () {
    try {
      final uk.ac.sanger.artemis.io.Location new_location =
        new uk.ac.sanger.artemis.io.Location (getRawRange ());
      
      if (getStrand ().isForwardStrand ()) {
        return new_location;
      } else {
        return new_location.getComplement ();
      }
    } catch (OutOfRangeException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    }
  }

  /**
   *  The start Marker created from the Strand object and the start position
   *  passed to the constructor.
   **/
  private Marker start = null;

  /**
   *  The end Marker created from the Strand object and the end position
   *  passed to the constructor.
   **/
  private Marker end = null;
}


