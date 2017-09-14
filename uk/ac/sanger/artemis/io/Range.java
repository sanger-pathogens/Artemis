/* Range.java
 *
 * created: Tue Oct  6 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/Range.java,v 1.1 2004-06-09 09:50:18 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.*;

import java.lang.Error;

/**
 *  This class represents an inclusive range of integers.
 *
 *  @author Kim Rutherford
 *  @version $Id: Range.java,v 1.1 2004-06-09 09:50:18 tjc Exp $
 **/

public class Range {
  /**
   *  Create a new Range object.  The start must be less than or equal to the
   *  end.
   *  @param start the start of the range
   *  @param end the end of the range
   *  @throws OutOfRangeException if the end value is less than the start
   *    value.
   **/
  public Range (int start, int end)
      throws OutOfRangeException {
    if (end < start) {
      throw new OutOfRangeException ("start: " + start + " > " + "end: " +
                                     end);
    }

    this.start = start;
    this.end = end;
  }

  /**
   *  Create a new Range object representing a single integer.
   *  @param number The single integer in the range.
   **/
  public Range (int number) {
    this.start = number;
    this.end = number;
  }

  /**
   *  Return the number of integers in the range, inclusive of the end points.
   **/
  public int getCount () {
    return end - start + 1;
  }

  /**
   *  Return a string representation of this object in the form "start..end",
   *  where start and end are the integers that were passed to the
   *  constructor.
   **/
  public String toString () {
    if (end > start) {
      return start + ".." + end;
    } else {
      return Integer.toString (start);
    }
  }

  /**
   *  Return the start position of this range, as passed to the constructor.
   **/
  public int getStart () {
    return start;
  }

  /**
   *  Return the start position of this range, as passed to the constructor.
   **/
  public int getEnd () {
    return end;
  }

  /**
   *  Return a copy of this object.
   **/
  public Range copy () {
    try {
      return new Range (getStart (), getEnd ());
    } catch (OutOfRangeException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    }
  }

  /**
   *  Returns true if and only if the test_range is completely contained within
   *  this Range.
   **/
  public boolean contains (final Range test_range) {
    if (this.overlaps (test_range)) {
      if (getStart () <= test_range.getStart () &&
          getEnd () >= test_range.getEnd ()) {
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }
  }

  /**
   *  Return a copy of this object with the start and end changed.
   **/
  public Range change (final int start, final int end)
      throws OutOfRangeException {
    return new Range (start, end);
  }

  /**
   *  This method translates the start and end of this Range into another
   *  coordinate system.  The Range will be truncated if necessary.
   *  @param constraint This contains the start and end base of the new
   *    coordinate system.  The position given by constraint.getStart () will
   *    be at postion/base 1 in the new coordinate system.
   *  @return the Range translated into the new coordinate system.  The Range
   *    will be truncated if it overlaps either end of the constraint.  null
   *    is returned if and only if the new Range lies outside of the new
   *    start and end.
   **/
  public Range truncate (final Range constraint) {
    if (constraint.getEnd () < getStart ()) {
      return null;
    } else {
      if (constraint.getStart () > getEnd ()) {
        return null;
      } else {
        final int new_start = getStart () - constraint.getStart () + 1;
        final int new_end = getEnd () - constraint.getStart () + 1;

        final Object start_object;

        if (new_start < 1) {
          start_object = new LowerInteger (1);
        } else {
          start_object = new Integer (new_start);
        }

        final Object end_object;

        if (new_end > constraint.getCount ()) {
          end_object = new UpperInteger (constraint.getCount ());
        } else {
          end_object = new Integer (new_end);
        }

        try {
          return FuzzyRange.makeRange (start_object, end_object);
        } catch (OutOfRangeException e) {
          throw new Error ("internal error - unexpected exception: " + e);
        }
      }
    }
  }

  /**
   *  Return true if and only if the argument equals this Range.
   **/
  public boolean equals (final Range test_range) {
    if (test_range instanceof Range) {
      if (test_range.getStart () == getStart () &&
          test_range.getEnd () == getEnd ()) {
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }
  }

  /**
   *  Return true if and only if there is an overlap between this Range and the
   *  given Range.
   **/
  public boolean overlaps (final Range arg_range) {
    if (getStart () < arg_range.getStart () &&
        getEnd () < arg_range.getStart ()) {
      return false;
    }
    if (getStart () > arg_range.getEnd () &&
        getEnd () > arg_range.getEnd ()) {
      return false;
    }
    return true;
  }
  
  
  /**
   *  Return true if and only if there is an overlap between this Range and the
   *  given Range.
   **/
  public boolean fuzzyOverlaps (final Range r2, int nbases) {
     
    if (getStart () < r2.getStart () - nbases &&
        getEnd () < r2.getStart () - nbases) {
      return false;
    }
    if (getStart () > r2.getEnd () + nbases &&
        getEnd () > r2.getEnd () + nbases) {
      return false;
    }
    return true;
  }
  

  /**
   *  The is the start value that was passed to the constructor
   */
  private int start;

  /**
   *  The is the end value that was passed to the constructor
   */
  private int end;
}


