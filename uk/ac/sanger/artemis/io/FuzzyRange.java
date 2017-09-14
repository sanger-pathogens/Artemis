/* FuzzyRange.java
 *
 * created: Wed Apr 28 1999
 *
 * This file is part of Artemis
 *
 * Copyright (C) 1999,2000,2001  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/FuzzyRange.java,v 1.1 2004-06-09 09:49:28 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.*;

/**
 *  This is a Range with a non-exact upper or lower bound.  It represents EMBL
 *  locations like (100.200)..300 or (10.20)..(40.50)
 *  An object of this class can be constructed with a Range object for the
 *  start or end or both.
 *
 *  @author Kim Rutherford
 *  @version $Id: FuzzyRange.java,v 1.1 2004-06-09 09:49:28 tjc Exp $
 **/

public class FuzzyRange extends Range {
  /**
   *  Create a new FuzzyRange object.
   *  @param start_object the start of the FuzzyRange
   *  @param end_object the end of the FuzzyRange
   **/
  private FuzzyRange (final Object start_object, final Object end_object)
      throws OutOfRangeException {
    super (getStartPosition (start_object),
           getEndPosition (end_object));
    this.start_object = start_object;
    this.end_object = end_object;
  }

  /**
   *  Create a new FuzzyRange object.
   *  @param position_object the position of the FuzzyRange
   **/
  private FuzzyRange (final Object position_object) {
    super (getPosition (position_object));
    this.start_object = position_object;
    this.end_object = null;
  }

  /**
   *  Make and return a new FuzzyRange or Range object.
   **/
  public static Range makeRange (final Object start_object,
                                 final Object end_object)
      throws OutOfRangeException {
    if (end_object == null) {
      return makeRange (start_object);
    }

    if (start_object instanceof Integer &&
        end_object instanceof Integer) {
      int start_pos = ((Integer) start_object).intValue ();
      int end_pos = ((Integer) end_object).intValue ();

      if (start_pos <= end_pos) {
        return new Range (start_pos, end_pos);
      } else {
        return new Range (end_pos, start_pos);
      }
    } else {
      return new FuzzyRange (start_object, end_object);
    }
  }

  /**
   *  Make and return a new FuzzyRange or Range object.
   **/
  public static Range makeRange (final Object start_object) {
    if (start_object instanceof Integer) {
      return new Range (((Integer) start_object).intValue ());
    } else {
      return new FuzzyRange (start_object);
    }
  }

  /**
   *  Return the integer position of the given object, which should be the
   *  start object of a range.
   **/
  private static int getStartPosition (final Object start_object) {
    if (start_object instanceof Integer) {
      return ((Integer) start_object).intValue ();
    } else {
      if (start_object instanceof LowerInteger) {
        return ((LowerInteger) start_object).getPosition ();
      } else {
        if (start_object instanceof Range) {
          return ((Range) start_object).getStart ();
        } else {
          throw new Error ("internal error: object is not recognised: " +
                           start_object);
        }
      }
    }
  }

  /**
   *  Return the integer position of the given object, which should be the
   *  end object of a range.
   **/
  private static int getEndPosition (final Object end_object) {
    if (end_object instanceof Integer) {
      return ((Integer) end_object).intValue ();
    } else {
      if (end_object instanceof UpperInteger) {
        return ((UpperInteger) end_object).getPosition ();
      } else {
        if (end_object instanceof Range) {
          return ((Range) end_object).getEnd ();
        } else {
          throw new Error ("internal error: object is not recognised: " +
                           end_object);
        }
      }
    }
  }

  /**
   *  Return the integer position of the given object.
   **/
  private static int getPosition (final Object position_object) {
    if (position_object instanceof Integer) {
      return ((Integer) position_object).intValue ();
    } else {
      if (position_object instanceof LowerInteger) {
        return ((LowerInteger) position_object).getPosition ();
      } else {
        if (position_object instanceof UpperInteger) {
          return ((UpperInteger) position_object).getPosition ();
        } else {
          if (position_object instanceof Range) {
            return ((Range) position_object).getStart ();
          } else {
            throw new Error ("internal error: object is not recognised: " +
                             position_object);
          }
        }
      }
    }
  }

  /**
   *  Return true if and only if the argument equals this Range.
   **/
  public boolean equals (final Range test_range) {
    if (test_range instanceof FuzzyRange) {

      final FuzzyRange fuzzy_range = (FuzzyRange) test_range;

      if (getStart () != fuzzy_range.getStart ()) {
        return false;
      }

      if (getEnd () != fuzzy_range.getEnd ()) {
        return false;
      }

      if (start_object instanceof Integer &&
          fuzzy_range.start_object instanceof Integer) {

        if (((Integer) start_object).intValue () !=
            ((Integer) fuzzy_range.start_object).intValue ()) {
          return false;
        }
      }

      if (start_object instanceof LowerInteger &&
          fuzzy_range.start_object instanceof LowerInteger) {

        if (((LowerInteger) start_object).getPosition () !=
            ((LowerInteger) fuzzy_range.start_object).getPosition ()) {
          return false;
        }
      }

      if (start_object instanceof Range &&
          fuzzy_range.start_object instanceof Range) {

        if (((Range) start_object).getStart () !=
            ((Range) fuzzy_range.start_object).getStart ()) {
          return false;
        }
      }

      if (end_object instanceof Integer &&
          fuzzy_range.end_object instanceof Integer) {

        if (((Integer) end_object).intValue () !=
            ((Integer) fuzzy_range.end_object).intValue ()) {
          return false;
        }
      }

      if (end_object instanceof UpperInteger &&
          fuzzy_range.end_object instanceof UpperInteger) {

        if (((UpperInteger) end_object).getPosition () !=
            ((UpperInteger) fuzzy_range.end_object).getPosition ()) {
          return false;
        }
      }

      if (end_object instanceof Range &&
          fuzzy_range.end_object instanceof Range) {

        if (((Range) end_object).getEnd () !=
            ((Range) fuzzy_range.end_object).getEnd ()) {
          return false;
        }
      }

      return true;
    }

    return false;
  }

  /**
   *  Return a String representation of this FuzzyRange object.
   **/
  public String toString () {
    final String start_string;

    if (start_object instanceof Range) {
      start_string ="(" +
        ((Range) start_object).getStart () + "." +
        ((Range) start_object).getEnd () +
         ")";
    } else {
      start_string = start_object.toString ();
    }

    if (end_object == null) {
      return start_string;
    } else {
      final String end_string;

      if (end_object instanceof Range) {
        end_string ="(" +
          ((Range) end_object).getStart () + "." +
          ((Range) end_object).getEnd () +
          ")";
      } else {
        end_string = end_object.toString ();
      }
      
      return start_string + ".." + end_string;
    }
  }

  /**
   *  Return a copy of this object with the start and end changed.  Currently
   *  if the start or end object is a Range then that object will be turned
   *  into an Integer.
   **/
  public Range change (final int start, final int end)
      throws OutOfRangeException {
    if ((start_object instanceof Integer ||
         start_object instanceof LowerInteger) &&
        (end_object instanceof Integer ||
         end_object instanceof UpperInteger)) {
      final Object new_start_object;
      
      if (start_object instanceof Integer) {
        new_start_object = new Integer (start);
      } else {
        new_start_object = new LowerInteger (new Integer (start));
      }

      final Object new_end_object;

      if (end_object instanceof Integer) {
        new_end_object = new Integer (end);
      } else {
        new_end_object = new UpperInteger (new Integer (end));
      }

      return new FuzzyRange (new_start_object, new_end_object);
    } else {
      return new Range (start, end);
    }
  }

  /**
   *  Return a copy of this object.
   **/
  public Range copy () {
    try {
      return makeRange (start_object, end_object);
    } catch (OutOfRangeException e) {
      throw new Error ("internal error - unexpected exception");
    }
  }
  
  public Object getStartObject()
  {
    return start_object;
  }
  
  public Object getEndObject()
  {
    return end_object;
  }

  /**
   *  The start Object of this FuzzyRange.
   **/
  private Object start_object = null;

  /**
   *  The end Range of this FuzzyRange.
   **/
  private Object end_object = null;
}

