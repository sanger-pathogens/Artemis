/* BetweenRange.java
 *
 * created: Thu Feb 11 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/BetweenRange.java,v 1.1 2004-06-09 09:48:50 tjc Exp $
 **/

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.*;

/**
 *  This is a Range that represents a position between two bases.
 *
 *  @author Kim Rutherford
 *  @version $Id: BetweenRange.java,v 1.1 2004-06-09 09:48:50 tjc Exp $
 **/

public class BetweenRange extends Range {
  /**
   *  Create a new BetweenRange object.
   *  @param start the start of the range
   *  @param end the end of the range
   **/
  public BetweenRange (int start, int end)
      throws OutOfRangeException {
    super (start, end);
  }

  /**
   *  Return a string representation of this object in the form "start^end",
   *  where start and end are the integers that were passed to the
   *  constructor.
   **/
  public String toString () {
    return getStart () + "^" + getEnd ();
  }

  /**
   *  Return true if and only if the argument equals this Range.
   **/
  public boolean equals (final Range test_range) {
    if (test_range instanceof BetweenRange) {
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
   *  Return a copy of this object.
   **/
  public Range copy () {
    try {
          return new BetweenRange (getStart (), getEnd ());
    } catch (OutOfRangeException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    }
  }
}


