/* OutOfRangeException.java
 *
 * created: Sun Jan 10 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/util/OutOfRangeException.java,v 1.1 2004-06-09 09:53:06 tjc Exp $
 */

package uk.ac.sanger.artemis.util;

/**
 *  This Exception is thrown when a value passed to a method is out of range
 *  for the given operation.  For example if a base position in a Location is
 *  less than 1 or greater than the length of the sequence.
 *
 *  @author Kim Rutherford
 *  @version $Id: OutOfRangeException.java,v 1.1 2004-06-09 09:53:06 tjc Exp $
 **/

public class OutOfRangeException extends Exception {
  /**
   *  Create a new OutOfRangeException object.
   **/
  public OutOfRangeException () {

  }

  /**
   *  Create a new OutOfRangeException object with the given message.
   **/
  public OutOfRangeException (final String message) {
    super (message);
  }
}
