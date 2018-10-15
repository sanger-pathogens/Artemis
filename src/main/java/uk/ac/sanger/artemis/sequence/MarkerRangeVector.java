/* MarkerRangeVector.java
 *
 * created: Sun Jan 24 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/sequence/MarkerRangeVector.java,v 1.1 2004-06-09 09:52:21 tjc Exp $
 */

package uk.ac.sanger.artemis.sequence;

import java.util.Vector;

/**
 *  This class implements a Vector of MarkerRange objects.
 *
 *  @author Kim Rutherford
 *  @version $Id: MarkerRangeVector.java,v 1.1 2004-06-09 09:52:21 tjc Exp $
 **/

public class MarkerRangeVector {
  /**
   *  Create a new (empty) MarkerRangeVector.
   **/
  public MarkerRangeVector () {

  }

  /**
   *  Performs the same function as Vector.addElement ()
   */
  public void add (final MarkerRange range) {
    vector.addElement (range);
  }

    /**
   *  Performs the same function as Vector.elementAt ()
   */
  public MarkerRange elementAt (final int index) {
    return (MarkerRange) vector.elementAt (index);
  }

  /**
   *  Return true if this object contains the given MarkerRange.
   **/
  public boolean contains (final MarkerRange range) {
    if (indexOf (range) == -1) {
      return false;
    } else {
      return true;
    }
  }

  /**
   *  Performs the same function as Vector.removeElement ()
   **/
  public int indexOf (final MarkerRange range) {
    return vector.indexOf (range);
  }

    /**
   *  Performs the same function as Vector.size ()
   */
  public int size () {
    return vector.size ();
  }

  /**
   *  Storage for MarkerRange objects.
   */
  private Vector vector = new Vector ();
}


