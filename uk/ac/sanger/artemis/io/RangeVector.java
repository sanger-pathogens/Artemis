/* RangeVector.java
 *
 * created: Thu Oct 29 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/RangeVector.java,v 1.1 2004-06-09 09:50:19 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import java.util.Vector;

/**
 *  This class implements a Vector of Range objects.
 *
 *  @author Kim Rutherford
 *  @version $Id: RangeVector.java,v 1.1 2004-06-09 09:50:19 tjc Exp $
 *
 **/

public class RangeVector {
  /**
   *  Create a new vector of Range objects.
   **/
  public RangeVector () {

  }

  /**
   *  Create a new vector of Range objects containing just the given Range.
   **/
  public RangeVector (final Range range) {
    add (range);
  }

  /**
   *  Performs the same function as Vector.addElement ()
   */
  public void addElement (Range node) {
    vector.addElement (node);
  }

  /**
   *  Performs the same function as Vector.insertElementAt ()
   **/
  public final void insertElementAt (Range range, int index) {
    vector.insertElementAt (range, index);
  }

  /**
   *  Performs the same function as Vector.addElement ()
   */
  public void add (Range node) {
    vector.addElement (node);
  }
  
  /**
   *  Performs the same function as Vector.elementAt ()
   */
  public Range elementAt (int index) {
    return (Range) vector.elementAt (index);
  }


  /**
   *  Performs the same function as Vector.removeElement ()
   **/
  public boolean removeElement (Range Range) {
    return vector.removeElement (Range);
  }
  

  /**
   *  Performs the same function as Vector.removeElement ()
   **/
  public int indexOf (Range Range) {
    return vector.indexOf (Range);
  }


  /**
   *  Performs the same function as Vector.size ()
   */
  public int size () {
    return vector.size ();
  }

  /**
   *  Reverse this RangeVector in place.
   **/
  public void reverse () {
    for (int i = 0 ; i < vector.size () / 2 ; ++i) {
      final int swap_position = vector.size () - i - 1;
      final Object tmp = vector.elementAt (i);
      vector.setElementAt (vector.elementAt (swap_position), i);
      vector.setElementAt (tmp, swap_position);
    }
  }
  
  /**
   *  Storage for Range objects.
   */
  final private Vector vector = new Vector ();
}


