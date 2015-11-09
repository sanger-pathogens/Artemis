/* AlignMatchVector.java
 *
 * created: Sat Jun 17 2000
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2000,2001  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/AlignMatchVector.java,v 1.3 2008-06-26 09:38:13 tjc Exp $
 */

package uk.ac.sanger.artemis;

import uk.ac.sanger.artemis.util.FastVector;

import java.util.Comparator;

/**
 *  This class is a Vector of AlignMatch objects.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: AlignMatchVector.java,v 1.3 2008-06-26 09:38:13 tjc Exp $
 **/

public class AlignMatchVector {
  /**
   *  Create a new (empty) AlignMatchVector object.
   **/
  public AlignMatchVector () {

  }

  /**
   *  Appends the given AlignMatch object to the vector if and only if it
   *  isn't already in the vector.
   **/
  public void addElement (AlignMatch item) {
    vector.add (item);
  }

  /**
   *  Appends the given AlignMatch object to the vector if and only if it
   *  isn't already in the vector.  (same as addElement ()).
   **/
  public void add (AlignMatch item) {
    addElement (item);
  }

  /**
   *  Performs the same function as Vector.elementAt ()
   **/
  public AlignMatch elementAt (int index) {
    return (AlignMatch) vector.get(index);
  }

  /**
   *  Performs the same function as Vector.removeElement ()
   **/
  public boolean remove (AlignMatch item) {
    return vector.remove (item);
  }

  /**
   *  Return true if this object contains the given AlignMatch.
   **/
  public boolean contains (AlignMatch item) {
    return vector.contains(item);
  }

  /**
   *  Performs the same function as Vector.removeAllElements ()
   **/
  public void removeAllElements () {
    vector.clear ();
  }

  /**
   *   Performs the same function as Vector.removeElement ()
   **/
  public int indexOf (AlignMatch item) {
    return vector.indexOf (item);
  }

  /**
   *  Performs the same function as Vector.size ()
   **/
  public int size () {
    return vector.size ();
  }

  /**
   *  Sort this vector.
   *  @param cmp The returned vector will be sorted with this Comparator.
   **/
  public void sort (final Comparator cmp) {
    vector = vector.mysort (cmp);
  }

  /**
   *  Create a new AlignMatchVector with the same contents as this one.
   **/
  public Object clone () {
    final AlignMatchVector return_vector = new AlignMatchVector ();
    return_vector.vector = (FastVector) vector.clone ();
    return return_vector;
  }

  /**
   *  Storage for AlignMatch objects.
   **/
  private FastVector vector = new FastVector ();
}
