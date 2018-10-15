/* FastVector.java
 *
 * created: Sun Feb 20 2000
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2000  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/util/FastVector.java,v 1.4 2006-08-09 16:35:31 tjc Exp $
 */

package uk.ac.sanger.artemis.util;

import java.util.ArrayList;
import java.util.Comparator;
import java.util.Collections;

/**
 *  This class implements a Vector of Objects with a fast version of
 *  contains().
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: FastVector.java,v 1.4 2006-08-09 16:35:31 tjc Exp $
 *
 **/

public class FastVector extends ArrayList
{
  /**
   *  Performs the same function as Vector.addElement()
   */
  public boolean add(Object object)
  {
    if(object == null)
      throw new Error("internal error - adding a null object");
    else if(contains(object))
      throw new Error("internal error - object added a second time");

    return super.add(object);
  }

  /**
   *  Performs the same function as Vector.lastElement()
   **/
  public Object lastElement()
  {
    return (Object)get(size() - 1);
  }

  /**
   *  Insert an Object after another.
   *  @param old_object The new_object will be inserted after this object
   *    or at the start if old_object isn't in the vector.
   *  @param new_object The new object to insert.
   **/
  public void insertElementAfter(Object old_object, Object new_object)
  {
    final int old_object_index = indexOf(old_object);

    if(old_object_index == -1)
      add(0, new_object);
    else
      add(old_object_index + 1, new_object);
  }

  /**
   *  Replace the Object at the given index. (Performs the same function as
   *  Vector.elementAt())
   **/
  public void setElementAt(final Object object, final int index)
  {
    remove(index);
    add(index, object);
  }

  /**
   *  Return a sorted copy of this vector.
   *  @param cmp The returned vector will be sorted with this Comparator.
   **/
  public FastVector mysort(final Comparator cmp)
  {
    final FastVector return_vector = (FastVector)clone();
    Collections.sort(return_vector, cmp);
    return return_vector;
  }

}
