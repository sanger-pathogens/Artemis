/* KeyVector.java
 *
 * created: Fri Apr 16 1999
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 1999  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/KeyVector.java,v 1.1 2004-06-09 09:49:45 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.FastVector;

/**
 *  This class implements a Vector of Key objects.
 *
 *  @author Kim Rutherford
 *  @version $Id: KeyVector.java,v 1.1 2004-06-09 09:49:45 tjc Exp $
 **/

public class KeyVector 
{
  /** Storage for Key objects. */
  private FastVector vector; 

  /**
   *  Create a new vector of Key objects.
   **/
  public KeyVector () 
  {
    vector = new FastVector ();
  }

  /**
   *  Create a new vector which contains only the given Key.
   **/
  public KeyVector (final Key new_key) 
  {
    vector = new FastVector ();
    add (new_key);
  }

  /**
   *  Performs the same function as Vector.addElement ()
   **/
  public void add (final Key node) 
  {
    vector.add (node);
  }

  /**
   *  Performs the same function as Vector.elementAt ()
   **/
  public Key elementAt (final int index) 
  {
    return (Key) vector.elementAt (index);
  }

  /**
   *  Performs the same function as Vector.setElementAt ()
   **/
  public void setElementAt (final Key key, final int index) 
  {
    vector.setElementAt (key, index);
  }

  /**
   *  Performs the same function as Vector.size ()
   **/
  public int size () 
  {
    return vector.size ();
  }

  /**
   *  Searches for the first occurence of the given argument, testing for
   *  equality using the equals method.
   *  @return the index of the first occurrence of the argument in this
   *    vector; returns -1 if the object is not found.
   **/
  public int indexOf (final Key key) 
  {
    return vector.indexOf (key);
  }

  /**
   *  Return true if this object contains the given Key, testing for
   *  equality using the equals method.
   **/
  public boolean contains (final Key key) 
  {
    return vector.contains (key);
  }

  /**
   *  Return a new copy of this object.
   **/
  public KeyVector copy () 
  {
    final KeyVector new_key_vector = new KeyVector ();

    new_key_vector.vector = (FastVector) vector.clone ();

    return new_key_vector;
  }

  /**
   * Sorts the elements of the vector using a simple O(n^2) selection
   * sort.
   */
  public void sort () 
  {
    int smallest;

    for (int i = 0; i < size (); ++i) 
    {
      //find smallest remaining element
      smallest = i;
      for(int j = i + 1 ; j < size () ; ++j) 
      {
        if(elementAt(j).compareTo (elementAt(smallest)) < 0) 
          smallest = j;
      }
      //exchange smallest and i
      if (smallest != i) 
      {
        final Key tmp = elementAt (i);
        setElementAt (elementAt(smallest), i);
        setElementAt (tmp, smallest);
      }
    }
  }

}

