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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/KeyVector.java,v 1.2 2006-08-09 16:35:31 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.FastVector;

/**
 *  This class implements a Vector of Key objects.
 *
 *  @author Kim Rutherford
 *  @version $Id: KeyVector.java,v 1.2 2006-08-09 16:35:31 tjc Exp $
 **/

public class KeyVector extends FastVector
{

  public KeyVector ()
  {
    super();
  }

  /**
   *  Create a new vector which contains only the given Key.
   **/
  public KeyVector (final Key new_key)
  {
    super();
    add (new_key);
  }

  /**
   *  Return a new copy of this object.
   **/
  public KeyVector copy ()
  {
    final KeyVector new_key_vector = (KeyVector)clone();

    return new_key_vector;
  }

  /**
   * Sorts the elements of the vector using a simple O(n^2) selection
   * sort.
   */
  public void mysort()
  {
    int smallest;

    for (int i = 0; i < size (); ++i)
    {
      //find smallest remaining element
      smallest = i;
      for(int j = i + 1 ; j < size () ; ++j)
      {
        if(((Key)get(j)).compareTo( (Key)get(smallest)) < 0)
          smallest = j;
      }
      //exchange smallest and i
      if (smallest != i)
      {
        final Key tmp = (Key)get(i);
        setElementAt (get(smallest), i);
        setElementAt (tmp, smallest);
      }
    }
  }

}

