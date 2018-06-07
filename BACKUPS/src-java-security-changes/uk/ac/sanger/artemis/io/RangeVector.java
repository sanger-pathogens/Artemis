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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/RangeVector.java,v 1.3 2007-10-25 19:25:00 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import java.util.Vector;

/**
 *  This class implements a Vector of Range objects.
 *
 *  @author Kim Rutherford
 *  @version $Id: RangeVector.java,v 1.3 2007-10-25 19:25:00 tjc Exp $
 *
 **/

public class RangeVector extends Vector<Range>
{
  private static final long serialVersionUID = 1L;

  /**
   *  Create a new vector of Range objects.
   **/
  public RangeVector() 
  {
  }

  /**
   *  Create a new vector of Range objects containing just the given Range.
   **/
  public RangeVector(final Range range) 
  {
    add (range);
  }


  /**
  *  Reverse this RangeVector in place.
  **/
  public void reverse()
  {
    for(int i = 0 ; i < size () / 2 ; ++i) 
    {
      final int swap_position = size () - i - 1;
      final Range tmp = elementAt(i);
      setElementAt(elementAt(swap_position), i);
      setElementAt(tmp, swap_position);
    }
  }
  
  public boolean containsRange(final Range r)
  {
    for(int i=0; i<size(); i++)
    {
      Range thisRange = elementAt(i);
      if(r.equals(thisRange))
        return true;
    }
    return false;
  }  
}

