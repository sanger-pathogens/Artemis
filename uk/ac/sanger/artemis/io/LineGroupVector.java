/* LineGroupVector.java
 *
 * created: Tue Oct 13 1998
 *
 * This file is part of Artemis
 * 
 * Copyright(C) 1998,1999,2000  Genome Research Limited
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or(at your option) any later version.
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/LineGroupVector.java,v 1.1 2004-06-09 09:49:47 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import java.util.Vector;

/**
 *  This class implements a Vector of LineGroup objects.
 *
 *  @author Kim Rutherford
 *  @version $Id: LineGroupVector.java,v 1.1 2004-06-09 09:49:47 tjc Exp $
 *
 */

public class LineGroupVector 
{

  /** Storage for LineGroup objects. */
  final private Vector vector = new Vector();

  /**
   *  Create a new vector of LineGroup objects.
   */
  public LineGroupVector() 
  {
  }

  /**
   *  Performs the same function as Vector.addElement()
   */
  public void addElement(LineGroup node) 
  {
    vector.addElement(node);
  }
  
  /**
   *  Performs the same function as Vector.insertElementAt()
   */
  public void insertElementAt(LineGroup node, int index) 
  {
    vector.insertElementAt(node, index);
  }
  
  /**
   *  Performs the same function as Vector.removeElement()
   **/
  public void removeElementAt(final int index) 
  {
    vector.removeElementAt(index);
  }
  
  /**
   *  Performs the same function as Vector.elementAt()
   */
  public LineGroup elementAt(int index) 
  {
    return (LineGroup)vector.elementAt(index);
  }

  /**
   *  Returns the last component of the vector.
   *  @return The LineGroup at index size() - 1. 
   **/
  public final LineGroup lastElement() 
  {
    return (LineGroup)vector.lastElement();
  }

  /**
   *  Performs the same function as Vector.size()
   */
  public int size() 
  {
    return vector.size();
  }
  
}


