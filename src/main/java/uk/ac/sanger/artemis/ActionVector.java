/* ActionVector.java
 *
 * created: Tue Sep 17 2002
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 2002  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/ActionVector.java,v 1.1 2004-06-09 09:44:05 tjc Exp $
 */

package uk.ac.sanger.artemis;

import java.util.Vector;

/**
 *  A Vector of Action objects.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: ActionVector.java,v 1.1 2004-06-09 09:44:05 tjc Exp $
 **/

public class ActionVector {
  /**
   *  Create a new, empty ActionVector.
   **/
  public ActionVector () {

  }

  /**
   *  Appends the given Action object to the vector.
   **/
  public void add (Action item) {
    if (vector.size () > 0 &&
        vector.size () >= Options.getOptions ().getUndoLevels ()) {
      vector.removeElementAt (0); 
    }
    vector.addElement (item);
  }
  
  /**
   *  Performs the same function as Vector.elementAt ()
   **/
  public Action elementAt (int index) {
    return (Action) vector.elementAt (index);
  }

  /**
   *  Remove and return the last Action in the Vector.
   **/
  public Action removeAndReturnLast ()
      throws ArrayIndexOutOfBoundsException
  {
    final Action return_action = (Action) vector.lastElement ();

    vector.removeElementAt (vector.size() - 1); 

    return return_action;
  }

  /**
   *  Return the size of the Vector.
   **/
  public int size () {
    return vector.size ();
  }

  /**
   *  Delegate.
   **/
  final Vector vector = new Vector ();
}
