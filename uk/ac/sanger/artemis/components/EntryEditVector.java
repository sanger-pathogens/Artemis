/* EntryEditVector.java
 *
 * created: Sat Oct 17 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/EntryEditVector.java,v 1.1 2004-06-09 09:46:27 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import java.util.Vector;

/**
 *  This class implements a Vector of EntryEdit objects.
 *
 *  @author Kim Rutherford
 *  @version $Id: EntryEditVector.java,v 1.1 2004-06-09 09:46:27 tjc Exp $
 *
 **/

public class EntryEditVector {

  /**
   *  Performs the same function as Vector.addElement ()
   */

  
  /**
   *  Performs the same function as Vector.elementAt ()
   */
  public EntryEdit elementAt (int index) {
    return vector.elementAt (index);
  }

  /**
   *  Performs the same function as Vector.removeElement ()
   **/
  public boolean removeElement (EntryEdit entry_edit) {
    return vector.removeElement (entry_edit);
  }

  /**
   *  Performs the same function as Vector.removeElement ()
   **/
  public int indexOf (EntryEdit entry_edit) {
    return vector.indexOf (entry_edit);
  }

  /**
   *  Performs the same function as Vector.size ()
   */
  public int size () {
    return vector.size ();
  }
  
  /**
   *  Storage for EntryEdit objects.
   */
  final private Vector<EntryEdit> vector = new Vector<EntryEdit> ();
}
