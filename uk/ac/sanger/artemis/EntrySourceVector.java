/* EntrySourceVector.java
 *
 * created: Wed Jun  7 2000
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/EntrySourceVector.java,v 1.1 2004-06-09 09:44:26 tjc Exp $
 */

package uk.ac.sanger.artemis;

import java.util.Vector;

/**
 *  This class is a Vector of EntrySource objects.
 *
 *  @author Kim Rutherford
 *  @version $Id: EntrySourceVector.java,v 1.1 2004-06-09 09:44:26 tjc Exp $
 *
 **/

public class EntrySourceVector {
  /**
   *  Create a new (empty) EntrySourceVector object.
   **/
  public EntrySourceVector () {

  }

  /**
   *  Appends the given EntrySource object to the vector if and only if it
   *  isn't already in the vector.
   **/
  public void addElement (EntrySource entry_source) {
    if (indexOf (entry_source) == -1) {
      vector.addElement (entry_source);
    }
  }
  
  /**
   *  Appends the given EntrySource object to the vector if and only if it
   *  isn't already in the vector.  (same as addElement ()).
   **/
  public void add (EntrySource entry_source) {
    addElement (entry_source);
  }
  
  /**
   *  Performs the same function as Vector.elementAt ()
   **/
  public EntrySource elementAt (int index) {
    return (EntrySource) vector.elementAt (index);
  }

  /**
   *  Performs the same function as Vector.removeElement ()
   **/
  public boolean removeElement (EntrySource entry_source) {
    return vector.removeElement (entry_source);
  }

  /**
   *  Return true if this object contains the given EntrySource.
   **/
  public boolean contains (EntrySource entry_source) {
    if (indexOf (entry_source) == -1) {
      return false;
    } else {
      return true;
    }
  }

  /**
   *  Performs the same function as Vector.removeAllElements ()
   **/
  public void removeAllElements () {
    vector.removeAllElements ();
  }

  /**
   *   Performs the same function as Vector.removeElement ()
   **/
  public int indexOf (EntrySource entry_source) {
    return vector.indexOf (entry_source);
  }

  /**
   *  Performs the same function as Vector.size ()
   **/
  public int size () {
    return vector.size ();
  }
  
  /**
   *  Create a new EntrySourceVector with the same contents as this one.
   **/
  public Object clone () {
    final EntrySourceVector return_vector = new EntrySourceVector ();
    return_vector.vector = (Vector) vector.clone ();
    return return_vector;
  }

  /**
   *  Storage for EntrySource objects.
   **/
  private Vector vector = new Vector ();
}
