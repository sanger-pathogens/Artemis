/* LocationParseNodeVector.java
 *
 * created: Wed Oct  7 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/LocationParseNodeVector.java,v 1.1 2004-06-09 09:49:54 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import java.util.Vector;

/**
 *  LocationParseNodeVector class is a Vector of objects of type
 *  LocationParseNode.
 *
 *  @author Kim Rutherford
 *  @version $Id: LocationParseNodeVector.java,v 1.1 2004-06-09 09:49:54 tjc Exp $
 * */

public class LocationParseNodeVector {
  /**
   *  Create a new LocationParseNodeVector.
   */
  public LocationParseNodeVector () {

  }

  /**
   *  Performs the same function as Vector.addElement ()
   */
  public void addElement (LocationParseNode node) {
    vector.addElement (node);
  }

  /**
   *  Performs the same function as Vector.elementAt ()
   */
  public LocationParseNode elementAt (int index) {
    return (LocationParseNode) vector.elementAt (index);
  }

  /**
   *  Return the last element of this vector.  This will throw an exception if
   *  the vector has nop elements.
   **/
  public LocationParseNode lastElement () {
    return elementAt (size () - 1);
  }

  /**
   *  Performs the same function as Vector.setElementAt ()
   **/
  void setElementAt (LocationParseNode node, int index) {
    vector.setElementAt (node, index);
  }

  /**
   *  Add the given node to the end of the vector.
   **/
  void addElementAtEnd (LocationParseNode node) {
    vector.insertElementAt (node, vector.size ());
  }

  /**
   *  Performs the same function as Vector.insertElementAt ().
   **/
  void insertElementAt (LocationParseNode node, int index) {
    vector.insertElementAt (node, index);
  }

  /**
   *  Performs the same function as Vector.removeElement ()
   **/
  boolean removeElement (LocationParseNode node) {
    return vector.removeElement (node);
  }

  /**
   *  Performs the same function as Vector.removeElementAt ().
   **/
  void removeElementAt (int index) {
    vector.removeElementAt (index);
  }

  /**
   *  Performs the same function as Vector.removeElement ()
   **/
  public int indexOf (Feature feature) {
    return vector.indexOf (feature);
  }

  /**
   *  Performs the same function as Vector.size ()
   */
  public int size () {
    return vector.size ();
  }


  /**
   *  Storage for LocationParseNode objects.
   */
  final private Vector vector = new Vector ();
}


