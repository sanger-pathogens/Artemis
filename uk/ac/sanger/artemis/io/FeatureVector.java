/* FeatureVector.java
 *
 * created: Tue Oct 13 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/FeatureVector.java,v 1.1 2004-06-09 09:49:27 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import java.util.Vector;

/**
 *  This class implements a Vector of Feature objects.
 *
 *  @author Kim Rutherford
 *  @version $Id: FeatureVector.java,v 1.1 2004-06-09 09:49:27 tjc Exp $
 *
 */

public class FeatureVector {
  /**
   *  Create a new vector of Feature objects with an initial capacity of 100.
   **/
  public FeatureVector () {
    vector = new Vector (100);
  }

  /**
   *  Create a new vector of Feature objects with the given initial capacity.
   **/
  public FeatureVector (final int inital_capacity) {
    vector = new Vector (inital_capacity);
  }

  /**
   *  Performs the same function as Vector.addElement ()
   */
  public void addElement (Feature node) {
    vector.addElement (node);
  }
  
  /**
   *  Performs the same function as addElement ().
   */
  public void add (Feature node) {
    addElement (node);
  }
  
  /**
   *  Add a feature to the end of the Vector.
   **/
  public final void addElementAtEnd (Feature feature) {
    vector.insertElementAt (feature, size ());
  }

  /**
   *  Performs the same function as Vector.elementAt ()
   **/
  public Feature elementAt (int index) {
    return (Feature) vector.elementAt (index);
  }

  /**
   *  Performs the same function as Vector.lastElement ()
   **/
  public Feature lastElement () {
    return (Feature) vector.lastElement ();
  }

  /**
   *  Performs the same function as Vector.removeElement ()
   **/
  public boolean removeElement (Feature feature) {
    return vector.removeElement (feature);
  }

  /**
   *  Performs the same function as Vector.removeElement ()
   **/
  public int indexOf (Feature feature) {
    return vector.indexOf (feature);
  }

  /**
   *  Performs the same function as Vector.size ()
   **/
  public int size () {
    return vector.size ();
  }

  /**
   *  Performs the same function as Vector.removeAllElement ()
   **/
  public void removeAllElements () {
    vector.removeAllElements ();
  }

  /**
   *  Performs the same function as Vector.insertElementAt ()
   **/
  public final void insertElementAt (Feature feature, int index) {
    vector.insertElementAt (feature, index);
  }

  /**
   *  Insert a Feature after another.
   *  @param old_feature The new_feature will be inserted after this feature
   *    or at the start if old_feature isn't in the vector.
   *  @param new_feature The new feature to insert.
   **/
  public void insertElementAfter (Feature old_feature, Feature new_feature) {
    final int old_feature_index = indexOf (old_feature);

    if (old_feature_index == -1) {
      insertElementAt (new_feature, 0);
    } else {
      insertElementAt (new_feature, old_feature_index + 1);
    }
  }

  /**
   *  Create a new FeatureVector with the same contents as this one.
   **/
  public Object clone () {
    final FeatureVector return_vector = new FeatureVector ();
    return_vector.vector = (Vector) vector.clone ();
    return return_vector;
  }

  /**
   *  Storage for Feature objects.
   **/
  private Vector vector = null;
}


