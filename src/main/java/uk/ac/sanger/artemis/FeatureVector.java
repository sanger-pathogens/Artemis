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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/FeatureVector.java,v 1.2 2006-08-09 16:35:31 tjc Exp $
 */

package uk.ac.sanger.artemis;

import uk.ac.sanger.artemis.io.IndexedGFFDocumentEntry;
import uk.ac.sanger.artemis.util.FastVector;

import java.util.*;

/**
 *  This class implements a Vector of Feature objects.
 *
 *  @author Kim Rutherford
 *  @version $Id: FeatureVector.java,v 1.2 2006-08-09 16:35:31 tjc Exp $
 *
 **/

public class FeatureVector {
  /**
   *  Create a new vector of Feature objects.
   **/
  public FeatureVector () {

  }

  /**
   *  Performs the same function as Vector.addElement ()
   */
  public void add (Feature feature) {
    vector.add (feature);
  }

  /**
   *  Add a feature to the end of the Vector.
   **/
  public final void addElementAtEnd (Feature feature) {
    vector.add (feature);
  }

  /**
   *  Performs the same function as Vector.elementAt ()
   */
  public Feature elementAt (int index) {
    return (Feature) vector.get(index);
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
  public boolean remove (Feature feature) {
    return vector.remove (feature);
  }

  /**
   *  Return true if this object contains the given Feature.
   **/
  public boolean contains (Feature feature) {
    if(feature.getEntry().getEMBLEntry() instanceof IndexedGFFDocumentEntry)
      return IndexedGFFDocumentEntry.contains(feature, this);

    return vector.contains (feature);
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
   *  Performs the same function as Vector.removeAllElement ()
   **/
  public void removeAllElements () {
    vector.clear ();
  }

  /**
   *  Performs the same function as Vector.removeElementAt ()
   **/
  public void removeElementAt(int index) {
    vector.remove (index);
  }

  /**
   *  Performs the same function as Vector.insertElementAt ()
   **/
  public final void insertElementAt (Feature feature, int index) {
    vector.add (index, feature);
  }

  /**
   *  Insert a Feature after another.
   *  @param old_feature The new_feature will be inserted after this feature
   *    or at the start if old_feature isn't in the vector.
   *  @param new_feature The new feature to insert.
   **/
  public void insertElementAfter (Feature old_feature, Feature new_feature) {
    vector.insertElementAfter (old_feature, new_feature);
  }

  /**
   *  Create a new FeatureVector with the same contents as this one.
   **/
  public Object clone () {
    final FeatureVector return_vector = new FeatureVector ();

    return_vector.vector = (FastVector) vector.clone ();

    return return_vector;
  }

  /**
   *  Return a sorted copy of this vector.
   *  @param cmp The returned vector will be sorted with this Comparator.
   **/
  public FeatureVector sort (final Comparator cmp) {
    final FeatureVector return_vector = (FeatureVector) clone ();

    return_vector.vector = return_vector.vector.mysort (cmp);

    return return_vector;
  }

  /**
   *  Storage for Feature objects.
   **/
  private FastVector vector = new FastVector ();
}
