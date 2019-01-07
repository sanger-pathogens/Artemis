/* FeaturePredicateVector.java
 *
 * created: Mon Oct 13 2003
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 2003  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/FeaturePredicateVector.java,v 1.1 2004-06-09 09:44:46 tjc Exp $
 */

package uk.ac.sanger.artemis;

import java.util.Vector;

/**
 *  A Vector of FeaturePredicate objects.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: FeaturePredicateVector.java,v 1.1 2004-06-09 09:44:46 tjc Exp $
 **/

public class FeaturePredicateVector {
  /**
   *  Appends the given FeaturePredicate object to the vector.
   **/
  public void add (FeaturePredicate item) {
    vector.addElement (item);
  }
  
  /**
   *  Performs the same function as Vector.elementAt ()
   **/
  public FeaturePredicate elementAt (int index) {
    return (FeaturePredicate) vector.elementAt (index);
  }

  /**
   *  Return the size of this Vector.
   **/
  public int size () {
    return vector.size ();
  }

  /**
   *  Return a new copy of this object.
   **/
  public FeaturePredicateVector copy () {
    final FeaturePredicateVector new_vector = new FeaturePredicateVector ();

    new_vector.vector = (Vector) vector.clone ();

    return new_vector;
  }

  /**
   *  Delegate.
   **/
  private Vector vector = new Vector ();
}
