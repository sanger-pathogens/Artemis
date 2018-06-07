/* FeatureFromVectorPredicate.java
 *
 * created: Wed Sep  8 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/FeatureFromVectorPredicate.java,v 1.1 2004-06-09 09:44:40 tjc Exp $
 */

package uk.ac.sanger.artemis;

/**
 *  Each object of this class can be used to test Feature objects to see if
 *  they are in the FeatureVector that is passed to the constructor.  See
 *  FeaturePredicate.
 *
 *  @author Kim Rutherford
 *  @version $Id: FeatureFromVectorPredicate.java,v 1.1 2004-06-09 09:44:40 tjc Exp $
 **/

public class FeatureFromVectorPredicate
  implements FeaturePredicate {
  /**
   *  Create a new FeatureFromVectorPredicate that tests as true (with
   *  testPredicate ()) only for those features in the given FeatureVector.
   **/
  public FeatureFromVectorPredicate (final FeatureVector feature_vector) {
    this.feature_vector = feature_vector;
  }

  /**
   *  Test a Feature against this FeatureKeyPredicate.
   *  @param feature The Feature to test the predicate against.
   *  @return Return true if and only if the given Feature is in the
   *    FeatureVector that was passed to the constructor.
   **/
  public boolean testPredicate (final Feature feature) {
    return feature_vector.contains (feature);
  }

  /**
   *  The FeatureVector that was passed to the constructor.
   **/
  private FeatureVector feature_vector;
}


