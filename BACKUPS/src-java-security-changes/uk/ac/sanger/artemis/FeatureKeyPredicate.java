/* FeatureKeyPredicate.java
 *
 * created: Tue Mar 30 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/FeatureKeyPredicate.java,v 1.1 2004-06-09 09:44:41 tjc Exp $
 */

package uk.ac.sanger.artemis;

import uk.ac.sanger.artemis.io.Key;

/**
 *  Each object of this class can be used to test Feature objects to see if
 *  they have the given Key.  See FeaturePredicate.
 *
 *  @author Kim Rutherford
 *  @version $Id: FeatureKeyPredicate.java,v 1.1 2004-06-09 09:44:41 tjc Exp $
 **/

public class FeatureKeyPredicate
    implements FeaturePredicate {
  /**
   *  Create a new FeatureKeyPredicate object.
   *  @param key The Key to test the Feature against.
   **/
  public FeatureKeyPredicate (final Key key) {
    this.key = key;
  }

  /**
   *  Test the given Feature against this FeatureKeyPredicate.
   *  @param feature The Feature to test the predicate against.
   *  @return Return true if and only if the given Feature has the same key
   *    as the one the was passed to the constructor.
   **/
  public boolean testPredicate (final Feature feature) {
    return feature.getKey ().equals (key);
  }
    
  /**
   *  The Key that was passed to the constructor.
   **/
  private Key key;
}
