/* FeatureComparator.java
 *
 * created: Tue Sep 14 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/FeatureComparator.java,v 1.2 2008-09-08 10:36:30 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import java.util.Comparator;

/**
 *  Objects of this class implement the Comparator interface for
 *  ComparableFeature objects.  The features are ordered by first base, then
 *  last base then using ComparableFeature.getNumericID ().  This is an
 *  example ordering: 1..100, 1..200, 50..100, 150..250.
 *
 *  @author Kim Rutherford
 *  @version $Id: FeatureComparator.java,v 1.2 2008-09-08 10:36:30 tjc Exp $
 **/

public class FeatureComparator implements Comparator {
  /**
   *  Compare two ComparableFeature Objects with respect to ordering.
   *
   *  @param fst first argument
   *  @param snd second argument
   *  @return a negative number if first is less than second; a
   *    positive number if first is greater than second; else 0
   **/
  public int compare (final Object first, final Object second) {
    if (first instanceof Integer && second instanceof Integer) {
      final int first_integer = ((Integer)first).intValue ();
      final int second_integer = ((Integer)second).intValue ();

      if (first_integer < second_integer) {
        return -1;
      } else {
        if (first_integer > second_integer) {
          return 1;
        } else {
          return 0;
        }
      }
    }

    // this nonsense exists so that FeatureTree.findByBase() will work
    if (first instanceof Integer) {
      final int first_integer = ((Integer)first).intValue ();

      final ComparableFeature second_feature = (ComparableFeature) second;
      final int second_feature_first_base = second_feature.getFirstBase (); 

      if (first_integer < second_feature_first_base) {
        return -1;
      } else {
        if (first_integer > second_feature_first_base) {
          return 1;
        } else {
          return 0;
        }
      }
    }

    if (second instanceof Integer) {
      final int second_integer = ((Integer)second).intValue ();

      final ComparableFeature first_feature = (ComparableFeature) first;
      final int first_feature_first_base = first_feature.getFirstBase ();

      if (first_feature_first_base < second_integer) {
        return -1;
      } else {
        if (first_feature_first_base > second_integer) {
          return 1;
        } else {
          return 0;
        }
      }
    }

    final ComparableFeature first_feature = (ComparableFeature) first;
    final ComparableFeature second_feature = (ComparableFeature) second;

    final int first_feature_first_base = first_feature.getFirstBase ();
    final int second_feature_first_base = second_feature.getFirstBase ();

    final int first_feature_last_base = first_feature.getLastBase ();
    final int second_feature_last_base = second_feature.getLastBase ();

    if (first_feature_first_base < second_feature_first_base) {
      return -1;
    }

    if (first_feature_first_base > second_feature_first_base) {
      return 1;
    }

    final int index_of_first =
      indexOfKey (first_feature.getKey ());
    final int index_of_second =
      indexOfKey (second_feature.getKey ());

    if (index_of_first != index_of_second) {
      if (index_of_first < index_of_second) {
        return -1;
      } else {
        if (index_of_first > index_of_second) {
          return 1;
        }
      }
    } else {
      // fall through
    }

    if (first_feature_last_base < second_feature_last_base) {
      return -1;
    }

    if (first_feature_last_base > second_feature_last_base) {
      return 1;
    }

    // we use the id counter as a tie breaker so that compare () will only
    // return 0 if the two ComparableFeature objects are the same
    if (first_feature.getNumericID () < second_feature.getNumericID ()) {
      return -1;
    } else {
      if (first_feature.getNumericID () > second_feature.getNumericID ()) {
        return 1;
      } else {
        return 0;
      }
    }
  }

  /**
   *  Return the index of the given Key in the keys_to_sort_first array.
   **/
  private int indexOfKey (final Key key) {
    for (int i = 0 ; i < keys_to_sort_first.length ; ++i) {
      if (keys_to_sort_first[i].equals (key.toString ())) {
        return i;
      }
    }

    return keys_to_sort_first.length;
  }

  final private String keys_to_sort_first[] = {
    "source",
    "gene",
    "CDS"
  };
}
