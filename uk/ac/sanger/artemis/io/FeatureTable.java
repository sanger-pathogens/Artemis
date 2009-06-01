/* FeatureTable.java
 *
 * created: Mon Oct 12 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/FeatureTable.java,v 1.2 2009-06-01 09:45:32 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

/**
 *  This object contains all the features in an entry.
 *
 *  @author Kim Rutherford
 *  @version $Id: FeatureTable.java,v 1.2 2009-06-01 09:45:32 tjc Exp $
 */

abstract class FeatureTable extends LineGroup {
  /**
   *  Create a new empty FeatureTable object.
   **/
  FeatureTable () {

  }

  /**
   *  Return the feature tree from this feature table.
   *  @return The features of this feature table.  The returned
   *    object is not a copy - changes to it will change the FeatureTable
   *    object itself
   */
  private FeatureTree getFeatures () {
    return features;
  }

  /**
   *  Return a count of the number of Feature obejcts in this FeatureTable.
   **/
  int getFeatureCount () {
    return getFeatures ().size ();
  }

  /**
   *  Return a vector containing the references of the Feature objects whose
   *  location overlap the given range.
   *  @param range Return features that overlap this range - ie the start of
   *    the feature is less than or equal to the end of the range and the end
   *    of the feature is greater than or equal to the start of the range.
   *  @return The features of this feature table the are within
   *    the given range.  The returned object is a copy - changes will not
   *    effect the FeatureTable object itself.
   **/
  FeatureVector getFeaturesInRange (Range range) {
    return getFeatures ().getFeaturesInRange (range);
  }

  /**
   *  Return a vector containing the references of the all Feature objects in
   *  this FeatureTable.
   *  @return The features of this feature table.  The returned
   *    object is a copy - changes will not effect the FeatureTable object
   *    itself.
   **/
  FeatureVector getAllFeatures () {
    final FeatureVector return_features = new FeatureVector ();

    final FeatureEnumeration enumerator = features ();

    while (enumerator.hasMoreFeatures ()) {
      final Feature this_feature = enumerator.nextFeature ();

      return_features.add (this_feature);
    }

    return return_features;
  }

  /**
   *  Add a feature to the feature table.  The features are ordered by first
   *  base and then last base.  This is an example ordering: 1..100, 1..200,
   *  50..100, 150..250.  This method resets cached_feature_index and
   *  cached_feature.
   *  @param new_feature The feature to add
   **/
  void add (final Feature new_feature) {
    cached_feature = null;
    cached_feature_index = -999;

    features.add (new_feature); 
  }

  /**
   *  Remove the given Feature from this FeatureTable.  This method resets
   *  cached_feature_index and cached_feature.
   **/
  Feature remove (final Feature feature) {
    cached_feature = null;
    cached_feature_index = -999;

    if (getFeatures ().contains (feature)) {
      getFeatures ().remove (feature);
      return feature;
    } else {
      return null;
    }
  }

  /**
   *  Return the ith Feature from this FeatureTable.  This Features are
   *  returned in a consistent order, sorted by the first base of each
   *  Feature.  This method uses and updates cached_feature_index and
   *  cached_feature.
   **/
  Feature getFeatureAtIndex (final int arg_index) {
    if (cached_feature != null) {
      if (arg_index == cached_feature_index) {
        return cached_feature;
      }

      if (arg_index == cached_feature_index + 1) {
        final Feature next_feature =
          getFeatures ().getNextFeature (cached_feature);

        if (next_feature == null) {
          throw new ArrayIndexOutOfBoundsException (arg_index);
        } else {
          cached_feature_index = arg_index;
          cached_feature = next_feature;

          return cached_feature;
        }
      }
    }

    final FeatureEnumeration enumerator = features ();

    int i = 0;

    while (enumerator.hasMoreFeatures ()) {
      final Feature this_feature = enumerator.nextFeature ();

      if (i == arg_index) {
        cached_feature = this_feature;
        cached_feature_index = i;

        return this_feature;
      }

      ++i;
    }

    return null;
  }

  /**
   *  Returns true if and only if this FeatureTable contains the given feature.
   **/
  public boolean contains (final Feature feature) {
    return getFeatures ().contains (feature);
  }

  /**
   *  Return the index of the given Feature.  This does the reverse of
   *  getFeatureAtIndex ().
   **/
  int indexOf (final Feature feature) {
    if (cached_feature == feature) {
      return cached_feature_index;
    }

    final FeatureEnumeration enumerator = features ();

    int i = 0;

    while (enumerator.hasMoreFeatures ()) {
      final Feature this_feature = enumerator.nextFeature ();

      if (this_feature == feature) {
        return i;
      }
      ++i;
    }

    return -1;
  }

  /**
   *  Returns an enumeration of the Feature objects in this FeatureTable.  The
   *  returned Enumeration object will generate all features in this object in
   *  turn. The first item generated is the item at index 0, then the item at
   *  index 1, and so on.
   **/
  public FeatureEnumeration features () {
    // return an enumeration from the FeatureTree
    return getFeatures ().features ();
  }

  /**
   *  This holds the last Feature that was returned from getFeatureAtIndex ()
   *  or the last Feature that was passed to indexOf ().
   **/
  private Feature cached_feature = null;

  /**
   *  This holds the last index that was passed to getFeatureAtIndex () or the
   *  last index that was returned from indexOf ().  This is used to improve
   *  the speed of those two methods.
   **/
  private int cached_feature_index = -999;

  /**
   *  This holds the features of this FeatureTable
   **/
  final private FeatureTree features =
    new FeatureTree (new FeatureComparator ());
}
