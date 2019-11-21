/* FeatureTree.java
 *
 * created: Thu Jan 28 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/FeatureTree.java,v 1.3 2005-11-28 16:46:38 tjc Exp $
 **/

package uk.ac.sanger.artemis.io;

import java.util.*;

/**
 *  A tree that stores StreamFeature objects ordered with a StreamFeatureComparator
 *  object.
 *
 *  @author Kim Rutherford
 **/

public class FeatureTree extends TreeSet {
  /**
   *  Create a new (empty) FeatureTree.
   **/
  public FeatureTree (final Comparator comparator) {
    super (comparator);

    this.comparator = comparator;
  }

  /**
   *  Wrapper for TreeSet.add () which sets size_of_largest_feature_seen.
   **/
  public synchronized boolean add (final Object element) {
    final Feature this_feature = (Feature) element;

    getBucket (this_feature).add (this_feature);

    return super.add (element);
  }

  /**
   *  Wrapper for TreeSet.() which removes the Feature from
   *  rbtree_buckets.
   **/
  public synchronized boolean remove (Object element) {
    final Feature this_feature = (Feature) element;

    getBucket (this_feature).remove (this_feature);

    return super.remove (element);
  }

  /**
   *  Find the sub-tree which has a first base that is >= the given base.
   **/
  private static SortedSet findByBase (final TreeSet tree, final int base) {
    final ComparableFeature test_feature = new ComparableFeature () {
      public void set (final Key k, final Location l,
                       final QualifierVector qv) {}
      public void setKey (final Key k) {}
      public void setLocation (final Location loc) {}
      public void setLocation (final Location loc, Entry entry) {}
      public void setQualifiers (final QualifierVector qualVec) {}
      public void setQualifier (final Qualifier qual) {}
      public void removeQualifierByName (final String qualName) {}
      private final Key dummy_key = new Key ("_dummy_key_");
      public Key getKey () {return dummy_key;}
      public Location getLocation () {return null;}
      public QualifierVector getQualifiers () {return null;}
      public Qualifier getQualifierByName (final String qualName) {return null;}
      public int getFirstBase () {return base;}
      public int getLastBase () {return base;}
      public long getNumericID () {return -1;}
      public Entry getEntry () {return null;}
      public Feature copy () {return null;}
      public void setUserData (final Object userData) {}
      public Object getUserData () {return null;}
      public boolean isReadOnly () {return false;}
    };

    return tree.tailSet (test_feature);
  }

  /**
   *  Add to features_in_range all the features in the given TreeSet which are
   *  within the given Range.
   *  @param max_feature_length the maximum length in bases of the features in
   *    the TreeSet.
   **/
  private void getFeaturesInRange (final TreeSet tree,
                                   final FeatureVector features_in_range,
                                   final Range range,
                                   final int max_feature_length) {
	
		
    // find the leftmost node in the range
    final SortedSet<Feature> tail_set =
      findByBase (tree, range.getStart () - max_feature_length);

    final int rangeEnd = range.getEnd ();

    // now loop over all the features that could possibly be in the range
    // (ie. those between range.getStart () - max_feature_length and
    // range.getEnd ())
    for (Feature this_feature : tail_set) {

      if (this_feature.getFirstBase () > rangeEnd) {
        return;
      }
      
      if (this_feature.getLocation ().getTotalRange ().overlaps (range)) {
        features_in_range.add (this_feature);
      }
    }
  }

  /**
   *  Return a vector containing the references of the Feature objects within
   *  the given range.
   *  @param range Return features that overlap this range - ie the start of
   *    the feature is less than or equal to the end of the range and the end
   *    of the feature is greater than or equal to the start of the range.
   *  @return The features that are within the given range.  The returned
   *    object is a copy - changes will not effect the FeatureTree object
   *    itself.
   **/
  public synchronized FeatureVector getFeaturesInRange (final Range range) {
    // this default size will cover many common cases
    final FeatureVector return_features = new FeatureVector();

    final int numBuckets = rbtree_buckets.size();
    
    for (int i = 0; i < numBuckets; ++i) {
    
      TreeSet bucket = (TreeSet)rbtree_buckets.elementAt (i);
    	
      if (bucket != null && bucket.size() > 0) {
    	  getFeaturesInRange (bucket,
                          return_features, range,
                          (int) Math.pow (BUCKET_MULTIPLIER, i + 1));
      }
    }
    
    return return_features;
  }

  /**
   *  Returns an enumeration of the Feature objects in this FeatureTree.  The
   *  returned Enumeration object will generate all features in this object in
   *  turn.
   **/
  public FeatureEnumeration features () {
    return new FeatureEnumerator ();
  }

  /**
   *  Return the reference of the Feature that is immediately after the
   *  argument Feature in the tree or null if there is no next Feature.
   *  The tree must contain this_feature.
   **/
  public Feature getNextFeature (final Feature this_feature) {
    final SortedSet tail_set = tailSet (this_feature);

    final Iterator tail_set_iterator = tail_set.iterator ();

    // must contain this_feature
    tail_set_iterator.next ();

    if (tail_set_iterator.hasNext ()) {
      return (Feature) tail_set_iterator.next ();
    } else {
      return null;
    }
  }

  /**
   *  An Enumeration of Feature objects.
   **/
  public class FeatureEnumerator implements FeatureEnumeration {
    /**
     *  Create a new FeatureEnumeration that will enumerate the enclosing
     *  DocumentEntry object.  The DocumentEntry object must not be changed
     *  while the enumeration is active.
     **/
    public FeatureEnumerator () {
      iterator = iterator ();
    }

    /**
     *  See the FeatureEnumeration interface for details.
     **/
    public boolean hasMoreFeatures () {
      return iterator.hasNext ();
    }

    /**
     *  See the FeatureEnumeration interface for details.
     **/
    public Feature nextFeature ()
        throws NoSuchElementException {
      return (Feature) iterator.next ();
    }

    /**
     *  This is the underlying enumeration that does all the work for us
     **/
    private Iterator iterator;
  }

  /**
   *  Features are stored in buckets.
   *  Bucket 0 will hold a reference to all features that are less than
   *  BUCKET_MULTIPLIER in length.  Bucket 1 will hold those from
   *  BUCKET_MULTIPLIER to BUCKET_MULTIPLIER*BUCKET_MULTIPLIER-1 (inclusive)
   *  in length.  Bucket 2 will hold those from BUCKET_MULTIPLIER^2 to
   *  BUCKET_MULTIPLIER^3-1 (inclusive) in length.  etc.
   **/
  private TreeSet getBucket (final Feature feature) {
    final int feature_length =
      feature.getLocation ().getTotalRange ().getCount ();

    final int feature_bucket;

    if (feature_length <= BUCKET_MULTIPLIER) {
      feature_bucket = 0;
    } else {
      // subtract 0.5 to feature_length to make sure that rounding errors go
      // in out favour
      feature_bucket =
        (int) (Math.log (feature_length + 0.5) / Math.log (BUCKET_MULTIPLIER));
    }

    // make the bucket we need and all smaller buckets
    while (rbtree_buckets.size () <= feature_bucket) {
      rbtree_buckets.addElement (new TreeSet (comparator));
    }

    return (TreeSet) rbtree_buckets.elementAt (feature_bucket);
  }

  /**
   *  The maximum number of buckets.
   **/
  private final int BUCKET_COUNT = 10;

  /**
   *  See comment on getBucket().
   **/
  private final int BUCKET_MULTIPLIER = 4;

  /**
   *  See comment above.
   **/
  private final Vector rbtree_buckets = new Vector (BUCKET_COUNT);

  /**
   *  The Comparator that was passed to the constructor.
   **/
  private Comparator comparator;
}
