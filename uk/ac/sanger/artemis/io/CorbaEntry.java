/* CorbaEntry.java
 *
 * created: Wed Dec 30 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/CorbaEntry.java,v 1.2 2008-06-11 15:12:20 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.*;

import nsdb.EmblSeq;
import nsdb.NucFeature;
import type.NoResult;

import java.util.Vector;
import java.io.IOException;
import java.util.NoSuchElementException;

/**
 *  This class extends the Entry class with the data for the entry coming from
 *  a Corba server.
 *
 *  @author Kim Rutherford
 *  @version $Id: CorbaEntry.java,v 1.2 2008-06-11 15:12:20 tjc Exp $
 **/

public class CorbaEntry extends ReadOnlyEntry
    implements Entry {
  /**
   *  Create a new CorbaEntry object from the given handle.
   *  @param entry_information The EntryInformation object of the new Entry.
   *  @param data This is the corba object that we will read from.
   *  @exception EntryInformationException Thrown if this
   *    Entry cannot contain the Key, Qualifier or Key/Qualifier combination of
   *    one of the features in the given Entry.
   **/
  public CorbaEntry (final EntryInformation entry_information,
                     final EmblSeq corba_handle)
      throws LocationParseException, InvalidKeyException, NoResult {
    this.corba_handle = corba_handle;
    this.sequence = new CorbaSequence (corba_handle);
    this.entry_information = entry_information;

    grabFeatures ();
  }

  /**
   *  Return the text of the EMBL header of this Entry or null if there is no
   *  header.
   **/
  public String getHeaderText () {
    return null;
  }

  /**
   *  Return the name of this Entry.
   **/
  public String getName () {
    return corba_handle.getBioSeqId ();
  }

  /**
   *  Return a count of the number of Feature objects in this Entry.
   **/
  public int getFeatureCount () {
    return features.size ();
  }

  /**
   *  Return the ith Feature from this Entry.  This Features are returned in a
   *  consistent order, sorted by the first base of each Feature.
   **/
  public Feature getFeatureAtIndex (int arg_index) {

    final FeatureEnumeration enumerator = features.features ();

    int i = 0;

    while (enumerator.hasMoreFeatures ()) {
      final Feature this_feature = enumerator.nextFeature ();
      
      if (i == arg_index) {
        return this_feature;
      }

      ++i;
    }

    return null;
  }

  /**
   *  Return the index of the given Feature.  This does the reverse of
   *  getFeatureAtIndex ().
   **/
  public int indexOf (Feature feature) {
    final FeatureEnumeration enumerator = features.features ();

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
   *  Returns true if and only if this Entry contains the given feature.
   **/
  public boolean contains (Feature feature) {
    return features.contains (feature);
  }

  /**
   *  Returns an enumeration of the Feature objects in this Entry.  The
   *  returned Enumeration object will generate all features in this object in
   *  turn. The first item generated is the item at index 0, then the item at
   *  index 1, and so on.
   **/
  public FeatureEnumeration features () {
    return features.features ();
  }

  /**
   *  Return a vector containing the references of the Feature objects within
   *  the given range.
   *  @param range Return features that overlap this range - ie the start of
   *    the feature is less than or equal to the end of the range and the end
   *    of the feature is greater than or equal to the start of the range.
   *  @return The features of this feature table the are within
   *    the given range.  The returned object is a copy - changes will not
   *    effect the FeatureTable object itself.
   **/
  public FeatureVector getFeaturesInRange (Range range) {
    return features.getFeaturesInRange (range);
  }

  /**
   *  Return a vector containing the references of all the Feature objects in
   *  this Entry.
   *  @return The features of this Entry.  The returned object
   *    is a copy - changes will not effect the Entry object itself.
   **/
  public FeatureVector getAllFeatures () {
    final FeatureVector return_features = new FeatureVector ();

    final FeatureEnumeration enumerator = features.features ();

    while (enumerator.hasMoreFeatures ()) {
      final Feature this_feature = enumerator.nextFeature ();
      
      return_features.add (this_feature);
    }

    return return_features;
  }

  /**
   *  Return the Sequence object from this entry or null if it does not
   *  contain one.
   *  @return a Sequence object for this Entry.  the returned object is
   *    not a copy - changes to it will change the Entry object itself
   **/
  public Sequence getSequence () {
    return sequence;
  }

  /**
   *  Return the EntryInformation object that was passed to the constructor.
   **/
  public EntryInformation getEntryInformation () {
    return entry_information; 
  }

  /**
   *  Read the features from the corba_handle and put them in a FeatureTree.
   **/
  private void grabFeatures ()
      throws LocationParseException, InvalidKeyException, NoResult {
    final NucFeature [] feature_handles = corba_handle.getNucFeatures ();

    for (int i = 0 ; i < feature_handles.length ; ++i) {
      try {
        final Feature new_feature = new CorbaFeature (feature_handles [i]);

        features.add (new_feature);
      } catch (InvalidRelationException e) {
        System.out.println ("exception while reading: " + e);
      }
    }
  }

  /**
   *  The EmblSeq object that was passed to the constructor.
   **/
  private EmblSeq corba_handle;

  /**
   *  This is created in the constructor 
   **/
  private CorbaSequence sequence = null;

  /**
   *  The FeatureTree containing the features of this Entry.
   **/
  private FeatureTree features = new FeatureTree (new FeatureComparator ());

  /**
   *  The EntryInformation object that was passed to the constructor.
   **/
  final private EntryInformation entry_information;

  public void dispose()
  {
    // TODO Auto-generated method stub
    
  }
}
