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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/RWCorbaEntry.java,v 1.3 2008-06-11 15:12:20 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.*;

import nsdb.EmblSeqWriter;
import nsdb.NucFeatureWriter;
import nsdb.NucFeature;
import type.NoResult;

import java.util.Vector;
import java.io.IOException;
import java.util.NoSuchElementException;

/**
 *  This class implements the Entry interface with the data for the entry
 *  coming from a Corba server.
 *
 *  @author Kim Rutherford
 *  @version $Id: RWCorbaEntry.java,v 1.3 2008-06-11 15:12:20 tjc Exp $
 **/

public class RWCorbaEntry extends EMBLObject
    implements Entry {
  /**
   *  Create a new CorbaEntry object from the given handle.
   *  @param entry_information The EntryInformation object of the new Entry.
   *  @param data This is the corba object that we will read from.
   *  @exception EntryInformationException Thrown if this
   *    Entry cannot contain the Key, Qualifier or Key/Qualifier combination of
   *    one of the features in the given Entry.
   **/
  public RWCorbaEntry (final EntryInformation entry_information,
                       final EmblSeqWriter corba_handle)
      throws LocationParseException, InvalidKeyException, NoResult,
             EntryInformationException {
    this.corba_handle = corba_handle;
    this.sequence = new CorbaSequence (corba_handle);
    this.entry_information = entry_information;

    if (this.sequence.length () == 0) {
      this.sequence = null;
    }
  }
 
 /**
   *  Returns true if and only if there have been some changes to this Entry
   *  since the last save.
   **/
  public boolean hasUnsavedChanges () {
    // XXX implement me (if necessary)

    return false;
  }

  /**
   *  Returns true if and only if this entry is read only.
   **/
  public boolean isReadOnly () {
    return false;
  }

  /**
   *  Return the text of the EMBL header of this Entry or null if there is no
   *  header.
   **/
  public String getHeaderText () {
    return null;
  }

  /**
   *  Set the header of this Entry to be the given text.
   *  @return true if and only if the header was successfully set.  Not all
   *    Entry objects can change their header, so it is up to the calling
   *    function to check the return value.
   *  @exception IOException thrown if there is a problem reading the header
   *    from the String - most likely ReadFormatException.
   **/
  public boolean setHeaderText (final String new_header) {
    return false;
  }

  /**
   *  Write this entry to the file/database/whatever it came from.  If this
   *  object is read only a ReadOnlyException will be thrown.
   **/
  public void save () throws IOException {
    try {
      corba_handle.commit ();
    } catch (nsdb.CommitFailed e) {
      throw new IOException ("save failed: " + e.reason);
    }
  }

  /**
   *  Return the name of this Entry.
   **/
  public String getName () {
    return "CORBA-" + corba_handle.getBioSeqId ();
  }

  /**
   *  Set the name of this Entry - if possible (the return value will let the
   *  caller know).
   *  @return true if and only if the name was successfully set.  Not all
   *    Entry objects can change there name, so it is up to the calling
   *    function to check the return value.
   **/
  public boolean setName (final String name) {
    return false;
  }

  /**
   *  Create a new Feature object of an appropriate type in this Entry.
   *  @param key The new feature key
   *  @param location The Location object for the new feature
   *  @param qualifiers The qualifiers for the new feature
   **/
  public Feature createFeature (Key key,
                                Location location,
                                QualifierVector qualifiers)
      throws InvalidRelationException, ReadOnlyException,
             EntryInformationException {
    final EmblSeqWriter embl_writer =
      nsdb.EmblSeqWriterHelper.narrow (corba_handle);
    if (embl_writer == null) {
      throw new ReadOnlyException ();
    } else {
      try {
        final nsdb.NucFeatureWriter nuc_feature_writer =
          embl_writer.createNucFeature (key.toString (), location.toString ());

        final RWCorbaFeature corba_feature =
          new RWCorbaFeature (getEntryInformation (), nuc_feature_writer);

        corba_feature.setRWCorbaEntry (this);

        corba_feature.setQualifiers (qualifiers);

        final NucFeatureHasher hasher =
          new NucFeatureHasher (nuc_feature_writer);

        feature_dictionary.put (hasher, corba_feature);

        return corba_feature;
      } catch (nsdb.LocationParse e) {
        throw new Error ("internal error - server doesn't understand " +
                         "this location: " + location.toString ());
      } catch (type.IndexOutOfRange e) {
        throw new Error ("internal error - location out of range: " + e);
      } catch (nsdb.InvalidKey e) {
        throw new Error ("internal error - server won't use " +
                         "this feature key: " + key.toString ());
      } catch (type.NoResult e) {
        throw new Error ("internal error - unexpected exception: " + e);
      } catch (nsdb.ReadOnlyException e) {
        throw new ReadOnlyException ();
      }
    }
  }

  /**
   *  Return a count of the number of Feature objects in this Entry.
   **/
  public int getFeatureCount () {
    return corba_handle.getNucFeatureCount ();
  }

  /**
   *  Add the given Feature to this Entry.  If the Feature is already in an
   *  Entry then Entry.remove () should be called on that Entry before calling
   *  Entry.add ().  An Error will be thrown otherwise.
   *  @exception ReadOnlyException If this entry is read only.
   *  @exception EntryInformationException Thrown if this Entry
   *    cannot contain the Key, Qualifier or Key/Qualifier combination of the
   *    given Feature.
   *  @return A reference that was passed to add (), if that Feature can be
   *    stored directly in this Entry, otherwise returns a reference to a new
   *    Feature, that is a copy of the argument.  The argument reference
   *    should not be used after the call to add (), unless the return
   *    reference happens to be the same as the argument.
   **/
  public Feature add (final Feature feature)
      throws EntryInformationException, ReadOnlyException {
    if (feature.getEntry () != null) {
      throw new Error ("internal error - a feature must have one owner");
    }

    try {
      return createFeature (feature.getKey (),
                            feature.getLocation (),
                            feature.getQualifiers ());
    } catch (InvalidRelationException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    }
  }

  /**
   *  Add the given Feature to this Entry.  If the Feature is already in an
   *  Entry then Entry.remove() should be called on that Entry before calling
   *  Entry.forcedAdd() (An Error will be thrown otherwise).  Invalid
   *  qualifiers will be quietly thrown away.  Features with invalid keys will
   *  not be added (and null will be returned).  "Invalid" means that the
   *  key/qualifier is non allowed to occur in an Entry of this type (probably
   *  determined by the EntryInformation object of this Entry).
   *  @exception ReadOnlyException If this entry is read only.
   *  @return A reference that was passed to add (), if that Feature can be
   *    stored directly in this Entry, otherwise returns a reference to a new
   *    Feature, that is a copy of the argument.  The argument reference
   *    should not be used after the call to add (), unless the return
   *    reference happens to be the same as the argument.  Returns null if and
   *    only if the new Feature has a key that is invalid for this Entry.
   **/
  public Feature forcedAdd (final Feature feature)
      throws ReadOnlyException {
    if (feature.getEntry () != null) {
      throw new Error ("internal error - a feature must have one owner");
    }

    try {
      return createFeature (feature.getKey (),
                            feature.getLocation (),
                            feature.getQualifiers ());
    } catch (InvalidRelationException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    } catch (EntryInformationException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    }
  }

  /**
   *  Remove the given Feature from this Entry.
   *  @return true if and only if the Feature was in this Entry.
   *  @exception ReadOnlyException If this entry is read only.
   **/
  public boolean remove (Feature feature)
      throws ReadOnlyException {
    if (feature instanceof RWCorbaFeature) {
      final RWCorbaFeature corba_feature = (RWCorbaFeature) feature;
      final NucFeatureWriter nuc_feature =
        corba_feature.getNucFeatureWriter ();
      if (contains (nuc_feature)) {
        try {
          corba_handle.remove (nuc_feature);
          corba_feature.setRWCorbaEntry (null);
          return true;
        } catch (nsdb.ReadOnlyException e) {
          throw new ReadOnlyException ();
        }
      } else {
        return false;
      }
    } else {
      return false;
    }
  }

  /**
   *  Return the ith Feature from this Entry.  This Features are returned in a
   *  consistent order, sorted by the first base of each Feature.
   **/
  public Feature getFeatureAtIndex (int index) {
    try {
      final NucFeature nuc_feature =
        corba_handle.getFeatureAtIndex (index);

      final NucFeatureWriter nuc_feature_writer =
        nsdb.NucFeatureWriterHelper.narrow (nuc_feature);

      return getFeatureOfNucFeature (nuc_feature_writer);
    } catch (type.IndexOutOfRange e) {
      throw new IndexOutOfBoundsException ();
    }
  }

  /**
   *  Return the index of the given Feature.  This does the reverse of
   *  getFeatureAtIndex ().
   **/
  public int indexOf (Feature feature) {
    if (!(feature instanceof RWCorbaFeature)) {
      return -1;
    }

    final NucFeature nuc_feature =
      ((RWCorbaFeature)feature).getNucFeatureWriter ();
      
    return corba_handle.indexOf (nuc_feature);
  }

  /**
   *  Returns true if and only if this Entry contains the given feature.
   **/
  public boolean contains (Feature feature) {
    if (feature instanceof RWCorbaFeature) {
      final NucFeatureWriter nuc_feature =
        ((RWCorbaFeature) feature).getNucFeatureWriter ();
      if (contains (nuc_feature)) {
        return true;
      } else {
        return false;
      }
    } else {
      return false;
    }
  }

  /**
   *  Returns an enumeration of the Feature objects in this Entry.  The
   *  returned Enumeration object will generate all features in this object in
   *  turn. The first item generated is the item at index 0, then the item at
   *  index 1, and so on.
   **/
  public FeatureEnumeration features () {
    return new FeatureEnumeration () {
        public boolean hasMoreFeatures () {
          return i < features.size ();
        }

        public Feature nextFeature () {
          return features.featureAt (i++);
        }

        private int i = 0;

        private FeatureVector features = getAllFeatures ();
      };
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
  public FeatureVector getFeaturesInRange (Range range)
      throws OutOfRangeException {
    final FeatureVector return_features = new FeatureVector ();

    try {
      final NucFeature [] feature_handles =
        corba_handle.getNucFeaturesInRange (range.getStart (),
                                            range.getEnd ());

      for (int i = 0 ; i < feature_handles.length ; ++i) {
        final NucFeatureWriter this_nuc_feature =
          nsdb.NucFeatureWriterHelper.narrow (feature_handles [i]);
        final RWCorbaFeature rw_corba_feature =
          getFeatureOfNucFeature (this_nuc_feature);
        return_features.add (rw_corba_feature);
      }
    } catch (NoResult e) {

    } catch (type.IndexOutOfRange e) {
      throw new OutOfRangeException (range.toString ());
    }

    return return_features;
  }

  /**
   *  Return a vector containing the references of all the Feature objects in
   *  this Entry.
   *  @return The features of this Entry.  The returned object
   *    is a copy - changes will not effect the Entry object itself.
   **/
  public FeatureVector getAllFeatures () {
    final FeatureVector return_features = new FeatureVector ();

    try {

      final NucFeature [] feature_handles =
        corba_handle.getNucFeatures ();

      for (int i = 0 ; i < feature_handles.length ; ++i) {
        final NucFeatureWriter this_nuc_feature =
          nsdb.NucFeatureWriterHelper.narrow (feature_handles [i]);
        
        final Feature this_feature = getFeatureOfNucFeature (this_nuc_feature);

        return_features.add (this_feature);
      }
    } catch (NoResult e) {

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
   *  Returns the RWCorbaFeature reference that corresponds to the given
   *  NucFeatureWriter reference.  If no RWCorbaFeature has yet been created
   *  for the given NucFeatureWriter, one will be created and returned.  This
   *  method uses and sets feature_dictionary.  This method will always return
   *  the same RWCorbaFeature reference for a given NucFeatureWriter.
   **/
  private RWCorbaFeature getFeatureOfNucFeature (final NucFeatureWriter
                                                   nuc_feature) {
    final NucFeatureHasher hasher = new NucFeatureHasher (nuc_feature);

    final RWCorbaFeature feature =
      (RWCorbaFeature) feature_dictionary.get (hasher);

    if (feature == null || feature.getEntry () == null) {
      // make a feature
      try {
        final RWCorbaFeature new_feature =
          new RWCorbaFeature (getEntryInformation (), nuc_feature);

//        System.out.println ("created: " + new_feature);

        feature_dictionary.put (hasher, new_feature);

        new_feature.setRWCorbaEntry (this);
        
        return new_feature;
      } catch (NoResult e) {
        throw new Error ("internal error - unexpected exception: " + e);
      } catch (EntryInformationException e) {
        throw new Error ("internal error - unexpected exception: " + e);
      }
    } else {
      return feature;
    }
  }

  /**
   *  Returns true if and only if the given NucFeatureWriter object is a
   *  feature in this entry.
   **/
  private boolean contains (final NucFeatureWriter nuc_feature) {
    final NucFeatureHasher hasher = new NucFeatureHasher (nuc_feature);

    final Object feature = feature_dictionary.get (hasher);

    if (feature == null) {
      return false;
    } else {
      return true;
    }
  }

  /**
   *  The EmblSeqWriter object that was passed to the constructor.
   **/
  private EmblSeqWriter corba_handle;

  /**
   *  This is created in the constructor
   **/
  private CorbaSequence sequence = null;

  /**
   *  The Dictionary is used to find the the embl.Feature reference that
   *  corresponds any given NucFeatureWriter object.  (see
   *  getFeatureOfNucFeature ()).
   **/
  private java.util.Hashtable feature_dictionary =
    new java.util.Hashtable ();

  /**
   *  The EntryInformation object that was passed to the constructor.
   **/
  final private EntryInformation entry_information;

  public void dispose()
  {
    // TODO Auto-generated method stub
    
  }
}


class NucFeatureHasher {
  NucFeatureHasher (final NucFeatureWriter nuc_feature) {
    this.nuc_feature = nuc_feature;
  }

  public int hashCode () {
    return nuc_feature._hash (Integer.MAX_VALUE);
  }

  public boolean equals (final Object object) {
    final NucFeatureWriter other_nuc_feature =
      ((NucFeatureHasher)object).nuc_feature;

//    return nuc_feature._is_equivalent (other_nuc_feature);

    // this works for JacORB, but may for work for other orbs too
    return nuc_feature.toString ().equals (other_nuc_feature.toString ());
  }

  final public NucFeatureWriter nuc_feature;
}
