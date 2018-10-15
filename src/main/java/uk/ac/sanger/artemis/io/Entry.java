/* Entry.java
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/Entry.java,v 1.2 2008-06-11 15:12:20 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.ReadOnlyException;
import uk.ac.sanger.artemis.util.OutOfRangeException;

import java.io.*;

/**
 *  Objects of this class represent one EMBL entry.
 *
 *  @author Kim Rutherford
 *  @version $Id: Entry.java,v 1.2 2008-06-11 15:12:20 tjc Exp $
 *
 **/

public interface Entry 
{
  /**
   *  Write this entry to the file/database/whatever it came from.  If this
   *  object is read only a ReadOnlyException will be thrown.
   **/
  void save() throws IOException;

  /**
   *  Returns true if and only if there have been some changes to this Entry
   *  since the last save.
   **/
  boolean hasUnsavedChanges();

  /**
   *  Returns true if and only if this entry is read only.
   **/
  boolean isReadOnly();

  /**
   *  Return the text of the EMBL header of this Entry or null if there is no
   *  header.
   **/
  String getHeaderText();

  /**
   *  Set the header of this Entry to be the given text.
   *  @return true if and only if the header was successfully set.  Not all
   *    Entry objects can change their header, so it is up to the calling
   *    function to check the return value.
   *  @exception IOException thrown if there is a problem reading the header
   *    from the String - most likely ReadFormatException.
   **/
  boolean setHeaderText (final String new_header)
      throws IOException;

  /**
   *  Return the name of this Entry.
   **/
  String getName();

  /**
   *  Set the name of this Entry - if possible (the return value will let the
   *  caller know).
   *  @return true if and only if the name was successfully set.  Not all
   *    Entry objects can change their name, so it is up to the calling
   *    function to check the return value.
   **/
  boolean setName(final String name);

  /**
   *  Create a new Feature object of an appropriate type in this Entry.
   *  @param key The new feature key
   *  @param location The Location object for the new feature
   *  @param qualifiers The qualifiers for the new feature (can be null if
   *    there are no qualifiers).
   *  @exception InvalidRelationException Thrown if this Feature cannot contain
   *    the given Qualifier.
   *  @exception EntryInformationException Thrown if a Feature in this Entry
   *    cannot contain the given Key, Qualifier or Key/Qualifier combination.
   **/
  Feature createFeature(Key key,
                        Location location,
                        QualifierVector qualifiers)
      throws EntryInformationException, ReadOnlyException, OutOfRangeException;

  /**
   *  Return a count of the number of Feature objects in this Entry.
   **/
  int getFeatureCount();

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
  Feature add(Feature feature)
      throws EntryInformationException, ReadOnlyException;

  /**
   *  Add the given Feature to this Entry.  If the Feature is already in an
   *  Entry then Entry.remove () should be called on that Entry before calling
   *  Entry.add ().  An Error will be thrown otherwise.  Invalid qualifiers
   *  will be quietly thrown away.  Features with invalid keys will not be
   *  added (and null will be returned).  "Invalid" means that the
   *  key/qualifier is non allowed to occur in an Entry of this type (probably
   *  determined by the EntryInformation object of this Entry).
   *  @exception ReadOnlyException If this entry is read only.
   *  @exception EntryInformationException Thrown if this Entry
   *    cannot contain the Key, Qualifier or Key/Qualifier combination of the
   *    given Feature
   *  @return A reference that was passed to add (), if that Feature can be
   *    stored directly in this Entry, otherwise returns a reference to a new
   *    Feature, that is a copy of the argument.  The argument reference
   *    should not be used after the call to add (), unless the return
   *    reference happens to be the same as the argument.  Returns null if and
   *    only if the new Feature has a key that is invalid for this Entry.
   **/
  Feature forcedAdd(Feature feature)
      throws ReadOnlyException;

  /**
   *  Remove the given Feature from this Entry.
   *  @return true if and only if the Feature was in this Entry.
   *  @exception ReadOnlyException If this entry is read only.
   **/
  boolean remove(Feature feature)
      throws ReadOnlyException;

  /**
   *  Return the ith Feature from this Entry.  The features are returned in a
   *  consistent order, sorted by the first base of each Feature.
   **/
  Feature getFeatureAtIndex (int i);

  /**
   *  Return the index of the given Feature.  This does the reverse of
   *  getFeatureAtIndex ().
   **/
  int indexOf(Feature feature);

  /**
   *  Returns true if and only if this Entry contains the given feature.
   **/
  boolean contains(Feature feature);

  /**
   *  Returns an enumeration of the Feature objects in this Entry.  The
   *  returned Enumeration object will generate all features in this object in
   *  turn. The first item generated is the item at index 0, then the item at
   *  index 1, and so on.
   **/
  FeatureEnumeration features();

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
  FeatureVector getFeaturesInRange(Range range)
      throws OutOfRangeException;

  /**
   *  Return a vector containing the references of all the Feature objects in
   *  this Entry.
   *  @return The features of this Entry.  The returned object
   *    is a copy - changes will not effect the Entry object itself.
   **/
  FeatureVector getAllFeatures();

  /**
   *  Return the Sequence object from this entry or null if it does not
   *  contain one.
   *  @return a Sequence object for this Entry.  the returned object is
   *    not a copy - changes to it will change the Entry object itself
   **/
  Sequence getSequence();

  /**
   *  Return the EntryInformation object for this Entry.
   **/
  EntryInformation getEntryInformation();
  
  /**
   * Dispose of entry objects
   */
  void dispose();
}
