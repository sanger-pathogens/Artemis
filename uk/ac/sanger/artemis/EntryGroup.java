/* EntryGroup.java
 *
 * created: Wed Nov 11 1998
 *
 * This file is part of Artemis
 *
 * Copyright(C) 1998,1999,2000,2001,2002  Genome Research Limited
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or(at your option) any later version.
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/EntryGroup.java,v 1.1 2004-06-09 09:44:21 tjc Exp $
 **/

package uk.ac.sanger.artemis;

import uk.ac.sanger.artemis.sequence.Bases;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.util.ReadOnlyException;
import uk.ac.sanger.artemis.util.OutOfRangeException;

/**
 *  Objects that implement this interface contain a vector of Entry object,
 *  with additional methods for querying and changing the feature tables of
 *  all the entries at once.
 *
 *  @author Kim Rutherford
 *  @version $Id: EntryGroup.java,v 1.1 2004-06-09 09:44:21 tjc Exp $
 **/

public interface EntryGroup
    extends FeatureChangeListener, EntryChangeListener
{
  /**
   *  Return the default Entry for this EntryGroup.  The "default" is the
   *  Entry where new features are created.
   **/
  Entry getDefaultEntry();

  /**
   *  Set the default Entry.  The "default" is the Entry where new features
   *  are created.
   *  @param entry The new default entry.  If this Entry is not active this
   *    method will return immediately.
   **/
  void setDefaultEntry(final Entry entry);

  /**
   *  Returns true if and only if there are any unsaved changes in any of the
   *  Entry objects in this EntryGroup.
   **/
  boolean hasUnsavedChanges();

  /**
   *  Return the index of a feature within this object.  This method treats
   *  all the features in all the active entries as if they were in one big
   *  array.  The first feature of the first entry will have index 1, the
   *  first from the second entry will have index 1 +(the number of features
   *  in the first entry), etc.
   *  @param feature The feature to find the index of.
   *  @return The index of the feature or -1 if the feature isn't in any of
   *    the entries.  The first index is 0 the last is the total number of
   *    features in all the entries of this object minus one.
   **/
  int indexOf(final Feature feature);

  /**
   *  Return the index of an Entry within this object.
   *  @return The index of the Entry or -1 if the Entry isn't in this Object.
   **/
  int indexOf(final Entry entry);

  /**
   *  Return true if any of the active entries in the group contains the given
   *  feature.
   **/
  boolean contains(final Feature feature);

  /**
   *  Return true if the given Entry is active(visible).  The Feature objects
   *  in an Entry that is not active will be ignored by the methods that deal
   *  will features: featureAt(), indexOf(), contains(), features(), etc.
   **/
  boolean isActive(final Entry entry);

  /**
   *  Set the "active" setting of the Entry at the given index.  If the index
   *  refers to the default entry and new_active is false, the default entry
   *  will be set to the active entry or null if there are no active entries.
   *  @param index The index of the Entry to change.
   *  @param new_active The new active setting.
   **/
  void setIsActive(final int index, final boolean new_active);

  /**
   *  Set the "active" setting of the given Entry.  The Entry is the default
   *  entry and new_active is false, the default entry will be set to the
   *  active entry or null if there are no active entries.
   *  @param entry The Entry to activate or deactivate.
   *  @param new_active The new active setting.
   **/
  void setIsActive(final Entry entry, final boolean new_active);

  /**
   *  Return the Entry from this EntryGroup that contains the sequence to view
   *  or return null if none of the entries contains a sequence.
   **/
  Entry getSequenceEntry();

  /**
   *  Returns the base length of the sequence of the first Entry in this group
   *  or 0 if this group is empty.
   **/
  int getSequenceLength();

  /**
   *  Returns the Bases object of the first Entry in this group or null if
   *  this group is empty.
   **/
  Bases getBases();

  /**
   *  Reverse and complement the sequence and all features in every Entry in
   *  this EntryGroup.
   **/
  void reverseComplement()
      throws ReadOnlyException;

  /**
   *  Return the number of entries in this EntryGroup.
   **/
  int size();

  /**
   *  Return the ith Entry in this EntryGroup.
   **/
  Entry elementAt(final int i);

  /**
   *  Return the Feature at the given index.  This method treats all the
   *  features in all the entries as if they were in one big array.  See
   *  the comment on indexOf().
   *  @param index The index of the required Feature.
   *  @return The Feature at the given index.  The first index is 0 the last
   *    is the total number of features in all the entries of this object minus
   *    one.  If the index is out of range then null will be returned.
   **/
  Feature featureAt(final int index);

  /**
   *  Return a vector containing the references of the Feature objects within
   *  the given range of indices.
   *  @param
   **/
  FeatureVector getFeaturesInIndexRange(final int start_index,
                                         final int end_index);

  /**
   *  Return a vector containing the references of the Feature objects within
   *  the given range for all the active entries in the EntryGroup.
   *  @param range Return features that overlap this range - ie the start of
   *    the feature is less than or equal to the end of the range and the end
   *    of the feature is greater than or equal to the start of the range.
   *  @return The non-source key features of this feature table the are within
   *    the given range.  The returned object is a copy - changes will not
   *    effect the FeatureTable object itself.
   **/
  FeatureVector getFeaturesInRange(final Range range)
      throws OutOfRangeException;

  /**
   *  Return a vector containing the references of the Feature objects from
   *  all the active entries in the EntryGroup.
   *  @return The non-source key features in active entries of this
   *    EntryGroup.  The returned object is a copy - changes will not effect
   *    the EntryGroup object itself.
   **/
  FeatureVector getAllFeatures();

  /**
   *  Return a count of the number of Feature objects from all the active
   *  entries in the EntryGroup.
   *  @return A count of the non-source key features in active entries of this
   *    EntryGroup.
   **/
  int getAllFeaturesCount();

  /**
   *  Add an Entry to this object and then emit the appropriate EntryChange
   *  events.
   **/
  void addElement(final Entry entry);

  /**
   *  A convenience method that does the same as addElement(Entry).
   **/
  void add(final Entry entry);

  /**
   *  Remove an Entry from this object and then emit the appropriate
   *  EntryGroupChange events.  The first entry in the group can only be
   *  removed if it is the only Entry because the first Entry contains the
   *  sequence.
   *  @return true if the removal succeeded, false if it fails(which can
   *    if the given Entry isn't in this EntryGroup or if the user tries to
   *    remove the first Entry).
   **/
  boolean removeElement(final Entry entry);

  /**
   *  A convenience method that does the same as removeElement(Entry).
   **/
  boolean remove(final Entry entry);

  /**
   *  Create a new(blank) Feature in the default Entry of this EntryGroup.
   *  See getDefaultEntry() and Entry.createFeature().
   *  @return The new Feature.
   **/
  Feature createFeature()
      throws ReadOnlyException;

  /**
   *  Create a new(empty) Entry in this EntryGroup.  See Entry.newEntry().
   *  @return The reference of the new Entry.
   **/
  Entry createEntry();

  /**
   *  Create a new(empty) Entry in this EntryGroup.  See Entry.newEntry().
   *  @param name The(file) name of the new Entry.
   *  @return The reference of the new Entry.
   **/
  Entry createEntry(final String name);

  /**
   *  Returns an enumeration of the Feature objects in this EntryGroup. The
   *  returned FeatureEnumeration object will generate all features in this
   *  object in turn. The first item generated is the item at index 0, then
   *  the item at index 1, and so on.
   **/
  FeatureEnumeration features();

  /**
   *  Adds the specified event listener to receive entry group change events
   *  from this object.
   *  @param l the event change listener.
   **/
  void addEntryGroupChangeListener(EntryGroupChangeListener l);

  /**
   *  Removes the specified event listener so that it no longer receives
   *  entry group change events from this object.
   *  @param l the event change listener.
   **/
  void removeEntryGroupChangeListener(EntryGroupChangeListener l);

  /**
   *  Adds the specified event listener to receive entry change events from
   *  this object.
   *  @param l the event change listener.
   **/
  void addEntryChangeListener(EntryChangeListener l);

  /**
   *  Removes the specified event listener so that it no longer receives
   *  entry change events from this object.
   *  @param l the event change listener.
   **/
  void removeEntryChangeListener(EntryChangeListener l);

  /**
   *  Adds the specified event listener to receive feature change events from
   *  this object.
   *  @param l the event change listener.
   **/
  void addFeatureChangeListener(FeatureChangeListener l);

  /**
   *  Removes the specified event listener so that it no longer receives
   *  feature change events from this object.
   *  @param l the event change listener.
   **/
  void removeFeatureChangeListener(FeatureChangeListener l);

  /**
   *  Return the reference of an EntryVector containing the active entries of
   *  this EntryGroup.
   **/
  EntryVector getActiveEntries();

  /**
   *  This method translates the start and end of every each Range in every
   *  Location into another coordinate system.  The Ranges will be truncated
   *  if necessary.
   *  @param constraint This contains the start and end base of the new
   *    coordinate system.  The position given by constraint.getStart() will
   *    be at postion/base 1 in the new coordinate system.
   *  @return a copy of the EntryGroup which has been translated into the new
   *    coordinate system.
   **/
  EntryGroup truncate(final Range constraint);

  /**
   *  Return the ActionController for this EntryGroup(for undo).
   **/
  ActionController getActionController();

  /**
   *  Return true if and only if one or more of the entries or features in
   *  this SimpleEntryGroup are read-only.
   **/
  boolean isReadOnly();

  /**
   *  Increment the reference count for this EntryGroup.
   **/
  void ref();

  /**
   *  Decrement the reference count for this EntryGroup.  When the reference
   *  count goes to zero a EntryGroupChangeEvent is sent to all
   *  EntryGroupChangeListeners with type EntryGroupChangeEvent.DONE_GONE.
   *  The listeners should then stop using the EntryGroup and release any
   *  associated resources.
   **/
  void unref();

  /**
   *  Return the current reference count.
   **/
  int refCount();
}
