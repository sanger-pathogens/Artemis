/* FilteredEntryGroup.java
 *
 * created: Fri Sep  3 1999
 *
 * This file is part of Artemis
 *
 * Copyright(C) 1999  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/FilteredEntryGroup.java,v 1.3 2006-10-26 12:40:30 tjc Exp $
 */

package uk.ac.sanger.artemis;

import uk.ac.sanger.artemis.sequence.*;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.util.ReadOnlyException;
import uk.ac.sanger.artemis.util.OutOfRangeException;

import java.util.NoSuchElementException;

/**
 *  An EntryGroup that filters out some features using a FeaturePredicate
 *  object.
 *
 *  @author Kim Rutherford
 *  @version $Id: FilteredEntryGroup.java,v 1.3 2006-10-26 12:40:30 tjc Exp $
 **/

public class FilteredEntryGroup implements EntryGroup
{
  /**
   *  This is the EntryGroup reference that was passed to the constructor.
   **/
  private EntryGroup entry_group;

  /**
   *  This is the FeaturePredicate reference that was passed to the
   *  constructor.
   **/
  private FeaturePredicate feature_predicate;

  /**
   *  This is a cache variable used by getFilteredFeatures().  It is set to
   *  null by featureChanged(), entryChanged() and entryGroupChanged().
   **/
  private FeatureVector filtered_features = null;

  /**
   *  The name of the filter use to create this object(a passed to the
   *  constructor).
   **/
  private String filter_name;

  /**
   *  Create a new FilteredEntryGroup.
   *  @param entry_group This is the EntryGroup to filter.
   *  @param feature_predicate This is the predicate to use for filtering.
   *  @param filter_name a short description of the filter(used for menus and
   *    title bars).  Can be null if there is no easy way to describe what
   *    the filter does.
   **/
  public FilteredEntryGroup(final EntryGroup entry_group,
                            final FeaturePredicate feature_predicate,
                            final String filter_name) 
  {
    this.entry_group = entry_group;
    this.feature_predicate = feature_predicate;
    this.filter_name = filter_name;

    getEntryGroup().addFeatureChangeListener(new FeatureChangeListener() 
    {
      public void featureChanged(FeatureChangeEvent event) 
      {
        filtered_features = null;
      }
    });

    getEntryGroup().addEntryChangeListener(new EntryChangeListener()
    {
      public void entryChanged(EntryChangeEvent event) 
      {
        switch(event.getType()) 
        {
          case EntryChangeEvent.FEATURE_DELETED:
            // fall through
          case EntryChangeEvent.FEATURE_ADDED:
            filtered_features = null;
            return;
        }
      }
    });

    final EntryGroupChangeListener entry_group_listener =
      new EntryGroupChangeListener() 
    {
        public void entryGroupChanged(EntryGroupChangeEvent event)
        {
          filtered_features = null;
        }
      };
    getEntryGroup().addEntryGroupChangeListener(entry_group_listener);
  }

  /**
   *  Create a new FilteredEntryGroup.
   *  @param entry_group This is the EntryGroup to filter.
   *  @param filtered_features This filtered features.
   *  @param filter_name a short description of the filter(used for menus and
   *    title bars).  Can be null if there is no easy way to describe what
   *    the filter does.
   **/
  public FilteredEntryGroup(final EntryGroup entry_group,
                            final FeatureVector filtered_features,
                            final String filter_name) 
  {
    this(entry_group, (FeaturePredicate)null, filter_name);
    this.filtered_features = filtered_features;
  }
  
  /**
   *  Return the default Entry for this EntryGroup.  The "default" is the
   *  Entry where new features are created.
   **/
  public Entry getDefaultEntry() 
  {
    return getEntryGroup().getDefaultEntry();
  }

  /**
   *  Set the default Entry.  The "default" is the Entry where new features
   *  are created.
   *  @param entry The new default entry.  If this Entry is not active this
   *    method will return immediately.
   **/
  public void setDefaultEntry(final Entry entry) 
  {
    getEntryGroup().setDefaultEntry(entry);
  }

  /**
   *  Return the name of the filter used to create this object(as passed to
   *  the constructor).  Can return null if there is no easy way to describe
   *  what the filter does.
   **/
  public String getFilterName() 
  {
    return filter_name;
  }

  /**
   *  Returns true if and only if there are any unsaved changes in any of the
   *  Entry objects in this EntryGroup.
   **/
  public boolean hasUnsavedChanges() 
  {
    return getEntryGroup().hasUnsavedChanges();
  }

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
  public int indexOf(final Feature feature) 
  {
    return getFilteredFeatures().indexOf(feature);
  }

  /**
   *  Return the index of an Entry within this object.
   *  @return The index of the Entry or -1 if the Entry isn't in this Object.
   **/
  public int indexOf(final Entry entry) 
  {
    return getEntryGroup().indexOf(entry);
  }

  /**
   *  Return true if any of the active entries in the group contains the given
   *  feature.
   **/
  public boolean contains(final Feature feature) 
  {
    if(getPredicate().testPredicate(feature)) 
      return getEntryGroup().contains(feature);
    else 
      return false;
  }

  /**
   *  Return true if the given Entry is active(visible).  The Feature objects
   *  in an Entry that is not active will be ignored by the methods that deal
   *  will features: featureAt(), indexOf(), contains(), features(), etc.
   **/
  public boolean isActive(final Entry entry) 
  {
    return getEntryGroup().isActive(entry);
  }

  /**
   *  Set the "active" setting of the Entry at the given index.  If the index
   *  refers to the default entry and new_active is false, the default entry
   *  will be set to the active entry or null if there are no active entries.
   *  @param index The index of the Entry to change.
   *  @param active The new active setting.
   **/
  public void setIsActive(final int index, final boolean new_active) 
  {
    getEntryGroup().setIsActive(index,new_active);
  }

  /**
   *  Set the "active" setting of the given Entry.  The Entry is the default
   *  entry and new_active is false, the default entry will be set to the
   *  active entry or null if there are no active entries.
   *  @param entry The Entry to activate or deactivate.
   *  @param new_active The new active setting.
   **/
  public void setIsActive(final Entry entry, final boolean new_active) 
  {
    setIsActive(indexOf(entry), new_active);
  }

  /**
   *  Return the Entry from this EntryGroup that contains the sequence to view
   *  or return null if none of the entries contains a sequence.
   **/
  public Entry getSequenceEntry() 
  {
    return getEntryGroup().getSequenceEntry();
  }

  /**
   *  Returns the base length of the sequence of the first Entry in this group
   *  or 0 if this group is empty.
   **/
  public int getSequenceLength() 
  {
    return getEntryGroup().getSequenceLength();
  }

  /**
   *  Returns the Bases object of the first Entry in this group or null if
   *  this group is empty.
   **/
  public Bases getBases() 
  {
    return getEntryGroup().getBases();
  }

  /**
   *  Reverse and complement the sequence and all features in every Entry in
   *  this EntryGroup.
   **/
  public void reverseComplement() throws ReadOnlyException
  {
    getEntryGroup().reverseComplement();
  }

  /**
   *  Return the number of entries in this EntryGroup.
   **/
  public int size() 
  {
    return getEntryGroup().size();
  }

  /**
   *  Return the ith Entry in this EntryGroup.
   **/
  public Entry elementAt(final int i) 
  {
    return getEntryGroup().elementAt(i);
  }

  /**
   *  Return the Feature at the given index.  This method treats all the
   *  features in all the entries as if they were in one big array.  See
   *  the comment on indexOf().
   *  @param index The index of the required Feature.
   *  @return The Feature at the given index.  The first index is 0 the last
   *    is the total number of features in all the entries of this object minus
   *    one.  If the index is out of range then null will be returned.
   **/
  public Feature featureAt(final int index) 
  {
    return getFilteredFeatures().elementAt(index);
  }

  /**
   *  Return a vector containing the references of the Feature objects within
   *  the given range of indices.
   *  @param
   **/
  public FeatureVector getFeaturesInIndexRange(final int start_index,
                                               final int end_index) 
  {
    final FeatureVector features = getFilteredFeatures();

    final FeatureVector return_vector = new FeatureVector();

    for(int i = start_index ; i <= end_index ; ++i) 
      return_vector.add(features.elementAt(i));

    return return_vector;
  }

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
  public FeatureVector getFeaturesInRange(final Range range)
      throws OutOfRangeException 
  {
    return filterFeatures(getEntryGroup().getFeaturesInRange(range));
  }

  /**
   *  Return a vector containing the references of the Feature objects from
   *  all the active entries in the EntryGroup.
   *  @return The non-source key features in active entries of this
   *    EntryGroup.  The returned object is a copy - changes will not effect
   *    the EntryGroup object itself.
   **/
  public FeatureVector getAllFeatures() 
  {
    return getFilteredFeatures();
  }

  /**
   *  Return a count of the number of Feature objects from all the active
   *  entries in the EntryGroup.
   *  @return A count of the non-source key features in active entries of this
   *    EntryGroup.
   **/
  public int getAllFeaturesCount() 
  {
    return getFilteredFeatures().size();
  }

  /**
   *  Add an Entry to this object and then emit the appropriate EntryChange
   *  events.
   **/
  public void addElement(final Entry entry) 
  {
    getEntryGroup().addElement(entry);
  }

  /**
   *  A convenience method that does the same as addElement(Entry).
   **/
  public void add(final Entry entry) 
  {
    getEntryGroup().add(entry);
  }

  /**
   *  Remove an Entry from this object and then emit the appropriate
   *  EntryGroupChange events.  The first entry in the group can only be
   *  removed if it is the only Entry because the first Entry contains the
   *  sequence.
   *  @return true if the removal succeeded, false if it fails(which can
   *    if the given Entry isn't in this EntryGroup or if the user tries to
   *    remove the first Entry).
   **/
  public boolean removeElement(final Entry entry) 
  {
    return getEntryGroup().removeElement(entry);
  }

  /**
   *  A convenience method that does the same as removeElement(Entry).
   **/
  public boolean remove(final Entry entry) 
  {
    return getEntryGroup().remove(entry);
  }

  /**
   *  Create a new(blank) Feature in the default Entry of this EntryGroup.
   *  See getDefaultEntry() and Entry.createFeature().
   *  @return The new Feature.
   **/
  public Feature createFeature() throws ReadOnlyException 
  {
    return getEntryGroup().createFeature();
  }

  /**
   *  Create a new(empty) Entry in this EntryGroup.  See Entry.newEntry().
   *  @return The reference of the new Entry.
   **/
  public Entry createEntry() 
  {
    return getEntryGroup().createEntry();
  }

  /**
   *  Create a new(empty) Entry in this EntryGroup.  See Entry.newEntry().
   *  @param name The(file) name of the new Entry.
   *  @return The reference of the new Entry.
   **/
  public Entry createEntry(final String name) 
  {
    return getEntryGroup().createEntry(name);
  }

  /**
   *  Returns an enumeration of the Feature objects in this EntryGroup. The
   *  returned FeatureEnumeration object will generate all features in this
   *  object in turn. The first item generated is the item at index 0, then
   *  the item at index 1, and so on.
   **/
  public FeatureEnumeration features() 
  {
    return new FeatureEnumerator();
  }

  /**
   *  An Enumeration of Feature objects.
   **/
  public class FeatureEnumerator implements FeatureEnumeration 
  {
    /**
     *  The EntryVector object that we are enumerating.  Set to null when there
     *  are no more Feature objects.
     **/
    private EntryVector active_entries;

    /** Index of the Entry that we will get the next Feature from. */
    private int entry_index = -1;

    /** Enumeration for the current entry. */
    private FeatureEnumeration feature_enumerator;

    /**
     *  The next Feature in the feature_enumerator.  This Feature reference
     *  will always contain a Feature that has pass the FeaturePredicate.
     **/
    private Feature next_feature;

    /**
     *  Create a new FeatureEnumeration that will enumerate the enclosing
     *  FilteredEntryGroup object.  The FilteredEntryGroup object must not be
     *  changed while the enumeration is active.
     **/
    public FeatureEnumerator() 
    {
      active_entries = getEntryGroup().getActiveEntries();

      entry_index = 0;

      if(active_entries.size() > 0) 
        feature_enumerator = active_entries.elementAt(entry_index).features();
      else 
        feature_enumerator = null;
    }

    /**
     *  See the FeatureEnumeration interface for details.
     **/
    public boolean hasMoreFeatures() 
    {
      if(feature_enumerator == null) 
        return false;

      while(true) 
      {
        if(feature_enumerator.hasMoreFeatures()) 
        {
          next_feature = feature_enumerator.nextFeature();

          if(getPredicate().testPredicate(next_feature)) 
            return true;
          else 
            continue;   // loop again
        }

        ++entry_index;

        if(entry_index == active_entries.size()) 
        {
          feature_enumerator = null;
          return false;
        } 
        else 
        {
          feature_enumerator =
            active_entries.elementAt(entry_index).features();
          // loop again
          continue;
        }
      }
    }

    /**
     *  See the FeatureEnumeration interface for details.
     **/
    public Feature nextFeature() throws NoSuchElementException
    {
      if(feature_enumerator == null) 
        throw new NoSuchElementException();

      if(next_feature == null) 
        throw new NoSuchElementException();
      else 
      {
        final Feature return_feature = next_feature;
        next_feature = null;
        return return_feature;
      }
    }

  }

  /**
   *  Implementation of the FeatureChangeListener interface.
   **/
  public void featureChanged(final FeatureChangeEvent event)
  {
    getEntryGroup().featureChanged(event);
  }

  /**
   *  Implementation of the EntryChangeListener interface.  We listen for
   *  changes from every entry in this group and pass the events though to all
   *  the object listening for EntryChangeEvents for the event from this
   *  EntryGroup.
   **/
  public void entryChanged(final EntryChangeEvent event) 
  {
    getEntryGroup().entryChanged(event);
  }

  /**
   *  Implementation of the EntryGroupChangeListener interface.  We listen for
   *  changed to the EntryGroup that this FilteredEntryGroup is filtering so
   *  that the filtered_features cache variable can be maintained.
   **/
  public void entryGroupChanged(final EntryGroupChangeEvent event) 
  {
    filtered_features = null;
  }

  /**
   *  Adds the specified event listener to receive entry group change events
   *  from this object.
   *  @param l the event change listener.
   **/
  public void addEntryGroupChangeListener(EntryGroupChangeListener l)
  {
    getEntryGroup().addEntryGroupChangeListener(l);
  }

  /**
   *  Removes the specified event listener so that it no longer receives
   *  entry group change events from this object.
   *  @param l the event change listener.
   **/
  public void removeEntryGroupChangeListener(EntryGroupChangeListener l) 
  {
    getEntryGroup().removeEntryGroupChangeListener(l);
  }

  /**
   *  Adds the specified event listener to receive entry change events from
   *  this object.
   *  @param l the event change listener.
   **/
  public void addEntryChangeListener(EntryChangeListener l) 
  {
    getEntryGroup().addEntryChangeListener(l);
  }

  /**
   *  Removes the specified event listener so that it no longer receives
   *  entry change events from this object.
   *  @param l the event change listener.
   **/
  public void removeEntryChangeListener(EntryChangeListener l) 
  {
    getEntryGroup().removeEntryChangeListener(l);
  }

  /**
   *  Adds the specified event listener to receive feature change events from
   *  this object.
   *  @param l the event change listener.
   **/
  public void addFeatureChangeListener(FeatureChangeListener l) 
  {
    getEntryGroup().addFeatureChangeListener(l);
  }

  /**
   *  Removes the specified event listener so that it no longer receives
   *  feature change events from this object.
   *  @param l the event change listener.
   **/
  public void removeFeatureChangeListener(FeatureChangeListener l) 
  {
    getEntryGroup().removeFeatureChangeListener(l);
  }

  /**
   *  Return the reference of an EntryVector containing the active entries of
   *  this EntryGroup.
   **/
  public EntryVector getActiveEntries() 
  {
    return getEntryGroup().getActiveEntries();
  }

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
  public EntryGroup truncate(final Range constraint) 
  {
    return getEntryGroup().truncate(constraint);
  }

  /**
   *  Return the ActionController for this EntryGroup(for undo).
   **/
  public ActionController getActionController() 
  {
    return getEntryGroup().getActionController();
  }

  /**
   *  Return true if and only if one or more of the entries or features in
   *  this SimpleEntryGroup are read-only.
   **/
  public boolean isReadOnly() 
  {
    return getEntryGroup().isReadOnly();
  }

  /**
   *  Increment the reference count for this EntryGroup.
   **/
  public void ref() 
  {
    getEntryGroup().ref();
  }

  /**
   *  Decrement the reference count for this EntryGroup.  When the reference
   *  count goes to zero a EntryGroupChangeEvent is sent to all
   *  EntryGroupChangeListeners with type EntryGroupChangeEvent.DONE_GONE.
   *  The listeners should then stop using the EntryGroup and release any
   *  associated resources.
   **/
  public void unref() 
  {
    getEntryGroup().unref();
  }

  /**
   *  Return the current reference count.
   **/
  public int refCount() 
  {
    return getEntryGroup().refCount();
  }

  /**
   *  Return those features in the given FeatureVector that pass
   *  feature_predicate.
   **/
  protected FeatureVector filterFeatures() 
  {
    final FeatureEnumeration test_enumerator = entry_group.features();
    final FeatureVector return_features = new FeatureVector();

    while(test_enumerator.hasMoreFeatures())
    {
      final Feature this_feature = test_enumerator.nextFeature();

      if (feature_predicate.testPredicate (this_feature))
            return_features.add (this_feature);
    }

    return return_features;
  }

  /**
   *  Return those features in the given FeatureVector that pass
   *  feature_predicate.
   **/
  protected FeatureVector filterFeatures(final FeatureVector features)
  {
    final FeatureVector return_features = new FeatureVector();

    for(int i = 0 ; i < features.size() ; ++i) 
    {
      final Feature this_feature = features.elementAt(i);

      if(getPredicate().testPredicate(this_feature)) 
        return_features.add(this_feature);
    }

    return return_features;
  }

  /**
   *  Return the features in the EntryGroup that pass the feature_predicate.
   **/
  private FeatureVector getFilteredFeatures() 
  {
    if(filtered_features == null) 
      filtered_features = filterFeatures();
    
    return filtered_features;
  }

  /**
   *  Return the EntryGroup reference that was passed to the constructor.
   **/
  private EntryGroup getEntryGroup() 
  {
    return entry_group;
  }

  /**
   *  Return the FeaturePredicate reference that was passed to the
   *  constructor.
   **/
  private FeaturePredicate getPredicate() 
  {
    return feature_predicate;
  }

}
