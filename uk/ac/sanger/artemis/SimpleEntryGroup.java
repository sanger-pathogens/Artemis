/* SimpleEntryGroup.java
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/SimpleEntryGroup.java,v 1.8 2008-06-11 15:15:23 tjc Exp $
 **/

package uk.ac.sanger.artemis;

import uk.ac.sanger.artemis.sequence.*;
import uk.ac.sanger.artemis.components.MessageDialog;
import uk.ac.sanger.artemis.io.GFFDocumentEntry;
import uk.ac.sanger.artemis.io.IndexedGFFDocumentEntry;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.io.StreamSequence;
import uk.ac.sanger.artemis.io.SimpleDocumentEntry;
import uk.ac.sanger.artemis.io.DatabaseDocumentEntry;
import uk.ac.sanger.artemis.io.DocumentEntry;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.ReadOnlyException;
import uk.ac.sanger.artemis.util.OutOfRangeException;

import java.util.Vector;
import java.util.NoSuchElementException;

/**
 *  This class implements a vector of Entry objects, with additional methods
 *  for querying and changing the feature tables of all the entries at
 *  once.  Objects of this class act a bit like single Entry objects.
 *
 *  @author Kim Rutherford
 *  @version $Id: SimpleEntryGroup.java,v 1.8 2008-06-11 15:15:23 tjc Exp $
 **/

public class SimpleEntryGroup extends EntryVector
                              implements EntryGroup
{

  /** vector of those objects listening for entry change events. */
  final private Vector entry_group_listener_list = new Vector();

  /** vector of those objects listening for entry change events. */
  final private Vector entry_listener_list = new Vector();

  /** vector of those objects listening for feature change events. */
  final private Vector feature_listener_list = new Vector();

  /** vector of entries that are currently active (visible). */
  final private EntryVector active_entries = new EntryVector();

  /**
   *  The default Entry for this SimpleEntryGroup.  The "default" is the Entry
   *  where new features are created.
   **/
  private Entry default_entry = null;

  /** Bases object that was passed to the constructor. */
  private Bases bases;

  /** Incremented by ref(), decremented by unref(). */
  private int reference_count = 0;

  /** The ActionController of this EntryGroup (used for undo). */
  final private ActionController action_controller = new ActionController();

  /** 
   *  Create a new empty SimpleEntryGroup object.
   **/
  public SimpleEntryGroup(final Bases bases) 
  {
    this.bases = bases;

    addFeatureChangeListener(getActionController());
    addEntryChangeListener(getActionController());
    getBases().addSequenceChangeListener(getActionController(),
                                         Bases.MIN_PRIORITY);
    addEntryGroupChangeListener(getActionController());
  }

  public SimpleEntryGroup() 
  {
    addFeatureChangeListener(getActionController());
    addEntryGroupChangeListener(getActionController());
  }
  
  /**
   *  Returns true if and only if there are any unsaved changes in any of the
   *  Entry objects in this EntryGroup.
   **/
  public boolean hasUnsavedChanges() 
  {
    final int my_size = size();
    for(int entry_index = 0; entry_index < my_size;
        ++entry_index) 
    {
      if(elementAt(entry_index).hasUnsavedChanges()) 
        return true;
    }

    return false;
  }

  /**
   *  Return the default Entry for this SimpleEntryGroup.  The "default" is the
   *  Entry where new features are created.
   **/
  public Entry getDefaultEntry() 
  {
    return default_entry;
  }

  /**
   *  Set the default Entry.  The "default" is the Entry where new features
   *  are created.
   *  @param entry The new default entry.  If this Entry is not active this
   *    method will return immediately.
   **/
  public void setDefaultEntry(Entry entry) 
  {
    if(entry != null && !isActive(entry))
      return;

    // do nothing
    if(default_entry == entry) 
      return;
    else 
      default_entry = entry;

    // now inform the listeners that a change has occured
    final EntryGroupChangeEvent event =
      new EntryGroupChangeEvent(this, getDefaultEntry(),
                                EntryGroupChangeEvent.NEW_DEFAULT_ENTRY);

    fireEvent(entry_group_listener_list, event);
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
  public int indexOf(Feature feature) 
  {
    int feature_count_of_previous_entries = 0;
    final int active_entries_size = active_entries.size();

    for(int entry_index = 0; entry_index < active_entries_size;
        ++entry_index) 
    {
      final Entry this_entry  = active_entries.elementAt(entry_index);
      final int feature_index = this_entry.indexOf(feature);

      if(feature_index != -1) 
        return feature_index + feature_count_of_previous_entries;

      feature_count_of_previous_entries += this_entry.getFeatureCount();
    }

    return -1;
  }

  /**
   *  Return true if any of the active entries in the group contains the given
   *  feature.
   **/
  public boolean contains(Feature feature) 
  {
    final int active_entries_size = active_entries.size();

    for(int i = 0; i < active_entries_size; ++i) 
    {
      final Entry current_entry = active_entries.elementAt(i);

      if(current_entry.contains(feature)) 
        return true;
    }

    return false;
  }

  /**
   *  Return true if the given Entry is active(visible).  The Feature objects
   *  in an Entry that is not active will be ignored by the methods that deal
   *  will features: featureAt(), indexOf(), contains(), features(), etc.
   **/
  public boolean isActive(Entry entry) 
  {
    if(active_entries.contains(entry)) 
      return true;
    else 
      return false;
  }

  /**
   *  Set the "active" setting of the Entry at the given index.  If the index
   *  refers to the default entry and new_active is false, the default entry
   *  will be set to the active entry or null if there are no active entries.
   *  @param index The index of the Entry to change.
   *  @param active The new active setting.
   **/
  public void setIsActive(int index, boolean new_active) 
  {
    final Entry entry = elementAt(index);

    if(new_active)
    {
      // no change
      if(isActive(entry)) 
        return;
      else 
      {
        // this is slow but it guarantees that the Entry references in the
        // active_entries vector are in the same order as in the
        // SimpleEntryGroup

        final EntryVector new_active_entries = new EntryVector();
        final int my_size = size();

        for(int i = 0; i < my_size; ++i) 
        {
          if(active_entries.contains(elementAt(i)) || index == i) 
            new_active_entries.add(elementAt(i));
        }

        active_entries.removeAllElements();

        final int new_active_entries_size = new_active_entries.size();

        for(int i = 0; i < new_active_entries_size; ++i) 
          active_entries.add(new_active_entries.elementAt(i));

        if(active_entries.size() >= 1 && getDefaultEntry() == null) 
        {
          // there was no default entry before calling addElement() so
          // make the first non-sequence entry the default entry
          if(active_entries.elementAt(0) == getSequenceEntry() &&
             active_entries.size() == 1)
          {
            // don't set the default entry to be the sequence entry unless
            // the user asks for it
            setDefaultEntry(null);
          }
          else 
          {
            if(active_entries.size() == 1) 
              setDefaultEntry(active_entries.elementAt(0));
            else 
              setDefaultEntry(active_entries.elementAt(1));
          }
        }
      }
    }
    else 
    {
      // no change
      if(!isActive(entry)) 
        return;
      else
      {
        active_entries.removeElement(entry);

        if(entry == getDefaultEntry()) 
        {
          if(active_entries.size() > 0) 
          {
            if(active_entries.elementAt(0) == getSequenceEntry()) 
            {
              // don't set the default entry to be the sequence entry unless
              // the user asks for it
              if(active_entries.size() > 1)
                setDefaultEntry(active_entries.elementAt(1));
              else
                setDefaultEntry(null);
            }
            else 
              setDefaultEntry(active_entries.elementAt(0));
          }
          else 
            setDefaultEntry(null);
        }
      }
    }

    // now inform the listeners that a change has occured
    final EntryGroupChangeEvent event;

    // change state
    if(new_active)        // become active
      event = new EntryGroupChangeEvent(this, entry,
                                        EntryGroupChangeEvent.ENTRY_ACTIVE);
    else                  // become inactive
      event = new EntryGroupChangeEvent(this, entry,
                                        EntryGroupChangeEvent.ENTRY_INACTIVE);

    fireEvent(entry_group_listener_list, event);
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
   *  Return the Entry from this SimpleEntryGroup that contains the sequence
   *  to view or return null if none of the entries contains a sequence.
   **/
  public Entry getSequenceEntry() 
  {
    if(size() == 0) 
      return null;
    else 
      return elementAt(0);
  }

  /**
   *  Returns the base length of the sequence of the first Entry in this group
   *  or 0 if this group is empty.
   **/
  public int getSequenceLength() 
  {
    return getBases().getLength();
  }

  /**
   *  Returns the Bases object of the first Entry in this group or null if
   *  this group is empty.
   **/
  public Bases getBases() 
  {
    return bases;
  }

  /**
   *  Reverse and complement the sequence and all features in every Entry in
   *  this SimpleEntryGroup.
   **/
  public void reverseComplement()
      throws ReadOnlyException 
  {
    if(isReadOnly()) 
      throw new ReadOnlyException();

    // reverse the sequence
    getBases().reverseComplement();
  }

  /**
   *  Return true if and only if one or more of the entries or features in
   *  this SimpleEntryGroup are read-only.
   **/
  public boolean isReadOnly() 
  {
    final int my_size = size();
    for(int i = 0; i < my_size; ++i) 
    {
      final Entry this_entry = elementAt(i);

      if(this_entry.isReadOnly()) 
        return true;

      final FeatureEnumeration feature_enum = this_entry.features();

      while(feature_enum.hasMoreFeatures()) 
      {
        if(feature_enum.nextFeature().isReadOnly()) 
          return true;
      }
    }

    return false;
  }

  /**
   *  Increment the reference count for this EntryGroup.
   **/
  public void ref() 
  {
    ++reference_count;
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
    --reference_count;

    if(reference_count == 0) 
    {
      // remove all the entries which will close any edit or view windows
      while(size() > 0) 
      {
        final Entry this_entry = elementAt(0);
        remove(this_entry);
        
        this_entry.getEMBLEntry().dispose();
      }

      // now inform the listeners that the EntryGroup is no more
      final EntryGroupChangeEvent event =
              new EntryGroupChangeEvent(this, null,
                                        EntryGroupChangeEvent.DONE_GONE);

      fireEvent(entry_group_listener_list, event);
    }
  }

  /**
   *  Return the current reference count.
   **/
  public int refCount() 
  {
    return reference_count;
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
  public Feature featureAt(int index) 
  {
    if(index < 0) 
      throw new Error("internal error - index out of range: " + index);

    final int active_entries_size = active_entries.size();

    for(int entry_index = 0; entry_index < active_entries_size; 
        ++entry_index)
    {
      final Entry this_entry = active_entries.elementAt(entry_index);

      if(index < this_entry.getFeatureCount()) 
        return this_entry.getFeature(index);

      index -= this_entry.getFeatureCount();
    }

    throw new Error("internal error - index out of range: " + index);
  }

  /**
   *  Return a vector containing the references of the Feature objects within
   *  the given range of indices.
   *  @param start_index The index of the first feature to return.
   *  @param end_index The index of the last feature to return.
   **/
  public FeatureVector getFeaturesInIndexRange(final int start_index,
                                               final int end_index) 
  {
    final FeatureVector return_vector = new FeatureVector();

    for(int i = start_index; i <= end_index; ++i) 
      return_vector.add(featureAt(i));

    return return_vector;
  }

  /**
   *  Return a vector containing the references of the Feature objects within
   *  the given range for all the active entries in the SimpleEntryGroup.
   *  @param range Return features that overlap this range - ie the start of
   *    the feature is less than or equal to the end of the range and the end
   *    of the feature is greater than or equal to the start of the range.
   *  @return The non-source key features of this feature table the are within
   *    the given range.  The returned object is a copy - changes will not
   *    effect the FeatureTable object itself.
   **/
  public FeatureVector getFeaturesInRange(Range range)
      throws OutOfRangeException 
  {
    final FeatureVector return_vector = new FeatureVector();
    final int my_size = size();

    for(int i = 0; i < my_size; ++i) 
    {
      final Entry this_entry = elementAt(i);

      if(isActive(this_entry)) 
      {
        final FeatureVector visible_entry_features =
          elementAt(i).getFeaturesInRange(range);

        final int visible_entry_features_size = visible_entry_features.size();

        for(int feature_index = 0; feature_index < visible_entry_features_size;
            ++feature_index) 
        {
          final Feature this_feature =
            visible_entry_features.elementAt(feature_index);
          return_vector.add(this_feature);
        }
      }
    }

    return return_vector;
  }

  /**
   *  Return a vector containing the references of the Feature objects from
   *  all the active entries in the SimpleEntryGroup.
   *  @return The non-source key features in active entries of this
   *    SimpleEntryGroup.  The returned object is a copy - changes will not
   *    effect the SimpleEntryGroup object itself.
   **/
  public FeatureVector getAllFeatures() 
  {
    final FeatureVector return_vector = new FeatureVector();
    final int my_size = size();

    for(int i = 0; i < my_size; ++i) 
    {
      final Entry this_entry = elementAt(i);

      if(isActive(this_entry)) 
      {
        final FeatureVector entry_features = elementAt(i).getAllFeatures();
        final int entry_features_size = entry_features.size();

        for(int feature_index = 0; feature_index < entry_features_size;
            ++feature_index) 
        {
          final Feature this_feature =
            entry_features.elementAt(feature_index);
          return_vector.add(this_feature);
        }
      }
    }

    return return_vector;
  }

  /**
   *  Return a count of the number of Feature objects from all the active
   *  entries in the SimpleEntryGroup.
   *  @return A count of the non-source key features in active entries of this
   *    SimpleEntryGroup.
   **/
  public int getAllFeaturesCount() 
  {
    int return_count = 0;
    final int my_size = size();

    for(int i = 0; i < my_size; ++i) 
    {
      final Entry this_entry = elementAt(i);

      if(isActive(this_entry)) 
        return_count += this_entry.getFeatureCount();
    }

    return return_count;
  }

  /**
   *  Add an Entry to this object and then emit the appropriate EntryChange
   *  events.
   **/
  public void addElement(Entry entry) 
  {
    super.addElement(entry);

    // set the default Entry to whichever Entry gets added first
    if(default_entry == null) 
      default_entry = entry;

    active_entries.add(entry);

    // now inform the listeners that an addition has occured
    final EntryGroupChangeEvent event =
      new EntryGroupChangeEvent(this, entry,
                                EntryGroupChangeEvent.ENTRY_ADDED);

    fireEvent(entry_group_listener_list, event);

    // make the new entry the default entry if and only if there was no entry
    // previously or there was only one entry previously and it contained
    // sequence but no features
    if(size() == 1) 
      setDefaultEntry(entry);
    else 
    {
      if(size() == 2) 
      {
        final Entry first_entry = elementAt(0);

        if(first_entry.getFeatureCount() == 0) 
        {
          final Bases first_entry_bases = first_entry.getBases();

          if(first_entry_bases != null &&
             first_entry_bases.getLength() > 0) 
            setDefaultEntry(entry);
        }
      }
    }

    entry.addEntryChangeListener(this);
    entry.addFeatureChangeListener(this);
  }

  /**
   *  A convenience method that does the same as addElement(Entry).
   **/
  public void add(final Entry entry) 
  {
    if(entry.getEMBLEntry() instanceof IndexedGFFDocumentEntry)
      ((IndexedGFFDocumentEntry)entry.getEMBLEntry()).setEntryGroup(this);
    else if(entry.getEMBLEntry() instanceof GFFDocumentEntry)
    {
      ((GFFDocumentEntry)entry.getEMBLEntry()).adjustCoordinates( getSequenceEntry() );
      if(!Options.isBlackBeltMode() && size() > 1 && 
          entry.getEMBLEntry().getSequence() != null )
      {
        new MessageDialog (null, "Warning", 
            "Overlaying a GFF with a sequence onto an entry with a sequence.",
            false);
      }
    }

    addElement(entry);
  }

  /**
   *  Remove an Entry from this object and then emit the appropriate
   *  EntryGroupChange events.  The first entry in the group can only be
   *  removed if it is the only Entry because the first Entry contains the
   *  sequence.
   *  @return true if the removal succeeded, false if it fails(which can if
   *    the given Entry isn't in this SimpleEntryGroup or if the user tries to
   *    remove the first Entry).
   **/
  public boolean removeElement(final Entry entry) 
  {
    // this call will sort out the default entry
    setIsActive(indexOf(entry), false);

    entry.dispose();

    final boolean remove_return = super.removeElement(entry);

    entry.removeEntryChangeListener(this);
    entry.removeFeatureChangeListener(this);

    active_entries.removeElement(entry);

    // now inform the listeners that a deletion has occured
    final EntryGroupChangeEvent event =
      new EntryGroupChangeEvent(this, entry,
                                EntryGroupChangeEvent.ENTRY_DELETED);

    fireEvent(entry_group_listener_list, event);

    return remove_return;
  }

  /**
   *  A convenience method that does the same as removeElement(Entry).
   **/
  public boolean remove(final Entry entry) 
  {
    return removeElement(entry);
  }

  /**
   *  Create a new(blank) Feature in the default Entry of this
   *  SimpleEntryGroup.  See getDefaultEntry() and Entry.createFeature().
   *  @return The new Feature.
   **/
  public Feature createFeature() throws ReadOnlyException
  {
    final Feature new_feature = getDefaultEntry().createFeature();
    return new_feature;
  }

  /**
   *  Create a new(empty) Entry in this SimpleEntryGroup.  See
   *  Entry.newEntry().
   *  @return The reference of the new Entry.
   **/
  public Entry createEntry() 
  {
    Entry new_entry = null;
    Entry default_entry = getDefaultEntry();
    if(default_entry != null &&
       default_entry.getEMBLEntry() != null &&
       default_entry.getEMBLEntry() instanceof DatabaseDocumentEntry)
    {
      DatabaseDocument doc =
        (DatabaseDocument)((DocumentEntry)default_entry.getEMBLEntry()).getDocument();
      DatabaseDocument new_doc = doc.createDatabaseDocument();
      
      try
      {
        DatabaseDocumentEntry new_doc_entry = 
          new DatabaseDocumentEntry();
        new_doc_entry.setDocument(new_doc);
        new_entry = new Entry(getBases(), new_doc_entry);
      }
      catch(Exception e)
      {
        e.printStackTrace();
      }
      
    }
    else
      new_entry = Entry.newEntry(getBases());
    
    add(new_entry);
    return new_entry;
  }

  /**
   *  Create a new(empty) Entry in this SimpleEntryGroup.  See
   *  Entry.newEntry().
   *  @param name The(file) name of the new Entry.
   *  @return The reference of the new Entry.
   **/
  public Entry createEntry(final String name) 
  {
    final Entry new_entry = createEntry();
    new_entry.setName(name);
    return new_entry;
  }

  /**
   *  Returns an enumeration of the Feature objects in this
   *  SimpleEntryGroup. The returned FeatureEnumeration object will generate
   *  all features in this object in turn. The first item generated is the
   *  item at index 0, then the item at index 1, and so on.
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

    /**
     *  The index of the Entry that we will get the next Feature from.
     **/
    private int entry_index = -1;

    /** Enumeration for the current entry */
    private FeatureEnumeration feature_enumerator;

    /**
     *  Create a new FeatureEnumeration that will enumerate the enclosing
     *  SimpleEntryGroup object.  The SimpleEntryGroup object must not be
     *  changed while the enumeration is active.
     **/
    public FeatureEnumerator() 
    {
      active_entries = getActiveEntries();

      entry_index = 0;

      if(active_entries.size() > 0) 
        feature_enumerator =
          active_entries.elementAt(entry_index).features();
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

      if(feature_enumerator.hasMoreFeatures()) 
        return true;

      ++entry_index;

      if(entry_index == active_entries.size()) 
        return false;
      else
      {
        feature_enumerator =
          active_entries.elementAt(entry_index).features();
        return hasMoreFeatures();
      }
    }

    /**
     *  See the FeatureEnumeration interface for details.
     **/
    public Feature nextFeature()
        throws NoSuchElementException 
    {
      if(feature_enumerator == null) 
        throw new NoSuchElementException();

      return feature_enumerator.nextFeature();
    }

  }

  /**
   *  Implementation of the FeatureChangeListener interface.  We listen for
   *  changes in every feature of every entry in this group.
   **/
  public void featureChanged(FeatureChangeEvent event) 
  {
    // pass the action straight through
    fireEvent(feature_listener_list, event);
  }

  /**
   *  Implementation of the EntryChangeListener interface.  We listen for
   *  changes from every entry in this group and pass the events though to all
   *  the object listening for EntryChangeEvents for the event from this
   *  SimpleEntryGroup.
   **/
  public void entryChanged(EntryChangeEvent event) 
  {
    // pass the action straight through
    fireEvent(entry_listener_list, event);
  }

  /**
   *  Send an event to those object listening for it.
   *  @param listeners A Vector of the objects that the event should be sent
   *    to.
   *  @param event The event to send
   **/
  private void fireEvent(Vector listeners, ChangeEvent event) 
  {
    final Vector targets;
    // copied from a book - synchronising the whole method might cause a
    // deadlock
    synchronized(this) 
    {
      targets = (Vector)listeners.clone();
    }

    //boolean seen_chado_manager = false;
    final int targets_size = targets.size();
    for(int i = 0; i < targets_size; ++i) 
    {
      ChangeListener target =(ChangeListener) targets.elementAt(i);

      if(event instanceof EntryGroupChangeEvent)
      {
        final EntryGroupChangeListener entry_group_change_listener =
         (EntryGroupChangeListener) target;
        final EntryGroupChangeEvent group_change_event =
         (EntryGroupChangeEvent) event;
        entry_group_change_listener.entryGroupChanged(group_change_event);
      } 
      else
      {
        if(event instanceof EntryChangeEvent) 
        {
          final EntryChangeListener entry_change_listener =
                                     (EntryChangeListener) target;

//        if(entry_change_listener instanceof ChadoTransactionManager)
//        {
            // just call this listener once
//          if(!seen_chado_manager)
//          {
//            entry_change_listener.entryChanged((EntryChangeEvent) event);
//            seen_chado_manager = true;
//          }
//        }
//        else  
            entry_change_listener.entryChanged((EntryChangeEvent) event);
        } 
        else 
        {
          final FeatureChangeListener feature_change_listener =
                                    (FeatureChangeListener) target;
          feature_change_listener.featureChanged((FeatureChangeEvent) event);
        }

      }
    }
  }

  /**
   *  Adds the specified event listener to receive entry group change events
   *  from this object.
   *  @param l the event change listener.
   **/
  public void addEntryGroupChangeListener(EntryGroupChangeListener l) 
  {
    entry_group_listener_list.addElement(l);
  }

  /**
   *  Removes the specified event listener so that it no longer receives
   *  entry group change events from this object.
   *  @param l the event change listener.
   **/
  public void removeEntryGroupChangeListener(EntryGroupChangeListener l) 
  {
    entry_group_listener_list.removeElement(l);
  }

  /**
   *  Adds the specified event listener to receive entry change events from
   *  this object.
   *  @param l the event change listener.
   **/
  public void addEntryChangeListener(EntryChangeListener l) 
  {
    entry_listener_list.addElement(l);
  }

  /**
   *  Removes the specified event listener so that it no longer receives
   *  entry change events from this object.
   *  @param l the event change listener.
   **/
  public void removeEntryChangeListener(EntryChangeListener l) 
  {
    entry_listener_list.removeElement(l);
  }

  /**
   *  Adds the specified event listener to receive feature change events from
   *  this object.
   *  @param l the event change listener.
   **/
  public void addFeatureChangeListener(FeatureChangeListener l) 
  {
    feature_listener_list.addElement(l);
  }

  /**
   *  Removes the specified event listener so that it no longer receives
   *  feature change events from this object.
   *  @param l the event change listener.
   **/
  public void removeFeatureChangeListener(FeatureChangeListener l) 
  {
    feature_listener_list.removeElement(l);
  }

  /**
   *  Return the reference of an EntryVector containing the active entries of
   *  this EntryGroup.
   **/
  public EntryVector getActiveEntries() 
  {
    return(EntryVector) active_entries.clone();
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
    final Bases new_bases = getBases().truncate(constraint);

    final EntryGroup new_entry_group = new SimpleEntryGroup(new_bases);
    final int my_size = size();

    for(int i = 0; i < my_size; ++i) 
    {
      final Entry this_entry = elementAt(i);
      final Entry new_entry = this_entry.truncate(new_bases, constraint);
      new_entry_group.add(new_entry);
    }

    if(size() > 0) 
    {
      final StreamSequence sequence =
       (StreamSequence) new_bases.getSequence();
      final SimpleDocumentEntry document_entry =
       (SimpleDocumentEntry) new_entry_group.elementAt(0).getEMBLEntry();
      document_entry.setSequence(sequence);
    }

    return new_entry_group;
  }

  /**
   *  Return the ActionController for this EntryGroup(for undo).
   **/
  public ActionController getActionController() 
  {
    return action_controller;
  }
  
}
