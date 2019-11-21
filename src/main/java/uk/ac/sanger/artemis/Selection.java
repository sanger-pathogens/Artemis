/* Selection.java
 *
 * created: Fri Nov  6 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/Selection.java,v 1.3 2008-06-26 09:36:50 tjc Exp $
 */

package uk.ac.sanger.artemis;

import uk.ac.sanger.artemis.sequence.*;

import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.io.RangeVector;

import java.io.*;
import java.util.Vector;
import java.awt.datatransfer.*;

/**
 *  Objects of this class hold the features/feature segments or entries that
 *  the user has selected.
 *
 *  @author Kim Rutherford
 *  @version $Id: Selection.java,v 1.3 2008-06-26 09:36:50 tjc Exp $
 *
 **/

public class Selection
    implements FeatureChangeListener, EntryChangeListener,
               Transferable, ClipboardOwner {
  /**
   *  Create a new Selection object.
   *  @param clipboard The system clipboard reference.  This may be null if
   *    this Selection should not set the system clipboard.
   **/
  public Selection (final Clipboard clipboard) {
    this.clipboard = clipboard;
  }

  private final DataFlavor flavors [] = {
    DataFlavor.stringFlavor,
    // deprecated at 1.3
    //   DataFlavor.plainTextFlavor
  };

  /**
   *  Returns the array of the data flavors that this object can provide.
   **/
  public synchronized DataFlavor[] getTransferDataFlavors () {
    return flavors;
  }

  /**
   *  Returns whether the requested flavor is supported by this object.
   *  @param flavor the requested flavor for the data
   **/
  public boolean isDataFlavorSupported (final DataFlavor flavor) {
    return
      flavor.equals (DataFlavor.stringFlavor)
// deprecated at 1.3
//      || flavor.equals (DataFlavor.plainTextFlavor)
      ;
  }

  /**
   *  If the data was requested in the "java.lang.String" flavor, return the
   *  String representing the selection, else throw an
   *  UnsupportedFlavorException.
   *  @param flavor the requested flavor for the data
   **/
  public synchronized Object getTransferData (final DataFlavor flavor)
      throws UnsupportedFlavorException, IOException {
    if (flavor.equals (DataFlavor.stringFlavor)) {
      return getSelectionText ();
// deprecated at 1.3
//    } else if (flavor.equals (DataFlavor.plainTextFlavor)) {
//      return new StringReader (getSelectionText ());
    } else {
      throw new UnsupportedFlavorException (flavor);
    }
  }

  /**
   *  Implementation of the ClipboardOwner interface.
   **/
  public void lostOwnership(final Clipboard clipboard,
                            final Transferable contents) {

  }

  /**
   *  Returns this Selection in a String.
   **/
  public String getSelectionText () {
    final FeatureVector selection_features = getAllFeatures ();

    final StringBuffer buffer = new StringBuffer ();

    final MarkerRange marker_range = getMarkerRange ();

    if (marker_range != null) {
      final String selection_bases = Strand.markerRangeBases (marker_range);

      buffer.append (selection_bases);
    }

    for (int i = 0 ; i < selection_features.size () && i < 50 ; ++i) {
      buffer.append (selection_features.elementAt (i).toString ());
    }

    return buffer.toString ();
  }

  /**
   *  Implementation of the FeatureChangeListener interface.  We need to
   *  listen to feature change events so a SelectionChangeEvent can be sent if
   *  a selected feature is changed.
   *  @param event The change event.
   **/
  public void featureChanged (FeatureChangeEvent event) {
    if (segments == null && features.contains (event.getFeature ()) ||
        getAllFeatures ().contains (event.getFeature ())) {
      if (event.getType () == FeatureChangeEvent.QUALIFIER_CHANGED ||
          event.getType () == FeatureChangeEvent.KEY_CHANGED) {
        // no need to reset the cache in this case
        /*final SelectionChangeEvent selection_change_event =
          new SelectionChangeEvent (this, SelectionChangeEvent.OBJECT_CHANGED);

        fireAction (selection_listener_list, selection_change_event);*/
      } else {
        changeSelection (SelectionChangeEvent.OBJECT_CHANGED);
      }
    }
  }

  /**
   *  Implementation of the EntryChangeListener interface.  We listen to
   *  EntryChange events so that we remove feature from the selection when
   *  they are removed from the Entry.
   **/
  public void entryChanged (final EntryChangeEvent event) {
    switch (event.getType ()) {
    case EntryChangeEvent.FEATURE_DELETED:
      
      features.remove (event.getFeature ());
      
      if (contains (event.getFeature ())) {
        // we have a segment of the feature in the Selection -
        // this will fire a SelectionChangeEvent
        removeSegmentsOf (event.getFeature ());
      } else {
        changeSelection (SelectionChangeEvent.SELECTION_CHANGED);
      }
      
      break;
    }
  }

  /**
   *  Adds the specified event listener to receive selection change events
   *  from this object.
   *  @param l the event change listener.
   **/
  public void addSelectionChangeListener (final SelectionChangeListener l) {
    selection_listener_list.addElement (l);
  }

  /**
   *  Removes the specified event listener so that it no longer receives
   *  selection change events from this object.
   *  @param l the event change listener.
   **/
  public void removeSelectionChangeListener (final SelectionChangeListener l) {
    selection_listener_list.removeElement (l);
  }

  /**
   *  Broadcast an event to notify SelectionChange event listeners that the
   *  selection has changed.
   *  @param type The type to use for the new event (see SelectionChangeEvent
   *    for details).
   **/
  private void changeSelection (final int type) {
    resetCache ();

    final SelectionChangeEvent event =
      new SelectionChangeEvent (this, type);

    fireAction (selection_listener_list, event);

    if (clipboard != null) {
      clipboard.setContents (this, this);
    }
  }

  /**
   *  Add a Feature to the selection if it isn't there already.
   **/
  public void add (final Feature feature) {
    if (addWithoutEvent (feature)) {
      changeSelection (SelectionChangeEvent.SELECTION_CHANGED);
    }
  }

  /**
   *  Add the given Feature objects to the selection, and then send an
   *  appropriate event to the SelectionChangeEvent listeners.
   **/
  public void add (final FeatureVector features) {
    for (int i = 0 ; i < features.size () ; ++i) {
      addWithoutEvent (features.elementAt (i));
    }
    changeSelection (SelectionChangeEvent.SELECTION_CHANGED);
  }

  /**
   *  Add a Feature to the selection if it isn't there already, but don't send
   *  an event to the SelectionChangeEvent listeners.
   *  @return true if and only if the segment with not already in the
   *    selection.
   **/
  private boolean addWithoutEvent (final Feature feature) {

    if (feature != null && features.contains (feature)) {
      return false;
      // do nothing
    } else {
      features.add (feature);
      return true;
    }
  }

  /**
   *  Add a FeatureSegment to the selection if it isn't there already, and
   *  send an event to the SelectionChangeEvent listeners..
   **/
  public void add (final FeatureSegment segment) {
    if (addWithoutEvent (segment)) {
      changeSelection (SelectionChangeEvent.SELECTION_CHANGED);
    }
  }

  /**
   *  Add a FeatureSegment to the selection if it isn't there already, but
   *  don't send an event to the SelectionChangeEvent listeners.
   *  @return true if and only if the segment with not already in the
   *    selection.
   **/
  private boolean addWithoutEvent (final FeatureSegment segment) {
    if (segments.contains (segment)) {
      return false;
      // do nothing
    } else {
      segments.addElement (segment);
      return true;
    }
  }

  /**
   *  Remove all the objects from the selection, add the given Feature object
   *  to the selection, and then send an appropriate event to the
   *  SelectionChangeEvent listeners.
   **/
  public void set (final Feature feature) {
    clearWithoutEvent ();
    addWithoutEvent (feature);
    changeSelection (SelectionChangeEvent.SELECTION_CHANGED);
  }

  /**
   *  Remove all the objects from the selection, add the given FeatureSegment
   *  object to the selection, and then send an appropriate event to the
   *  SelectionChangeEvent listeners.
   **/
  public void set (final FeatureSegment feature_segment) {
    clearWithoutEvent ();
    addWithoutEvent (feature_segment);
    changeSelection (SelectionChangeEvent.SELECTION_CHANGED);
  }

  /**
   *  Remove all the objects from the selection, add the given Feature objects
   *  to the selection, and then send an appropriate event to the
   *  SelectionChangeEvent listeners.
   **/
  public void set (final FeatureVector features) {
    clearWithoutEvent ();
    for (int i = 0 ; i < features.size () ; ++i) {
      addWithoutEvent (features.elementAt (i));
    }
    changeSelection (SelectionChangeEvent.SELECTION_CHANGED);
  }

  /**
   *  Set the MarkerRange object that this Selection holds and remove all
   *  other objects.
   **/
  public void setMarkerRange (final MarkerRange marker_range) {
    if (this.marker_range != marker_range) {
      final MarkerRange old_marker_range = this.marker_range;
      clearWithoutEvent ();
      this.marker_range = marker_range;
      if (old_marker_range == null) {
        changeSelection (SelectionChangeEvent.SELECTION_CHANGED);
      } else {
        changeSelection (SelectionChangeEvent.OBJECT_CHANGED);
      }
    }
  }


  /**
   *  Return true if and only if there is nothing in the selection.
   **/
  public boolean isEmpty () {
    if (features.size () == 0 &&
        segments.size () == 0 &&
        marker_range == null) {
      return true;
    } else {
      return false;
    }
  }

  /**
   *  Remove all the objects from the selection, and send an event to the
   *  SelectionChangeEvent listeners.
   **/
  public void clear () {
    if (!isEmpty ()) {
      clearWithoutEvent ();
      changeSelection (SelectionChangeEvent.SELECTION_CHANGED);
    }
  }

  /**
   *  Remove all the objects from the selection, but don't send an event to
   *  the SelectionChangeEvent listeners.
   **/
  private void clearWithoutEvent () {
    features.removeAllElements ();
    segments.removeAllElements ();
    marker_range = null;
  }

  /**
   *  Remove a Feature from the selection.
   **/
  public void remove (final Feature feature) {
    if (features.remove (feature)) {
      changeSelection (SelectionChangeEvent.SELECTION_CHANGED);
    }
  }

  /**
   *  Remove the FeatureSegments of the given Feature from the selection.
   **/
  public void removeSegmentsOf (final Feature feature) {
    final FeatureSegmentVector segments = feature.getSegments ();

    for (int segment_index = 0 ;
         segment_index < segments.size () ;
         ++segment_index) {
      final FeatureSegment this_segment = segments.elementAt (segment_index);

      remove (this_segment);
    }
  }

  /**
   *  Remove a FeatureSegment from the selection.
   **/
  public void remove (final FeatureSegment segment) {
    if (segments.removeElement (segment)) {
      changeSelection (SelectionChangeEvent.SELECTION_CHANGED);
    }
  }

  /**
   *  Return true if this selection contains the given Feature.
   **/
  public boolean contains (final Feature feature) {
    if (getSelectedFeatures ().contains (feature)) {
      return true;
    }

    for (int i = 0 ; i < segments.size () ; ++i) {
      final FeatureSegment this_segment = segments.elementAt (i);
      final Feature this_feature = this_segment.getFeature ();

      if (feature == this_feature) {
        return true;
      }
    }

    return false;
  }

  /**
   *  Return true if this selection contains the given Feature.
   **/
  public boolean contains (final FeatureSegment feature_segment) {
    if (getAllSegments ().indexOf (feature_segment) == -1) {
      return false;
    } else {
      return true;
    }
  }

  /**
   *  Return a Range object that exactly covers the current selection.
   *  Returns null if and only if nothing is selected.
   **/
  public Range getSelectionRange () {
    final Marker lowest_base_marker = getLowestBaseOfSelection ();
    final Marker highest_base_marker = getHighestBaseOfSelection ();

    if (lowest_base_marker == null || highest_base_marker == null) {
      return null;
    } else {
      final Strand strand = lowest_base_marker.getStrand ();

      final int first_postion = lowest_base_marker.getRawPosition ();
      final int last_postion = highest_base_marker.getRawPosition ();

      try {
        return new Range (first_postion, last_postion);
      } catch (OutOfRangeException e) {
        throw new Error ("internal error - unexpected exception: " + e);
      }
    }
  }

  /**
   *  Return a Range object for each of the selected objects.
   **/
  public RangeVector getSelectionRanges () {
    final RangeVector return_ranges = new RangeVector ();

    if (getSelectedFeatures ().size () > 0 ||
        getSelectedSegments ().size () > 0) {
      final FeatureSegmentVector segments = getAllSegments ();

      for (int i = 0 ; i < segments.size () ; ++i) {
        final FeatureSegment this_segment = segments.elementAt (i);

        return_ranges.add (this_segment.getRawRange ());
      }
    } else {
      final MarkerRange marker_range = getMarkerRange ();

      if (marker_range != null) {
        return_ranges.add (marker_range.getRawRange ());
      }
    }

    return return_ranges;
  }

  /**
   *  Return a Marker for the base in this selection that is closest to the
   *  first base of the sequence, or null if nothing is selected.
   **/
  public Marker getLowestBaseOfSelection () {
    if (lowest_base_marker != null) {
      return lowest_base_marker;
    }

    if (getSelectedFeatures ().size () > 0 ||
        getSelectedSegments ().size () > 0) {
      final FeatureSegmentVector segments = getAllSegments ();

      int current_min = 99999999;
      Marker current_min_marker = null;

      for (int i = 0 ; i < segments.size () ; ++i) {
        final FeatureSegment this_segment = segments.elementAt (i);

        final Marker start_marker = this_segment.getStart ();
        if (start_marker.getRawPosition () < current_min) {
          current_min = start_marker.getRawPosition ();
          current_min_marker = start_marker;
        }

        final Marker end_marker = this_segment.getEnd ();
        if (end_marker.getRawPosition () < current_min) {
          current_min = end_marker.getRawPosition ();
          current_min_marker = end_marker;
        }
      }

      lowest_base_marker = current_min_marker;
    } else {
      final MarkerRange range = getMarkerRange ();

      if (range == null) {
        lowest_base_marker = null;
      } else {
        if (range.isForwardMarker ()) {
          lowest_base_marker = range.getStart ();
        } else {
          lowest_base_marker = range.getEnd ();
        }
      }
    }

    return lowest_base_marker;
  }

  /**
   *  Return a Marker for the base in this selection that is closest to the
   *  last base of the sequence, or null if nothing is selected.
   **/
  public Marker getHighestBaseOfSelection () {
    if (highest_base_marker != null) {
      return highest_base_marker;
    }

    if (getSelectedFeatures ().size () > 0 ||
        getSelectedSegments ().size () > 0) {
      final FeatureSegmentVector segments = getAllSegments ();

      int current_max = -1;
      Marker current_max_marker = null;

      for (int i = 0 ; i < segments.size () ; ++i) {
        final FeatureSegment this_segment = segments.elementAt (i);

        final Marker start_marker = this_segment.getStart ();
        if (start_marker.getRawPosition () > current_max) {
          current_max = start_marker.getRawPosition ();
          current_max_marker = start_marker;
        }

        final Marker end_marker = this_segment.getEnd ();
        if (end_marker.getRawPosition () > current_max) {
          current_max = end_marker.getRawPosition ();
          current_max_marker = end_marker;
        }
      }

      highest_base_marker = current_max_marker;
    } else {
      final MarkerRange range = getMarkerRange ();

      if (range == null) {
        highest_base_marker = null;
      } else {
        if (range.isForwardMarker ()) {
          highest_base_marker = range.getEnd ();
        } else {
          highest_base_marker = range.getStart ();
        }
      }
    }

    return highest_base_marker;
  }


  /**
   *  Return the first base of the MarkerRange (if it is set), or the first
   *  base of the selected features (if there are any selected features), or
   *  the first base of the selected FeatureSegment objects (if there are any
   *  selected FeatureSegment objects), otherwise null.  If the selection
   *  contains features or segments from both strands, then the result of
   *  calling this method is the same as calling getLowestBaseOfSelection ().
   **/
  public Marker getStartBaseOfSelection () {
    if (start_base_marker != null) {
      return start_base_marker;
    }

    if (getSelectedFeatures ().size () > 0 ||
        getSelectedSegments ().size () > 0) {
      final FeatureSegmentVector segments = getAllSegments ();

      boolean seen_forward_feature = false;
      boolean seen_reverse_feature = false;

      int current_min = 99999999;
      FeatureSegment current_min_segment = null;

      for (int i = 0 ; i < segments.size () ; ++i) {
        final FeatureSegment this_segment = segments.elementAt (i);

        if (this_segment.getFeature ().isForwardFeature ()) {
          seen_forward_feature = true;
        } else {
          seen_reverse_feature = true;
        }

        if (seen_forward_feature && seen_reverse_feature) {
          return getLowestBaseOfSelection ();
        }

        final Marker start_marker = this_segment.getStart ();

        if (current_min_segment == null ||
            start_marker.getPosition () < current_min) {
          current_min = start_marker.getPosition ();
          current_min_segment = this_segment;
        }
      }

      if (current_min_segment == null) {
        // there are no segments with real ranges
        start_base_marker = null;
      } else {
        start_base_marker = current_min_segment.getStart ();
      }
    } else {
      final MarkerRange range = getMarkerRange ();

      if (range == null) {
        start_base_marker = null;
      } else {
        start_base_marker = range.getStart ();
      }
    }

    return start_base_marker;
  }

  /**
   *  Return the last base of the MarkerRange (if it is set), or the last base
   *  of the selected features (if there are any selected features), or the
   *  last base of the selected FeatureSegment objects (if there are any
   *  selected FeatureSegment objects), otherwise null.  If the selection
   *  contains features or segments from both strands, then the result of
   *  calling this method is the same as calling getLowestBaseOfSelection ().
   **/
  public Marker getEndBaseOfSelection () {
    if (end_base_marker != null) {
      return end_base_marker;
    }

    if (getSelectedFeatures ().size () > 0 ||
        getSelectedSegments ().size () > 0) {
      final FeatureSegmentVector segments = getAllSegments ();

      boolean seen_forward_feature = false;
      boolean seen_reverse_feature = false;

      int current_max = -1;
      FeatureSegment current_max_segment = null;

      for (int i = 0 ; i < segments.size () ; ++i) {
        final FeatureSegment this_segment = segments.elementAt (i);

        if (this_segment.getFeature ().isForwardFeature ()) {
          seen_forward_feature = true;
        } else {
          seen_reverse_feature = true;
        }

        if (seen_forward_feature && seen_reverse_feature) {
          return getHighestBaseOfSelection ();
        }

        final Marker start_marker = this_segment.getStart ();

        if (current_max_segment == null ||
            start_marker.getPosition () > current_max) {
          current_max = start_marker.getPosition ();
          current_max_segment = this_segment;
        }
      }

      if (current_max_segment == null) {
        // there are no segments with real ranges
        end_base_marker = null;
      } else {
        end_base_marker = current_max_segment.getEnd ();
      }
    } else {
      final MarkerRange range = getMarkerRange ();

      if (range == null) {
        end_base_marker = null;
      } else {
        end_base_marker = range.getEnd ();
      }
    }

    return end_base_marker;
  }

  /**
   *  Return the bases of the selected features or the selected marker range.
   **/
  public String getSelectedBases () {
    if (getMarkerRange () == null) {
      final StringBuffer buffer = new StringBuffer ();

      final FeatureVector all_features = getAllFeatures ();

      for (int i = 0 ; i < all_features.size () ; ++i) {
        final Feature this_feature = all_features.elementAt (i);
        buffer.append (this_feature.getBases ());
      }

      return buffer.toString ();
    } else {
      return Strand.markerRangeBases (getMarkerRange ());
    }
  }

  /**
   *  Return the object that was set with the setMarkerRange (MarkerRange)
   *  method.
   **/
  public MarkerRange getMarkerRange () {
    return marker_range;
  }

  /**
   *  Return a vector of the Feature objects of this selection.
   **/
  public FeatureVector getSelectedFeatures () {
    return features;
  }

  /**
   *  Return a vector of the FeatureSegment objects of this selection.
   **/
  public FeatureSegmentVector getSelectedSegments () {
    return segments;
  }

  /**
   *  Return a vector of all the FeatureSegment objects of this selection.
   *  This will include the FeatureSegment objects that are owned by Feature
   *  objects in the selection.
   **/
  public FeatureSegmentVector getAllSegments () {
    final FeatureSegmentVector return_segments = new FeatureSegmentVector ();

    for (int i = 0 ; i < features.size () ; ++i) {
      final Feature selection_feature = features.elementAt (i);
      final FeatureSegmentVector segments = selection_feature.getSegments ();

      for (int segment_index = 0 ;
           segment_index < segments.size () ;
           ++segment_index) {
        final FeatureSegment this_segment = segments.elementAt (segment_index);

        if (!return_segments.contains (this_segment)) {
          return_segments.addElement (this_segment);
        }
      }
    }

    for (int i = 0 ; i < segments.size () ; ++i) {
      final FeatureSegment this_segment = segments.elementAt (i);

      if (!return_segments.contains (this_segment)) {
        return_segments.addElement (this_segment);
      }
    }

    return return_segments;
  }

  /**
   *  Return a vector of all the Feature objects of this selection.  This will
   *  include the Feature objects that own the FeatureSegment objects in the
   *  selection.
   **/
  public FeatureVector getAllFeatures () {
    if (all_features == null) {
      all_features = (FeatureVector) features.clone ();
      
      for (int i = 0 ; i < segments.size () ; ++i) {
        final FeatureSegment this_segment = segments.elementAt (i);
        final Feature this_feature = this_segment.getFeature ();
        
        if (!all_features.contains (this_feature)) {
          all_features.add (this_feature);
        }
      }
    }

    return all_features;
  }

  /**
   *  Send an event to those object listening for it.
   *  @param listeners A Vector of the objects that the event should be sent
   *    to.
   *  @param event The event to send
   **/
  private void fireAction (final Vector listeners, final ChangeEvent event) {
    final Vector targets;
    // copied from a book - synchronising the whole method might cause a
    // deadlock
    synchronized (this) {
      targets = (Vector) listeners.clone ();
    }

    for ( int i = 0 ; i < targets.size () ; ++i ) {
      ChangeListener target = (ChangeListener) targets.elementAt (i);

      if (event instanceof SelectionChangeEvent) {
        final SelectionChangeListener selection_change_listener =
          (SelectionChangeListener) target;
        selection_change_listener.selectionChanged ((SelectionChangeEvent) event);
      } else {
        throw new Error ("Selection.fireAction () - unknown event");
      }
    }
  }

  /**
   *  Reset all cached values in the Selection.
   **/
  private void resetCache () {
    lowest_base_marker = null;
    highest_base_marker = null;
    start_base_marker = null;
    end_base_marker = null;
    all_features = null;
  }

  /**
   *  A Vector containing the Feature objects that this selection currently
   *  holds.
   **/
  private FeatureVector features = new FeatureVector ();

  /**
   *  A Vector containing the FeatureSegment objects that this selection
   *  currently holds.
   **/
  private FeatureSegmentVector segments = new FeatureSegmentVector ();

  /**
   *  This is a cache used by getAllFeatures ().  This will be set to null
   *  anytime the selection changes.
   **/
  private FeatureVector all_features = null;

  /**
   *  Each Selection object can hold one MarkerRange.
   **/
  private MarkerRange marker_range = null;

  /**
   *  A vector of those objects listening for selection change events.
   **/
  final private Vector selection_listener_list = new Vector ();

  /**
   *  The system clipboard as passed to the constructor.
   **/
  final private Clipboard clipboard;

  /**
   *  The cached value that getLowestBaseOfSelection () returns/sets.
   **/
  private Marker lowest_base_marker = null;

  /**
   *  The cached value that getHighestBaseOfSelection () returns/sets.
   **/
  private Marker highest_base_marker = null;

  /**
   *  The cached value that getStartBaseOfSelection () returns/sets.
   **/
  private Marker start_base_marker = null;

  /**
   *  The cached value that getEndBaseOfSelection () returns/sets.
   **/
  private Marker end_base_marker = null;
}
