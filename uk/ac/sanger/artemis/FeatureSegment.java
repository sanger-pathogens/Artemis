/* FeatureSegment.java
 *
 * created: Sun Oct 11 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/FeatureSegment.java,v 1.4 2005-11-28 16:46:38 tjc Exp $
 */

package uk.ac.sanger.artemis;

import uk.ac.sanger.artemis.util.*;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.io.FuzzyRange;
import uk.ac.sanger.artemis.io.Location;
import uk.ac.sanger.artemis.sequence.*;

import java.util.Vector;

/**
 *  Objects of this class will generally represent one exon.  We use the name
 *  "segment" because EMBL does not specify that a feature with a join or order
 *  must be a coding sequence.
 *
 *  @author Kim Rutherford
 *  @version $Id: FeatureSegment.java,v 1.4 2005-11-28 16:46:38 tjc Exp $
 *
 **/

public class FeatureSegment
    implements Selectable, SequenceChangeListener, MarkerChangeListener 
{

  private boolean complement;

  /**
   *  Create a new FeatureSegment object.
   *  @param feature The reference of the Feature that contains this segment.
   *  @param range The Range of this segment.
   **/
  public FeatureSegment (final Feature feature,
                         final Range range) 
  {
    this.feature = feature;

    if(feature == null)
      throw new Error();

    setRange(range);

    // check which strand this is on (used to test for trans spliced)
    final Location location = feature.getLocation();
    this.complement = location.isComplement(range);
  }

  /**
   *  Set the Range of this FeatureSegment.  Updates all listeners and
   *  Markers.
   **/
  public void setRange (final Range range) {
    this.range = range;

    final int start_position;
    final int end_position;

    if (getFeature ().isForwardFeature ()) {
      start_position = range.getStart ();
      end_position = range.getEnd ();
    } else {
      // markers are reversed on the reverse strand
      start_position = range.getEnd ();
      end_position = range.getStart ();
    }

    final Strand strand = getFeature ().getStrand ();

    try {
      if (start == null) {
        start = strand.makeMarkerFromRawPosition (start_position);
      } else {
        start.removeMarkerChangeListener (this);
        start.setRawPosition (start_position);
      }
      if (end == null) {
        end = strand.makeMarkerFromRawPosition (end_position);
      } else {
        end.removeMarkerChangeListener (this);
        end.setRawPosition (end_position);
      }
    } catch (OutOfRangeException e) {
      throw new Error (getFeature().getIDString()+ " internal error - " +
                       "unexpected OutOfRangeException for position: " +
                       start_position);
    }

    start.addMarkerChangeListener (this);
    end.addMarkerChangeListener (this);
  }

  /**
   *  Return value of getFrameID ().
   **/
  public final static int NO_FRAME         = -1;

  /**
   *  Return value of getFrameID ().
   **/
  public final static int FORWARD_FRAME_1  =  0;

  /**
   *  Return value of getFrameID ().
   **/
  public final static int FORWARD_FRAME_2  =  1;

  /**
   *  Return value of getFrameID ().
   **/
  public final static int FORWARD_FRAME_3  =  2;

  /**
   *  Return value of getFrameID ().
   **/
  public final static int FORWARD_STRAND   =  3;

  /**
   *  Return value of getFrameID ().
   **/
  public final static int REVERSE_STRAND  =  4;

  /**
   *  Return value of getFrameID ().
   **/
  public final static int REVERSE_FRAME_3 =  5;

  /**
   *  Return value of getFrameID ().
   **/
  public final static int REVERSE_FRAME_2 =  6;

  /**
   *  Return value of getFrameID ().
   **/
  public final static int REVERSE_FRAME_1 =  7;

  /**
   *  Return value of getFrameID ().
   **/
  public final static int SCALE_LINE       =  8;

  /**
   *  Return the owning Feature of this segment (as passed to the constructor).
   **/
  public Feature getFeature () {
    return feature;
  }

  /**
   *  Return the start position Marker of this object.
   **/
  public Marker getStart () {
    return start;
  }

  /**
   *  Return the end position Marker of this object.
   **/
  public Marker getEnd () {
    return end;
  }

  /**
   *  Set the position of start marker of this FeatureSegment.
   **/
  void setStartPosition (final int position)
      throws OutOfRangeException {
    getStart ().setPosition (position);
    updateRange ();
  }

  /**
   *  Set the position of end marker of this FeatureSegment.
   **/
  void setEndPosition (final int position)
      throws OutOfRangeException {
    getEnd ().setPosition (position);
    updateRange ();
  }

  /**
   *  Return the range (in bases of this segment) that was passed to the
   *  constructor.
   **/
  public Range getRawRange () {
    return range;
  }

  /**
   *  Return a MarkerRange that exactly covers this segment.
   **/
  public MarkerRange getMarkerRange () {
    try {
      return new MarkerRange (getFeature ().getStrand (),
                              getStart ().getPosition (),
                              getEnd ().getPosition ());
    } catch (OutOfRangeException e) {
      throw new Error ("internal error - " + e);
    }
  }

  /**
   *  Remove this FeatureSegment from the Feature that contains it.
   **/
  public void removeFromFeature ()
      throws ReadOnlyException, LastSegmentException{
    getFeature ().removeSegment (this);
    feature = null;
  }

  /**
   *  Return the number of bases in this feature (total of all segments).
   **/
  public int getBaseCount () {
    return getEnd ().getPosition () - getStart ().getPosition () + 1;
  }

  /**
   *  Return the bases in this FeatureSegment.
   **/
  public String getBases () {
    final Strand strand = getFeature ().getStrand ();
    try {
      return strand.getSubSequence (new Range (getStart ().getPosition (),
                                               getEnd ().getPosition ()));
    } catch (OutOfRangeException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    }
  }

  /**
   *  Return true if and only if the Marker objects at the ends of this
   *  segment can be changed.
   **/
  public boolean canDirectEdit () {
    return ! (range instanceof FuzzyRange);
  }

  /**
   *  Implementation of the SequenceChangeListener interface.  We listen for
   *  this event so that we can update the range of this segment.
   **/
  public void sequenceChanged (final SequenceChangeEvent event) 
  {
    // the location of the feature itself will change later
    if(event.getType () != SequenceChangeEvent.REVERSE_COMPLEMENT &&
       event.getType () != SequenceChangeEvent.CONTIG_REVERSE_COMPLEMENT &&
       event.getType () != SequenceChangeEvent.CONTIG_REORDER) 
    {
      updateRange ();
    }
  }

  /**
   *  Add this FeatureSegment to the SequenceChangeListener list.
   **/
  void startListening () {
    final Bases bases = getFeature ().getEntry ().getBases ();

    // it doesn't matter what the priority is as long as it is lower than
    // the Marker priority, so that all Marker objects get updated before
    // we try to update the FeatureSegment
    final int PRIORITY = Marker.LISTENER_PRIORITY - 1;
    bases.addSequenceChangeListener (this, PRIORITY);

    this.bases = bases;
  }

  /**
   *  Remove this FeatureSegment from all listener lists.
   **/
  void stopListening () {
    bases.removeSequenceChangeListener (this);
  }

  /**
   *  Adjust the range of this FeatureSegment to match the start and end
   *  markers of this object.
   **/
  private void updateRange () {
    final Range new_range;

    try {
      if (getFeature ().isForwardFeature ()) {
        new_range = range.change (getStart ().getRawPosition (),
                                  getEnd ().getRawPosition ());
      } else {
        new_range = range.change (getEnd ().getRawPosition (),
                                  getStart ().getRawPosition ());
      }

      range = new_range;
    } catch (OutOfRangeException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    }
  }

  /**
   *  Adds the specified event listener to receive marker change events from
   *  this object.  An event will be sent when either the start or end Marker
   *  of this segment changes.
   *  @param l the change event listener.
   **/
  public void addMarkerChangeListener (MarkerChangeListener l) {
    marker_listener_list.addElement (l);
    if (marker_listener_list.size () == 1) {
      // we need to start listening
      start.addMarkerChangeListener (this);
      end.addMarkerChangeListener (this);
    }
  }

  /**
   *  Removes the specified event listener so that it no longer receives
   *  marker change events from this object.
   *  @param l the change event listener.
   **/
  public void removeMarkerChangeListener (MarkerChangeListener l) {
    marker_listener_list.removeElement (l);
    if (marker_listener_list.size () == 0) {
      // stop listening
      start.removeMarkerChangeListener (this);
      end.removeMarkerChangeListener (this);
    }
  }

  /**
   *  This method fixes up the internals of this Feature when a Marker
   *  changes by calling initialise ().
   **/
  public void markerChanged (final MarkerChangeEvent event) {
    updateRange ();
    fireEvent (event);
  }

  /**
   *  Send a MarkerChangeEvent to each object that is listening for it.
   **/
  private void fireEvent (MarkerChangeEvent event) {
    final Vector targets;
    // copied from a book - synchronising the whole method might cause a
    // deadlock
    synchronized (this) {
      targets = (Vector) marker_listener_list.clone ();
    }

    for (int i = 0 ; i < targets.size () ; ++i) {
      final MarkerChangeListener target =
        (MarkerChangeListener) targets.elementAt (i);

      target.markerChanged (event);
    }
  }

  /**
   *  Returns 0, 1 or 2 depending on which translation frame this segment is
   *  in.  A frame shift of zero means that the bases should be translated
   *  starting at the start position of this segment, 1 means start
   *  translating one base ahead of the start position and 2 means start
   *  translating two bases ahead of the start position.
   **/
  private int getFrameShift() 
  {
    // find the number of bases in the segments before this one
    int base_count = 0;

    final FeatureSegmentVector segments =
      getFeature ().getSegments ();

    int direction = 0;
    for(int i = 0; i < segments.size(); ++i) 
    {
      int this_direction;
      if(segments.elementAt(i).isForwardSegment())
        this_direction = 1;
      else
        this_direction = -1;

      if(segments.elementAt (i) == this) 
      {
        if(i != 0 && this_direction != direction)
          base_count = 0;

        break;
      }
      else 
      {
        if(i == 0)
          direction = this_direction;
        else if(this_direction != direction)
          base_count = 0;

        base_count += segments.elementAt (i).getBaseCount ();
      }
    }

    int mod_value = (base_count + 3
                       - (getFeature ().getCodonStart () - 1)) % 3;

    if (mod_value == 1) {
      return 2;
    } else {
      if (mod_value == 2) {
        return 1;
      } else {
        return 0;
      }
    }
  }

  /**
   *  Return the frame ID of the segment.  There are three forward and three
   *  reverse frame lines as well a forward strand and a backward strand line.
   *  This method returns an ID of the line to place this segment on.
   **/
  public int getFrameID () {
    // this will be 0, 1 or 2 depending on which frame the segment is in
    final int start_base_modulo =
      (getStart ().getPosition () - 1 + getFrameShift ()) % 3;

//  if (getFeature ().getStrand ().isForwardStrand ()) 
    if(isForwardSegment())
    {
      switch (start_base_modulo) {
      case 0:
        return FORWARD_FRAME_1;
      case 1:
        return FORWARD_FRAME_2;
      case 2:
        return FORWARD_FRAME_3;
      }
    } else {
      switch (start_base_modulo) {
      case 0:
        return REVERSE_FRAME_1;
      case 1:
        return REVERSE_FRAME_2;
      case 2:
        return REVERSE_FRAME_3;
      }
    }

    return NO_FRAME;
  }

  public boolean isForwardSegment()
  {
    return !complement;
  }

  /**
   *  The reference of the Feature that contains this segment, as passed to
   *  the constructor.
   **/
  private Feature feature;

  /**
   *  The start position of the segment, set in getStart ().
   **/
  private Marker start = null;

  /**
   *  The end position of the segment, set in getEnd ().
   **/
  private Marker end = null;

  /**
   *  This is the Range reference that was passed to the constructor.
   **/
  private Range range = null;

  /**
   *  A vector containing the references of those objects listening for
   *  marker change events.
   **/
  private final Vector marker_listener_list = new Vector ();

  /**
   *  The object used by startListening () and stopListening ().
   **/
  private Bases bases;
}
