/* Marker.java
 *
 * created: Wed Oct 28 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/sequence/Marker.java,v 1.4 2007-07-05 11:58:16 tjc Exp $
 */

package uk.ac.sanger.artemis.sequence;

import uk.ac.sanger.artemis.io.PartialSequence;
import uk.ac.sanger.artemis.util.*;
import java.util.Vector;

/**
 *  Objects of this class mark postions on one strand of DNA.
 *
 *  @author Kim Rutherford
 *  @version $Id: Marker.java,v 1.4 2007-07-05 11:58:16 tjc Exp $
 *
 **/

public class Marker {

  /**
   *  Create a new strand marker.  Methods Outside this package should use
   *  Strand.makeMarker () to make a new Marker.
   *  @param strand The strand that this marker is associated with - a marker
   *    only makes sense in the context of a Strand.
   *  @param position The position on the strand that this marker points to.
   *  @exception OutOfRangeException Thrown if the position is less than 1 or
   *    greater than the length of the sequence.
   **/
  Marker (final Strand strand, final int position)
      throws OutOfRangeException {
    markerinternal = new MarkerInternal (strand, position);

    strand.getBases ().addSequenceChangeListener (markerinternal,
                                                  LISTENER_PRIORITY);
  }

  /**
   *  The priority value use when adding a MarkerInternal object as a
   *  SequenceChangeListener.
   **/
  public static final int LISTENER_PRIORITY = Bases.MEDIUM_PRIORITY;

  /**
   *  Return the strand that this marker is associated with.
   **/
  public Strand getStrand () {
    return getInternal ().getStrand ();
  }

  /**
   *  Return the position on the strand that this marker points to.
   **/
  public int getPosition () {
    return getInternal ().getPosition ();
  }

  /**
   *  Set the position of the Marker, without changing the Strand.
   *  @exception OutOfRangeException Thrown if the new position is less than 1
   *    or greater than the length of the sequence.
   **/
  public void setPosition (final int position)
      throws OutOfRangeException {
    getInternal ().setPosition (position);
  }

  /**
   *  Return postion of this Marker on the underlying Bases object.
   **/
  public int getRawPosition () {
    return getStrand ().getRawPosition (getPosition ());
  }

  /**
   *  Set the raw position of the Marker, without changing the Strand.
   *  @exception OutOfRangeException Thrown if the new position is less than 1
   *    or greater than the length of the sequence.
   **/
  public void setRawPosition (final int position)
      throws OutOfRangeException {
    getInternal ().setPosition (getStrand ().getRawPosition (position));
  }

  /**
   *  Returns true if and only this Marker and argument Marker refer to the
   *  same position on the same Strand.
   **/
  public boolean equals (final Marker test_marker) {
    if (test_marker.getPosition () == getPosition () &&
        test_marker.getStrand () == getStrand ()) {
      return true;
    } else {
      return false;
    }
  }

  /**
   *  Adds the specified event listener to receive marker change events from
   *  this object.
   *  @param l the change event listener.
   **/
  public void addMarkerChangeListener (MarkerChangeListener l) {
    if (marker_listener_list.size () == 0) {
      getInternal ().parent = this;
    }
    marker_listener_list.addElement (l);
  }

  /**
   *  Removes the specified event listener so that it no longer receives
   *  marker change events from this object.
   *  @param l the change event listener.
   **/
  public void removeMarkerChangeListener (MarkerChangeListener l) {
    marker_listener_list.removeElement (l);
    if (marker_listener_list.size () == 0) {
      getInternal ().parent = null;
    }
  }

  /**
   *  Return the argument which has the lowest raw position or the first if
   *  the have equal raw positions.
   **/
  public static Marker getRawLowest (final Marker marker_1,
                                     final Marker marker_2) {
    if (marker_1.getRawPosition () <= marker_2.getRawPosition ()) {
      return marker_1;
    } else {
      return marker_2;
    }
  }

  /**
   *  Send a MarkerChangeEvent to those object listening for it.
   *  @param event The event to send
   **/
  void fireEvent (MarkerChangeEvent event) {
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
   *  Return a new Marker object that is on the same Strand as this Marker but
   *  has it's position moved by the given amount.  For example if this Marker
   *  has position 100 then calling moveBy (-50) will return a Marker at
   *  position 50.
   *  @exception OutOfRangeException Thrown if the new position is less than 1
   *    or greater than the length of the sequence.
   **/
  public Marker moveBy (final int offset)
      throws OutOfRangeException {
    return new Marker (getStrand (), getPosition () + offset);
  }

  /**
   *  Calls Strand.removeSequenceChangeListener () for the MarkerInternal
   *  object of this Marker.
   **/
  public void finalize () {
    getStrand ().getBases ().removeSequenceChangeListener (markerinternal);
  }

  /**
   *  Return the MarkerInternal that was created by the constructor and which
   *  contains the guts of the Marker.
   **/
  private MarkerInternal getInternal () {
    return markerinternal;
  }

  /**
   *  This is created by the constructor.  See the documentation for the
   *  MarkerInternal class.
   **/
  private /*final*/ MarkerInternal markerinternal;

  /**
   *  A vector containing the references of those objects listening for
   *  marker change events.
   **/
  private final Vector marker_listener_list = new Vector ();
}

/**
 *  This is an internal class used by Marker to make sure everything is
 *  garbage collected correctly.  The constructor of the Marker class makes a
 *  new MarkerInternal object and calls
 *  Strand.addSequenceChangeListener(new_markerinternal).  The finalizer of the
 *  Marker class calls Strand.removeSequenceChangeListener().  This ensures
 *  that when there are no explicit references to the Marker it won't be keep
 *  alive by having a reference on the sequence change listener list.
 *
 *  @author Kim Rutherford
 *  @version $Id: Marker.java,v 1.4 2007-07-05 11:58:16 tjc Exp $
 *
 **/
class MarkerInternal
    implements SequenceChangeListener {
  /**
   *  Create a new strand marker.  Methods Outside this package should use
   *  Strand.makeMarker () to make a new Marker.
   *  @param strand The strand that this marker is associated with - a marker
   *    only makes sense in the context of a Strand.
   *  @param position The position on the strand that this marker points to.
   *  @exception OutOfRangeException Thrown if the position is less than 1 or
   *    greater than the length of the sequence.
   **/
  public MarkerInternal (final Strand strand, final int position)
      throws OutOfRangeException {
    this.strand   = strand;
    this.position = position;

    checkPosition (position);
  }

  /**
   *  Return the strand that this marker is associated with.
   **/
  public Strand getStrand () {
    return strand;
  }

  /**
   *  Return the position on the strand that this marker points to.
   **/
  public int getPosition () {
    return position;
  }

  /**
   *  Set the position of the Marker, without changing the Strand.
   *  @exception OutOfRangeException Thrown if the new position is less than 1
   *    or greater than the length of the sequence.
   **/
  public void setPosition (final int position)
      throws OutOfRangeException {
    final int old_position = this.position;

    checkPosition (position);

    this.position = position;

    if (parent != null) {
      parent.fireEvent (new MarkerChangeEvent (parent,
                                               getStrand (),
                                               old_position));
    }
  }

  /**
   *  Implementation of the SequenceChangeListener interface.  We listen here
   *  so that the Marker can be updated when the sequence changes.
   **/
  public void sequenceChanged (final SequenceChangeEvent event) {
    if(event.getType () == SequenceChangeEvent.REVERSE_COMPLEMENT ||
       event.getType () == SequenceChangeEvent.CONTIG_REVERSE_COMPLEMENT ||
       event.getType () == SequenceChangeEvent.CONTIG_REORDER) 
      return;

    final int sub_sequence_length = event.getSubSequence ().length ();

    final Bases bases = getStrand ().getBases ();

    final int event_start_base;

    if (getStrand ().isForwardStrand ()) {
      event_start_base = event.getPosition ();
    } else {
      if (event.getType () == SequenceChangeEvent.INSERTION) {
        // this is the length of sequence before the insertion
        final int new_sequence_length =
          bases.getLength () - sub_sequence_length;
        event_start_base = new_sequence_length - event.getPosition ();
      } else {
        event_start_base = bases.getLength () - event.getPosition ();
      }
    }
      
    if (event.getType () == SequenceChangeEvent.DELETION) {
      if (position <= event_start_base) {
        return;
      } else {
        if (event_start_base + sub_sequence_length > position) {
          position = event_start_base;
        } else {
          position -= sub_sequence_length;
        }
      }
    } else {
      if (position < event_start_base) {
        return;
      } else {
        position += sub_sequence_length;
      }
    }
  }

  /**
   *  Check that the given position is within the range of the Strand.
   *  @exception OutOfRangeException Thrown if the position is less than 1
   *    or greater than the length of the sequence.
   **/
  private void checkPosition (final int position)
      throws OutOfRangeException {
    if (position < 1 || position > strand.getSequenceLength ()) 
    {
      if(!(strand.getBases().getSequence() instanceof PartialSequence))
        throw new OutOfRangeException ("position: " + position);
    }
  }

  /**
   *  The strand that this marker is associated with.
   **/
  private Strand strand;

  /**
   *  The position on the strand that this marker points to.
   **/
  private int position;

  /**
   *  This will be null when the Marker object has know registered listeners
   *  and will contain the reference of the Marker that created this
   *  MarkerInternal object when the Marker has any listeners.
   **/
  Marker parent = null;
}
