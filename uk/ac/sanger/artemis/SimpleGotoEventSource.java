/* SimpleGotoEventSource.java
 *
 * created: Sat Jun 17 2000
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 2000  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/SimpleGotoEventSource.java,v 1.1 2004-06-09 09:45:10 tjc Exp $
 */

package uk.ac.sanger.artemis;

import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.sequence.*;

import uk.ac.sanger.artemis.util.*;

import java.util.Vector;
import java.util.EventObject;

/**
 *  An simple implementation of GotoEventSource.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: SimpleGotoEventSource.java,v 1.1 2004-06-09 09:45:10 tjc Exp $
 **/

public class SimpleGotoEventSource
    implements GotoEventSource {
  /**
   *  Create a new SimpleGotoEventSource that operates on the given
   *  EntryGroup.  We need an EntryGroup reference so we know where the last
   *  base is.
   **/
  public SimpleGotoEventSource (final EntryGroup entry_group) {
    this.entry_group = entry_group;
  }

  /**
   *  This method sends an GotoEvent to all the GotoEvent listeners that will
   *  make the given base visible.
   *  @param base_marker The base to make visible.
   **/
  public void gotoBase (final Marker base_marker) {
    final GotoEvent new_event = new GotoEvent (this, base_marker);

    sendGotoEvent (new_event);
  }

  /**
   *  This method sends a GotoEvent to all the GotoEvent listeners that will
   *  make the first base of the sequence visible.
   **/
  public void gotoFirstBase () {
    gotoBase (1);
  }

  /**
   *  This method sends a GotoEvent to all the GotoEvent listeners that will
   *  make the last base of the sequence visible.
   **/
  public void gotoLastBase () {
    gotoBase (getEntryGroup ().getSequenceLength ());
  }

  /**
   *  This method sends an GotoEvent to all the GotoEvent listeners that will
   *  make the given base visible.
   *  @param destination_base The base to make visible.
   *  @return The reference of a MarkerRange object for the given base or null
   *    if the call fails.  It will fail only if the destination_base is less
   *    than 1 or greater than the length of the sequence.
   **/
  public MarkerRange gotoBase (final int destination_base) {
    try {
      final Strand forward_strand =
        getEntryGroup ().getBases ().getForwardStrand ();
      final Marker destination_marker =
        forward_strand.makeMarker (destination_base);
      gotoBase (destination_marker);
      return new MarkerRange (destination_marker.getStrand (),
                              destination_base, destination_base);
    } catch (OutOfRangeException e) {
      return null;
    }
  }

  /**
   *  Send the given event to all the GotoListeners.
   **/
  public void sendGotoEvent (final GotoEvent goto_event) {
    // don't send if there is nowhere to go to
    if (goto_event.getMarker () != null) {
      fireAction (goto_listener_list, goto_event);
    }
  }

  /**
   *  Send an event to those object listening for it.
   *  @param listeners A Vector of the objects that the event should be sent
   *    to.
   *  @param event The event to send
   **/
  private void fireAction (final Vector listeners, final EventObject event) {
    final Vector targets;
    // copied from a book - synchronising the whole method might cause a
    // deadlock
    synchronized (this) {
      targets = (Vector) listeners.clone ();
    }

    for (int i = 0 ; i < targets.size () ; ++i) {
      GotoListener target = (GotoListener) targets.elementAt (i);

      if (event instanceof GotoEvent) {
        final GotoListener goto_listener = (GotoListener) target;
        goto_listener.performGoto ((GotoEvent) event);
      } else {
        throw new Error ("EntryEdit.fireAction () - unknown event");
      }
    }
  }

  /**
   *  Adds the specified event listener to receive Goto events from this
   *  object.
   *  @param l the GotoEvent listener.
   **/
  public void addGotoListener (final GotoListener l) {
    goto_listener_list.addElement (l);
  }

  /**
   *  Removes the specified event listener so that it no longer receives Goto
   *  events from this object.
   *  @param l the GotoEvent listener.
   **/
  public void removeGotoListener (final GotoListener l) {
    goto_listener_list.removeElement (l);
  }

  /**
   *  Return the EntryGroup object that was passed to the constructor.
   **/
  public EntryGroup getEntryGroup () {
    return entry_group;
  }

  /**
   *  The EntryGroup object that was passed to the constructor.
   **/
  private EntryGroup entry_group;

  /**
   *  A vector of those objects listening for Goto events.
   **/
  final private Vector goto_listener_list = new Vector ();
}
