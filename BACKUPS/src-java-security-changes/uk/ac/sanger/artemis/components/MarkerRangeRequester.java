/* MarkerRangeRequester.java
 *
 * created: Mon Jul 10 2000
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/MarkerRangeRequester.java,v 1.1 2004-06-09 09:47:04 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import java.util.Vector;

/**
 *  A requester that gets a MarkerRange or a Range from the user.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: MarkerRangeRequester.java,v 1.1 2004-06-09 09:47:04 tjc Exp $
 **/

public class MarkerRangeRequester extends TextRequester {
  /**
   *  Create a new MarkerRangeRequester component with the given prompt. Other
   *  components can listen for MarkerRangeRequesterEvent object.
   *  @param prompt A message that is displayed in the component beside the
   *    TextArea that the user types into.  This String is also used as the
   *    JFrame title.
   *  @param width The width of the TextRequester in the new requester.
   *  @param initial_text The initial text to put in the TextRequester.
   **/
  public MarkerRangeRequester (final String prompt,
                               final int width,
                               final String initial_text) {
    super (prompt, width, initial_text);
  }

  /**
   *  Add the given object as a listen for MarkerRangeRequester events from
   *  this MarkerRangeRequester.
   **/
  public void
    addMarkerRangeRequesterListener (final MarkerRangeRequesterListener l) {
    listeners.addElement (l);
  }

  /**
   *  Send a MarkerRangeRequesterEvent of type OK to all the listeners.
   **/
  protected void performOK () {
    final MarkerRangeRequesterEvent new_event =
      new MarkerRangeRequesterEvent (this, getText (),
                                     MarkerRangeRequesterEvent.OK);

    sendEvent (new_event);

    super.performOK ();
  }

  /**
   *  Send a MarkerRangeRequesterEvent of type CANCEL to all the listeners.
   **/
  protected void performCancel () {
    final MarkerRangeRequesterEvent new_event =
      new MarkerRangeRequesterEvent (this, getText (),
                                     MarkerRangeRequesterEvent.CANCEL);

    sendEvent (new_event);

    super.performCancel ();
  }

  /**
   *  Send the given MarkerRangeRequesterEvent to all the object that are
   *  listening for the event.
   **/
  private void sendEvent (final MarkerRangeRequesterEvent event) {
    for (int i = 0 ; i < listeners.size () ; ++i) {
      final MarkerRangeRequesterListener listener =
        ((MarkerRangeRequesterListener) listeners.elementAt (i));

      listener.actionPerformed (event);
    }
  }

  /**
   *  This contains the objects that are listening for MarkerRangeRequester
   *  events from this MarkerRangeRequester.
   **/
  private Vector listeners = new Vector ();
}
