/* GotoEventSource.java
 *
 * created: Sat Jan  9 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/GotoEventSource.java,v 1.1 2004-06-09 09:44:52 tjc Exp $
 */

package uk.ac.sanger.artemis;

import uk.ac.sanger.artemis.sequence.*;

/**
 *  This interface should be implemented by those classes whose objects need
 *  to send GotoEvents and can register GotoListeners.
 *
 *  @author Kim Rutherford
 *  @version $Id: GotoEventSource.java,v 1.1 2004-06-09 09:44:52 tjc Exp $
 **/

public interface GotoEventSource {
  /**
   *  This method sends a GotoEvent to all the GotoEvent listeners that will
   *  make the given base visible.
   *  @param base_marker The Marker of the base to make visible.
   **/
  void gotoBase (final Marker base_marker);

  /**
   *  This method sends a GotoEvent to all the GotoEvent listeners that will
   *  make the given base visible.  This just calls makeBaseVisible ().
   *  @param destination_base The base on the forward strand to make visible.
   *  @return The reference of a MarkerRange object for the given base or null
   *    if the call fails.  It will fail only if the destination_base is less
   *    than 1 or greater than the length of the sequence.
   **/
  MarkerRange gotoBase (final int destination_base);

  /**
   *  This method sends a GotoEvent to all the GotoEvent listeners that will
   *  make the first base of the sequence visible.
   **/
  void gotoFirstBase ();

  /**
   *  This method sends a GotoEvent to all the GotoEvent listeners that will
   *  make the last base of the sequence visible.
   **/
  void gotoLastBase ();

  /**
   *  Send the given event to all the GotoListeners.
   **/
  void sendGotoEvent (final GotoEvent goto_event);

  /**
   *  Adds the specified event listener to receive Goto events from this
   *  object.
   *  @param l the GotoEvent listener.
   **/
  void addGotoListener (final GotoListener l);

  /**
   *  Removes the specified event listener so that it no longer receives Goto
   *  events from this object.
   *  @param l the GotoEvent listener.
   **/
  void removeGotoListener (final GotoListener l) ;
  
}


