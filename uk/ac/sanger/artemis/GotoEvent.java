/* GotoEvent.java
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/GotoEvent.java,v 1.1 2004-06-09 09:44:51 tjc Exp $
 */

package uk.ac.sanger.artemis;

import uk.ac.sanger.artemis.sequence.Marker;

/**
 *  This event is sent when a view component (eg FeatureDisplay) should
 *  scroll to centre on a new base.
 *
 *  @author Kim Rutherford
 *  @version $Id: GotoEvent.java,v 1.1 2004-06-09 09:44:51 tjc Exp $
 **/

public class GotoEvent extends java.util.EventObject {
  /**
   *  Create a new GotoEvent object from the given Marker.
   *  @param source The Object that generated the event - probably a component.
   *  @param goto_position This is the position that the listeners should goto.
   **/
  public GotoEvent (final Object source, final Marker goto_position) {
    super (source);
    this.goto_position = goto_position;
  }

  /**
   *  Return the Marker that was passed to the constructor.
   **/
  public Marker getMarker () {
    return goto_position;
  }

  /**
   *  This is the Marker that was passed to the constructor.
   **/
  final private Marker goto_position;
}


