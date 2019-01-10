/* SelectionChangeEvent.java
 *
 * created: Tue Oct 27 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/SelectionChangeEvent.java,v 1.1 2004-06-09 09:45:05 tjc Exp $
 **/

package uk.ac.sanger.artemis;

/**
 *  This event is sent when the currently selected object changes.
 *
 *  @author Kim Rutherford
 *  @version $Id: SelectionChangeEvent.java,v 1.1 2004-06-09 09:45:05 tjc Exp $
 **/

public class SelectionChangeEvent extends ChangeEvent {
  /**
   *  Create a new SelectionChangeEvent object.
   *  @param source The Object that generated the event - probably a component.
   *  @param type The type of event.  Should be one of OBJECT_ADDED,
   *    OBJECT_REMOVED or OBJECT_CHANGED.
   **/
  public SelectionChangeEvent (final Object source,
                               final int type) {
    super (source);
    this.type = type;
  }

  /**
   *  Return the type that was passed to the constructor.
   **/
  public int getType () {
    return type;
  }

  /**
   *  Flag used when the event is caused by the addition or removal of an
   *  object from the selection.
   **/
  public final static int SELECTION_CHANGED = 1;

  /**
   *  Flag used for all other events.  Examples include the case when a
   *  qualifier or location changes for a selected feature.
   **/
  public final static int OBJECT_CHANGED = 3;
  
  /**
   *  The event type that was passed to the constructor.
   **/
  private int type;
}


