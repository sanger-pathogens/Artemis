/* EntryGroupChangeEvent.java
 *
 * created: Mon Dec  7 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/EntryGroupChangeEvent.java,v 1.1 2004-06-09 09:44:22 tjc Exp $
 */

package uk.ac.sanger.artemis;

/**
 *  This event is sent when a change occurs in an EntryGroup.  eg. an Entry
 *  is deleted or added.
 *
 *  @author Kim Rutherford
 *  @version $Id: EntryGroupChangeEvent.java,v 1.1 2004-06-09 09:44:22 tjc Exp $
 **/

public class EntryGroupChangeEvent extends ChangeEvent {
  /**
   *  Event type - Entry removed.
   **/
  public static final int ENTRY_DELETED = 1;

  /**
   *  Event type - Entry added.
   **/
  public static final int ENTRY_ADDED = 2;

  /**
   *  Event type - Entry has been made active.
   **/
  public static final int ENTRY_ACTIVE = 3;

  /**
   *  Event type - Entry has been made inactive.
   **/
  public static final int ENTRY_INACTIVE = 4;

  /**
   *  Event type - There is now a different default Entry.
   **/
  public static final int NEW_DEFAULT_ENTRY = 5;

  /**
   *  Event type - This event means that the entry group is now not used by
   *  anything important.
   **/
  public static final int DONE_GONE = 6;

  /**
   *  Create a new EntryChangeEvent object.
   *  @param entry_group This EntryGroup object that this event refers to.
   *  @param entry This Entry object that this event refers to (if any).
   *  @param type This type of the event.
   **/
  public EntryGroupChangeEvent (final EntryGroup entry_group,
                                final Entry entry,
                                final int type) {
    super (entry_group);
    this.entry = entry;
    this.type = type;
  }

  /**
   *  Return the type of this event ie. the type passed to the constructor.
   **/
  public int getType () {
    return type;
  }

  /**
   *  Return the target Entry object for this event.  Will return null for
   *  DONE_GONE events.
   **/
  public Entry getEntry () {
    return entry;
  }

  /**
   *  This is a convenience method that returns the EntryGroup the generated
   *  this event.
   **/
  public EntryGroup getEntryGroup () {
    return (EntryGroup) getSource ();
  }

  /**
   *  The Entry object that was passed to the constructor.
   **/
  private Entry entry;

  /**
   *  This is the type of this event (eg ENTRY_ADDED, ENTRY_DELETED, etc), as
   *  passed to the constructor
   **/
  private int type;
}
