/* EntryChangeEvent.java
 *
 * created: Sat Oct 17 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/EntryChangeEvent.java,v 1.2 2006-05-31 10:38:48 tjc Exp $
 */

package uk.ac.sanger.artemis;

/**
 *  This event is sent when a change occurs in an entry.  eg. a feature is
 *  deleted.
 *
 *  @author Kim Rutherford
 *  @version $Id: EntryChangeEvent.java,v 1.2 2006-05-31 10:38:48 tjc Exp $
 *
 */

public class EntryChangeEvent extends ChangeEvent {
  /**
   *  Event type - feature removed.
   **/
  public static final int FEATURE_DELETED = 1;

  /**
   *  Event type - feature added.
   **/
  public static final int FEATURE_ADDED = 2;

  /**
   *  Event type - The name of an entry has changed.
   **/
  public static final int NAME_CHANGED = 3;

  /**
   *  Event type - The header of an entry has changed.
   **/
  public static final int HEADER_CHANGED = 4;

  /**
   *  The Entry that was passed to the constructor (if any).
   **/
  private Entry entry;
  
  /**
   *  The Feature object that was passed to the constructor.
   **/
  private Feature feature;
    
  /**
   *  
   **/
  private boolean duplicate; 
  
  /**
   *  This is the type of this event (eg FEATURE_DELETED, FEATURE_ADDED, etc),
   *  as passed to the constructor
   **/
  private int type;
  
  /**
   *  Create a new EntryChangeEvent object.
   *  @param entry This Entry object that this event refers to.
   *  @param feature This Feature object that this event refers to.
   *  @param type This type of the event.
   **/
  public EntryChangeEvent (Entry entry,
                           Feature feature,
                           int type) 
  {
    super (entry);
    this.entry = entry;
    this.feature = feature;
    this.type = type;
  }

  /**
   *  Create a new EntryChangeEvent object.
   *  @param entry This Entry object that this event refers to.
   *  @param feature This Feature object that this event refers to.
   *  @param type This type of the event.
   **/
  public EntryChangeEvent (Entry entry,
                           Feature feature,
                           boolean duplicate,
                           int type) 
  {
    super (entry);
    this.entry = entry;
    this.feature = feature;
    this.duplicate = duplicate;
    this.type = type;
  }
  
  /**
   *  Create a new EntryChangeEvent object.  This constructor is used for
   *  HEADER_CHANGED and NAME_CHANGED events.
   *  @param source The source of this event.
   *  @param entry This Entry object that this event refers to.
   *  @param type This type of the event.
   **/
  public EntryChangeEvent (Object source,
                           Entry entry,
                           int type)
  {
    super (source);
    this.entry = entry;
    this.feature = null;
    this.type = type;
  }

  /**
   *  Return the type of this event, ie the type passed to the
   *  constructor.
   **/
  public int getType () 
  {
    return type;
  }

  /**
   *  Return the target Feature object for this event.
   **/
  public Feature getFeature () 
  {
    return feature;
  }

  
  public boolean isDuplicate()
  {
    return duplicate;
  }

  /**
   *  This is a convenience method that returns the Entry the generated this
   *  event.
   **/
  public Entry getEntry () 
  {
    return entry;
  }
}
