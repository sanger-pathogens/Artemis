/* Action.java
 *
 * created: Tue Sep 17 2002
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/Action.java,v 1.2 2008-01-28 16:27:54 tjc Exp $
 */

package uk.ac.sanger.artemis;


/**
 *  An Action is anything done by the user that causes the state of the data
 *  in Artemis to change.  It is implemented as Vector of uk.ac.sanger.artemis.ChangeEvent
 *  objects.  For example an Action might be the result of the user doing
 *  "Trim to Met".  In that case there will be one Action that contains a
 *  FeatureChangeEvent for each Feature that was trimmed.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: Action.java,v 1.2 2008-01-28 16:27:54 tjc Exp $
 **/

public class Action 
{
  /**
   *  Delegate.
   **/
  private final ChangeEventVector change_vector = new ChangeEventVector ();
  
  /**
   *  Create a new, empty Action.
   **/
  public Action () 
  {
  }

  /**
   *  Returns true if and only if nothing has happened during this action.
   **/
  public boolean isEmpty () 
  {
    if (change_vector.size () == 0)
      return true;
    else
      return false;
  }

  /**
   *  Add a new ChangeEvent to this Action.
   **/
  public void addChangeEvent (final ChangeEvent change_event) 
  {
    if (change_event instanceof FeatureChangeEvent) 
      addFeatureChangeEvent ((FeatureChangeEvent) change_event);
    else if (change_event instanceof EntryChangeEvent) 
      addEntryChangeEvent ((EntryChangeEvent) change_event);
    else 
      throw new Error ("internal error - unknown event type: " +
                       change_event);
  }

  /**
   *  Add a new FeatureChangeEvent to this Action.
   **/
  public void addFeatureChangeEvent (final FeatureChangeEvent
                                       feature_change_event) 
  {
    if (feature_change_event.featureHasChanged ())
      change_vector.add (feature_change_event);
  }
 
  /**
   *  Add a new EntryChangeEvent to this Action.
   **/
  public void addEntryChangeEvent (final EntryChangeEvent entry_change_event)
  {
    change_vector.add (entry_change_event);
  }

  /**
   *  Return a ChangeEventVector containing all the ChangeEvents that occured
   *  during this Action.
   **/
  public ChangeEventVector getChangeEvents () 
  {
    return change_vector;
  }
}
