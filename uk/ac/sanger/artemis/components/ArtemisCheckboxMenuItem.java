/* CheckboxMenuItem.java
 *
 * created: Mon Aug 19 2002
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 2001,2002  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/ArtemisCheckboxMenuItem.java,v 1.1 2004-06-09 09:46:02 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import java.awt.*;
import java.awt.event.*;

import java.util.Vector;

/**
 *  This component is a replacement for the CheckboxMenuItem component.  It is
 *  needed on GNU/Linux with versions 1.2, 1.3 and 1.4 or the VM.  See this
 *  bug for more:
 *    http://developer.java.sun.com/developer/bugParade/bugs/4533641.html
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: ArtemisCheckboxMenuItem.java,v 1.1 2004-06-09 09:46:02 tjc Exp $
 **/

public class ArtemisCheckboxMenuItem
    extends MenuItem
    implements ActionListener, ItemSelectable {
  /**
   *  Make a new ArtemisCheckboxMenuItem with the given label.
   **/
  public ArtemisCheckboxMenuItem (final String label) {
    super (label);

    this.orig_label = label;

    addActionListener (this);
  }

  /**
   *  Implementation of ActionListener.
   **/
  public void actionPerformed (ActionEvent _) {
    final ItemEvent e;
    
    if (getState ()) {
      setState (false);
      e = new ItemEvent (this, ItemEvent.ITEM_STATE_CHANGED,
                         this, ItemEvent.DESELECTED);
    } else {
      setState (true);
      e = new ItemEvent (this, ItemEvent.ITEM_STATE_CHANGED,
                         this, ItemEvent.SELECTED);
    }

    for (int i = 0 ; i < listeners.size () ; ++i) {
      ((ItemListener) listeners.elementAt (i)).itemStateChanged (e) ;
    }
  }

  /**
   *  Returns null (implementation of ItemSelectable).
   **/
  public Object[] getSelectedObjects() {
    return null;
  }

  /**
   *  Sets this check box menu item to the specifed state. The boolean value
   *  true indicates "on" while false indicates "off."
   *  @param state true if the check box menu item is on, otherwise false
   **/
  public void setState (final boolean state) {
    this.state = state;
    if (state) {
      setLabel ("Disable " + orig_label + " (current on)");
    } else {
      setLabel ("Enable " + orig_label + " (current off)");
    }
  }

  /**
   *  Determines whether the state of this check box menu item is "on" or
   *  "off."
   *  @return the state of this check box menu item, where true indicates "on"
   *    and false indicates "off" 
   **/
  public boolean getState () {
    return state;
  }

  /**
   *  Add the given listener.
   **/
  public void addItemListener (final ItemListener listener) {
    listeners.addElement (listener);
  }

  /**
   *  remove the given listener.
   **/
  public void removeItemListener (final ItemListener listener) {
    listeners.removeElement (listener);
  }

  private String orig_label;
  private boolean state = false;

  private Vector listeners = new Vector ();
}
