/* KeyChoice.java
 *
 * created: Mon Sep  6 1999
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 1999  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/KeyChoice.java,v 1.4 2008-05-29 15:16:28 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;

import uk.ac.sanger.artemis.util.StringVector;
import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.io.KeyVector;
import uk.ac.sanger.artemis.io.EntryInformation;

import java.awt.*;
import java.awt.event.*;

import javax.swing.*;

/**
 *  This component is a Choice component that shows the possible feature keys.
 *
 *  @author Kim Rutherford
 *  @version $Id: KeyChoice.java,v 1.4 2008-05-29 15:16:28 tjc Exp $
 **/
public class KeyChoice extends JPanel {
  /**
   *  Create a new KeyChoice component with CDS as the default key.
   *  @param entry_information The object to get the list of possible
   *    keys from.
   **/
  public KeyChoice (final EntryInformation entry_information) {
    this (entry_information, Key.CDS);
  }

  /**
   *  Create a new KeyChoice component with the given key as the default.
   *  @param entry_information The object to get the list of possible
   *    keys from.
   **/
  public KeyChoice (final EntryInformation entry_information,
                    final Key default_key) {
    this.entry_information = entry_information;
    this.default_key = default_key;

    final Font font = Options.getOptions ().getFont ();
    setFont (font);

    key_chooser = new JComboBox ();
    key_chooser.setEditable(true);

    final int MAX_VISIBLE_ROWS = 30;
    
    key_chooser.setMaximumRowCount (MAX_VISIBLE_ROWS);

    updateChoice ();
  }

  /**
   *  Return the currently selected key.
   **/
  public Key getSelectedItem () {
    final Key key;

    key = new Key ((String) key_chooser.getSelectedItem ());

    return key;
  }

  /**
   *  Set the selected Key.
   **/
  public void setKey (final Key new_key) {
    final int key_index = keyIndex (new_key);

    if (key_index == -1) {
      // add the key
      key_chooser.addItem (new_key.toString ());
      key_chooser.setSelectedItem (new_key.toString ());
    } else {
      key_chooser.setSelectedIndex (key_index);
    }
  }

  /**
   *  Adds the specified item listener to receive item events from the Choice
   *  component of this KeyChoice.
   *  @param l The item listener.
   **/
  public void addItemListener(ItemListener l) {
    key_chooser.addItemListener (l);
  }


  /**
   *  Removes the specified item listener so that it no longer receives item
   *  events from the Choice component of this KeyChoice.
   *  @param l The item listener.
   **/
  public void removeItemListener(ItemListener l) {
    key_chooser.removeItemListener (l);
  }

  /**
   *  Update the key_chooser.
   **/
  private void updateChoice () {
    removeAll ();

    final Font font = Options.getOptions ().getFont ();
    key_chooser.setFont (font);

    add (key_chooser);

    KeyVector keys = getSortedKeys ();

    if (keys == null) {
      keys = new KeyVector ();
      keys.add (Key.CDS);
    }

    for (int i = 0 ; i < keys.size () ; ++i) {
      key_chooser.addItem ( ((Key)keys.get(i)).toString ());
    }

    if (keyIndex (default_key) != -1) {
      setKey (default_key);
    }

    // XXX change to revalidate().
    validate ();

    if (getParent () != null) {
      // XXX change to revalidate().
      getParent ().validate ();
      if (getParent ().getParent () != null) {
        // XXX change to revalidate().
        getParent ().getParent ().validate ();
      }
    }
  }

  /**
   *  Return a alphanumerically sorted vector containing the String
   *  representations of the common keys (those listed in the common_keys
   *  property of the options file).  If there is no common_keys option then
   *  all the legal keys are returned.
   **/
  private KeyVector getSortedKeys () {
    return entry_information.getSortedValidKeys ();
  }

  /**
   *  Return the index in the key_chooser component of the given Key.
   **/
  private int keyIndex (final Key key) {
    for (int i = 0 ; i < key_chooser.getItemCount () ; ++i) {
      if (key.toString ().equals ((String)key_chooser.getItemAt (i))) {
        return i;
      }
    }
    return -1;
  }

  public void setEnabled(final boolean isEnabled)
  {
    key_chooser.setEnabled(isEnabled);
    repaint();
  }
  
  /**
   *  The Key that was passed to the constructor.
   **/
  private Key default_key = null;

  /**
   *  The JComboBox component that will show the feature keys.
   **/
  private JComboBox key_chooser = null;

  /**
   *  This toggle sets whether to show common keys or uncommon keys.
   **/
  private JCheckBox common_keys_checkbox = null;

  /**
   *  The EntryInformation object that was passed to the constructor.
   **/
  private EntryInformation entry_information = null;
}
