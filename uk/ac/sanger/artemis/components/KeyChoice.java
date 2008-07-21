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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/KeyChoice.java,v 1.5 2008-07-21 15:43:16 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.io.KeyVector;
import uk.ac.sanger.artemis.io.EntryInformation;

import java.awt.Font;
import java.awt.event.ItemListener;

import javax.swing.JComboBox;
import javax.swing.JPanel;

/**
 *  This component is a Choice component that shows the possible feature keys.
 *
 *  @author Kim Rutherford
 *  @version $Id: KeyChoice.java,v 1.5 2008-07-21 15:43:16 tjc Exp $
 **/
public class KeyChoice extends JPanel 
{
  private static final long serialVersionUID = 1L;

  /** The JComboBox component that will show the feature keys. */
  private JComboBox key_chooser = null;
  
  /**
   *  Create a new KeyChoice component with CDS as the default key.
   *  @param entry_information The object to get the list of possible
   *    keys from.
   **/
  public KeyChoice (final EntryInformation entry_information) 
  {
    this (entry_information, Key.CDS);
  }

  /**
   *  Create a new KeyChoice component with the given key as the default.
   *  @param entry_information The object to get the list of possible
   *    keys from.
   **/
  public KeyChoice (final EntryInformation entry_information,
                    final Key default_key) 
  {
    final Font font = Options.getOptions ().getFont ();
    setFont (font);

    key_chooser = new JComboBox ();
    key_chooser.setEditable(true);

    final int MAX_VISIBLE_ROWS = 30;
    key_chooser.setMaximumRowCount (MAX_VISIBLE_ROWS);

    addChoice (entry_information, default_key);
  }

  /**
   *  Return the currently selected key.
   **/
  public Key getSelectedItem () 
  {
    return new Key ((String) key_chooser.getSelectedItem ());
  }

  /**
   *  Set the selected Key.
   **/
  public void setKey (final Key new_key) 
  {
    final int key_index = keyIndex (new_key);

    if (key_index == -1) 
    {
      // add the key
      key_chooser.addItem (new_key.toString ());
      key_chooser.setSelectedItem (new_key.toString ());
    } 
    else
      key_chooser.setSelectedIndex (key_index);
  }

  /**
   *  Adds the specified item listener to receive item events from the Choice
   *  component of this KeyChoice.
   *  @param l The item listener.
   **/
  public void addItemListener(ItemListener l) 
  {
    key_chooser.addItemListener (l);
  }


  /**
   *  Removes the specified item listener so that it no longer receives item
   *  events from the Choice component of this KeyChoice.
   *  @param l The item listener.
   **/
  public void removeItemListener(ItemListener l)
  {
    key_chooser.removeItemListener (l);
  }

  /**
   *  Add the key_chooser.
   **/
  private void addChoice (final EntryInformation entry_information,
                          final Key default_key) 
  {
    final Font font = Options.getOptions ().getFont ();
    key_chooser.setFont (font);

    KeyVector keys = entry_information.getSortedValidKeys();

    if (keys == null) 
    {
      keys = new KeyVector ();
      keys.add (Key.CDS);
    }

    for (int i = 0 ; i < keys.size () ; ++i)
      key_chooser.addItem ( ((Key)keys.get(i)).toString ());

    if (keyIndex (default_key) != -1)
      setKey (default_key);
    add (key_chooser);
  }

  /**
   *  Return the index in the key_chooser component of the given Key.
   **/
  private int keyIndex (final Key key) 
  {
    for (int i = 0 ; i < key_chooser.getItemCount () ; ++i) 
    {
      if (key.toString ().equals ((String)key_chooser.getItemAt (i))) 
        return i;
    }
    return -1;
  }

  public void setEnabled(final boolean isEnabled)
  {
    key_chooser.setEnabled(isEnabled);
    repaint();
  }
}