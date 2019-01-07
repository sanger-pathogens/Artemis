/* EntryGroupMenu.java
 *
 * created: Fri Aug 27 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/EntryGroupMenu.java,v 1.2 2004-12-03 18:11:28 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;

import java.util.*;
import java.awt.*;
import java.awt.event.*;

import javax.swing.*;

/**
 *  A menu containing the current entries in an EntryGroup.
 *
 *  @author Kim Rutherford
 *  @version $Id: EntryGroupMenu.java,v 1.2 2004-12-03 18:11:28 tjc Exp $
 **/

public class EntryGroupMenu extends JMenu
    implements EntryGroupChangeListener, EntryChangeListener {
  /**
   *  Create a new EntryGroupMenu object that will show the entries in the
   *  given EntryGroup.
   *  @param frame The JFrame that owns this JMenu.
   *  @param entry_group The EntryGroup object to display.
   *  @param menu_name The name of the new menu.
   **/
  public EntryGroupMenu (final JFrame frame,
                         final EntryGroup entry_group,
                         final String menu_name) {
    super (menu_name);

    this.frame = frame;
    this.entry_group = entry_group;

    entry_group.addEntryGroupChangeListener (this);
    entry_group.addEntryChangeListener (this);

    refreshMenu ();
  }

  /**
   *  Create a new EntryGroupMenu object that will show the entries in the
   *  given EntryGroup.
   *  @param frame The JFrame that owns this JMenu.
   *  @param entry_group The EntryGroup object to display.
   **/
  public EntryGroupMenu (final JFrame frame,
                         final EntryGroup entry_group) {
    this (frame, entry_group, "Entries");
  }

  /**
   *  Implementation of the EntryGroupChangeListener interface.  We listen to
   *  EntryGroupChange events so that we can update the display if entries
   *  are added or deleted.
   **/
  public void entryGroupChanged (final EntryGroupChangeEvent event) {
    switch (event.getType ()) {
    case EntryGroupChangeEvent.ENTRY_ADDED:
    case EntryGroupChangeEvent.ENTRY_DELETED:
    case EntryGroupChangeEvent.ENTRY_INACTIVE:
    case EntryGroupChangeEvent.ENTRY_ACTIVE:
    case EntryGroupChangeEvent.NEW_DEFAULT_ENTRY:
      refreshMenu ();
      break;
    }
  }

  /**
   *  Implementation of the EntryChangeListener interface.
   **/
  public void entryChanged (final EntryChangeEvent event) {
    if (event.getType () == EntryChangeEvent.NAME_CHANGED) {
      refreshMenu ();
    }
  }

  /**
   *  Update the menus to the reflect the current contents of the EntryGroup.
   **/
  private void refreshMenu () {
    removeAll ();

    if (entry_group == null || entry_group.size () == 0) {
      add (new JMenuItem ("(No Entries Currently)"));
      return;
    }

    final JMenu set_entry_name_menu = new JMenu ("Set Name Of Entry");

    add (set_entry_name_menu);

    final JMenu set_default_menu = new JMenu ("Set Default Entry");
    final JMenu delete_entry_menu = new JMenu ("Remove An Entry");

    addSeparator ();

    add (set_default_menu);
    add (delete_entry_menu);

    final JMenuItem delete_active_entries_menu =
      new JMenuItem ("Remove Active Entries");

    delete_active_entries_menu.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent event) {
        deleteActiveEntries ();
      }
    });
    
    add (delete_active_entries_menu);

    final JMenuItem deactivate_all = new JMenuItem ("Deactivate All Entries");

    deactivate_all.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent event) {
        for (int i = 0 ; i < entry_group.size () ; ++i) {
          final Entry this_entry = entry_group.elementAt (i);

          entry_group.setIsActive (this_entry, false);
        }
      }
    });

    add (deactivate_all);

    entry_components = new Vector ();

    addSeparator ();

    for (int i = 0 ; i < entry_group.size () ; ++i) {
      final Entry this_entry = entry_group.elementAt (i);

      String entry_name = this_entry.getName ();

      if (entry_name == null) {
        entry_name = "no name";
      }

      final JMenuItem set_entry_name_item = new JMenuItem (entry_name);
      set_entry_name_item.addActionListener (new ActionListener () {
        public void actionPerformed (ActionEvent event) {
         setEntryName (this_entry);
        }
      });

      set_entry_name_menu.add (set_entry_name_item);

      final JMenuItem delete_entry_item = new JMenuItem (entry_name);
      delete_entry_item.addActionListener (new ActionListener () {
        public void actionPerformed (ActionEvent event) {
          if (this_entry.hasUnsavedChanges ()) {
            final String message;

            if (this_entry.getName () != null) {
              message = "there are unsaved changes in " +
                this_entry.getName () + " - really remove?";
            } else {
              message = "there are unsaved changes in entry #" +
                entry_group.indexOf (this_entry) + " - really remove?";
            }

            final YesNoDialog yes_no_dialog =
              new YesNoDialog (frame, message);

            if (!yes_no_dialog.getResult ()) {
              return;
            }
          }
          entry_group.remove (this_entry);
        }
      });

      delete_entry_menu.add (delete_entry_item);

      final JMenuItem set_default_entry_item = new JMenuItem (entry_name);
      set_default_entry_item.addActionListener (new ActionListener () {
        public void actionPerformed (ActionEvent event) {
          entry_group.setDefaultEntry (this_entry);
        }
      });

      set_default_menu.add (set_default_entry_item);

      add (entry_group.elementAt (i));
    }
  }

  /**
   *  Add a Label or Checkbox for the given Entry to this component.
   **/
  private void add (final Entry entry) {
    String entry_name = entry.getName ();

    if (entry_name == null) {
      entry_name = "no name";
    }

    if (entry_group.getDefaultEntry () == entry) {
      entry_name = entry_name + "  (default entry)";
    }

    final JCheckBoxMenuItem new_component =
      new JCheckBoxMenuItem (entry_name, entry_group.isActive (entry));

    new_component.addItemListener (new ItemListener () {
      public void itemStateChanged (ItemEvent event) {
        final int button_index =
          entry_components.indexOf (event.getSource ());

        if (event.getStateChange () == ItemEvent.SELECTED) {
          entry_group.setIsActive (button_index, true);
        } else {
          entry_group.setIsActive (button_index, false);
        }
      }
    });

    entry_components.addElement (new_component);
    add (new_component);
  }

  /**
   *  Create a new TextRequester component and set the name of the default
   *  Entry to whatever the user types.
   **/
  private void setEntryName (final Entry entry) {
    final TextRequester text_requester =
      new TextRequester ("New name for the entry?", 18, "");

    text_requester.addTextRequesterListener (new TextRequesterListener () {
      public void actionPerformed (final TextRequesterEvent event) {
        if (event.getType () == TextRequesterEvent.CANCEL) {
          return;
        }

        final String requester_text = event.getRequesterText ().trim ();
        if (requester_text.length () > 0) {
          if (entry.setName (requester_text)) {
            // it worked
          } else {
            new MessageDialog (frame,
                               "could not set the name of the default entry");
          }
        }
      }
    });

    text_requester.setVisible(true);
  }

  /**
   *  Delete the active entries after asking the user.
   **/
  private void deleteActiveEntries () {
    if (Options.getOptions ().isNoddyMode ()) {
      final YesNoDialog dialog =
        new YesNoDialog (frame,
                         "Are you sure you want to remove the " +
                         "active entries?");
          
      if (!dialog.getResult ()) {
        return;
      }
    }
        
    for (int i = entry_group.size () - 1 ; i >= 0 ; --i) {
      final Entry this_entry = entry_group.elementAt (i);
          
      if (this_entry.hasUnsavedChanges ()) {
        final String message;
            
        if (this_entry.getName () != null) {
          message = "there are unsaved changes in " +
            this_entry.getName () + " - really remove?";
        } else {
          message = "there are unsaved changes in entry #" +
            (i + 1) + " - really remove?";
        }
            
        final YesNoDialog yes_no_dialog =
          new YesNoDialog (frame, message);
            
        if (!yes_no_dialog.getResult ()) {
          continue;
        }
      }
          
      if (entry_group.isActive (this_entry)) {
        entry_group.remove (this_entry);
      }
    }
  }

  /**
   *  A vector containing the Entry objects that this EntryEdit object knows
   *  about.  This reference is obtained from owning_component.
   **/
  private EntryGroup entry_group;

  /**
   *  A vector containing one Checkbox or Label for each Entry in the
   *  EntryGroup object.
   **/
  private Vector entry_components = new Vector ();

  /**
   *  The JFrame reference that was passed to the constructor.
   **/
  private JFrame frame = null;
}

