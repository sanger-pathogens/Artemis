/* EntryHeaderEdit.java
 *
 * created: Wed Jun 16 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/EntryHeaderEdit.java,v 1.2 2007-03-23 14:15:54 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;

import java.awt.event.*;
import javax.swing.*;
import javax.swing.event.*;

/**
 *  Objects of this class are used to edit the header of an Entry.
 *
 *  @author Kim Rutherford
 *  @version $Id: EntryHeaderEdit.java,v 1.2 2007-03-23 14:15:54 tjc Exp $
 **/

public class EntryHeaderEdit
    implements EntryChangeListener, EntryGroupChangeListener,
               DocumentListener {
  /**
   *  Create a new EntryHeaderEdit object for the given Entry.
   **/
  public EntryHeaderEdit (final EntryGroup entry_group,
                          final Entry edit_entry) {
    
    this.edit_entry = edit_entry;
    this.entry_group = entry_group;
    
    file_viewer = new FileViewer ("Artemis Entry Header Editor: " +
                                  edit_entry.getName () == null ?
                                  "" : edit_entry.getName ());

    file_viewer.getContentPane ().add (error_text, "North");

    readHeader ();

    file_viewer.getTextPane().setEditable (true);
    file_viewer.getTextPane().getDocument ().addDocumentListener (this);

    getEntry ().addEntryChangeListener (this);
    entry_group.addEntryGroupChangeListener (this);

    file_viewer.addWindowListener (new WindowAdapter () {
      public void windowClosed (WindowEvent event) {
        stopListening ();
      }
    });
  }


  /**
   *  Remove this object as a feature and entry change listener.
   **/
  public void stopListening () {
    getEntry ().removeEntryChangeListener (this);
    entry_group.removeEntryGroupChangeListener (this);
  }

  /**
   *  Implementation of the EntryChangeListener interface.  We listen to
   *  EntryChange events so we can update this component if the header changes
   *  for another reason.
   **/
  public void entryChanged (EntryChangeEvent event) {
    switch (event.getType ()) {
    case EntryChangeEvent.HEADER_CHANGED:
      
      if (event.getSource () == current_text) {
        // don't bother with events from us
        return;
      }
    
      // re-read the information from the entry
      readHeader ();
      break;
    default:
      // do nothing
      break;
    }
  }

  /**
   *  Implementation of the EntryGroupChangeListener interface.  We listen to
   *  EntryGroupChange events so we can notify the user if of this component
   *  if the entry gets deleted.
   **/
  public void entryGroupChanged (final EntryGroupChangeEvent event) {
    switch (event.getType ()) {
    case EntryGroupChangeEvent.ENTRY_DELETED:
      if (event.getEntry () == edit_entry) {
        stopListening ();
        file_viewer.dispose ();
      }
      break;
    default:
      // do nothing
      break;
    }
  }

  public void changedUpdate(DocumentEvent e) {
    textValueChanged ();
  }
  public void insertUpdate(DocumentEvent e) {
    textValueChanged ();
  }
  public void removeUpdate(DocumentEvent e) {
    textValueChanged ();
  }
  
  /**
   *  Implementation of the TextListener interface.  When the text changes we
   *  update the Feature object that we are showing.
   **/
  public void textValueChanged () {
    if (current_text.equals (file_viewer.getText ())) {
      // the text hasn't really changed
      return;
    }
    
    current_text = file_viewer.getText ();

    try {
      // ignore text change events while reading
      edit_entry.setHeaderText (current_text);
      error_text.setText ("");
    } catch (uk.ac.sanger.artemis.io.ReadFormatException e) {
      error_text.setText (e + (e.getLineNumber () > 1 ?
                               " at line: " + e.getLineNumber () :
                               ""));
    } catch (java.io.IOException e) {
      error_text.setText (e.toString ());
    }
  }

  /**
   *  Read the header of edit_entry into this component.
   **/
  public void readHeader () {
    final String header = edit_entry.getHeaderText ();

    if (header != null) {
      file_viewer.setText (header);

      current_text = file_viewer.getText ();
    }
  }

  /**
   *  Return the Entry that this object is displaying.
   **/
  private Entry getEntry () {
    return edit_entry;
  }

  /**
   *  The EntryGroup that passed to the constructor.
   **/
  private EntryGroup entry_group;

  /**
   *  The Entry that contains the Feature this object is displaying.
   **/
  private Entry edit_entry;

  /**
   *  The FileViewer object that is displaying the feature.
   **/
  private FileViewer file_viewer;

  /**
   *  The text version of the Feature we are currently displaying.
   **/
  private String current_text = "";

  /**
   *  A Label for showing errors and messages.
   **/
  private JLabel error_text = new JLabel ("");
}

