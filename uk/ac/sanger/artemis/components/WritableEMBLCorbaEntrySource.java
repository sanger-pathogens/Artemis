/* WritableEMBLCorbaEntrySource.java
 *
 * created: Mon Jun 12 2000
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2000,2002x  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/WritableEMBLCorbaEntrySource.java,v 1.1 2004-06-09 09:47:59 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.sequence.*;

import uk.ac.sanger.artemis.util.*;
import uk.ac.sanger.artemis.io.LocationParseException;
import uk.ac.sanger.artemis.io.InvalidKeyException;

import nsdb.ServerInfo;
import nsdb.EntryStats;
import nsdb.EmblSeq;
import nsdb.EmblPackage.Superceded;
import type.NoResult;

import java.io.*;
import java.awt.*;
import java.net.*;

import javax.swing.*;

/**
 *  This is an EntrySource that can get writable Entry objects from an EMBL
 *  CORBA server.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: WritableEMBLCorbaEntrySource.java,v 1.1 2004-06-09 09:47:59 tjc Exp $
 **/

public class WritableEMBLCorbaEntrySource extends EMBLCorbaEntrySource {
  /**
   *  Create a new WritableEMBLCorbaEntrySource from the given String.
   *  @param frame The component that created this EntrySource.  (Used for
   *    requesters.)
   *  @param ior_url_string A String containing the URL of the IOR for the
   *    server.
   **/
  public WritableEMBLCorbaEntrySource (final JFrame frame,
                                       final String ior_url_string)
      throws MalformedURLException {
    super (frame, ior_url_string);
  }

  /**
   *  Get an Entry object from the Ensembl CORBA server.
   *  @param bases The Bases object to pass to the Entry constructor.
   *  @exception OutOfRangeException Thrown if one of the features in
   *    the Entry is out of range of the Bases object.
   *  @param show_progress If true a InputStreamProgressDialog will be shown
   *    while loading.
   *  @return null if and only if the user cancels the read or if the read
   *    fails.
   **/
  public Entry getEntry (final Bases bases, final boolean show_progress)
      throws OutOfRangeException, IOException {
    return makeCorbaDialog (bases, false);
  }

  /**
   *  Get an Entry object from the Ensembl CORBA server.
   *  @exception OutOfRangeException Thrown if one of the features in
   *    the Entry is out of range of the Bases object.
   *  @exception NoSequenceException Thrown if the entry that we read has no
   *    sequence.
   *  @param show_progress If true a InputStreamProgressDialog will be shown
   *    while loading.
   *  @return null if and only if the user cancels the read or if the read
   *    fails.
   **/
  public Entry getEntry (final boolean show_progress)
      throws OutOfRangeException, NoSequenceException, IOException {
    return makeCorbaDialog (null, false);
  }

  /**
   *  Create a TextRequester, wait for the user to type an accession number
   *  and then read that entry from the EMBL CORBA server.
   *  @param bases If this is null a new Bases object will be created for the
   *    Entry once it has been read from the server.  If not null then it will
   *    be passed to the Entry constructor.
   *  @param read_only true if and only if a read-only Entry should be created
   *    (some are always read only).
   *  @return null if and only if the user cancels the read or if the read
   *    fails.
   *  @exception OutOfRangeException Thrown if one of the features in
   *    the Entry is out of range of the Bases object.
   **/
  protected Entry makeCorbaDialog (final Bases bases,
                                   final boolean read_only)
      throws OutOfRangeException, IOException {
    final org.omg.CORBA.ORB orb =
      org.omg.CORBA.ORB.init (new String [0], new java.util.Properties());

    final nsdb.Embl corba_handle = getServerHandle ();

    final nsdb.EmblWriter embl_writer =
      nsdb.EmblWriterHelper.narrow (corba_handle);

    if (embl_writer == null) {
      final String message = "Server reference is not an EmblWriter: " +
        corba_handle;
      new MessageDialog (getFrame (), message);
      return null;
    }

    final ListDialog list_dialog =
      new ListDialog (getFrame (), "Select an entry");

    final ServerInfo stats = embl_writer.getServerInfo ();

    final EntryStats [] entry_stats_list = stats.entry_stats_list;
    final EntryStats [] file_stats_list = stats.file_stats_list;

    list_dialog.getList ().setModel (new AbstractListModel () {
      public int getSize () {
        return entry_stats_list.length + file_stats_list.length;
      }
      public Object getElementAt (int i) {
        if (i < entry_stats_list.length) {
          return makeListString (entry_stats_list[i], false);
        } else {
          return makeListString (file_stats_list[i-entry_stats_list.length],
                                 true);
        }
      }
    });

    // returns when user hits ok or cancel
    final String selected_entry_string =
      (String) list_dialog.getSelectedValue ();

    if (selected_entry_string == null) {
      // user hit cancel
      return null;
    }

    final int end_of_entry_name = selected_entry_string.indexOf ("   ");

    final String corba_id =
      selected_entry_string.substring (0, end_of_entry_name);

    if (corba_id.length () > 0) {
      final MessageDialog message_frame =
        new MessageDialog (getFrame (),
                           "reading entry - please wait", false);

      try {
        return makeFromCorbaID (bases, corba_id, read_only);
      } finally {
        message_frame.dispose ();
      }
    } else {
      return null;
    }
  }

  /**
   *  Make a ListDialog then let the user choose an entry.
   **/
  private String makeListString (final EntryStats this_entry_stats,
                                 final boolean is_file_entry) {
    final StringBuffer buffer = new StringBuffer ();

    buffer.append (this_entry_stats.name).append ("   ");

    for (int j = 0 ; j < (30 - this_entry_stats.name.length ()) ; ++j) {
      buffer.append (" ");
    }

    final long last_change_time = this_entry_stats.last_change_time;

    if (is_file_entry) {
      buffer.append ("(file)");
    } else {
      if (last_change_time > 0) {
        buffer.append (new java.util.Date (last_change_time * 1000L));
      } else {
        buffer.append ("(not saved)");
      }
    }

    return buffer.toString ();
  }

  /**
   *  Return the name of this source (for display to the user in menus).
   **/
  public String getSourceName () {
    return "Database";
  }

  /**
   *  Returns true if and only if this EntrySource always returns "full"
   *  entries.  ie. entries that contain features and sequence.  Entries that
   *  are read from a file may contain just features so in this class this
   *  method returns false.
   **/
  public boolean isFullEntrySource () {
    return false;
  }
}
