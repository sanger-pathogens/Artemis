/* EMBLCorbaEntrySource.java
 *
 * created: Wed Jun  7 2000
 *
 * This file is part of Artemis
 *
 * Copyright (C) 1999,2000  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/EMBLCorbaEntrySource.java,v 1.1 2004-06-09 09:46:17 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.sequence.*;

import uk.ac.sanger.artemis.util.*;
import uk.ac.sanger.artemis.io.LocationParseException;
import uk.ac.sanger.artemis.io.InvalidKeyException;
import uk.ac.sanger.artemis.io.EntryInformationException;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.SimpleEntryInformation;

import nsdb.EmblSeq;
import nsdb.EmblPackage.Superceded;
import type.NoResult;

import java.io.*;
import java.awt.*;
import java.net.*;

import javax.swing.*;

/**
 *  This is an EntrySource that reads Entry objects from the EMBL CORBA
 *  server.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: EMBLCorbaEntrySource.java,v 1.1 2004-06-09 09:46:17 tjc Exp $
 **/

public class EMBLCorbaEntrySource
    extends CorbaEntrySource
    implements EntrySource 
{
  /**
   *  Create a new EMBLCorbaEntrySource from the given String.
   *  @param frame The component that created this EntrySource.  (Used for
   *    requesters.)
   *  @param ior_url_string A String containing the URL of the IOR for the
   *    server.
   **/
  public EMBLCorbaEntrySource (final JFrame frame,
                               final String ior_url_string)
      throws MalformedURLException//, IOException
  {
    super (frame, ior_url_string);
  }

  /**
   *  Get an Entry object from the Ensembl CORBA server.
   *  @param bases The Bases object to pass to the Entry constructor.
   *  @param progress_thread Progress thread to monitor entry reading.
   *  @param show_progress If true a InputStreamProgressDialog will be shown
   *    while loading.  (Not implemented)
   *  @exception OutOfRangeException Thrown if one of the features in
   *    the Entry is out of range of the Bases object.
   *  @return null if and only if the user cancels the read or if the read
   *    fails.
   **/
  public Entry getEntry (final Bases bases, final ProgressThread progress_thread,
                         final boolean show_progress)
      throws OutOfRangeException, IOException
  {
    return makeCorbaDialog (bases, false);
  }

  /**
   *  Get an Entry object from the Ensembl CORBA server.
   *  @param bases The Bases object to pass to the Entry constructor.
   *  @param show_progress If true a InputStreamProgressDialog will be shown
   *    while loading.  (Not implemented)
   *  @exception OutOfRangeException Thrown if one of the features in
   *    the Entry is out of range of the Bases object.
   *  @return null if and only if the user cancels the read or if the read
   *    fails.
   **/
  public Entry getEntry (final Bases bases, final boolean show_progress)
      throws OutOfRangeException, IOException 
  {
    return makeCorbaDialog (bases, false);
  }

   /**
   *  Get an Entry object from the Ensembl CORBA server.
   *  @param show_progress If true a InputStreamProgressDialog will be shown
   *    while loading.  (Not implemented)
   *  @param progress_thread Progress thread to monitor entry reading.
   *  @exception OutOfRangeException Thrown if one of the features in
   *    the Entry is out of range of the Bases object.
   *  @return null if and only if the user cancels the read or if the read
   *    fails.
   **/
  public Entry getEntry (final boolean show_progress, 
                         final ProgressThread progress_thread)
      throws OutOfRangeException, IOException 
  {
    return makeCorbaDialog (null, false);
  }


  /**
   *  Get an Entry object from the Ensembl CORBA server.
   *  @param show_progress If true a InputStreamProgressDialog will be shown
   *    while loading.  (Not implemented)
   *  @exception OutOfRangeException Thrown if one of the features in
   *    the Entry is out of range of the Bases object.
   *  @exception NoSequenceException Thrown if the entry that we read has no
   *    sequence.
   *  @return null if and only if the user cancels the read or if the read
   *    fails.
   **/
  public Entry getEntry (final boolean show_progress)
      throws OutOfRangeException, NoSequenceException, IOException {
    return makeCorbaDialog (null, false);
  }

  /**
   *  Return the name of this source (for display to the user in menus).
   **/
  public String getSourceName () {
    return "EBI - CORBA";
  }

  /**
   *  Returns true if and only if this EntrySource always returns "full"
   *  entries.  ie. entries that contain features and sequence.  Entries that
   *  are read from EMBL always contain sequence so in this class this method
   *  returns false.
   **/
  public boolean isFullEntrySource () {
    return true;
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

    final TextDialog text_dialog =
      new TextDialog (getFrame (), "Enter an accession number:", 10, "");

    final String text = text_dialog.getText ();

    if (text == null) {
      // user cancel
      return null;
    } else {
      final String corba_id = text.trim ();

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
  }

  /**
   *  Given an accession number and the handle of an EMBL corba server, this
   *  method will ask the user (using a TextRequester) for the id of a entry
   *  in the server and will then attempt to get it.
   *  @param bases If this is null a new Bases object will be created for the
   *    Entry once it has been read from the server.  If not null then it will
   *    be passed to the Entry constructor.
   *  @param corba_handle The handle of the nsdb.Embl object from which we
   *    will read the entry.
   *  @param corba_id The id of the entry in the database
   *  @param read_only true if and only if a read-only Entry should be created
   *    (some are always read only).
   *  @exception OutOfRangeException Thrown if one of the features in
   *    the Entry is out of range of the Bases object.
   **/
  protected Entry makeFromCorbaID (final Bases bases,
                                   final String corba_id,
                                   final boolean read_only)
      throws OutOfRangeException, IOException {
    try {
      final nsdb.Embl corba_handle = getServerHandle ();

      final uk.ac.sanger.artemis.io.Entry new_embl_entry;

      final nsdb.EmblWriter embl_writer =
        nsdb.EmblWriterHelper.narrow (corba_handle);

      final EntryInformation entry_information =
        new SimpleEntryInformation (Options.getArtemisEntryInformation ());

      if (read_only || embl_writer == null) {
   // first try to make a plain EmblSeq object
        final EmblSeq embl_seq = corba_handle.getEmblSeq (corba_id);

        new_embl_entry =
          new uk.ac.sanger.artemis.io.CorbaEntry (entry_information,
                                                      embl_seq);
      } else {
        // make a read-write object
        final nsdb.EmblSeqWriter embl_write_seq =
          embl_writer.getEmblSeqWriter (corba_id);

        new_embl_entry =
          new uk.ac.sanger.artemis.io.RWCorbaEntry (entry_information,
                                                        embl_write_seq);
      }

      final Bases real_bases;

      if (bases == null) {
        if (new_embl_entry.getSequence () == null) {
          final String message =
            "the entry contains no sequence: " + corba_id;
          new MessageDialog (getFrame (), message);
          return null; 
       }

        real_bases = new Bases (new_embl_entry.getSequence ());
      } else {
        real_bases = bases;
      }

      return new Entry (real_bases, new_embl_entry);
    } catch (NoResult e) {
      final String message =
        "Database query failed (no result) while getting id: " + corba_id;
      new MessageDialog (getFrame (), message);
    } catch (Superceded e) {
      //  Superceded is thrown by getEmblSeq method if accession number
      //  doesn't exist anymore because it was merged or split
      final String message =
        "This accession number has been superceded: " + corba_id;
      new MessageDialog (getFrame (), message);
    } catch (LocationParseException e) {
      final String message =
        "Unexpected error while accessing " + corba_id + ": " + e;
      new MessageDialog (getFrame (), message);
    } catch (InvalidKeyException e) {
      final String message =
        "Unexpected error while accessing " + corba_id + ": " + e;
      new MessageDialog (getFrame (), message);
    } catch (org.omg.CORBA.OBJECT_NOT_EXIST e) {
      final String message =
        "the object you requested (" + corba_id + ") does not exist";
      new MessageDialog (getFrame (), message);
    } catch (org.omg.CORBA.COMM_FAILURE e) {
      final String message =
        "Failed to get an object from Corba: " + e;
      new MessageDialog (getFrame (), message);
    } catch (EntryInformationException e) {
      final String message =
        "Failed to get an object from Corba: " + e;
      new MessageDialog (getFrame (), message);
    }

    return null;
  }

  /**
   *  Return the handle of the EMBL database.  The database IOR is found by
   *  calling getIOR ().
   **/
  protected nsdb.Embl getServerHandle ()
      throws org.omg.CORBA.COMM_FAILURE, IOException {
    final org.omg.CORBA.ORB orb =
      org.omg.CORBA.ORB.init (new String [0], new java.util.Properties());

    final org.omg.CORBA.Object obj;

    obj = orb.string_to_object (getIOR ());

    return nsdb.EmblHelper.narrow (obj);
  }
}
