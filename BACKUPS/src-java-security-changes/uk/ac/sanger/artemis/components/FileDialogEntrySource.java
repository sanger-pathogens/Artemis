/* FileDialogEntrySource.java
 *
 * created: Thu Jun  8 2000
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/FileDialogEntrySource.java,v 1.2 2004-12-03 17:47:04 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.sequence.*;
import uk.ac.sanger.artemis.util.*;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.SimpleEntryInformation;

import javax.swing.JFrame;

/**
 *  This is an EntrySource that reads Entry objects from the local filesystem.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: FileDialogEntrySource.java,v 1.2 2004-12-03 17:47:04 tjc Exp $
 **/

public class FileDialogEntrySource
    implements EntrySource 
{
   
   /** 
   *  The component that created this FileEntrySource.  (Used for requesters.)
   **/
  final JFrame frame;

  /** InputStreamProgressEvents are sent to this object. */
  private final InputStreamProgressListener stream_progress_listener;

  /**
   *  Create a new FileDialogEntrySource.
   *  @param frame The component that created this FileDialogEntrySource.
   *    (Used for requesters.)
   *  @param listener InputStreamProgressEvent objects will be sent to this
   *    listener as progress on reading is made.
   **/
  public FileDialogEntrySource (final JFrame frame,
                                final InputStreamProgressListener listener) 
  {
    this.frame = frame;
    this.stream_progress_listener = listener;
  }

  /**
   *  Get an Entry object from this source (by reading from a file).
   *  @param bases The Bases object to pass to the Entry constructor.
   *  @exception OutOfRangeException Thrown if one of the features in
   *    the Entry is out of range of the Bases object.
   *  @param show_progress If true a InputStreamProgressDialog will be shown
   *    while loading.
   *  @return null if and only if the user cancels the read.
   **/
  public Entry getEntry(final Bases bases,
                        final boolean show_progress)
      throws OutOfRangeException
  {
    try
    {
      return getEntryInternal(bases, show_progress);
    }
    catch (NoSequenceException e)
    {
      throw new Error ("internal error - unexpected exception: " + e);
    }
  }

  /**
   *  Get an Entry object from this source (by reading from a file).
   *  @exception OutOfRangeException Thrown if one of the features in
   *    the Entry is out of range of the Bases object.
   *  @exception NoSequenceException Thrown if the entry that we read has no
   *    sequence.
   *  @param show_progress If true a InputStreamProgressDialog will be shown
   *    while loading.
   *  @return null if and only if the user cancels the read.
   **/
  public Entry getEntry (final boolean show_progress)
      throws OutOfRangeException, NoSequenceException 
  {
    return getEntryInternal (null, show_progress);
  }

  /**
   *  Returns true if and only if this EntrySource always returns "full"
   *  entries.  ie. entries that contain features and sequence.  Entries that
   *  are read from a file may contain just features so in this class this
   *  method returns false.
   **/
  public boolean isFullEntrySource () 
  {
    return false;
  }
  

////  change ArtemisMain.readArgsAndOptions () to use this?:

//    /**
//     *  Get an Entry object from this source by name (by reading from a file).
//     *  @exception OutOfRangeException Thrown if one of the features in
//     *    embl_entry is out of range of the Bases object.
//     *  @exception NoSequenceException Thrown if the entry that we read has no
//     *    sequence.
//     *  @return null if and only if there is no Entry with that name.
//     **/
//    public Entry getEntryByName (final String entry_file_name)
//        throws OutOfRangeException, NoSequenceException, IOException {
//      final Document new_document =
//        new FileProgressDocument (new File (entry_file_name),
//                                  getInputStreamProgressListener ());

//      final EntryInformation new_entry_information =
//        new SimpleEntryInformation (Options.getArtemisEntryInformation ());

//      boolean seen_error = false;

//      while (true) {
//        try {
//          final uk.ac.sanger.artemis.io.Entry new_entry =
//            DocumentEntryFactory.makeDocumentEntry (new_entry_information,
//                                                    new_document);

//          return new Entry (new_entry);
//        } catch (EntryInformationException e) {

//          if (!seen_error) {
//            final String message =
//              "warning while reading " + entry_file_name + " - " +
//              e.getMessage ();

//            System.err.println (message);

//            new MessageDialog (frame, message);

//            seen_error = true;
//          }

//          EntryFileDialog.handleOpenException (new_entry_information, e);

//          // go around the loop again
//        }
//      }
//    }

  /**
   *  Make a new Entry.
   *  @param bases The Bases object to pass to the Entry constructor.  If null
   *    a new Bases object will be created.
   *  @exception OutOfRangeException Thrown if one of the features in
   *    the Entry is out of range of the Bases object.
   *  @exception NoSequenceException Thrown if bases is null and the Entry that
   *    we read has no sequence.
   **/
  private Entry makeEntry (final Bases bases,
                           final uk.ac.sanger.artemis.io.Entry embl_entry)
      throws OutOfRangeException, NoSequenceException 
  {
    if (bases == null) 
      return new Entry (embl_entry);
    else 
      return new Entry (bases, embl_entry);
  }

  /**
   *  Return the InputStreamProgressListener that was passed to the
   *  constructor.
   **/
  public InputStreamProgressListener getInputStreamProgressListener () 
  {
    return stream_progress_listener;
  }

  /**
   *  Return the name of this source (for display to the user in menus).
   **/
  public String getSourceName () 
  {
    return "Filesystem";
  }

  /**
   *  Implementation of getEntry ().
   *  @param bases The Bases object to pass to the Entry constructor.
   *  @exception OutOfRangeException Thrown if one of the features in
   *    Entry is out of range of the Bases object.
   *  @exception NoSequenceException Thrown if bases is null and the entry that
   *    we read has no sequence.
   *  @param show_progress If true a InputStreamProgressDialog will be shown
   *    while loading.
   *  @return null if and only if the user cancels the read.
   **/
  private Entry getEntryInternal(final Bases bases,
                                 final boolean show_progress)
      throws OutOfRangeException, NoSequenceException 
  {
    final EntryInformation new_entry_information =
      new SimpleEntryInformation(Options.getArtemisEntryInformation());

    EntryFileDialog dialog;

    if(bases == null) 
      dialog = new EntryFileDialog(frame, true);
    else 
      dialog = new EntryFileDialog(frame, false);


    uk.ac.sanger.artemis.io.Entry new_embl_entry =
      dialog.getEntry(new_entry_information, stream_progress_listener,
                      show_progress);

    if(new_embl_entry == null) 
      return null;

    return makeEntry(bases, new_embl_entry);
  }

}

