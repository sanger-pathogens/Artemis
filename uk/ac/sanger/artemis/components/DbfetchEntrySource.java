/* DbfetchEntrySource.java
 *
 * created: Fri Nov 28 2003
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2003  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/DbfetchEntrySource.java,v 1.3 2004-12-16 10:44:33 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.sequence.*;

import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.util.Document;
import uk.ac.sanger.artemis.util.DocumentFactory;
import uk.ac.sanger.artemis.io.DocumentEntryFactory;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.EntryInformationException;
import uk.ac.sanger.artemis.io.SimpleEntryInformation;

import java.io.*;
import java.util.regex.Pattern;

import javax.swing.*;

/**
 *  This is an EntrySource that reads Entry objects from the EMBL Dbfetch
 *  server.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: DbfetchEntrySource.java,v 1.3 2004-12-16 10:44:33 tjc Exp $
 **/

public class DbfetchEntrySource
    implements EntrySource 
{
  private static Pattern REFSEQ_PATTERN = Pattern.compile("^[a-zA-Z]{2}?_\\w.*$");
  /**
   *  Create a new DbfetchEntrySource.
   *  @param frame The component that created this EntrySource.  (Used for
   *    requesters.)
   **/
  public DbfetchEntrySource (final JFrame frame) 
  {
  }

  /**
   *  Get an Entry object from the Ensembl Dbfetch server.
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
    final TextDialog text_dialog =
      new TextDialog (getFrame (), "Enter an accession number:", 10, "");

    final String text = text_dialog.getText ();

    if (text == null) {
      // user cancel
      return null;
    }

    final String embl_id = text.trim ();

    if(embl_id.length () > 0) 
    {
      final LogReadListener read_event_logger = new LogReadListener (embl_id);

      final EntryInformation entry_information =
        new SimpleEntryInformation(Options.getArtemisEntryInformation ());

//    final MessageDialog message_frame =
//      new MessageDialog (getFrame (),
//                         "reading entry - please wait", false);

      String db = "EMBL";
      if(REFSEQ_PATTERN.matcher(embl_id).matches())
        db = "refseq";
      
      final String url_string =
    		  "https://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db="+db+"&id=" + embl_id +"&format=fasta&style=raw";

      final Document url_document =
        DocumentFactory.makeDocument(url_string);

      try
      {
        final uk.ac.sanger.artemis.io.Entry new_embl_entry =
          DocumentEntryFactory.makeDocumentEntry(entry_information,
                                                 url_document,
                                                 read_event_logger);

        if (read_event_logger.seenMessage()) {
          final YesNoDialog yes_no_dialog =
            new YesNoDialog (frame,
                             "there were warnings while reading - view now?");

          if (yes_no_dialog.getResult ()) {
            Splash.showLog ();
          }
        }

        final Bases real_bases;

        if (bases == null) {
          if (new_embl_entry.getSequence () == null) {
            final String message =
              "the entry contains no sequence: " + embl_id;
            new MessageDialog (getFrame (), message);
            return null;
          }

          real_bases = new Bases (new_embl_entry.getSequence ());
        } else {
          real_bases = bases;
        }

        return new Entry (real_bases, new_embl_entry);
      }
      catch (EntryInformationException e) 
      {
        throw new Error ("internal error - unexpected exception: " + e);
      } 
//    finally
//    {
//      message_frame.dispose ();
//    }
    }

    return null;
  }

  /**
   *  Get an Entry object from the Ensembl Dbfetch server.
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
      throws OutOfRangeException, NoSequenceException, IOException 
  {
    return getEntry (null, show_progress);
  }

  /**
   *  Return the name of this source (for display to the user in menus).
   **/
  public String getSourceName () {
    return "EBI - Dbfetch";
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
   *  Return the JFrame that was passed to the constructor.
   **/
  public JFrame getFrame () {
    return frame;
  }

  /**
   *  The JFrame that was passed to the constructor.
   **/
  private JFrame frame = null;
}
