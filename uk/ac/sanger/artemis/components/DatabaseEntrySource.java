/* DatabaseEntrySource.java
 *
 * created: Mar 2005
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2001  Genome Research Limited
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
 */

package uk.ac.sanger.artemis.components;

import javax.swing.JOptionPane;
import java.net.*;
import java.io.*;
import java.util.*;

import uk.ac.sanger.artemis.util.*;
import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.sequence.*;
import uk.ac.sanger.artemis.io.DatabaseDocumentEntry;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.SimpleEntryInformation;
import uk.ac.sanger.artemis.io.InvalidKeyException;
import uk.ac.sanger.artemis.io.EntryInformationException;

/**
 *
 *  This is an EntrySource that reads Entry objects from a relational
 *  database.
 *
 **/

public class DatabaseEntrySource implements EntrySource 
{
  /**
   *  Create a new DatabaseEntrySource.
   **/
  public DatabaseEntrySource() 
  {
  }

  /**
   *  @param bases The Bases object to pass to the Entry constructor.
   *  @param show_progress If true a InputStreamProgressDialog will be shown
   *    while loading.  (Not implemented)
   *  @exception OutOfRangeException Thrown if one of the features in
   *    the Entry is out of range of the Bases object.
   *  @return null if and only if the read is cancelled by the user or if the
   *    read fails.
   **/
  public Entry getEntry(final Bases bases, final boolean show_progress)
      throws OutOfRangeException, IOException
  {
    return null;
  }

  public Entry getEntry(final boolean show_progress)
      throws OutOfRangeException, IOException
  {
    return null;
  }

  /**
   *  Get an Entry object from this source.
   *  @param id Feature ID to read in.
   *  @exception OutOfRangeException Thrown if one of the features in
   *    the Entry is out of range of the Bases object.
   *  @exception NoSequenceException Thrown if the entry that we read has no
   *    sequence.
   *  @return null if and only if the user cancels the read or if the read
   *    fails.
   **/
  public Entry getEntry(String id, InputStreamProgressListener progress_listener)
      throws OutOfRangeException, NoSequenceException, IOException
  {
    return makeFromID(null,id,false,progress_listener);
  }

  /**
   *  Returns true if and only if this EntrySource always returns "full"
   *  entries.  ie. entries that contain features and sequence.
   **/
  public boolean isFullEntrySource()
  {
    return true;
  }

  /**
   *  Return the name of this source (for display to the user in menus).
   **/
  public String getSourceName()
  {
    return "Database";
  }


  protected Hashtable getDatabaseEntries()
  {
    DatabaseDocument doc = new DatabaseDocument("jdbc:postgresql://pcs3:13001/chadoCVS?user=es2");
    return doc.getDatabaseEntries();
  }


  /**
   *
   *  Given an accession number and the handle of an EMBL corba server, this
   *  method will ask the user (using a TextRequester) for the id of a entry
   *  in the server and will then attempt to get it.
   *  @param bases If this is null a new Bases object will be created for the
   *    Entry once it has been read from the server.  If not null then it will
   *    be passed to the Entry constructor.
   *  @param corba_id The id of the entry in the database
   *  @param read_only true if and only if a read-only Entry should be created
   *    (some are always read only).
   *  @exception OutOfRangeException Thrown if one of the features in
   *    the Entry is out of range of the Bases object.
   **/
  protected Entry makeFromID(final Bases bases,
                             final String id,
                             final boolean read_only,
                             InputStreamProgressListener progress_listener)
      throws OutOfRangeException, IOException 
  {
    try 
    {
      DatabaseDocumentEntry db_entry = null;

      final EntryInformation entry_information =
        new SimpleEntryInformation(Options.getArtemisEntryInformation());

      if(read_only) 
      {
      } 
      else 
      {
        DatabaseDocument doc = new DatabaseDocument("jdbc:postgresql://pcs3:13001/chadoCVS?user=es2",
                                                    id, progress_listener);
        db_entry = new DatabaseDocumentEntry(entry_information, doc);
      }

      final Bases real_bases;

      if(bases == null)
      {
        if(db_entry.getSequence() == null)
        {
          JOptionPane.showMessageDialog(null,
                       "The selected entry contains no sequence: " + id,
                       "No Sequence",
                       JOptionPane.ERROR_MESSAGE);

          return null; 
        }

        real_bases = new Bases(db_entry.getSequence());
      } 
      else
        real_bases = bases;

      return new Entry(real_bases, db_entry);
    } 
    catch(InvalidKeyException e) 
    {
      JOptionPane.showMessageDialog(null,
                       "Unexpected error while accessing " + id + ": " + e,
                       "Invalid Key",
                       JOptionPane.ERROR_MESSAGE);
    }
    catch(EntryInformationException e) 
    {
      JOptionPane.showMessageDialog(null,
                       "Failed to get entry: " + e,
                       "Entry Information Exception",
                       JOptionPane.ERROR_MESSAGE);
    }

    return null;
  }

}
