/* BioJavaEntrySource.java
 *
 * created: Tue Apr 10 2001
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/BioJavaEntrySource.java,v 1.1 2004-06-09 09:46:06 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import java.io.*;

import org.biojava.bio.seq.io.EmblLikeFormat;

import uk.ac.sanger.artemis.util.*;
import uk.ac.sanger.artemis.io.BioJavaEntry;
import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.sequence.*;

/**
 *  This is an EntrySource that reads Entry objects from a BioJava Sequence
 *  object.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 **/

public class BioJavaEntrySource implements EntrySource 
{
  /**
   *  Create a new BioJavaEntrySource.
   **/
  public BioJavaEntrySource () {

  }

  /**
   *  Get an Entry object from this source (by reading from a file or whatever).
   *  @param bases The Bases object to pass to the Entry constructor.
   *  @param show_progress If true a InputStreamProgressDialog will be shown
   *    while loading.  (Not implemented)
   *  @exception OutOfRangeException Thrown if one of the features in
   *    the Entry is out of range of the Bases object.
   *  @return null if and only if the read is cancelled by the user or if the
   *    read fails.
   **/
    public Entry getEntry (final Bases bases, final ProgressThread progress_thread,
                           final boolean show_progress)
        throws OutOfRangeException, IOException
    {
      final String fileName = "/nfs/team81/kmr/pow/java2/AB000095.embl";
      final FileDocument document = new FileDocument (new File (fileName));

      final BioJavaEntry emblEntry =
        new BioJavaEntry (document, new EmblLikeFormat ());

      return new Entry (bases, emblEntry);
    }

  /**
   *  Get an Entry object from this source (by reading from a file or whatever).
   *  @param bases The Bases object to pass to the Entry constructor.
   *  @param show_progress If true a InputStreamProgressDialog will be shown
   *    while loading.  (Not implemented)
   *  @exception OutOfRangeException Thrown if one of the features in
   *    the Entry is out of range of the Bases object.
   *  @return null if and only if the read is cancelled by the user or if the
   *    read fails.
   **/
    public Entry getEntry (final Bases bases, final boolean show_progress)
        throws OutOfRangeException, IOException
    {
      return getEntry(bases, null, show_progress);
    }

  /**
   *  Get an Entry object from this source (by reading from a file or whatever).  
   *  A Bases object will be created for the sequence of the new Entry.
   *  @param show_progress If true a InputStreamProgressDialog will be shown
   *    while loading.  (Not implemented)
   *  @exception OutOfRangeException Thrown if one of the features in
   *    the Entry is out of range of the Bases object.
   *  @exception NoSequenceException Thrown if the entry that we read has no
   *    sequence.
   *  @return null if and only if the user cancels the read or if the read
   *    fails.
   **/
    public Entry getEntry(final boolean show_progress, 
                          final ProgressThread progress_thread)
        throws OutOfRangeException, NoSequenceException, IOException
    {
      return getEntry(show_progress);
    }

  /**
   *  Get an Entry object from this source (by reading from a file or whatever).  
   *  A Bases object will be created for the sequence of the new Entry.
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
      final String fileName = "/nfs/team81/kmr/pow/java2/AE002734.game";
      final FileDocument document = new FileDocument (new File (fileName));
      
      final BioJavaEntry emblEntry =
        new BioJavaEntry (document, new EmblLikeFormat ());
      
      return new Entry (emblEntry);
    }

  /**
   *  Returns true if and only if this EntrySource always returns "full"
   *  entries.  ie. entries that contain features and sequence.
   **/
  public boolean isFullEntrySource () {
    return true;
  }

  /**
   *  Return the name of this source (for display to the user in menus).
   **/
  public String getSourceName () {
    return "BioJava";
  }
}
