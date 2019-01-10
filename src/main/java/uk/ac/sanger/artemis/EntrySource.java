/* EntrySource.java
 *
 * created: Wed Jun  7 2000
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/EntrySource.java,v 1.2 2004-12-03 17:47:04 tjc Exp $
 */

package uk.ac.sanger.artemis;

import uk.ac.sanger.artemis.sequence.*;
import uk.ac.sanger.artemis.util.*;

import java.io.*;

/**
 *  This interface is implemented by those objects that can produce an Entry
 *  object.  Examples include objects that can read an Entry from a file.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 **/

public interface EntrySource 
{
  /**
   *  Get an Entry object from this source (by reading from a file or whatever).
   *  @param bases The Bases object to pass to the Entry constructor.
   *  @exception OutOfRangeException Thrown if one of the features in 
   *    the Entry is out of range of the Bases object.
   *  @param show_progress If true show a Dialog showing the progress while
   *    loading (may be ignored).
   *  @return null if and only if the read is cancelled by the user or if the
   *    read fails.
   **/
  Entry getEntry(final Bases bases,
                 final boolean show_progress)
      throws OutOfRangeException, IOException;

  /**
   *  Get an Entry object from this source (by reading from a file or whatever).  
   *  A Bases object will be created for the sequence of the new Entry.
   *  @exception OutOfRangeException Thrown if one of the features in
   *    the Entry is out of range of the Bases object.
   *  @exception NoSequenceException Thrown if the entry that we read has no
   *    sequence.
   *  @param show_progress If true show a Dialog showing the progress while
   *    loading (may be ignored).
   *  @return null if and only if the user cancels the read or if the read
   *    fails.
   **/
  Entry getEntry(final boolean show_progress)
      throws OutOfRangeException, NoSequenceException, IOException;

  /**
   *  Returns true if and only if this EntrySource always returns "full"
   *  entries.  ie. entries that contain features and sequence.
   **/
  boolean isFullEntrySource();

  /**
   *  Return the name of this source (to display to the user in menus).
   **/
  String getSourceName ();
}
