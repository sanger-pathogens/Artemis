/* ReadOnlyEntry.java
 *
 * created: Tue Feb 15 2000
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/ReadOnlyEntry.java,v 1.1 2004-06-09 09:50:26 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.*;

import java.io.*;

/**
 *  Base class for those Entry classes that have only read-only methods.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: ReadOnlyEntry.java,v 1.1 2004-06-09 09:50:26 tjc Exp $
 **/

abstract public class ReadOnlyEntry {
  /**
   *  Always returns true for objects of this class.
   **/
  public boolean isReadOnly () {
    return true;
  }

  /**
   *  This method always returns false for objects of this class.
   **/
  public boolean hasUnsavedChanges () {
    return false;
  }

  /**
   *  Attempt to set the header of this Entry to be the given text - always
   *  fails (returns false).
   **/
  public boolean setHeaderText (final String new_header) {
    return false;
  }

  /**
   *  Always throws a ReadOnlyException exception for objects of this class.
   **/
  public void save () throws IOException {
    throw new ReadOnlyException ("Save is not implemented for this entry");
  }

  /**
   *  Set the name of this Entry - always fails (returns false).
   **/
  public boolean setName (final String name) {
    return false;
  }

  /**
   *  Always throws a ReadOnlyException exception for objects of this class.
   **/
  public Feature createFeature (Key key,
                                Location location,
                                QualifierVector qualifiers)
      throws ReadOnlyException {
    throw new ReadOnlyException ();
  }

  /**
   *  Always throws ReadOnlyException for this type of Entry.
   **/
  public Feature add (final Feature feature)
      throws ReadOnlyException {
    throw new ReadOnlyException ();
  }

  /**
   *  Always throws ReadOnlyException for this type of Entry.
   **/
  public Feature forcedAdd (final Feature feature)
      throws ReadOnlyException {
    throw new ReadOnlyException ();
  }

  /**
   *  Always throws a ReadOnlyException exception for objects of this class.
   **/
  public boolean remove (Feature feature)
      throws ReadOnlyException {
    throw new ReadOnlyException ();
  }
}
