/* GenbankDocumentEntry.java
 *
 * created: Sun Sep 12 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/GenbankDocumentEntry.java,v 1.1 2004-06-09 09:49:35 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.*;
import java.io.IOException;

/**
 *  A DocumentEntry that can read a GENBANK entry from a Document.
 *
 *  @author Kim Rutherford
 *  @version $Id: GenbankDocumentEntry.java,v 1.1 2004-06-09 09:49:35 tjc Exp $
 **/

public class GenbankDocumentEntry extends PublicDBDocumentEntry {
  /**
   *  Create a new GenbankDocumentEntry object associated with the given
   *  Document.
   *  @param entry_information The EntryInformation object for this Entry.
   *  @param document This is the file that we will read from.  This is also
   *    used for saving the entry back to the file it came from and to give
   *    the new object a name.
   *  @param listener The object that will listen for ReadEvents.
   *  @exception IOException thrown if there is a problem reading the entry -
   *    most likely ReadFormatException.
   *  @exception EntryInformationException Thrown if force is false and if this
   *    Entry cannot contain the Key, Qualifier or Key/Qualifier combination of
   *    one of the features in the given Entry.
   **/
  public GenbankDocumentEntry (final EntryInformation entry_information,
                               final Document document,
                               final ReadListener listener)
      throws IOException, EntryInformationException {
    super (entry_information, document, listener);
  }

  /**
   *  Create a new GenbankDocumentEntry that will be a copy of the given
   *  Entry and has no Document associated with it.  The new
   *  GenbankDocumentEntry cannot be saved to a file with save () unless save
   *  (Document) has been called first.  The save (Document) method will
   *  assign a Document..
   *  @param entry_information The EntryInformation object for this Entry.
   *  @param force If true then invalid qualifiers and any features with
   *    invalid keys in the new Entry will be quietly thrown away.  "Invalid"
   *    means that the key/qualifier is not allowed to occur in an Entry of
   *    this type (probably determined by the EntryInformation object of this
   *    Entry).  If false an EntryInformationException will be thrown for
   *    invalid keys or qualifiers.
   *  @exception EntryInformationException Thrown if force is false and if this
   *    Entry cannot contain the Key, Qualifier or Key/Qualifier combination of
   *    one of the features in the given Entry.
   **/
  public GenbankDocumentEntry (final EntryInformation entry_information,
                               final Entry new_entry, final boolean force)
      throws EntryInformationException {
    super (entry_information, new_entry, force);
  }

  /**
   *  Create a new GenbankDocumentEntry that will be a copy of the given
   *  Entry and has no Document associated with it.  The new
   *  GenbankDocumentEntry cannot be saved to a file with save () unless save
   *  (Document) has been called first.  The save (Document) method will
   *  assign a Document.  The new DocumentEntry will use a copy of the
   *  EntryInformation object from the given Entry. 
   *  @exception EntryInformationException Thrown if this Entry cannot contain
   *    the Key, Qualifier or Key/Qualifier combination of one of the features
   *    in the given Entry.
   **/
  public GenbankDocumentEntry (final Entry new_entry)
      throws EntryInformationException {
    super (new SimpleEntryInformation (new_entry.getEntryInformation ()),
           new_entry, false);
  }

  /**
   *  Create a new empty GenbankDocumentEntry object that has no Document
   *  associated with it.  The new GenbankDocumentEntry cannot be saved to a
   *  file with save () unless save (Document) has been called first.  The
   *  save (Document) method will assign a Document.
   *  @param entry_information The EntryInformation object for this Entry.
   **/
  public GenbankDocumentEntry (final EntryInformation entry_information) {
    super (entry_information);
  }
}
