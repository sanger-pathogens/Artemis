/* PublicDBDocumentEntry.java
 *
 * created: Sat Sep 11 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/PublicDBDocumentEntry.java,v 1.1 2004-06-09 09:50:03 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.*;

import java.util.Vector;
import java.io.Reader;
import java.io.StringReader;
import java.io.BufferedReader;
import java.io.Writer;
import java.io.FileWriter;
import java.io.File;
import java.io.IOException;

/**
 *  This class extends the Entry class with the data for the entry coming from
 *  a Document object.  The Document must contain an EMBL entry or a GENBANK
 *  entry.
 *
 *  @author Kim Rutherford
 *  @version $Id: PublicDBDocumentEntry.java,v 1.1 2004-06-09 09:50:03 tjc Exp $
 **/

public class PublicDBDocumentEntry extends SimpleDocumentEntry
    implements DocumentEntry 
{
  /**
   *  Create a new PublicDBDocumentEntry object associated with the given
   *  Document.
   *  @param entry_information The EntryInformation object of the new Entry.
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
  PublicDBDocumentEntry(final EntryInformation entry_information,
                        final Document document, final ReadListener listener)
      throws IOException, EntryInformationException 
  {
    super(entry_information, document, listener);
  }

  /**
   *  Create a new PublicDBDocumentEntry that will be a copy of the given
   *  Entry and has no Document associated with it.  The new
   *  PublicDBDocumentEntry cannot be saved to a file with save () unless save
   *  (Document) has been called first.
   *  @param entry_information The EntryInformation object of the new Entry.
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
  public PublicDBDocumentEntry(final EntryInformation entry_information,
                               final Entry new_entry, final boolean force)
      throws EntryInformationException 
  {
    super(entry_information, new_entry, force);
  }

  /**
   *  Create a new empty PublicDBDocumentEntry object that has no Document
   *  associated with it.  The new PublicDBDocumentEntry cannot be saved to a
   *  file with save () unless save (Document) has been called first.  The
   *  save (Document) method will assign a Document.
   *  @param entry_information The EntryInformation object of the Entry that
   *    will contain this Feature.
   **/
  public PublicDBDocumentEntry(final EntryInformation entry_information) 
  {
    super(entry_information);
  }

  /**
   *  If the given feature can be added directly to this Entry, then return
   *  it, otherwise create and return a new feature of the appropriate type.
   *  @param copy if true then always new a new copy of the Feature.
   **/
  protected SimpleDocumentFeature makeNativeFeature(final Feature feature,
                                                    final boolean copy) 
  {
    if (!copy && (feature instanceof EmblStreamFeature &&
                  this instanceof EmblDocumentEntry ||
                  feature instanceof GenbankStreamFeature &&
                  this instanceof GenbankDocumentEntry)) 
    {
      return (PublicDBStreamFeature) feature;
    } 
    else 
    {
      final PublicDBStreamFeature feature_copy;

      if (this instanceof EmblDocumentEntry) 
        feature_copy = new EmblStreamFeature (feature);
      else 
        feature_copy = new GenbankStreamFeature (feature);
      
      return feature_copy;
    }
  }

  /**
   *  If the given Sequence can be added directly to this Entry, then return a
   *  copy of it, otherwise create and return a new feature of the appropriate
   *  type for this Entry.
   **/
  protected StreamSequence makeNativeSequence (final Sequence sequence) 
  {
    if(this instanceof EmblDocumentEntry) 
      return new EmblStreamSequence (sequence);
    else 
      return new GenbankStreamSequence (sequence);
  }
}
