/* DocumentEntry.java
 *
 * created: Wed Dec 30 1998
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 1998,1999,2000  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/DocumentEntry.java,v 1.2 2005-02-03 15:17:50 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.*;

import java.io.*;
import java.util.Date;

/**
 *  This class extends the Entry class with the data for the entry coming from
 *  a Document object.
 *
 *  @author Kim Rutherford
 *  @version $Id: DocumentEntry.java,v 1.2 2005-02-03 15:17:50 tjc Exp $
 **/

public interface DocumentEntry extends Entry {
  /**
   *  Write this Entry to the Document it came from.  This method uses the
   *  current Document reference (as given by getDocument () or the
   *  constructor).  This method will throw a NullPointerException if the is
   *  no current Document ie. if getDocument () return null.  Use 
   *  save (Document) to save and set the current Document.
   *  @exception IOException thrown if there is a problem saving the entry.
   **/
  void save ()
      throws IOException;

  /**
   *  Write this Entry to the given Document.
   *  @param document This is the file that we will write to.
   *  @exception IOException thrown if there is a problem saving the entry.
   **/
  void save (final Document document)
      throws IOException;

  /**
   *  Write this Entry to the given stream.
   *  @param writer The stream to write to.
   *  @exception IOException thrown if there is a problem writing the entry.
   **/
  void writeToStream (final Writer writer)
      throws IOException;
 

  /**
   *  Arrange for hasUnsavedChanges () to return true until the next save.
   **/
  void setDirtyFlag ();

  /**
   *  Return the Date when this Entry last changed or null if this Entry
   *  hasn't changed since the last save.
   **/
  Date getLastChangeTime ();

  /**
   *  Return the File reference that was passed to the constructor or null if
   *  none was passed.
   **/
  Document getDocument ();
}


