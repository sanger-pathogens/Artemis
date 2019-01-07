/* MiscLineGroup.java
 *
 * created: Sun Sep 26 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/MiscLineGroup.java,v 1.1 2004-06-09 09:50:01 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.LinePushBackReader;

import java.io.IOException;
import java.io.Writer;

/**
 *  Class for objects that contain lines that aren't handle by a more specific
 *  class.
 *
 *  @author Kim Rutherford
 *  @version $Id: MiscLineGroup.java,v 1.1 2004-06-09 09:50:01 tjc Exp $
 **/

abstract public class MiscLineGroup extends LineGroup {
  /**
   *  Create a new MiscLineGroup object.  One line is read from the in_stream
   *  and stored.
   **/
  public MiscLineGroup (LinePushBackReader in_stream)
      throws IOException {
    line = in_stream.readLine ();
  }

  /**
   *  Create a new MiscLineGroup object that consists of the given String
   *  (which should not be terminated with a newline).
   **/
  public MiscLineGroup (String line) {
    this.line = line;
  }

  /**
   *  Return the line that was read from the input stream.
   **/
  public String getLine () {
    return line;
  }

  /**
   *  Write the line stored in this EmblMisc object to the given stream.
   *  @param writer The stream to write to.
   **/
  public void writeToStream (final Writer writer)
      throws IOException {
    
    writer.write (toString ());
  }

  /**
   *  Return this LineGroup object as a String.
   **/
  public String toString () {
    return getLine () + "\n";
  }

  /**
   *  Set the text of this MiscLineGroup.
   **/
  protected void setLine (final String line) {
    this.line = line;
  }

  /**
   *  The line that was read or passed to the constructor.
   **/
  private String line;
}
