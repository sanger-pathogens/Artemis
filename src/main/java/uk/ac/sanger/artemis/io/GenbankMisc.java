/* GenbankMisc.java
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/GenbankMisc.java,v 1.1 2004-06-09 09:49:36 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.LinePushBackReader;

import java.io.IOException;
import java.io.Writer;

/**
 *  This class is used to store EMBL entry lines that are not handled by other
 *  classes.
 *
 *  @author Kim Rutherford
 *  @version $Id: GenbankMisc.java,v 1.1 2004-06-09 09:49:36 tjc Exp $
 **/

public class GenbankMisc extends MiscLineGroup {
  /**
   *  Create a new GenbankMisc object.  One or more lines are read from the
   *  in_stream and stored.
   **/
  public GenbankMisc (LinePushBackReader in_stream)
      throws IOException {
    super (in_stream.readLine ());

    if (getLine ().startsWith ("FEATURES ")) {
      // special case (hack) for "FEATURES" line - just read one line.
      return;
    }
    while (true) {
      final String temp_line = in_stream.readLine ();

      if (temp_line != null && temp_line.startsWith (" ")) {
        setLine (getLine () + "\n" + temp_line);
      } else {
        // end of line group
        in_stream.pushBack (temp_line);
        return;
      }
    }
  }

  /**
   *  Write the line stored in this GenbankMisc object to the given stream.
   *  @param writer The stream to write to.
   **/
  public void writeToStream (final Writer writer)
      throws IOException {
    
    writer.write (toString ());
  }
}
