/* FeatureHeader.java
 *
 * created: Thu Jul 20 2000
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/FeatureHeader.java,v 1.1 2004-06-09 09:49:24 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.LinePushBackReader;

import java.io.IOException;

/**
 *  Class used to store FH lines.
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: FeatureHeader.java,v 1.1 2004-06-09 09:49:24 tjc Exp $
 **/

public class FeatureHeader extends EmblMisc {
  /**
   *  Create a new FeatureHeader object by reading the current "FH" line and
   *  the ones directly after it.
   **/
  public FeatureHeader (final LinePushBackReader in_stream)
      throws IOException {
    super (getLines (in_stream));
  }

  /**
   *  Read the next "FH" lines from the stream.
   **/
  private static String getLines (final LinePushBackReader in_stream)
      throws IOException {
    final StringBuilder buffer = new StringBuilder ();
    while (true) {
      final String this_line = in_stream.readLine ();

      if (this_line == null) {
        break;
      }

      if (getLineType (this_line) == EMBL_FEATURE_HEADER) {
        if (buffer.length () != 0) {
          buffer.append ("\n");
        }
        buffer.append (this_line);
      } else {
        in_stream.pushBack (this_line);
        break;
      }
    }

    return buffer.toString ();
  }
}
