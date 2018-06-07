/* EmblMisc.java
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/EmblMisc.java,v 1.1 2004-06-09 09:49:10 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.LinePushBackReader;

import java.io.IOException;

/**
 *  This class is used to store EMBL entry lines that are not handled by other
 *  classes.
 *
 *  @author Kim Rutherford
 *  @version $Id: EmblMisc.java,v 1.1 2004-06-09 09:49:10 tjc Exp $
 **/

public class EmblMisc extends MiscLineGroup {
  /**
   *  Create a new EmblMisc object.  One line is read from the in_stream
   *  and stored.
   **/
  public EmblMisc (final LinePushBackReader in_stream)
      throws IOException {
    super (in_stream);
  }

  /**
   *  Create a new EmblMisc object that consists of the given String
   *  (which should not be terminated with a newline).
   **/
  public EmblMisc (final String line) {
    super (line);
  }
}
