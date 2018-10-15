/* GFFMisc.java
 *
 * created: Thu Feb 17 2000
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/GFFMisc.java,v 1.1 2004-06-09 09:49:32 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.LinePushBackReader;

import java.io.IOException;

/**
 *  GFFMisc class
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: GFFMisc.java,v 1.1 2004-06-09 09:49:32 tjc Exp $
 **/

public class GFFMisc extends MiscLineGroup {
  /**
   *  Create a new MiscLineGroup object.  One line is read from the in_stream
   *  and stored.
   **/
  public GFFMisc (LinePushBackReader in_stream)
      throws IOException {
    super (in_stream);
  }
}
