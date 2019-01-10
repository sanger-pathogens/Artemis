/* Logger.java
 *
 * created: Wed Aug 30 2000
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/Logger.java,v 1.1 2004-06-09 09:44:55 tjc Exp $
 */

package uk.ac.sanger.artemis;

/**
 *  An interface for simple logging of warnings/errors/messages.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: Logger.java,v 1.1 2004-06-09 09:44:55 tjc Exp $
 **/

import java.io.*;

public interface Logger 
{
  /** Send the given String to the log. */
  void log (final String message);

  /** Read from the given Reader and send it to the log. */
  void log (final Reader reader)
      throws IOException;
}
