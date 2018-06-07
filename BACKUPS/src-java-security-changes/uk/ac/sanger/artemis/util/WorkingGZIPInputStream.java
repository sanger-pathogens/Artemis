/* WorkingGZIPInputStream.java
 *
 * created: Thu May 30 2002
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2001  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/util/WorkingGZIPInputStream.java,v 1.1 2004-06-09 09:53:17 tjc Exp $
 */

package uk.ac.sanger.artemis.util;

import java.io.*;

/**
 *  A wrapper class to work around a bug in GZIPInputStream.  The read()
 *  method sometimes throws a IOException complaining about a Corrupt GZIP
 *  trailer on the jdk 1.1.8-13 on the Alphas.  This wrapper catches and
 *  ignores that exception.
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: WorkingGZIPInputStream.java,v 1.1 2004-06-09 09:53:17 tjc Exp $
 **/

public class WorkingGZIPInputStream extends java.util.zip.GZIPInputStream {
  /**
   *  Creates a new input stream with the specified buffer size.
   *  @param in the input stream
   *  @param size the input buffer size
   *  @exception IOException if an I/O error has occurred
   **/

  public WorkingGZIPInputStream (InputStream in, int size)
      throws IOException {
    super (in, size);
  }

  /**
   *  Creates a new input stream with a default buffer size.
   *  @param in the input stream
   *  @exception IOException if an I/O error has occurred
   **/
  public WorkingGZIPInputStream (InputStream in)
      throws IOException {
    super (in);
  }

  /**
   *  Calls super.read() and then catch and ignore any IOExceptions that
   *  mention "Corrupt GZIP trailer".
   **/
  public int read (byte buf[], int off, int len)
      throws IOException {
    try {
      return super.read (buf, off, len);
    } catch (IOException e) {
      if (e.getMessage ().indexOf ("Corrupt GZIP trailer") != -1) {
        return -1;
      } else {
        throw e;
      }
    }
  }
}
