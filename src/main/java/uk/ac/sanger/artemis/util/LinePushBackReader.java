/* LinePushBackReader.java
 *
 * created: Mon Oct 12 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/util/LinePushBackReader.java,v 1.1 2004-06-09 09:53:05 tjc Exp $
 */

package uk.ac.sanger.artemis.util;

import java.io.Reader;
import java.io.LineNumberReader;
import java.io.IOException;

/**
 *  This class buffers exactly one line of input, which can be pushed back if
 *  needed.  All methods in Reader can be called on objects of this class, but
 *  only readLine () understands about pushBack ().  This class also counts
 *  lines using LineNumberReader (see the getLineNumber () method).
 *
 *  @author Kim Rutherford
 *  @version $Id: LinePushBackReader.java,v 1.1 2004-06-09 09:53:05 tjc Exp $
 */

public class LinePushBackReader extends Reader {
  /**
   *  Create a new LinePushBackReader object.
   *  @param in_stream the stream to create the object from
   */
  public LinePushBackReader (Reader in_stream) {
    buffered_reader = new LineNumberReader (in_stream);
  }


  /**
   *  Does the same thing as Reader.readLine (), but understands about
   *  pushBack ()
   *  @return the next line from the input stream or null if we are at the end
   *    of file
   */
  public String readLine () throws IOException {
    if (line_buffer != null) {
      String tmp_string = line_buffer;
      line_buffer = null;
      return tmp_string;
    } else {
      String line = buffered_reader.readLine ();
      return line;
    }
  }


  /**
   *  Does the same thing as LineNumberReader.getLineNumber ()
   */
  public int getLineNumber () {
    return buffered_reader.getLineNumber ();
  }

  
  /**
   *  Does the same thing as Reader.read ().
   */
  public int read(char [] cbuf, int off, int len)
       throws IOException {
    return buffered_reader.read (cbuf, off, len);
  }


  /**
   *  Does the same thing as Reader.close ().
   */
  public void close ()
       throws IOException {
    buffered_reader.close ();
  }

  /**
   *  Push a String back into the LinePushBackReader stream.  This can be
   *  called once between calls to readLine ().  Note that you don't have to
   *  push back the any String can pushed back.
   *  @param line The string to push back
   *  @exception PushBackException Thrown if this method is called twice before
   *   a call to readLine ()
   */
  public void pushBack (String line)
       throws PushBackException {

    if (line_buffer != null) {
      throw new PushBackException ("LinePushBackReader.pushBack () " +
                                   "called twice before calling readLine ()");
    } else {
      line_buffer = line;
    }
  }

  /**
   *  A single line buffer that is used only when the user invokes
   *  pushBack () or calls readLine () after pushBack ().
   */
  private String line_buffer;

  /**
   *  This is the stream that was passed to the constructor.
   */
  private LineNumberReader buffered_reader;
}


