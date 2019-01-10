/* ProgressInputStream.java
 *
 * created: Mon Sep 20 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/util/ProgressInputStream.java,v 1.1 2004-06-09 09:53:08 tjc Exp $
 */

package uk.ac.sanger.artemis.util;

import java.io.*;

/**
 *  This is a FileInputStream which has the readLine () and close () methods
 *  overridden to show the progress in a Label component.
 *
 *  @author Kim Rutherford
 *  @version $Id: ProgressInputStream.java,v 1.1 2004-06-09 09:53:08 tjc Exp $
 **/

public class ProgressInputStream extends FilterInputStream {
  /**
   *  Creates an InputStream to read from the specified InputStream object and
   *  show the progress in a Label as it goes.
   *  @param input_stream the InputStream to be read from
   *  @param listener InputStreamProgressEvent objects will be sent to this
   *    listener as progress on reading is made.
   *  @exception IOException  if an I/O error occurs.
   **/
  public
    ProgressInputStream (final InputStream input_stream,
                         final InputStreamProgressListenerVector listeners) {
    super (input_stream);

    this.listeners = listeners;
  }

  /**
   * Reads up to <code>b.length</code> bytes of data from this input
   * stream into an array of bytes. This method blocks until some input
   * is available.
   *
   * @param      b   the buffer into which the data is read.
   * @return     the total number of bytes read into the buffer, or
   *             <code>-1</code> if there is no more data because the end of
   *             the file has been reached.
   * @exception  IOException  if an I/O error occurs.
   */
  public int read(byte b[]) throws IOException {
    final int result = super.read (b);

    if (result > 0) {
      byte_count += b.length;
      maybeFireEvent ();
    } else {
      if (result == -1) {
        final InputStreamProgressEvent event =
          new InputStreamProgressEvent (InputStreamProgressEvent.EOF);
        fireEvent (event);
      }
    }

    return result;
  }

  /**
   * Reads up to <code>len</code> bytes of data from this input stream
   * into an array of bytes. This method blocks until some input is
   * available.
   *
   * @param      b     the buffer into which the data is read.
   * @param      off   the start offset of the data.
   * @param      len   the maximum number of bytes read.
   * @return     the total number of bytes read into the buffer, or
   *             <code>-1</code> if there is no more data because the end of
   *             the file has been reached.
   * @exception  IOException  if an I/O error occurs.
   */
  public int read(byte b[], int off, int len) throws IOException {
    final int result = super.read (b, off, len);

    if (result > 0) {
      byte_count += b.length;
      maybeFireEvent ();
    } else {
      if (result == -1) {
        final InputStreamProgressEvent event =
          new InputStreamProgressEvent (InputStreamProgressEvent.EOF);
        fireEvent (event);
      } 
    }

    return result;
  }

  /**
   *  Close this file input stream and release any system resources associated
   *  with the stream.
   *
   *  @exception  IOException  if an I/O error occurs.
   **/
  public void close ()
      throws IOException {
    super.close ();

    final InputStreamProgressEvent event =
      new InputStreamProgressEvent (InputStreamProgressEvent.EOF);
    fireEvent (event);
  }

  /**
   *  If we haven't sent a InputStreamProgressEvent to the listener recently,
   *  then send one.
   **/
  private void maybeFireEvent () {
    if (byte_count > 50000) {
      // do some rounding (we divide by 2 because chars are twice the size
      // of bytes) 
      final int char_count = (byte_count / 10000) * 10000;

      if (last_count + 10000 < char_count) {
        fireEvent (new InputStreamProgressEvent (char_count));
        last_count = char_count;
      }
    }
  }

  /**
   *  Send the event to all the listeners.
   **/
  private void fireEvent (final InputStreamProgressEvent event) {
    if (listeners != null) {
      for (int i = 0 ; i < listeners.size () ; ++i) {
        listeners.elementAt (i).progressMade (event);
      }
    }
  }

  /**
   *  InputStreamProgressEvents are sent to these object.
   **/
  private final InputStreamProgressListenerVector listeners;

  /**
   *  The number of bytes that have been read so far.
   **/
  private int byte_count = 0;

  /**
   *  This is the char count that was sent to the listener last time.
   **/
  private int last_count = 0;
}
