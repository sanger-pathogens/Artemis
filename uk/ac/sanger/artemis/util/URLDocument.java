/* URLDocument.java
 *
 * created: Fri Dec 18 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/util/URLDocument.java,v 1.1 2004-06-09 09:53:16 tjc Exp $
 */

package uk.ac.sanger.artemis.util;

import java.net.*;
import java.io.*;

/**
 *  Objects of this class are Documents created from a URL.
 *
 *  @author Kim Rutherford
 *  @version $Id: URLDocument.java,v 1.1 2004-06-09 09:53:16 tjc Exp $
 **/

public class URLDocument extends Document {
  /**
   *  Create a new Document from a URL.
   *  @param location This should be a URL string giving the location on the
   *    server where a file or directory is to be found (for example
   *    getDocumentBase ()).
   **/
  public URLDocument (URL location) {
    super (location);
  }

  /**
   *  Append a String to the Document location with the correct separator.
   *  @param name The name to append.
   **/
  public Document append (String name) throws IOException {
    return new URLDocument (new URL (getURL ().toString () + '/' + name));
  }

  /**
   *  Return the name of this Document (the last element of the Document
   *  location).
   **/
  public String getName () {
    final String url_string = getURL ().toString ();
    final String file_url_string =
      url_string.substring (url_string.lastIndexOf ('/') + 1);

    return file_url_string;
  }

  /**
   *  Return a Document with the last element stripped off.
   **/
  public Document getParent () {
    final String url_string = getURL ().toString ();
    final String directory_url_string =
      url_string.substring (0, url_string.lastIndexOf ('/'));

    try {
      return new URLDocument (new URL (directory_url_string));
    } catch (IOException e) {
      return null;
    }
  }
  
  /**
   *  Return true if and only if the Document refered to by this object exists
   *  and is readable.  Always returns true.
   **/
  public boolean readable () {
    // XXX
    return true;
  }


  /**
   *  Return true if and only if the Document refered to by this object exists
   *  and can be written to.  Always returns false.
   **/
  public boolean writable () {
    return false;
  }


  /**
   *  Create a new InputStream object from this Document.  The contents of the
   *  Document can be read from the InputStream.
   *  @exception IOException Thrown if the Document can't be read from
   *    (for example if it doesn't exist).
   **/
  public InputStream getInputStream () throws IOException {
    // location must be a URL object
    final URL url = (URL) getLocation ();
    final URLConnection connection = url.openConnection ();

    connection.connect ();

    final InputStream in_stream =
      new ProgressInputStream (connection.getInputStream (),
                               getProgressListeners ());

    if (getURL ().toString ().endsWith (".gz")) {
      // assume this file is gzipped
      System.out.println (getName ());
      return new java.util.zip.GZIPInputStream (in_stream);
    } else {
      return in_stream;
    }
  }

  /**
   *  Create a new OutputStream object from this Document.  The contents of the
   *  Document can be written from the stream.
   *  @exception IOException Thrown if the Document can't be written.
   **/
  public OutputStream getOutputStream () throws IOException {
    throw new ReadOnlyException ("this Document can not be written to");
  }

  /**
   *  Return the URL object that this URLDocument is encapsulating.
   **/
  private URL getURL () {
    return (URL) getLocation ();
  }
}


