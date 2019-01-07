/* Document.java
 *
 * created: Thu Dec 17 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/util/Document.java,v 1.3 2007-03-01 15:42:00 tjc Exp $
 */

package uk.ac.sanger.artemis.util;

import java.io.*;

/**
 *  Each object of this class represents a file/directory on the server (when
 *  the program is running as an applet) or a file/directory in the current
 *  directory (when the program is running stand alone).  Each object
 *  encapsulates a File or URL object.
 *
 *  @author Kim Rutherford
 *  @version $Id: Document.java,v 1.3 2007-03-01 15:42:00 tjc Exp $
 **/

public abstract class Document 
{
  /**
   *  Create a new Document from a (directory) location and a file name.
   *  @param location A directory Document which will be appended to the file
   *    argument
   *  @param name The name of the document to access.
   *  @exception IOException If the directory_document is a URL Document and
   *    the name argument isn't well formed then a MalformedURLException
   *    exception is thrown.
   **/
  public Document(Document directory_document, String name)
      throws IOException 
  {
    location = directory_document.append (name).getLocation ();
  }

  /**
   *  Append a String to the Document location with the correct separator.
   *  @param name The name to append.
   **/
  public abstract Document append(String name) throws IOException;

  /**
   *  Return the name of this Document (the last element of the Document
   *  location).
   **/
  public abstract String getName();

  /**
   *  Return a Document with the last element stripped off.
   **/
  public abstract Document getParent();

  /**
   *  Return true if and only if the Document refered to by this object exists
   *  and is readable.
   **/
  public abstract boolean readable();

  /**
   *  Return true if and only if the Document refered to by this object exists
   *  and can be written to.
   **/
  public abstract boolean writable();

  /**
   *  Create a new InputStream object from this Document.  The contents of the
   *  Document can be read from the stream.
   *  @exception IOException Thrown if the Document can't be read from
   *    (for example if it doesn't exist).
   **/
  public abstract InputStream getInputStream() throws IOException;

  /**
   *  Create a new Reader object from this Document.
   *  @exception IOException An exception is thrown is the Document can't be
   *    read from (for example if it doesn't exist).
   **/
  public Reader getReader() throws IOException 
  {
    return new InputStreamReader(getInputStream ());
  }

  /**
   *  Create a new LinePushBackReader object from this Document.
   *  @exception IOException An exception is thrown is the Document can't be
   *    read from (for example if it doesn't exist).
   *  @return A LinePushBackReader object.  One object is created when
   *    getLinePushBackReader() is first called, subsequent calls will get the
   *    same object.
   **/
  public LinePushBackReader getLinePushBackReader() throws IOException 
  {
    if (line_push_push_reader == null) 
    {
      line_push_push_reader =
        new LinePushBackReader(new InputStreamReader(getInputStream()));
    }

    return line_push_push_reader;
  }

  /**
   *  Add listener that will be passed to the ProgressInputStream constructor
   *  when getInputStream () is called.
   **/
  public void addInputStreamProgressListener(final InputStreamProgressListener listener) 
  {
    listeners.add (listener);
  }

  /**
   *  Return all the InputStreamProgressListener objects for this Document.
   **/
  protected InputStreamProgressListenerVector getProgressListeners () 
  {
    return listeners;
  }

  /**
   *  Create a new OutputStream object from this Document.  The contents of the
   *  Document can be written from the stream.
   *  @exception IOException Thrown if the Document can't be written.
   **/
  public abstract OutputStream getOutputStream() throws IOException;

  /**
   *  Create a new Writer object from this Document.  The Document can then be
   *  written to using the new Writer object.  The old centents of the
   *  Document will be lost.
   *  @exception ReadOnlyException is thrown if the Document is read only (a
   *    URLDocument will always be read only).
   **/
  public Writer getWriter() throws IOException 
  {
    final int BUFFER_SIZE = 100000;
    return new BufferedWriter (new OutputStreamWriter (getOutputStream ()),
                                                       BUFFER_SIZE);
  }

  /**
   *  Return a String representation of this Document.
   **/
  public String toString() 
  {
    return getLocation().toString();
  }

  /**
   *  Create a new Document with the given location.
   **/
  protected Document (Object location) 
  {
    if(location == null) 
      throw new Error ("internal error - created a null Document");
   
    this.location = location;
  }

  /**
   *  Return the raw location object from this Document - eg a File or a URL.
   **/
  public Object getLocation() 
  {
    return location;
  }

  protected void setLocation(Object location)
  {
    this.location = location; 
  }
  
  /**
   *  The actual location object representing this Document - eg a File or URL
   *  object.
   **/
  private Object location = null;

  /**
   *  Returned by getLinePushBackReader().
   **/
  private LinePushBackReader line_push_push_reader = null;

  /**
   *  InputStreamProgressEvents are sent to this objects.
   **/
  private final InputStreamProgressListenerVector listeners =
    new InputStreamProgressListenerVector ();
}
