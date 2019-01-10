/* FileDocument.java
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/util/FileDocument.java,v 1.1 2004-06-09 09:53:01 tjc Exp $
 */

package uk.ac.sanger.artemis.util;

import java.io.*;

/**
 *  Objects of this class are Documents created from a file.
 *
 *  @author Kim Rutherford
 *  @version $Id: FileDocument.java,v 1.1 2004-06-09 09:53:01 tjc Exp $
 **/

public class FileDocument extends Document {
  /**
   *  Create a new FileDocument from a File.
   *  @param location This should be a file or directory name.
   **/
  public FileDocument (File location) {
    super (location);
  }

  /**
   *  Append a String to the Document location with the correct separator.
   *  @param name The name to append.
   **/
  public Document append (String name) throws IOException {
    return new FileDocument (new File (getFile (), name));
  }

  /**
   *  Return the name of this Document (the last element of the Document
   *  location).
   **/
  public String getName () {
    return getFile ().getName ();
  }

  /**
   *  Return a Document with the last element stripped off.
   **/
  public Document getParent () {
    try {
      final File canonical_file = new File (getFile ().getCanonicalPath ());
      return new FileDocument (new File (canonical_file.getParent ()));
    } catch (IOException e) {
      return null;
    }
  }
  
  /**
   *  Return true if and only if the Document refered to by this object exists
   *  and is readable.
   **/
  public boolean readable () {
    if (getFile ().exists () && getFile ().canRead ()) {
      return true;
    } else {
      return false;
    }
  }


  /**
   *  Return true if and only if the Document refered to by this object exists
   *  and can be written to.
   **/
  public boolean writable () {
    if (getFile ().exists () && getFile ().canWrite ()) {
      return true;
    } else {
      return false;
    }
  }

  /**
   *  Create a new InputStream object from this Document.  The contents of the
   *  Document can be read from the InputStream.
   *  @exception IOException Thrown if the Document can't be read from
   *    (for example if it doesn't exist).
   **/
  public InputStream getInputStream ()
      throws IOException {
    final File read_file = (File) getLocation ();

    final InputStream file_input_stream =
      new ProgressInputStream (new FileInputStream (read_file),
                               getProgressListeners ());
    
    if (read_file.getName ().endsWith (".gz")) {
      final BufferedInputStream ins = new BufferedInputStream(file_input_stream);
      
      if (! System.getProperty("java.version").startsWith("1.5.")) 
      {
        if(DocumentBlockCompressed.isValidFile(ins)) // BGZIP
          return DocumentBlockCompressed.getBlockCompressedInputStream(ins);
      }
      
      ins.close();
      file_input_stream.close();

      // assume GZIP
      return new WorkingGZIPInputStream (
          new ProgressInputStream (new FileInputStream (read_file),
          getProgressListeners ()));
    } else
      return file_input_stream;      
  }

  /**
   *  Create a new OutputStream object from this Document.  The Document can
   *  then be written to using the new object.  The old centents of the
   *  Document will be lost.
   *  @exception ReadOnlyException is thrown if the Document is read only.
   **/
  public OutputStream getOutputStream () throws IOException {
    final File write_file = (File) getLocation ();

    final FileOutputStream file_output_stream =
      new FileOutputStream (write_file);

    if (write_file.getName ().endsWith (".gz")) {
      // assume this file should be gzipped
      return new java.util.zip.GZIPOutputStream (file_output_stream);
    } else {
      return file_output_stream;
    }
  }

  /**
   *  Return the File object that this FileDocument is encapsulating.
   **/
  public File getFile () {
    return (File) getLocation ();
  }
  
}

class DocumentBlockCompressed
{
  protected static boolean isValidFile(BufferedInputStream ins) throws IOException
  {
    return htsjdk.samtools.util.BlockCompressedInputStream.isValidFile(ins);
  }
  
  protected static InputStream getBlockCompressedInputStream(BufferedInputStream ins) throws IOException
  {
    return new htsjdk.samtools.util.BlockCompressedInputStream(ins);
  }
}
