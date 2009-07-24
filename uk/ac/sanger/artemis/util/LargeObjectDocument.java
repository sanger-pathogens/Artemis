/* LargeObjectDocument.java
 *
 * created: 2009
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 2009  Genome Research Limited
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
 */

package uk.ac.sanger.artemis.util;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.sql.SQLException;

import org.postgresql.largeobject.LargeObject;

/**
 *  Objects of this class are Documents created from a database
 *  and for reading blobs as a LargeObject.
 **/
public class LargeObjectDocument extends Document 
{
  private String name;
  private LargeObject obj;
  
  /**
   *  Create a new FileDocument from a File.
   *  @param location This should be a file or directory name.
   **/
  public LargeObjectDocument (String location, String name, LargeObject obj)
  {
    super (location);
    this.obj = obj;
    this.name = name;
  }

  /**
   *  Return the name of this Document (the last element of the Document
   *  location).
   **/
  public String getName ()
  {
    return name;
  }

  /**
   *  Create a new InputStream object from this Document.  The contents of the
   *  Document can be read from the InputStream.
   *  @exception IOException Thrown if the Document can't be read from
   *    (for example if it doesn't exist).
   **/
  public InputStream getInputStream ()
      throws IOException 
  {
    try
    {
      try
      {
        // assume zipped 
        return new WorkingGZIPInputStream (obj.getInputStream());
      }
      catch(IOException e)
      {
        obj.seek(0);
        return obj.getInputStream();
      }
    }
    catch (SQLException e)
    {
      e.printStackTrace();
      throw new IOException( e.getMessage() );
    }
  }

  /**
   *  Create a new OutputStream object from this Document.  The Document can
   *  then be written to using the new object.  The old centents of the
   *  Document will be lost.
   *  @exception ReadOnlyException is thrown if the Document is read only.
   **/
  public OutputStream getOutputStream () throws IOException 
  {
    return null;
  }

  public Document append(String name) throws IOException
  {
    return null;
  }

  public Document getParent()
  {
    return null;
  }

  public boolean readable()
  {
    return false;
  }

  public boolean writable()
  {
    return false;
  }
}
