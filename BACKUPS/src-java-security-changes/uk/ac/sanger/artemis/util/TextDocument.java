/* TextDocument.java
 *
 * created: 2006
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 2006  Genome Research Limited
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
 */

package uk.ac.sanger.artemis.util;

import javax.swing.JOptionPane;
import java.awt.Toolkit;
import java.awt.datatransfer.Clipboard;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.Transferable;
import java.awt.datatransfer.UnsupportedFlavorException;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

/**
 *  Objects of this class are Documents created from a file.
 **/
public class TextDocument extends Document 
{
  /**
  *  Create a new TextDocument from a String.
  *  @param location This should be a file or directory name. 
  *  @param remote_file File on server
  **/
  public TextDocument() 
  {
    super("Copy and Paste");
  }

  /**
  *
  *  Append a String to the Document location with the correct separator.
  *  @param name The name to append.
  *
  **/
  public Document append(String name) throws IOException 
  {
    return null;
  }

  /**
  *
  *  Return the name of this Document (the last element of the Document
  *  location).
  *
  **/
  public String getName() 
  {
    return null;
  }

  /**
   *  Return a Document with the last element stripped off.
   **/
  public Document getParent()
  {
    return null;
  }
  
  /**
   *  Return true if and only if the Document refered to by this object exists
   *  and is readable.
   **/
  public boolean readable() 
  {
    return true;
  }


  /**
   *  Return true if and only if the Document refered to by this object exists
   *  and can be written to.
   **/
  public boolean writable() 
  {
    return true;
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
    Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
    //odd: the Object param of getContents is not currently used
    Transferable contents = clipboard.getContents(null);

    boolean hasTransferableText = (contents != null) &&
                                  contents.isDataFlavorSupported(DataFlavor.stringFlavor);
    byte[] bytes = null;

    if(hasTransferableText) 
    {
      try
      {
        bytes = ((String)contents.getTransferData(DataFlavor.stringFlavor)).getBytes();
      }
      catch(UnsupportedFlavorException ex)
      {
        //highly unlikely since we are using a standard DataFlavor
        System.out.println(ex);
      }
      catch(IOException ex) 
      {
        System.out.println(ex);
      }
    }

    final InputStream file_input_stream =
      new ProgressInputStream(new ByteArrayInputStream(bytes),
                              getProgressListeners());
    
    return file_input_stream;      
  }

  /**
   *  Create a new OutputStream object from this Document.  The Document can
   *  then be written to using the new object.  The old centents of the
   *  Document will be lost.
   *  @exception ReadOnlyException is thrown if the Document is read only.
   **/
  public OutputStream getOutputStream() throws IOException 
  {
    final File write_file = new File(getName());

    if(write_file.exists())
    {
      int n = JOptionPane.showConfirmDialog(null,
                 "The file :\n"+getName()+
                 "\nalready exists on the local disk.\nOverwrite?",
                 "Overwrite "+getName(),
                 JOptionPane.YES_NO_OPTION);
      if(n == JOptionPane.NO_OPTION)
        return null;
    }

    final FileOutputStream file_output_stream =
      new FileOutputStream(write_file);

    if(write_file.getName().endsWith(".gz")) 
      return new java.util.zip.GZIPOutputStream(file_output_stream);
    else
      return file_output_stream;
  }

}
