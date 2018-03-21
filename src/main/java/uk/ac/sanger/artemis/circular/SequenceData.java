/*
 * Copyright (C) 2008  Genome Research Limited
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
 *  @author: Tim Carver
 */

package uk.ac.sanger.artemis.circular;

import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.Transferable;
import java.awt.datatransfer.UnsupportedFlavorException;
import java.io.*;

/**
*
* Object to represent the content of each row in the DragJTable
* of a SequenceList.
*
*/
public class SequenceData implements Transferable, Serializable
{
  /** sequence file/database  */
  public String s_name;       
  /** sequence start */
  public String s_beg;  
  /** sequence end */
  public String s_end;   
  /** use as the default */
  public Boolean s_default;
  /** file on remote file system */
  public Boolean s_remote;    
  /** sequence list file */
  public Boolean s_listFile;

  /** dataflavor for drag and drop transfers */
  public static DataFlavor SEQUENCEDATA =
           new DataFlavor(SequenceData.class, "Sequence data");

  /** available dataflavors */
  static DataFlavor flavors[] = { SEQUENCEDATA };

  public SequenceData()
  {
    s_name = new String();
    s_beg = new String();
    s_end = new String();
    s_default = new Boolean(false);
    s_remote = new Boolean(false);
    s_listFile = new Boolean(false);
  }

  public SequenceData(String name, String beg, String end,
                      Boolean lis, Boolean def, Boolean remote)
  {
    s_name = name;
    s_beg = beg;
    s_end = end;
    s_default = def;
    s_listFile = lis;
    s_remote = remote;
  }

/**
*
* Returns the sequence name
* @return s_name sequence name
*
*/
  public String getSequenceName()
  {
    return s_name;
  }

// Transferable
  public DataFlavor[] getTransferDataFlavors()
  {
    return flavors;
  }

  public boolean isDataFlavorSupported(DataFlavor f)
  {
    if(f.equals(SEQUENCEDATA))
      return true;
    return false;
  }

  public Object getTransferData(DataFlavor d)
      throws UnsupportedFlavorException, IOException
  {
    if(d.equals(SEQUENCEDATA))
      return this;
    else throw new UnsupportedFlavorException(d);
  }

// Serializable
   private void writeObject(java.io.ObjectOutputStream out) throws IOException
   {
     out.defaultWriteObject();
   }

   private void readObject(java.io.ObjectInputStream in)
     throws IOException, ClassNotFoundException
   {
     in.defaultReadObject();
   }

}

