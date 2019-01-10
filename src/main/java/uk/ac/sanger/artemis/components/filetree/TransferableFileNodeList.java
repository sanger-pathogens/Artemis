/********************************************************************
*
*  This library is free software; you can redistribute it and/or
*  modify it under the terms of the GNU Library General Public
*  License as published by the Free Software Foundation; either
*  version 2 of the License, or (at your option) any later version.
*
*  This library is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
*  Library General Public License for more details.
*
*  You should have received a copy of the GNU Library General Public
*  License along with this library; if not, write to the
*  Free Software Foundation, Inc., 59 Temple Place - Suite 330,
*  Boston, MA  02111-1307, USA.
*
*  @author: Copyright (C) Tim Carver
*
********************************************************************/

package uk.ac.sanger.artemis.components.filetree;

import java.awt.datatransfer.*;
import java.io.*;
import java.util.*;

/**
*
* File node for remote file tree manager
*
*/
public class TransferableFileNodeList extends Vector
                    implements Transferable, Serializable
{
  final public static DataFlavor TRANSFERABLEFILENODELIST = 
         new DataFlavor(TransferableFileNodeList.class, "file list");
  static DataFlavor remoteFlavors[] = { TRANSFERABLEFILENODELIST }; 

  public TransferableFileNodeList(int n)
  {
    super(n);
  }

// Transferable
  public DataFlavor[] getTransferDataFlavors()
  {
    return remoteFlavors;
  }

  public boolean isDataFlavorSupported(DataFlavor f)
  {
    return f.equals(TRANSFERABLEFILENODELIST);
  }

  public Object getTransferData(DataFlavor d) 
      throws UnsupportedFlavorException, IOException
  {
    if(d.equals(TRANSFERABLEFILENODELIST))
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

