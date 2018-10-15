/* TransferableContig.java
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2005  Genome Research Limited
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or(at your option) any later version.
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

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.Feature;

import java.awt.datatransfer.*;
import java.io.*;
import java.util.*;

/**
*
* File node for remote file tree manager
*
*/
public class TransferableContig 
                    implements Transferable, Serializable
{
  final public static DataFlavor TRANSFERABLECONTIG = 
         new DataFlavor(TransferableContig.class, "contig feature");
  static DataFlavor contig_flavors[] = { TRANSFERABLECONTIG }; 

  private Feature feature;

  public TransferableContig(final Feature feature)
  {
    this.feature = feature;
  }

// Transferable
  public DataFlavor[] getTransferDataFlavors()
  {
    return contig_flavors;
  }

  public boolean isDataFlavorSupported(DataFlavor f)
  {
    return f.equals(TRANSFERABLECONTIG);
  }

  public Object getTransferData(DataFlavor d) 
      throws UnsupportedFlavorException, IOException
  {
    if(d.equals(TRANSFERABLECONTIG))
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

