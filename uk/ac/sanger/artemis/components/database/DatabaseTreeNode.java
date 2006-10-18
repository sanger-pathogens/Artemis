/* DatabaseTreeNode.java
 *
 * created: October 2006
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2006  Genome Research Limited
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

package uk.ac.sanger.artemis.components.database;

import java.awt.datatransfer.Transferable;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.UnsupportedFlavorException;
import javax.swing.tree.DefaultMutableTreeNode;

import java.io.Serializable;
import java.io.IOException;

/**
*
* File node for local file tree manager
*
*/
public class DatabaseTreeNode extends DefaultMutableTreeNode 
                 implements Transferable, Serializable
{
  /** data flavour of a database node */
  //public static DataFlavor DATABASETREENODE =
  //       new DataFlavor(DatabaseTreeNode.class, "Database Node");
  /** flavours database node and string */
  //static DataFlavor flavors[] = { DATABASETREENODE, DataFlavor.stringFlavor };
  
  public static final DataFlavor DATABASETREENODE = 
          new DataFlavor(DatabaseTreeNode.class, "Work Package");
  public static final DataFlavor STRING_DATA_FLAVOUR = 
          new DataFlavor(String.class, "text/plain");
  private static final DataFlavor flavors[] = { DATABASETREENODE};
  
  
  private DatabaseEntrySource entrySource;
  private String featureId;
  private String schema;
    
  public DatabaseTreeNode(final String name, 
                          final DatabaseEntrySource entrySource)
  { 
    super(name);
    this.entrySource = entrySource;
  }

  public DatabaseTreeNode(final String name, 
                          final DatabaseEntrySource entrySource,
                          final String featureId,
                          final String schema)
  { 
    super(name);
    this.entrySource = entrySource;
    this.featureId   = featureId;
    this.schema      = schema;
  }
    
// Transferable
  public DataFlavor[] getTransferDataFlavors()
  {
    return flavors;
  }

  public boolean isDataFlavorSupported(DataFlavor f)
  {
    if(f.equals(DATABASETREENODE) || f.equals(DataFlavor.stringFlavor))
      return true;
    return false;
  }

  public Object getTransferData(DataFlavor d)
        throws UnsupportedFlavorException, IOException
  {
    if(d.equals(DATABASETREENODE))
      return this;
    else throw new UnsupportedFlavorException(d);
  }


  public DatabaseEntrySource getEntrySource()
  {
    return entrySource;
  }

  public String getFeatureId()
  {
    return featureId;
  }

  public String getSchema()
  {
    return schema;
  }
  
//Serializable
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
