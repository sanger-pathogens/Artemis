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
*  Copyright (C) Genome Research Limited
*
********************************************************************/

package uk.ac.sanger.artemis.components;

import java.awt.datatransfer.*;
import javax.swing.tree.*;
import java.io.*;
import java.util.*;

/**
*
* File node for local file tree manager
*
*/
public class FileNode extends DefaultMutableTreeNode 
                 implements Transferable, Serializable
{
    /** true if explored */
    private boolean explored = false;
    /** data flavour of a file node */
    public static DataFlavor FILENODE =
           new DataFlavor(FileNode.class, "Local file");
    /** flavours file node and string */
    static DataFlavor flavors[] = { FILENODE, DataFlavor.stringFlavor };

    /**
    *
    * @param file	file node file
    *
    */
    public FileNode(File file)
    { 
      setUserObject(file); 
    }

    /** Determine if this is a directory */
    public boolean getAllowsChildren() { return isDirectory(); }
    /** Determine if this is a file */
    public boolean isLeaf() { return !isDirectory(); }
    /** Get the File this node represents */
    public File getFile() { return (File)getUserObject(); }
    /** Determine if this node has been explored */
    public boolean isExplored() { return explored; }
    /** Determine if this is a directory */
    public boolean isDirectory() 
    {
      File file = getFile();
      return file.isDirectory();
    }

    /**
    *
    * Returns the name of the file 
    *
    */
    public String toString() 
    {
      File file = (File)getUserObject();
      String filename = file.toString();
      int index = filename.lastIndexOf(File.separator);

      return (index != -1 && index != filename.length()-1) ? 
                          filename.substring(index+1) : 
                                            filename;
    }

    /**
    *
    * Explores a directory adding a FileNode for each
    * child
    *
    */
    public void explore(FileFilter filter) 
    {
      if(!isDirectory())
        return;

      if(!isExplored()) 
      {
        File file = getFile();
        explored = true;
        File[] children;
// filter files
        children = file.listFiles(filter);
        
// sort into alphabetic order
        java.util.Arrays.sort(children);
        for(int i=0; i < children.length; ++i)
          add(new FileNode(children[i]));
      }
    }

    /**
    *
    * Forces the directory to be re-explored
    *
    */
    public void reExplore(FileFilter filter)
    {
      explored = false;
      removeAllChildren();
      explore(filter);
    }
 
// Transferable
    public DataFlavor[] getTransferDataFlavors()
    {
      return flavors;
    }

    public boolean isDataFlavorSupported(DataFlavor f)
    {
      if(f.equals(FILENODE) || f.equals(DataFlavor.stringFlavor))
        return true;
      return false;
    }

    public Object getTransferData(DataFlavor d)
        throws UnsupportedFlavorException, IOException
    {
      if(d.equals(FILENODE))
        return this;
      else if(d.equals(DataFlavor.stringFlavor))
        return getFile().getAbsolutePath();
      else throw new UnsupportedFlavorException(d);
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
