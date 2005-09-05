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

import uk.ac.sanger.artemis.j2ssh.FTProgress;
import java.awt.datatransfer.*;
import javax.swing.tree.*;
import java.io.*;
import java.util.*;

/**
*
* File node for remote file tree manager
*
*/
public class RemoteFileNode extends DefaultMutableTreeNode 
                    implements Transferable, Serializable
{
  /** true if node is explored */
  private boolean explored = false;
  /** true if node is a directory */
  private boolean isDir = false;
  /** full name of node */
  private String fullname;
  /** path to the file on the server */
  private String serverPathToFile;
  /** root directory */
  private String rootdir; 

  /** parent directory listing */
  private transient FileList parentList;        // make transient for
  /** remote server file roots */
  private transient String froots;
  /** file separator for server files */
  private String fs = "/";

  final public static DataFlavor REMOTEFILENODE = 
         new DataFlavor(RemoteFileNode.class, "Remote file");
  static DataFlavor remoteFlavors[] = { REMOTEFILENODE, 
                             DataFlavor.stringFlavor };

  public RemoteFileNode(boolean isDir)
  {
    this.isDir = isDir;
  }

  /**
  *
  * @param froots		remote server file roots
  * @param file		file for this node
  * @param parentList		parent directory listing
  * @param parent		parent to this node
  *
  */
  public RemoteFileNode(String froots, String file,
                        FileList parentList, String parent)
  {
    this(froots, file, parentList, parent, false);
  }

  /**
  *
  * @param froots             remote server file roots
  * @param file               file for this node
  * @param parentList         parent directory listing
  * @param parent             parent to this node
  * @param ldir		true if the node is a directory
  *
  */
  public RemoteFileNode(String froots, String file,
                        FileList parentList, String parent,
                        boolean ldir)
  { 
    this.froots = froots;
    this.parentList = parentList;
    isDir = ldir;
    rootdir = froots;
    serverPathToFile = froots;

    if(parent != null)
    {
      if(parent.endsWith("/."))
        parent = parent.substring(0,parent.length()-1);
      else if(parent.endsWith(fs))
        parent = parent.substring(0,parent.length());

      if(parent.equals("."))
        fullname = file;
      else
      {
        fullname = parent + fs + file;
        if(serverPathToFile.endsWith(fs))
          serverPathToFile = serverPathToFile.concat(parent);
        else
          serverPathToFile = serverPathToFile.concat(fs+parent);
      }
    }
    
    if(parentList != null)
    {
      if(parentList.isDirectory(file))
        isDir = true;
    }
    else if(parent == null)
      fullname = ".";

    setUserObject(file); 
  }
   
  /** @return		true if node is a directory */
  public boolean getAllowsChildren() { return isDir; }
  /** @return         true if node is a file */
  public boolean isLeaf() { return !isDir; }
  /** @return         true if node is a directory */
  public boolean isDirectory() { return isDir; }
  /** @return         the node name */
  public String getFile() { return (String)getUserObject(); }
  /** @return         root directory */
  public String getRootDir() { return rootdir; }
  /** @return         full name of the node */
  public String getFullName() { return fullname; }
  /** @return         path on server */
  public String getPathName() { return serverPathToFile; }
  /** @return         true if explored */
  public boolean isExplored() { return explored; }
   
  /**
  *
  * Get the server name 
  * @return 	server name
  *
  */
  public String getServerName() 
  { 
    String prefix = serverPathToFile;
    if(!prefix.endsWith(fs))
      prefix = prefix.concat(fs);

    if(fullname.equals("."))
      return prefix;

    return prefix + (String)getUserObject();
  }

  /**
  *
  * Explore the node and add new child nodes 
  *
  */
  public void explore() 
  {
    if(!isDir)
      return;

    if(!explored)
    {
      FileList flist = new FileList();

      String dir;
      if(getRootDir().equals(""))
        dir = new String("~/"+getFullName());
      else
        dir = new String(getRootDir()+"/"+getFullName());
      dir = dir.trim();

      flist.getDirList(dir);
      Vector children = flist.fileVector();
      for(int i=0;i<children.size();i++)
        add(new RemoteFileNode(froots,(String)children.get(i),
                               flist,fullname));
    }
      explored = true;
  }

  protected boolean delete()
  {
    FileList flist = new FileList();
    return flist.delete(getRootDir()+"/"+getFullName());
  }

  protected boolean mkdir(String dir)
  {
    FileList flist = new FileList();
    return flist.mkdir(dir);
  }

  protected boolean rename(String new_file)
  {
    FileList flist = new FileList();
    return flist.rename(getRootDir()+"/"+getFullName(), new_file);
  }

  public boolean put(File local_file, FTProgress monitor)
  {
    FileList flist = new FileList();
    String dir;
    if(!isDirectory())
    {
      dir = getRootDir()+"/"+getFullName();
      int index = dir.lastIndexOf("/");
      dir = dir.substring(0,index);
    }
    else
      dir = getRootDir()+"/"+getFullName();

    return flist.put(dir, local_file, monitor);
  }


  public byte[] getFileContents(FTProgress monitor)
  {
    FileList flist = new FileList();
    return flist.getFileContents(getRootDir()+"/"+getFullName(), monitor);
  }

// Transferable
  public DataFlavor[] getTransferDataFlavors()
  {
    return remoteFlavors;
  }

  public boolean isDataFlavorSupported(DataFlavor f)
  {
    if(f.equals(REMOTEFILENODE) || f.equals(DataFlavor.stringFlavor))
      return true;
    return false;
  }

  public Object getTransferData(DataFlavor d) 
      throws UnsupportedFlavorException, IOException
  {
    if(d.equals(REMOTEFILENODE))
      return this;
    else if(d.equals(DataFlavor.stringFlavor))
      return getServerName();
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

