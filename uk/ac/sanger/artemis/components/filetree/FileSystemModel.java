/*
 *
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

package uk.ac.sanger.artemis.components.filetree;

import java.io.File;
import java.io.FileFilter;
import java.util.Date;
import java.text.SimpleDateFormat;
import java.awt.Cursor;

import javax.swing.JFrame;
import javax.swing.tree.TreeNode;
import javax.swing.tree.DefaultTreeModel;
import uk.ac.sanger.artemis.components.filetree.FileNode;
import uk.ac.sanger.artemis.components.filetree.RemoteFileNode;

/**
 *
 * FileSystemModel is a TreeTableModel representing a hierarchical file 
 * system. Nodes in the FileSystemModel are FileNodes which, when they 
 * are directory nodes, cache their children to avoid repeatedly querying 
 * the real file system. 
 * 
 */
public class FileSystemModel extends DefaultTreeModel 
{

  // Names of the columns.
  static protected String[]  cNames = {"Name", "Size", "Modified"};

  // Types of the columns.
  static protected Class[]  cTypes = {TreeTableModel.class, Integer.class, String.class};
  /** busy cursor */
  private Cursor cbusy = new Cursor(Cursor.WAIT_CURSOR);
  /** done cursor */
  private Cursor cdone = new Cursor(Cursor.DEFAULT_CURSOR);

  private FileFilter filter;
  private JFrame frame;
  SimpleDateFormat formatter  = new SimpleDateFormat("MMM dd HH:mm yyyy");

  public FileSystemModel(final JFrame frame) 
  {
    super(new FileNode(new File("")));
    this.frame = frame;

    filter = new FileFilter()
    {
      public boolean accept(File pathname)
      {
        if(pathname.isDirectory() ||
           !pathname.getName().startsWith("."))
          return true;
        return false;
      }
    };

    FileNode rootNode = (FileNode)getRoot();
    rootNode.setDirectory(true);

    rootNode.add( new FileNode(new File( System.getProperty("user.dir") )) );
    rootNode.add( new FileNode(new File( System.getProperty("user.home") )) );
  }

  public FileSystemModel(File rt[], FileFilter filter, JFrame frame)
  {
    super(new FileNode(new File("")));

    this.frame = frame;
    this.filter = filter;
    FileNode rootNode = (FileNode)getRoot();
    rootNode.setDirectory(true);

    File homeDir = new File(System.getProperty("user.home"));
    rootNode.add(new FileNode(homeDir));

    for(int i=0; i<rt.length; i++)
    {
      if(rt[i].compareTo(homeDir) != 0)
        rootNode.add( new FileNode(rt[i]) );
    }
  }

  public FileSystemModel(String froots[], final JFrame frame)
  {
    super(new RemoteFileNode(true));
    this.frame = frame;

    RemoteFileNode rootNode = (RemoteFileNode)getRoot();
    for(int i=0; i<froots.length; i++)
    {
      File f = new File(froots[i]);
      RemoteFileNode node = new RemoteFileNode(froots[i], f.getName(),
                                              null, null, true);
      rootNode.add(node);
    }
  }

  public FileFilter getFilter()
  {
    return filter;
  }


  public void setFilter(FileFilter filter)
  {
    this.filter = filter;
  }

  //
  // Some convenience methods. 
  //

  protected Object[] getChildren(Object node)
  {
    if(node instanceof FileNode)
    {
      FileNode fileNode = ((FileNode)node); 

      if(fileNode.getFile().getName().equals(""))
      {
        FileNode[] children = new FileNode[fileNode.getChildCount()];
        for(int i=0; i<children.length; i++)
          children[i] = (FileNode)fileNode.getChildAt(i);
        return children;
      }

      return fileNode.getChildren(filter); 
    }
    else
    {
      RemoteFileNode fileNode = ((RemoteFileNode)node);
      return fileNode.getChildren();
    }
  }

  //
  // The TreeModel interface
  //
  public int getChildCount(Object node) 
  { 
    if(node instanceof RemoteFileNode)
      frame.setCursor(cbusy);
    Object[] children = getChildren(node); 

    if(node instanceof RemoteFileNode)
      frame.setCursor(cdone);

    return (children == null) ? 0 : children.length;
  }

  public Object getChild(Object node, int i) 
  { 
    return getChildren(node)[i]; 
  }

  // The superclass's implementation would work, but this is more efficient. 
  public boolean isLeaf(Object node) 
  {
    return ((TreeNode)node).isLeaf(); 
  }

  //
  //  The TreeTableNode interface. 
  //

  public int getColumnCount()
  {
    return cNames.length;
  }

  public String getColumnName(int column) 
  {
    return cNames[column];
  }

  public Class getColumnClass(int column) 
  {
    return cTypes[column];
  }
 
  public Object getValueAt(Object node, int column)
  {
    if(node instanceof FileNode)
    {
      File file = ((FileNode)node).getFile(); 
      try
      {
        switch(column)
        {
          case 0:
            return file.getName();
          case 1:
            return file.isFile() ? new Integer((int)file.length()) : null;
//   case 2:
//return file.isFile() ?  "File" : "Directory";
          case 2:
            return file.exists() ? 
                   formatter.format(new Date(file.lastModified())) : null;
//                 DateFormat.getDateTimeInstance().format(new Date(file.lastModified())) : null;
        }
      }
      catch(SecurityException se) { }
    }
    else
    {
      try
      {
        switch(column)
        {
          case 0:
            return ((RemoteFileNode)node).getFile();
          case 1:
            return !((RemoteFileNode)node).isDirectory() ? 
                    ((RemoteFileNode)node).length() : null;
//   case 2:
//return file.isFile() ?  "File" : "Directory";
          case 2:
            return ((RemoteFileNode)node).getModifiedTime() != null ?
               formatter.format(((RemoteFileNode)node).getModifiedTime()) : null;
//             DateFormat.getDateTimeInstance().format(((RemoteFileNode)node).getModifiedTime()) : null;
        }
      }
      catch(SecurityException se) { }
    }
    
    return null; 
  }
}

