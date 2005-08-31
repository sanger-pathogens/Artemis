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
import uk.ac.sanger.artemis.j2ssh.SshFileManager;
import java.util.Vector;
import java.util.Collections;
import java.io.IOException;
import java.io.File;

public class FileList 
{
  /** vector containing directories */
  private Vector vdir;
  private Vector vfile;
  private static SshFileManager ssh_client = new SshFileManager();
 
  public FileList()
  {
  }

  /**
  *
  * Routine to obtain a directory listing distinguishing between
  * files and directories.
  * @param dir	remote directory to list
  *
  */
  protected void getDirList(String dir)
  {
    try
    {
      ssh_client.remoteList(dir);
    }
    catch(IOException ioe)
    {
      ioe.printStackTrace();
    }
    vdir  = ssh_client.getDirList(); 
    vfile = ssh_client.getFileList();
    Collections.sort(vfile);
  }

  /**
  *
  * Delete the given file / directory
  *
  */
  protected boolean delete(String file)
  {
    return ssh_client.delete(file);
  }

  /**
  *
  * Make a directory
  *
  */
  protected boolean mkdir(String dir)
  {
    return ssh_client.mkdir(dir);
  }


  /**
  *
  * Print working directory
  *
  */
  protected String pwd()
  {
    return ssh_client.pwd();
  }

  /**
  *
  * Rename a file
  *
  */
  protected boolean rename(String old_file, String new_file)
  {
    return ssh_client.rename(old_file, new_file);
  }


  /**
  *
  * Put a file
  *
  */
  protected boolean put(String dir, File local_file, FTProgress monitor)
  {
    return ssh_client.put(dir, local_file, monitor);
  }


  /**
  *
  * Get the file contents
  *
  */
  protected byte[] getFileContents(String file, FTProgress monitor)
  {
    return ssh_client.getFileContents(file, monitor);
  }


  /**
  *
  * Gets the list of files as a Vector
  * @return     list of files as a Vector
  *
  */
  public Vector fileVector() 
  {
    final Vector vfile_filtered = new Vector();
    for(int i=0; i<vfile.size(); i++)
    {
      String sfile = (String)vfile.get(i);
      if(!sfile.startsWith("."))
        vfile_filtered.add(sfile);
    }

    return vfile_filtered;
  }


  /**
  *
  * Gets whether this name is a directory
  * @return     true if it is a directory
  *
  */
  public boolean isDirectory(String d) 
  {
    return vdir.contains(d);
  }

  protected boolean isConnected()
  {
    return ssh_client.isConnected();
  }
}
