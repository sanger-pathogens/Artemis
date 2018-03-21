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
import uk.ac.sanger.artemis.j2ssh.SshPSUClient;

import java.util.Hashtable;
import java.util.zip.GZIPOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.File;

import com.sshtools.j2ssh.sftp.FileAttributes;

public class FileList 
{
  /** vector containing directories */
  private Hashtable vdir;
  protected static SshFileManager ssh_client;
 
  public FileList()
  {
    if(ssh_client == null)
      ssh_client = new SshFileManager();
  }

  public FileList(SshFileManager ssh_client)
  {
    FileList.ssh_client = ssh_client;
  }

  /**
  *
  * Routine to obtain a directory listing distinguishing between
  * files and directories.
  * @param dir	remote directory to list
  *
  */
  protected Hashtable getDirList(String dir)
  {
    try
    {
      ssh_client.remoteList(dir);
    }
    catch(IOException ioe)
    {
      return null;
    }
    vdir  = ssh_client.getDirList(); 
    return ssh_client.getFileList();
//  Collections.sort(vfile);
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
  protected boolean put(String dir, File local_file, 
                        FTProgress monitor, boolean force)
  {
    return ssh_client.put(dir, local_file, monitor, force);
  }

  /**
  *
  * @param name of file to get status for
  *
  */ 
  public FileAttributes stat(String filename)
  {
    return ssh_client.stat(filename);
  }


  /**
  *
  * Get the file contents
  *
  */
  public byte[] getFileContents(String file, FTProgress monitor)
  {
    return ssh_client.getFileContents(file, monitor);
  }
  
  /**
   * Transfer an entry from a zip file on the remote file system and
   * create a copy on the local file system.
   * @param zipFile
   * @param entry
   * @param dir_name
   * @return
   */
  public File getZipEntryContents(String zipFile, String entry, File dir_name)
  {  
    String stdOut = getZipEntry(zipFile, entry, "unzip -p "+zipFile+" "+entry+".gz | gunzip");
    if(stdOut == null || stdOut.equals(""))
      stdOut = getZipEntry(zipFile, entry, "unzip -p "+zipFile+" "+entry);
  
    if(stdOut == null || stdOut.equals(""))
      return null;
    // write to the local file system
    try
    {
      File localfn = new File(dir_name.getAbsoluteFile(), entry+".gz");
      localfn.getParentFile().mkdirs();
      if(localfn.createNewFile())
      {
        GZIPOutputStream out = new GZIPOutputStream(new FileOutputStream(localfn));
        out.write(stdOut.getBytes());
        out.close();
        return localfn;
      }
    }
    catch (IOException e)
    {
      System.err.println(e.getMessage());
    }
    return null;
  }
  
  private String getZipEntry(String zipFile, String entry, String cmd)
  {
    SshPSUClient sshClient = new SshPSUClient(cmd);
    sshClient.start();
    while(sshClient.isAlive() || sshClient.isAlive())
      try
      {
        Thread.sleep(2);
      }
      catch (InterruptedException e1){}

    //String stdErr = sshClient.getStdErr();
    return sshClient.getStdOut();
  }


  /**
  *
  * Gets whether this name is a directory
  * @return     true if it is a directory
  *
  */
  public boolean isDirectory(String d) 
  {
    return vdir.containsKey(d);
  }

  public static boolean isConnected()
  {
    if(ssh_client == null)
      return false;
    return ssh_client.isConnected();
  }
}
