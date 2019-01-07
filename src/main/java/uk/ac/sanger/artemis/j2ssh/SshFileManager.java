/* ExternalProgram.java
 *
 * created: Aug 2005
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
 **/

package uk.ac.sanger.artemis.j2ssh;

import javax.swing.JOptionPane;

import java.io.File;
import java.io.ByteArrayOutputStream;
import java.io.IOException;

import java.util.Hashtable;

import com.sshtools.j2ssh.SshClient;
import com.sshtools.j2ssh.sftp.FileAttributes;
import com.sshtools.j2ssh.sftp.SftpFile;
import com.sshtools.j2ssh.SftpClient;
import com.sshtools.j2ssh.SshException;

/**
*
* Client to use ssh connection to server to run blast/fasta
* remotely. 
*
*/
public class SshFileManager
{
  private Hashtable dir_list;
  private Hashtable file_list;
  private SshClient ssh;

  public SshFileManager()
  {
    SshLogin sshLogin = new SshLogin();
    ssh = sshLogin.getSshClient();
  }

  public SshFileManager(SshLogin sshLogin)
  {
    ssh = sshLogin.getSshClient(true);
    pwd();
  }

  private void rescue()
  {
    try
    {
      ssh.disconnect();
      SshLogin sshLogin = new SshLogin();
      ssh = sshLogin.getSshClient();
    }
    catch(Exception exp)
    {  
      System.out.println("SshFileManager.rescue()");
      exp.printStackTrace();
    }
  }


  /**
  *
  * Return an active SftpClient object
  *
  */
  private SftpClient getSftpClient()
             throws IOException
  {
    SftpClient sftp;  
    try
    {
      if(!ssh.hasActiveSftpClient())
        return ssh.openSftpClient(); 
 
      sftp = ssh.getActiveSftpClient();
      return sftp; 
    }
    catch(IOException ioe)
    {} 
    return ssh.openSftpClient();
  }

  /**
  *
  * Remote directory listing
  *
  */
  public boolean remoteList(String remoteRootDir)
                    throws IOException
  {
    SftpClient sftp = getSftpClient();
    Object list[] = null;

    try
    {
      list = sftp.ls(remoteRootDir).toArray();     
    }
    catch(java.io.FileNotFoundException fnf)
    {
      return false;
    }

    dir_list  = new Hashtable();
    file_list = new Hashtable(); 

    for(int i=0; i < list.length; i++)
    {
      SftpFile sfile = (SftpFile)list[i];
      FileAttributes fat = sfile.getAttributes();

      if(sfile.isDirectory() || sfile.isLink())
        dir_list.put(sfile.getFilename(), fat);

      file_list.put(sfile.getFilename(), fat);
    }
     
    return true;
  }

  /**
  *
  * Delete file or directory
  *
  */
  public boolean delete(String filename)
  {
    try
    {
      SftpClient sftp = getSftpClient();
      sftp.rm(filename);
//    sftp.quit();
    }
    catch(IOException ioe)
    {
      rescue();
      ioe.printStackTrace();
      return false;
    }

    return true;
  }

  /**
  *
  * Make directory
  *
  */
  public boolean mkdir(String dir)
  {
    try
    {
      SftpClient sftp = getSftpClient();
      sftp.mkdir(dir);
//    sftp.quit();
    }
    catch(IOException ioe)
    {
      rescue();
      ioe.printStackTrace();
      return false;
    }

    return true;
  }

  /**
  *
  * Get working directory
  *
  */
  public String pwd()
  {
    String pwd = null;
    try
    {
      SftpClient sftp = getSftpClient();
      pwd = sftp.pwd();
//    sftp.quit();
    }
    catch(SshException exp)
    {
      JOptionPane.showMessageDialog(null,
               "Cannot start SSH session.\n",
               "SSH Error",
               JOptionPane.ERROR_MESSAGE);
    }
    catch(IOException ioe)
    {
      rescue();
      ioe.printStackTrace();
    }

    return pwd;
  }

  /**
  *
  *
  *
  */
  public boolean rename(String old_file, String new_file)
  {
    try
    {
      SftpClient sftp = getSftpClient();
      sftp.rename(old_file,new_file);
//    sftp.quit();
    }
    catch(IOException ioe)
    {
      rescue();
//    ioe.printStackTrace();
      return false;
    }

    return true;
  }

  /**
  *
  * @param name of file to get status for
  *
  */
  public FileAttributes stat(String filename) 
  {
    SftpClient sftp = null;
    try
    {
      sftp = getSftpClient();
      FileAttributes fat = sftp.stat(filename);
      return fat;
    }
    catch(SshException sshe)
    {
      rescue();
      return stat(filename);
    }  
    catch(IOException ioe)
    {
      return null;
      // remote file doesn't exist
    }

  }

  /**
  *
  * @param dir name of directory to write the file to
  * @param local_file file to copy to the server
  *
  */
  public boolean put(final String dir, final File local_file,
                     final FTProgress monitor, final boolean force)
  { 
    SftpClient sftp = null;

    try
    {
      sftp = getSftpClient();
      FileAttributes attr = sftp.stat(dir+"/"+local_file.getName());

      if(attr.isDirectory())
      {
        JOptionPane.showMessageDialog(null,
               "Cannot overwrite the directory\n"+
               dir+"/"+local_file.getName(),"Cannot Overwrite",
               JOptionPane.ERROR_MESSAGE);
        return false;
      }
      else if(!force)
      {
        int n = JOptionPane.showConfirmDialog(null,
               "Overwrite\n"+
               dir+"/"+local_file.getName() + "\n?",
               "Confirm the sequence entry",
               JOptionPane.YES_NO_OPTION);
        if(n == JOptionPane.NO_OPTION)
          return false;
      }      
    }
    catch(SshException sshe)
    {
      if(System.getProperty("debug") != null)
      {
        System.out.println("put() ");
        sshe.printStackTrace();
      }

      rescue();
      return put(dir, local_file, monitor, force);
    }
    catch(IOException ioe)
    {
      // remote file doesn't exist
    }

    if(sftp == null)
      return false;

     
    try
    {
      sftp.put(local_file.getCanonicalPath(), 
         dir+"/"+local_file.getName(), monitor);
      return true;
    }
    catch(SshException sshe)
    {
      if(System.getProperty("debug") != null)
      {
        System.out.println("put() 2");
        sshe.printStackTrace();
      }
      rescue();
      return put(dir, local_file, monitor, true);
    }
    catch(IOException ioe)
    {
      rescue();
      ioe.printStackTrace();
      return false;
    }

  }


  /**
  *
  * Return the file contents as a byte array
  *
  */
  public byte[] getFileContents(String file, final FTProgress monitor)
  {
    // sometimes SftpClient.get() hangs so put in 
    // separate thread

    SshGet get = new SshGet(file, monitor);
    get.start();

    try
    {
      int count = 0;
      while(get.isAlive())
      {
        Thread.currentThread().sleep(50);
        count++;
        if(count > 100 && monitor.getProgress() < 1.)
        {
          get.destroy();
          SshLogin sshLogin = new SshLogin();
          ssh = sshLogin.getSshClient();
          return getFileContents(file,monitor);
        }
      }
    }
    catch(InterruptedException iex){}

    return get.getByteArray();
  }


  public Hashtable getFileList()
  {
    return file_list;
  }

  public Hashtable getDirList()
  {
    return dir_list;
  }
 
  public boolean isConnected()
  {
    if(ssh == null || !ssh.isConnected())
      return false;

    return true;
  }

  public class SshGet extends Thread
  {
    private String file;
    private FTProgress monitor;
    private ByteArrayOutputStream os;
    private boolean done = false;

    public SshGet(final String file, final FTProgress monitor)
    {
      this.file = file;
      this.monitor = monitor; 
    }
  
    public void run()
    {
      while(!done)
      {
        try
        {
          SftpClient sftp = getSftpClient();

          os = new ByteArrayOutputStream();
          sftp.get(file, os, monitor);
          done = true;
          os.close();
        }
        catch(SshException se2)
        {
          rescue();
        }
        catch(IOException ioe2)
        {
          done = true; // file doesn't exist?
        }
      }
    }

    public void destroy()
    {
      done = true;
    }

    public byte[] getByteArray()
    {
      return os.toByteArray();
    }
  }
}
