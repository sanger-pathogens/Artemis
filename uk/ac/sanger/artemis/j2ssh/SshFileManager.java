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

import uk.ac.sanger.artemis.components.SwingWorker;

import javax.swing.JOptionPane;
import javax.swing.JFileChooser;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.JPasswordField;
import javax.swing.SwingConstants;

import java.awt.GridLayout;
import java.io.File;
import java.io.FileReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.IOException;

import java.util.Date;
import java.util.Hashtable;
import java.util.Vector;
import java.util.Properties;
import java.util.logging.FileHandler;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;

import com.sshtools.j2ssh.io.UnsignedInteger32;
import com.sshtools.j2ssh.SshClient;
import com.sshtools.j2ssh.authentication.AuthenticationProtocolState;
import com.sshtools.j2ssh.authentication.PasswordAuthenticationClient;
import com.sshtools.j2ssh.session.SessionChannelClient;
import com.sshtools.j2ssh.sftp.FileAttributes;
import com.sshtools.j2ssh.sftp.SftpFile;
import com.sshtools.j2ssh.sftp.SftpFileOutputStream;
import com.sshtools.j2ssh.SftpClient;
import com.sshtools.j2ssh.configuration.ConfigurationLoader;
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
      sftp = ssh.getActiveSftpClient();
      return sftp;
    }
    catch(IOException ioe){}
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

    try
    {
      sftp.cd(remoteRootDir);
    }
    catch(java.io.FileNotFoundException fnf)
    {
      return false;
    }

    Object list[] = sftp.ls().toArray();
    
    dir_list  = new Hashtable();
    file_list = new Hashtable(); 

    for(int i=0; i < list.length; i++)
    {
      SftpFile sfile = (SftpFile)list[i];
      FileAttributes fat = sfile.getAttributes();
//    String modTime = fat.getModTimeString();
      long modTime = fat.getModifiedTime().longValue()*1000;

      if(sfile.isDirectory() || sfile.isLink())
        dir_list.put(sfile.getFilename(), new Date(modTime));

      file_list.put(sfile.getFilename(), new Date(modTime));
    }
     
//  sftp.quit();
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

  private boolean putTransfer(final SftpClient sftp, 
                              final String dir, final File local_file)
  {
    
    SwingWorker progressWorker = new SwingWorker()
    {
      public Object construct()
      {
        try
        {
          FileTransferProgressMonitor monitor =
              new FileTransferProgressMonitor(null);
          FTProgress progress = monitor.add(local_file.getName());

          sftp.put(local_file.getCanonicalPath(),
             dir+"/"+local_file.getName(), progress);
          return new Boolean(true);
        }
        catch(IOException ioe)
        {
          rescue();
          ioe.printStackTrace();
          return new Boolean(false);
        }

      }
    };
    progressWorker.start();

 
    return true;
//  return ((Boolean)progressWorker.get()).booleanValue();
  }

  /**
  *
  * Return the file contents as a byte array
  *
  */
  public byte[] getFileContents(String file, final FTProgress monitor)
  {
    try 
    {
      SftpClient sftp = getSftpClient();

      ByteArrayOutputStream os = new ByteArrayOutputStream();
      sftp.get(file, os, monitor);
//    sftp.quit();
      return os.toByteArray();
    }
    catch(SshException se)
    {
      rescue();
      return getFileContents(file, monitor);
    }
    catch(IOException ioe)
    {
      rescue();
//    ioe.printStackTrace();
      return null;
    }
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
}
