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

import java.util.Vector;
import java.util.Properties;
import java.util.logging.FileHandler;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;

import com.sshtools.j2ssh.SshClient;
import com.sshtools.j2ssh.authentication.AuthenticationProtocolState;
import com.sshtools.j2ssh.authentication.PasswordAuthenticationClient;
import com.sshtools.j2ssh.session.SessionChannelClient;
import com.sshtools.j2ssh.sftp.FileAttributes;
import com.sshtools.j2ssh.sftp.SftpFile;
import com.sshtools.j2ssh.sftp.SftpFileOutputStream;
import com.sshtools.j2ssh.SftpClient;
import com.sshtools.j2ssh.configuration.ConfigurationLoader;


/**
*
* Client to use ssh connection to server to run blast/fasta
* remotely. 
*
*/
public class SshFileManager
{

  private Vector dir_list;
  private Vector file_list;

  private SshClient ssh;

  public SshFileManager()
  {
    SshLogin sshLogin = new SshLogin();
    ssh = sshLogin.getSshClient();
  }

  private void rescue()
  {
    ssh.disconnect();
    SshLogin sshLogin = new SshLogin();
    ssh = sshLogin.getSshClient();
  }

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
    }

    Object list[] = sftp.ls().toArray();
    
    dir_list  = new Vector();
    file_list = new Vector();

    for(int i=0; i < list.length; i++)
    {
      SftpFile sfile = (SftpFile)list[i];
      if(sfile.isDirectory())
        dir_list.add(sfile.getFilename());

      file_list.add(sfile.getFilename());
    }
     
//  sftp.quit();
    return true;
  }

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

  public String pwd()
  {
    String pwd = null;
    try
    {
      SftpClient sftp = getSftpClient();
      pwd = sftp.pwd();
//    sftp.quit();
    }
    catch(IOException ioe)
    {
      rescue();
      ioe.printStackTrace();
    }

    return pwd;
  }

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
  public boolean put(String dir, File local_file)
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
      else
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
    catch(IOException ioe)
    {
      // remote file doesn't exist
    }

    if(sftp == null)
      return false;

    try
    {
      sftp.put(local_file.getCanonicalPath(), 
               dir+"/"+local_file.getName());
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
  * Return the file contents as a byte array
  *
  */
  public byte[] getFileContents(String file)
  {
    try 
    {
      SftpClient sftp = getSftpClient();

      ByteArrayOutputStream os = new ByteArrayOutputStream();
      sftp.get(file, os);
//    sftp.quit();
      return os.toByteArray();
    }
    catch(IOException ioe)
    {
      rescue();
      ioe.printStackTrace();
      return null;
    }
  }


  public Vector getFileList()
  {
    return file_list;
  }

  public Vector getDirList()
  {
    return dir_list;
  }
 
}
