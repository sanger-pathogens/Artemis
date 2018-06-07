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

import org.apache.log4j.Logger;

import java.io.File;
import java.io.ByteArrayOutputStream;
import java.io.IOException;

import java.util.Hashtable;

import com.sshtools.ssh.ChannelOpenException;
import com.sshtools.ssh.SshClient;
import com.sshtools.sftp.SftpFileAttributes;
import com.sshtools.sftp.SftpStatusException;
import com.sshtools.sftp.TransferCancelledException;
import com.sshtools.sftp.SftpFile;
import com.sshtools.sftp.SftpClient;
import com.sshtools.ssh.SshException;

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
  
  private int retryIntervalSecs = 5;
  private int numRetries = 10;
  
  /** Logging instance. */
  private static org.apache.log4j.Logger logger = Logger.getLogger(SshFileManager.class);

  /**
   * Constructor.
   * Creates a new SshLogin.
   */
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
  
  /**
   * Exit the given sftp client and just log any errors.
   * @param sftpClient SftpClient
   */
  public void closeSftpClient(SftpClient sftpClient)
  {
	  try 
	  {
		  if (sftpClient != null) 
		  {
			  sftpClient.exit();
		  }
	  }
	  catch (Exception e)
	  {
		  logger.error("Failed to cleanly close down sftp client: " + e.getMessage());
	  }
  }

  /**
   * Attempt to clean-up after a connection drop and then
   * reconnect.
   */
  private void rescue()
  {
    try
    {
      ssh.disconnect();
    }
    catch(Exception exp)
    {  
      logger.warn("rescue(): Disconnect failed", exp);
    }
    
    try
    {
      SshLogin sshLogin = new SshLogin();
      ssh = sshLogin.getSshClient();
    }
    catch(Exception exp)
    {  
      logger.error("rescue(): Failed to create new SSH connection", exp);
    }
  }

  /** 
   * Determine if the given status exception indicates that a file
   * does not exists.
   * @param e SftpStatusException
   * @return true if file does not exist
   */
  private boolean fileDoesNotExist(SftpStatusException e)
  {
	return (
		e.getStatus() == SftpStatusException.SSH_FX_NO_SUCH_FILE ||
		e.getStatus() == SftpStatusException.SSH_FX_NO_SUCH_PATH ||
		e.getStatus() == SftpStatusException.SSH_FX_FILE_IS_A_DIRECTORY);  
  }
  
  /** 
   * Determine if the given status exception indicates that a directory
   * does not exists.
   * @param e SftpStatusException
   * @return true if directory does not exist
   */
  private boolean directoryDoesNotExist(SftpStatusException e)
  {
	return (
		e.getStatus() == SftpStatusException.SSH_FX_NO_SUCH_PATH ||
		e.getStatus() == SftpStatusException.SSH_FX_FILE_IS_A_DIRECTORY);  
  }
  
  /**
  *
  * Return an active SftpClient object
  * @return SftpClient
  */
  private SftpClient getSftpClient()
             throws SshException, SftpStatusException, ChannelOpenException
  {
	try
    {
        return new SftpClient(ssh);
    }
    catch(Exception e)
    {
    		// connection may already be corrupted, so we try one more time...
    		logger.error("Failed to create sftp client connection. Retrying.", e);
    }
    
    return new SftpClient(ssh);
  }
  
  /**
   * Convert the given exception to an IOException.
   * Added to minimise changes for legacy client code.
   * @param e Exception
   * @return IOException
   */
  private IOException convertToIoException(Exception e)
  {
	  return new IOException(e.getMessage(), e);
  }

  /**
   * Remote directory listing
   * @param remoteRootDir String
   */
  public boolean remoteList(String remoteRootDir)
                    throws IOException
  {
    
    Object list[] = null;
    SftpClient sftp = null;

    try
    {
    	  sftp = getSftpClient();
      list = sftp.ls(remoteRootDir);
      
      dir_list  = new Hashtable();
      file_list = new Hashtable(); 

      for(int i=0; i < list.length; i++)
      {
        SftpFile sfile = (SftpFile)list[i];
        SftpFileAttributes fat = sfile.getAttributes();

        if(sfile.isDirectory() || sfile.isLink())
          dir_list.put(sfile.getFilename(), fat);

        file_list.put(sfile.getFilename(), fat);
      }
    }
    catch(SftpStatusException sse)
    {
    	  logger.warn("Caught SftpStatusException for remote directory listing, with status " + sse.getStatus() + ": " + sse.getMessage());
	  if (directoryDoesNotExist(sse))
	  {
	    return false;
	  }
	  else
	  {
		throw convertToIoException(sse);
	  }
    }
    catch (SshException | ChannelOpenException e)
    {
    	  logger.error("Caught exception while attempting to list remote directory contents" , e);
    	  rescue();
      throw convertToIoException(e);
    }
    finally
    {
    	  // cleanup
  	  closeSftpClient(sftp);
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
	SftpClient sftp = null;
			  
    try
    {
      sftp = getSftpClient();
      sftp.rm(filename);
    }
    catch(ChannelOpenException | SshException | SftpStatusException e)
    {
    	  logger.error("Caught exception while attempting to delete remote file", e);
      rescue();
      return false;
    }
    finally
    {
    	  // cleanup
  	  closeSftpClient(sftp);
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
	SftpClient sftp = null;
	  
    try
    {
      sftp = getSftpClient();
      sftp.mkdir(dir);
    }
    catch(ChannelOpenException | SshException | SftpStatusException e)
    {
    	  logger.error("Caught exception while attempting to create new remote directory", e);
      rescue();
      return false;
    }
    finally
    {
    	  // cleanup
  	  closeSftpClient(sftp);
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
    SftpClient sftp = null;
    
    try
    {
      sftp = getSftpClient();
      pwd = sftp.pwd();
    }
    catch(ChannelOpenException | SshException | SftpStatusException e)
    {
    	  logger.error("Caught exception while attempting to execute a remote pwd" , e);
    	
      JOptionPane.showMessageDialog(null,
               "Cannot start SSH session.\n",
               "SSH Error",
               JOptionPane.ERROR_MESSAGE);
      
      rescue();
    }
    finally
    {
    	  // cleanup
  	  closeSftpClient(sftp);
    }
    

    return pwd;
  }

  /**
  *
  * Rename a given file.
  *
  */
  public boolean rename(String old_file, String new_file)
  {
	SftpClient sftp = null;
	
    try
    {
      sftp = getSftpClient();
      sftp.rename(old_file,new_file);
    }
    catch(ChannelOpenException | SshException | SftpStatusException e)
    {
    	  logger.error("Caught exception while attempting to rename a remote file", e);
      rescue();
      return false;
    }
    finally
    {
    	  // cleanup
  	  closeSftpClient(sftp);
    }

    return true;
  }

  /**
  *
  * @param name of file to get status for
  *
  */
  public SftpFileAttributes stat(String filename) 
  {
	  // TODO KEV - has dodgy recursion
	  
    SftpClient sftp = null;
    try
    {
      sftp = getSftpClient();
      SftpFileAttributes fat = sftp.stat(filename);
      return fat;
    }
    catch(ChannelOpenException | SshException e)
    {
    	  logger.error("Caught exception while calling stat", e);
      rescue();
      // TODO FIX dodgy recursion
      //return stat(filename);
    }  
    catch(SftpStatusException sse)
    {
    	  if (fileDoesNotExist(sse))
    	  {
    		  return null;
    	  }
    	  else
    	  {
    		  rescue();
    		  // TODO FIX dodgy recursion
    	      //return stat(filename);
    	  }
    }
    finally
    {
    	  // cleanup
  	  closeSftpClient(sftp);
    }

      // TODO KEV added
	  return null;
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
	// TODO KEV dodgy recursion
	  
    SftpClient sftp = null;

    try
    {
	    try
	    {
	      sftp = getSftpClient();
	      SftpFileAttributes attr = sftp.stat(dir+"/"+local_file.getName());
	
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
	    catch(ChannelOpenException | SshException e)
	    {
	      logger.error("Attempt to connect to remote directory failed. Reconnecting. Error message: " + e.getMessage());
	      rescue();
	      
	      // TODO FIX THIS RECURSION
	      //return put(dir, local_file, monitor, force);
	    }
	    catch(SftpStatusException sse)
	    {
	      if (!fileDoesNotExist(sse))
	      {
	    	  	if(System.getProperty("debug") != null)
	        {
	            sse.printStackTrace();
	        }
	    	  	else
	    	  	{
	    	  	  rescue();
	    	  	  // TODO FIX THIS RECURSION
	          //return put(dir, local_file, monitor, force);
	    	  	}
	      }
	    }
	
	    if(sftp == null)
	      return false;
	
	     
	    try
	    {
	      sftp.put(local_file.getCanonicalPath(), 
	         dir+"/"+local_file.getName(), monitor);
	      return true;
	    }
	    catch(SshException | TransferCancelledException e)
	    {
	    	  logger.error("Attempt to put file to remote host failed.  Error message: " + e.getMessage());
	      rescue();
	      // TODO FIX THIS RECURSION
	      //return put(dir, local_file, monitor, true);
	    }
	    catch(SftpStatusException sse)
	    {
	      rescue();
	      sse.printStackTrace();
	    }
	    catch(IOException ioe)
	    {
	      rescue();
	      ioe.printStackTrace();
	    }
	    
	    return false;
	    
    }
	finally
	{
	   // cleanup
	   closeSftpClient(sftp);
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
    	  SftpClient sftp = null;
    	
      while(!done)
      {
        try
        {
          sftp = getSftpClient();

          os = new ByteArrayOutputStream();
          sftp.get(file, os, monitor);
          done = true;
          os.close();
        }
        catch(ChannelOpenException | SshException | TransferCancelledException e)
        {
        	  logger.error("Failed to connect to remote host. Error message: " + e.getMessage());
          rescue();
        }
        catch(SftpStatusException ste)
        {
        	  done = true; // e.g. file doesn't exist
        }
        catch (IOException ioe)
        {
        	  done = true;
        }
        finally 
        {
        	  // cleanup
        	  closeSftpClient(sftp);
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
