/* SshPSUClient.java
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

import uk.ac.sanger.artemis.components.MessageDialog;


import java.io.File;
import java.io.FileReader;
import java.io.InputStream;
import java.io.BufferedReader;
import java.io.IOException;

import java.util.Vector;
import java.util.Properties;

import com.sshtools.j2ssh.SshException;
import com.sshtools.j2ssh.SshClient;
import com.sshtools.j2ssh.session.SessionChannelClient;
import com.sshtools.j2ssh.SftpClient;
import com.sshtools.j2ssh.configuration.ConfigurationLoader;


/**
*
* Client to use ssh connection to server to run blast/fasta
* remotely. 
*
*/
public class SshPSUClient extends Thread
{
  public static org.apache.log4j.Logger logger4j = 
      org.apache.log4j.Logger.getLogger(SshPSUClient.class);
  
  private String cmd = null;
  private SshClient ssh;
  private boolean zipResults = false;
  
  //
  StdoutStdErrHandler stdouth;
  StdoutStdErrHandler stderrh;

  
  public SshPSUClient(final String cmd)
  {
    this.cmd = cmd;
    SshLogin sshLogin = new SshLogin();
    ssh = sshLogin.getSshClient();
  }

  private SshClient rescue()
  {
    try
    {
      ssh.disconnect();
      SshLogin sshLogin = new SshLogin();
      ssh = sshLogin.getSshClient();
    }
    catch(Exception exp)
    {
      logger4j.warn("SshPSUClient.rescue()");
      exp.printStackTrace();
    }
    return ssh;
  }

  public void run()
  {
    String program = cmd;
    boolean completed = false;
    try
    {
      ConfigurationLoader.initialize(false);

      if(ssh == null)
        return;

      runProgram();

      // Quit
      //ssh.disconnect();
    }
    catch(IOException ioe){}
    finally
    {
      if(completed)
        new MessageDialog(null,
            "Finished \n" + program, 
            "Process Finished",
            false);
    }
  }


  /**
  *
  * Read the sequence filenames in a list file 
  * @param String file	list filename
  * @return the sequence filename collection
  * 
  */
  private Vector readListFile(String file)
  {
    Vector seqfiles = new Vector();
    try
    {
      String line;
      BufferedReader in = new BufferedReader(new FileReader(file));
      while((line = in.readLine()) != null )
      {
        File seq = new File(line);
        if(seq.exists())
        {
          seqfiles.add(seq.getAbsolutePath());
        }
      }
      in.close();
    }
    catch (IOException e)
    {
      logger4j.warn("Problem reading list file");
    }
    return seqfiles;
  }

  /**
  *
  * Get the properties from the j2ssh.properties file.
  *
  */
  private Properties getProperties()
  {
    Properties settings = new Properties();
    ClassLoader cl = this.getClass().getClassLoader();
    // try out of the classpath
    try
    {
      settings.load(cl.getResourceAsStream("j2ssh.properties"));
    }
    catch (Exception e)
    {
    }

    if(settings.getProperty("zip") != null)
    {
      String zipValue = settings.getProperty("zip");
      zipResults = Boolean.parseBoolean(zipValue);
      
      logger4j.debug("zip results :: "+zipResults);
    }
    
    return settings;
  }

  /**
  *
  * Run fasta or blast on the server ssh'ed into
  *
  */
  private boolean runProgram()
                    throws IOException
  {
    SessionChannelClient session = null;

    try 
    {
      if(!ssh.isConnected())
        rescue();

      session = ssh.openSessionChannel();
    }
    catch(IOException exp)
    {
      logger4j.warn("NOT STARTED runProgram()");
      if(System.getProperty("debug") != null)
      {
        exp.printStackTrace();
      }
      rescue();
    }
    
    // run the application

    try
    {
      session.executeCommand(cmd);
    }
    catch(IOException exp)
    {
      logger4j.warn("session   : "+cmd);
      logger4j.warn("runProgram: "+exp.getMessage());
      if(System.getProperty("debug") != null)
      {
        exp.printStackTrace();
      }
    }

    logger4j.debug("STARTED session "+cmd);

    // Reading from the session InputStream
    stdouth = new StdoutStdErrHandler(session, true);
    stderrh = new StdoutStdErrHandler(session, false);
    
    stdouth.start();
    stderrh.start();

    try
    {
      // make sure we hang around for stdout
      while(stdouth.isAlive() || stderrh.isAlive())
        Thread.sleep(2);
    }
    catch(InterruptedException ie)
    {
      ie.printStackTrace();
    }
    
    // stdout & stderr
    if(stderrh.getBufferSize() > 0)
      logger4j.debug("STDERR :"+stderrh.getOutput());

    session.close();
    return true;
  }

  
  public StringBuffer getStdOutBuffer()
  {
    return stdouth.getStdOutBuffer();
  }
  
  public String getStdOut()
  {
    return stdouth.getOutput();
  }

  public String getStdErr()
  {
    return stderrh.getOutput();
  }
  
  /**
  *
  * Return an active SftpClient object
  *
  */
  private synchronized SftpClient getSftpClient()
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
    catch(SshException sshe)
    {
      try
      {
        sftp = ssh.openSftpClient();
        return sftp;
      }
      catch(IOException ioe)
      {
        logger4j.debug("getSftpClient() -- 2");
        if(System.getProperty("debug") != null)
        {
          ioe.printStackTrace();
        }
        rescue();
      }
    }
    catch(IOException ioe2)
    { 
      logger4j.debug("getSftpClient()");
      if(System.getProperty("debug") != null)
      {
        ioe2.printStackTrace();
      }
      rescue();
    }
    return ssh.openSftpClient();
  }

  /**
  *
  * Thread to handle stdout/stderr reading without blocking.
  *
  */
  class StdoutStdErrHandler extends Thread
  {
    private SessionChannelClient session;
    private boolean isStdout;
    private StringBuffer buff = new StringBuffer();

    protected StdoutStdErrHandler(final SessionChannelClient session, 
                                  final boolean isStdout)
    {
      this.session  = session;
      this.isStdout = isStdout;
    }

    public void run()
    {
      try
      {
        final InputStream in;
        if(isStdout)
          in = session.getInputStream();
        else
          in = session.getStderrInputStream();

        final byte buffer[] = new byte[400];
        int read;
        while((read = in.read(buffer)) > 0)
          buff.append(new String(buffer, 0, read));

        in.close();
      }
      catch(Exception ioe)
      {
        ioe.printStackTrace();
      }
    }
    
    public synchronized StringBuffer getStdOutBuffer()
    {
      return buff;
    }

    public synchronized String getOutput()
    {
      return buff.toString();
    }
    
    public int getBufferSize()
    {
      return buff.length();
    }
  }

}
