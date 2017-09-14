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

import javax.swing.JFileChooser;

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
import com.sshtools.j2ssh.sftp.SftpFile;
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
  // defaults
  private String listfilepath = null;
  private String cmd      = null;
  private String bsub     = null;
  private String logfile  = null;
  private String db       = null;
  private String wdir     = null;
  private boolean justProg = false;
  
  //
  private SshClient ssh;
  private String user;
  private boolean keep = false;
  private boolean zipResults = false;
  
  //
  StdoutStdErrHandler stdouth;
  StdoutStdErrHandler stderrh;

  public SshPSUClient(String args[])
  {
    // process arguments
    if(args != null && args.length > 0)
    {
      for(int i=0; i<args.length; i++)
      {
        if(args[i].equals("-f") && i < args.length-1)
          listfilepath = args[i+1];
        else if(args[i].equals("-cmd") && i < args.length-1)
          cmd = args[i+1];
        else if(args[i].equals("-bsub") && i < args.length-1)
          bsub = args[i+1];
        else if(args[i].equals("-l") && i < args.length-1)
          logfile = args[i+1];
        else if(args[i].equals("-d") && i < args.length-1)
          db = args[i+1];
        else if(args[i].equals("-wdir") && i < args.length-1)
          wdir = args[i+1];
        else if(args[i].equals("-keep"))
          keep = true;
      }
    }

    SshLogin sshLogin = new SshLogin();
    ssh = sshLogin.getSshClient();
    user = sshLogin.getUser();
  }
  
  public SshPSUClient(final String cmd)
  {
    this.cmd = cmd;
    SshLogin sshLogin = new SshLogin();
    ssh = sshLogin.getSshClient();
    justProg = true;
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

      if(justProg)
        runProgram();
      else
        completed = runBlastOrFasta(program);

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
  *  Wait until a file appears on the server.  
  *
  */
  private boolean waitUntilFileAppears(String file)
                   throws InterruptedException, IOException
  {
    for(int i=0; i < 500; i++)
    {
      logger4j.debug("waitUntilFileAppears() "+file);
      Thread.sleep(1000);
      try
      {
        if(fileExists(getSftpClient() , wdir, file))
          return true;
      }
      catch(SshException sshe)
      {
        if(System.getProperty("debug") != null)
        {
          logger4j.warn("waitUntilFileAppears()");
          sshe.printStackTrace();
        }
        try
        {
          rescue();
          continue;
        } catch(Exception exp) {}
      }
    }

    return false;
  }

  private boolean fileExists(SftpClient sftp, String wdir, String file)
             throws SshException, IOException
  {
    Object list[] = null;
    try
    {
      list = sftp.ls(wdir).toArray();
    }
    catch(SshException sshe)
    {
      sftp = getSftpClient();
      list = sftp.ls(wdir).toArray();
    }

    for(int j=0; j<list.length;j++)
    {
      if( ((SftpFile)list[j]).getFilename().equals(file) )
        return true;
    }
    return false;
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

    if(bsub == null && settings.getProperty("bsub") != null)
      bsub = settings.getProperty("bsub");
    if(db == null)
    {
      if(settings.getProperty("default_db") != null)
        db = settings.getProperty("default_db");
      else
        db = "%uniprot";
    } 

    if(settings.getProperty("zip") != null)
    {
      String zipValue = settings.getProperty("zip");
      zipResults = Boolean.parseBoolean(zipValue);
      
      logger4j.debug("zip results :: "+zipResults);
    }
    
    if(wdir == null && settings.getProperty("wdir") != null)
      wdir = settings.getProperty("wdir");

    if(cmd != null)
    {
      if(cmd.equals("blastp") && settings.getProperty("blastp") != null)
        cmd = settings.getProperty("blastp");
      else if(cmd.equals("blastn") && settings.getProperty("blastn") != null)
        cmd = settings.getProperty("blastn");
      else if(cmd.equals("blastx") && settings.getProperty("blastx") != null)
        cmd = settings.getProperty("blastx");
      else if(cmd.equals("tblastx") && settings.getProperty("tblastx") != null)
        cmd = settings.getProperty("tblastx"); 
      else if(cmd.equals("fasta") && settings.getProperty("fasta") != null) 
        cmd = settings.getProperty("fasta");
      else if(cmd.equals("fastx") && settings.getProperty("fastx") != null)
        cmd = settings.getProperty("fastx");
    }

    return settings;
  }

 
  /**
  *
  * Run fasta or blast on the server ssh'ed into
  *
  */
  private boolean runBlastOrFasta(String program)
                    throws IOException
  {
    Properties settings = getProperties();

    // prompt for local listfile
    if(listfilepath == null)
    {
      JFileChooser chooser = new JFileChooser();
      int returnVal = chooser.showOpenDialog(null);
      if(returnVal == JFileChooser.APPROVE_OPTION) 
        listfilepath = chooser.getSelectedFile().getAbsolutePath();
      else
        return false;
    }

    SftpClient sftp = getSftpClient();

    // loop over sequence files in the listfile
    Vector seqfile = readListFile(listfilepath);
    for(int i=0; i<seqfile.size();i++)
    {
      String filepath = (String)seqfile.get(i);
      int index = filepath.lastIndexOf(System.getProperty("file.separator"));
      String filename = filepath;
      if(index > -1)
        filename = filename.substring(index+1);

      if(i == 0)
      {
        try
        {
          if(wdir.endsWith("scratch118") || wdir.endsWith("scratch118/"))
          {
            if(fileExists(sftp , wdir+"/bacteria/", user))
              wdir = wdir+"/bacteria/";
            else if(fileExists(sftp , wdir+"/parasites/", user))
              wdir = wdir+"/parasites/";
            else if(fileExists(sftp , wdir+"/pathogen/", user))
              wdir = wdir+"/pathogen/";
            else if(fileExists(sftp , wdir+"/viruses/", user))
              wdir = wdir+"/viruses/";
          }
          
          if(!keep)
            wdir = wdir + "/" + user;
          sftp.mkdir(wdir);
          wdir = wdir + "/" + program + "/";

          sftp.mkdir(wdir);
          logger4j.debug("mkdir() " + wdir);
          // sftp.put(filepath, wdir+filename);
        }
        catch(SshException sshe)
        {
          logger4j.debug("runBlastOrFasta()");
          if(System.getProperty("debug") != null)
          {
            sshe.printStackTrace();
          }
          rescue();
          sftp = getSftpClient();
          if(!wdir.endsWith(program + "/"))
            wdir = wdir + "/" + program + "/";
        }
        catch(IOException ioe)
        {
          // directory already created
        }
      }
 

      try
      {
        sftp.put(filepath, wdir+filename);

        logger4j.debug("PUT SUCCESS "+wdir+filename);
      }
      catch(SshException ioe)
      {
        logger4j.debug("runBlastOrFasta() - 2");
        if(System.getProperty("debug") != null)
        {
          ioe.printStackTrace();
        }
        rescue();
        sftp = getSftpClient();
        sftp.put(filepath, wdir+filename);
      }

      logger4j.debug("STARTING session");

      SessionChannelClient session = null;

      try 
      {
        if(!ssh.isConnected())
          rescue();

        session = ssh.openSessionChannel();
      }
      catch(IOException exp)
      {
        logger4j.debug("NOT STARTED runBlastOrFasta() ----- 3 "+filename);
        if(System.getProperty("debug") != null)
        {
          exp.printStackTrace();
        }
        rescue();
      }

      String outputfile = wdir+filename+".out";
      final String actualCMD;
     
      if(bsub == null)
      {
        if( ((cmd.indexOf("fasta3") > -1) || (cmd.indexOf("fastx3") > -1))
            && settings.getProperty(db) != null)
          db = settings.getProperty(db);
        else if(db.startsWith("%"))
          db = db.substring(1,db.length());

        if( (cmd.indexOf("fasta3") > -1) ||
            (cmd.indexOf("fastx3") > -1) )
          actualCMD = cmd+" "+wdir+filename+" "+db+" > "+outputfile;
        else
          actualCMD = cmd+" -d "+db+" -i "+wdir+filename+" -o "+outputfile;
      }
      else
      {
        if( (cmd.indexOf("fasta3") > -1) ||
            (cmd.indexOf("fastx3") > -1) )
        {
          if(settings.getProperty(db) != null)
            db = settings.getProperty(db);
          actualCMD = bsub+" -o "+ outputfile +" -e "+ outputfile + ".err " +
                       cmd+" "+wdir+filename+" "+db;
        }
        else
          actualCMD = bsub+" -o "+ outputfile +" -e "+ outputfile + ".err " +
                       cmd+" "+db+" "+wdir+filename;
      }

      // run the application
      logger4j.debug(actualCMD);

      try
      {
        session.executeCommand(actualCMD);
      }
      catch(IOException exp)
      {
        logger4j.debug("runBlastOrFasta() - 3");
        if(System.getProperty("debug") != null)
        {
          exp.printStackTrace();
        }
      }

      logger4j.debug("STARTED session "+filename);

      // Reading from the session InputStream
      StdoutStdErrHandler stdouth = new StdoutStdErrHandler(session, true);
      StdoutStdErrHandler stderrh = new StdoutStdErrHandler(session, false);
    
      stdouth.start();
      stderrh.start();

      boolean isFile = false;
      try
      {
        // make sure we hang around for stdout
        while(stdouth.isAlive() || stderrh.isAlive())
          Thread.sleep(10);

        isFile = waitUntilFileAppears(filename+".out");
      }
      catch(InterruptedException ie)
      {
        ie.printStackTrace();
      }
       
      // stdout & stderr
      logger4j.debug("STDOUT "+filename+"\n"+stdouth.getOutput());
      logger4j.debug("STDERR "+filename+"\n"+stderrh.getOutput());


      try
      {
        sftp = getSftpClient();
        sftp.get(outputfile, filepath+".out");
      }
      catch(Exception ioe)
      {
        logger4j.debug("runBlastOrFasta() - 3");
        if(System.getProperty("debug") != null)
        {
          ioe.printStackTrace();
        }
        rescue();
        sftp = getSftpClient();
        sftp.get(outputfile, filepath+".out");
      }

      logger4j.debug("GET SUCCESS "+filepath+".out");
      sftp.rm(wdir+filename);
      
      if(!keep)
      {
        sftp.rm(outputfile);
      }
      else if(zipResults)
      {
        cmd = "gzip "+outputfile+"; zip -j "+wdir+program+".zip "+
                                    outputfile+".gz; rm -f "+outputfile+".gz";
        logger4j.debug(cmd);
        SshPSUClient sshClient = new SshPSUClient(cmd);
        sshClient.start();
      }
      sftp = getSftpClient();
      sftp.rm(outputfile+".err");
      session.close();
    }

    return true;
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


  /**
   * The main program for the PasswordConnect class
   *
   * @param args The command line arguments
   */
  public static void main(String args[]) 
  {
    new SshPSUClient(args);
  }
}
