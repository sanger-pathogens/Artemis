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
public class SshPSUClient extends Thread
{
  // defaults
  private String listfilepath = null;
  private String cmd      = null;
  private String bsub     = null;
  private String logfile  = null;
  private String db       = null;
  private String wdir     = "/nfs/pathscratch1/scratch";

  // login variables
  private String hostname = null;
  private String user     = null;
  private int port        = -1;
  private static JPasswordField pfield = new JPasswordField(16);
  private static JTextField portfield  = new JTextField(16);
  private static JTextField hostfield  = new JTextField(16);
  private static JTextField ufield  = new JTextField(16);

  public SshPSUClient(String args[])
  {
    // process arguments
    if(args != null && args.length > 0)
    {
      for(int i=0; i<args.length; i++)
      {
        if(args[i].equals("-h") && i < args.length-1)
          hostname = args[i+1];
        else if(args[i].equals("-u") && i < args.length-1)
          user = args[i+1];
        else if(args[i].equals("-f") && i < args.length-1)
          listfilepath = args[i+1];
        else if(args[i].equals("-cmd") && i < args.length-1)
          cmd = args[i+1];
        else if(args[i].equals("-bsub") && i < args.length-1)
          bsub = args[i+1];
        else if(args[i].equals("-l") && i < args.length-1)
          logfile = args[i+1];
        else if(args[i].equals("-p") && i < args.length-1)
          port = Integer.parseInt(args[i+1]);
        else if(args[i].equals("-d") && i < args.length-1)
          db = args[i+1];
        else if(args[i].equals("-wdir") && i < args.length-1)
          wdir = args[i+1];
      }
    }

    try
    {
      // Setup a logfile
      if(logfile != null)
      {
        Handler fh = new FileHandler(logfile);
        fh.setFormatter(new SimpleFormatter());
        Logger.getLogger("com.sshtools").setUseParentHandlers(false);
        Logger.getLogger("com.sshtools").addHandler(fh);
        Logger.getLogger("com.sshtools").setLevel(Level.ALL);
      }
      else
        Logger.getLogger("com.sshtools").setLevel(Level.OFF);
    }
    catch(IOException ioe){}
  }

  public void run()
  {
    String program = cmd;
    // get properties from j2ssh.properties
    Properties settings = getProperties();

    boolean completed = false;
    try
    {
      ConfigurationLoader.initialize(false);

      SshClient ssh = login();
      if(ssh == null)
        return;

      completed = runBlastOrFasta(ssh, program, settings);

      // Quit
      ssh.disconnect();
    }
    catch(IOException ioe){}
    finally
    {
      if(completed)
        JOptionPane.showMessageDialog(null,
            "Finished \n" + program,
            "Process Finished",
            JOptionPane.INFORMATION_MESSAGE);
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
    }
    catch (IOException e)
    {
      System.out.println("Problem reading list file");
    }
    return seqfiles;
  }

  /**
  *
  *  Wait until a file appears on the server.  
  *
  */
  private boolean waitUntilFileAppears(SftpClient sftp, String file)
                   throws InterruptedException, IOException
  {
    for(int i=0; i < 100; i++)
    {
      Thread.currentThread().sleep(5000);
      Object list[] = sftp.ls(wdir).toArray();

      for(int j=0; j<list.length;j++)
      {
        if( ((SftpFile)list[j]).getFilename().equals(file) )
          return true;
      }
    }

    return false;
  }

  /**
  *
  * Log the user in.
  *
  */
  private SshClient login()
            throws IOException
  {
    SshClient ssh = null;
    int result = AuthenticationProtocolState.FAILED;

    while(result != AuthenticationProtocolState.COMPLETE)
    {
      if(!setLogin())
        return null;

      // Create a password authentication instance
      PasswordAuthenticationClient pwd = new PasswordAuthenticationClient();
      user = ufield.getText().trim();
      pwd.setUsername(user);
      pwd.setPassword(new String(pfield.getPassword()));

      // Make a client connection
      ssh = new SshClient();
      hostname = hostfield.getText().trim();
      if(portfield.getText().trim().equals(""))
        port = -1;
      else
        port = Integer.parseInt(portfield.getText().trim());

      if(port < 0)
        ssh.connect(hostname);
      else
        ssh.connect(hostname,port);

      // Try the authentication
      result = ssh.authenticate(pwd);
    }
    return ssh;
  }

  /**
  *
  * Set the login information.
  *
  */
  private boolean setLogin()
  {
    JPanel promptPanel = new JPanel(new GridLayout(4,2));

    if(hostname != null && hostfield.getText().equals(""))
      hostfield.setText(hostname);
   
    if(port >-1 && portfield.getText().equals(""))
      portfield.setText(Integer.toString(port));

    if(user != null && ufield.getText().equals(""))
      ufield.setText(user);

    JLabel hostlab = new JLabel(" Hostname:", SwingConstants.LEFT);
    JLabel portlab = new JLabel("     Port:", SwingConstants.LEFT);

    JLabel ulab = new JLabel(" Username:", SwingConstants.LEFT);
    JLabel plab = new JLabel(" Password:", SwingConstants.LEFT);
    //add labels etc
    promptPanel.add(hostlab);
    promptPanel.add(hostfield);

    promptPanel.add(portlab);
    promptPanel.add(portfield);

    promptPanel.add(ulab);
    promptPanel.add(ufield);

    promptPanel.add(plab);
    promptPanel.add(pfield);
  
    Object[] options = { "CANCEL", "LOGIN AND RUN"};

    int select = JOptionPane.showOptionDialog(null, promptPanel,
                          "LOGIN",
                           JOptionPane.YES_NO_CANCEL_OPTION,
                           JOptionPane.QUESTION_MESSAGE,
                           null,
                           options,
                           options[1]);

    if(select == 0)
      return false;

    return true;
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

    if(hostname == null && settings.getProperty("host") != null)
      hostname = settings.getProperty("host");
    if(port < 0 && settings.getProperty("port") != null)
      port = Integer.parseInt(settings.getProperty("port"));
    if(bsub == null && settings.getProperty("bsub") != null)
      bsub = settings.getProperty("bsub");
    if(user == null)
      user = System.getProperty("user.name");
    if(db == null)
    {
      if(settings.getProperty("default_db") != null)
        db = settings.getProperty("default_db");
      else
        db = "%uniprot";
    } 

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
  private boolean runBlastOrFasta(SshClient ssh, String program, Properties settings)
                    throws IOException
  {
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

    SftpClient sftp = ssh.openSftpClient();

    // loop over sequence files in the listfile
    Vector seqfile = readListFile(listfilepath);
    for(int i=0; i<seqfile.size();i++)
    {
      String filepath = (String)seqfile.get(i);
      int index = filepath.lastIndexOf(System.getProperty("file.separator"));
      String filename = filepath;
      if(index > -1)
        filename = filename.substring(index+1);

      try
      {
        wdir = wdir+"/"+user;
        sftp.mkdir(wdir);
        wdir = wdir+"/"+program+"/";
        sftp.mkdir(wdir);

        sftp.put(filepath, wdir+filename);
      }
      catch(IOException ioe)
      {}
  
      SessionChannelClient session = ssh.openSessionChannel();

      String outputfile = wdir+filename+".out";
      final String actualCMD;

      if( (cmd.indexOf("fasta33") > -1) ||
          (cmd.indexOf("fastx33") > -1) )
      {
        if(settings.getProperty(db) != null)
          db = settings.getProperty(db);
        actualCMD = bsub+" -o "+ outputfile +" -e "+ outputfile + ".err " +
                       cmd+" "+wdir+filename+" "+db;
      }
      else
        actualCMD = bsub+" -o "+ outputfile +" -e "+ outputfile + ".err " +
                       cmd+" "+db+" "+wdir+filename;

      // run the application
      if(System.getProperty("debug") != null)
        System.out.println(actualCMD);
      session.executeCommand(actualCMD);

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
          Thread.currentThread().sleep(10);

        isFile = waitUntilFileAppears(sftp, filename+".out");
      }
      catch(InterruptedException ie)
      {
        ie.printStackTrace();
      }
        
      if(System.getProperty("debug") != null)
      {
        // stdout & stderr
        System.out.println(stdouth.getOutput());
        System.out.println(stderrh.getOutput());
      }

//    ByteArrayOutputStream os = new ByteArrayOutputStream();
//    sftp.get(outputfile, os);
//    System.out.println(os.toString());

      sftp.get(outputfile, filepath+".out");
      session.close();
    }
    sftp.quit();

    return true;
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

    protected StdoutStdErrHandler(SessionChannelClient session, 
                                  boolean isStdout)
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

        byte buffer[] = new byte[100];
        int read;
        while((read = in.read(buffer)) > 0)
          buff.append(new String(buffer, 0, read));
      }
      catch(IOException ioe){}
    }

    public String getOutput()
    {
      return buff.toString();
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
