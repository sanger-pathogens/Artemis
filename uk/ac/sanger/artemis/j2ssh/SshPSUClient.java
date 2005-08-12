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


import javax.swing.*;
import java.awt.*;
import java.io.*;
import com.sshtools.j2ssh.SshClient;
import com.sshtools.j2ssh.authentication.AuthenticationProtocolState;
import com.sshtools.j2ssh.authentication.PasswordAuthenticationClient;
import com.sshtools.j2ssh.io.UnsignedInteger32;
import com.sshtools.j2ssh.session.SessionChannelClient;
import com.sshtools.j2ssh.sftp.FileAttributes;
import com.sshtools.j2ssh.sftp.SftpFile;
import com.sshtools.j2ssh.sftp.SftpFileOutputStream;
import com.sshtools.j2ssh.SftpClient;
import com.sshtools.j2ssh.configuration.ConfigurationLoader;

import java.util.Vector;
import java.util.Properties;
import java.util.logging.FileHandler;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;

/**
 *
 */
public class SshPSUClient extends Thread
{

  private String hostname = null;
  private String user     = null;
  private String listfilepath = null;
  private String cmd      = null;
  private String bsub     = null;
  private String logfile  = null;
  private String db       = null;
  private String wdir     = "/nfs/pathscratch1/scratch";
  private int port        = -1;

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
  }

  public void run()
  {
    String program = cmd;
    // get properties from j2ssh.properties
    Properties settings = getProperties();
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
    if(cmd.equals("blastp") && settings.getProperty("blastp") != null)
      cmd = settings.getProperty("blastp");
    if(cmd.equals("fasta") && settings.getProperty("fasta") != null)
      cmd = settings.getProperty("fasta");

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

      ConfigurationLoader.initialize(false);

//    if(hostname == null)
//      hostname = (String)JOptionPane.showInputDialog(
//                           null, "Name of server machine:",
//                           "Server hostname",
//                           JOptionPane.PLAIN_MESSAGE, null,
//                           null, null);
     
      BufferedReader reader =
          new BufferedReader(new InputStreamReader(System.in));

      // Create a password authentication instance
      PasswordAuthenticationClient pwd = new PasswordAuthenticationClient();
      // Get the users name

      JPanel promptPanel = new JPanel(new GridLayout(4,2));

      JTextField hostfield  = new JTextField(16);
      if(hostname != null)
        hostfield.setText(hostname);
      
      JTextField portfield  = new JTextField(16);
      if(port >-1)
        portfield.setText(Integer.toString(port));

      JTextField ufield  = new JTextField(16);
      if(user != null)
        ufield.setText(user);
      JPasswordField pfield = new JPasswordField(16);

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
     
      Object[] options = { "CANCEL", "LOGIN"};

      int select = JOptionPane.showOptionDialog(null, promptPanel,
                              "LOGIN",
                               JOptionPane.YES_NO_CANCEL_OPTION,
                               JOptionPane.QUESTION_MESSAGE,
                               null,
                               options,
                               options[1]);

      if(select == 0)
        return;

      // Make a client connection
      SshClient ssh = new SshClient();
      hostname = hostfield.getText().trim();
      if(portfield.getText().trim().equals(""))
        port = -1;
      else
        port = Integer.parseInt(portfield.getText().trim());

      // Connect to the host
      if(port < 0)
        ssh.connect(hostname);
      else
        ssh.connect(hostname,port);

      user = ufield.getText().trim();
      pwd.setUsername(user);
      pwd.setPassword(new String(pfield.getPassword()));

      // Try the authentication
      int result = ssh.authenticate(pwd);
   
      // Evaluate the result
      if(result == AuthenticationProtocolState.COMPLETE)
      {
        if(listfilepath == null)
        {
          JFileChooser chooser = new JFileChooser();
          int returnVal = chooser.showOpenDialog(null);
          if(returnVal == JFileChooser.APPROVE_OPTION) 
            listfilepath = chooser.getSelectedFile().getAbsolutePath();
          else
            return;
        }

        SftpClient sftp = ssh.openSftpClient();
        Vector seqfile = readListFile(listfilepath);
        for(int i=0; i<seqfile.size();i++)
        {
          String filepath = (String)seqfile.get(i);
          int index = filepath.lastIndexOf(System.getProperty("file.separator"));
          String filename = filepath;
          if(index > -1)
            filename = filename.substring(index+1);

          wdir = wdir+"/"+user+"/";
          try
          {
            sftp.mkdir(wdir);
            sftp.put(filepath, wdir+filename);
          }
          catch(IOException ioe)
          {}
      
          SessionChannelClient session = ssh.openSessionChannel();

          String outputfile = wdir+filename+".out";
          final String actualCMD;

          if(cmd.indexOf("fasta33") > -1)
          {
            if(settings.getProperty(db) != null)
              db = settings.getProperty(db);
            actualCMD = bsub+" -o "+ outputfile +" -e "+ outputfile + ".err " +
                           cmd+" "+wdir+filename+" "+db;
          }
          else
            actualCMD = bsub+" -o "+ outputfile +" -e "+ outputfile + ".err " +
                           cmd+" "+db+" "+wdir+filename;

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

            int count = 0;

            while(!isFile && count < 100)
            {
              Thread.currentThread().sleep(5000);
              Object list[] = sftp.ls(wdir).toArray();

              for(int j=0; j<list.length;j++)
              {
                if(((SftpFile)list[j]).getFilename().equals(filename+".out"))
                  isFile = true;
              } 
              count++;
            }
          }
          catch(InterruptedException ie)
          {
            ie.printStackTrace();
          }
          catch(java.io.IOException ioe)
          {
            ioe.printStackTrace();
          }
          
          // stdout
          System.out.println(stdouth.getOutput());
          System.out.println(stderrh.getOutput());

//        ByteArrayOutputStream os = new ByteArrayOutputStream();
//        sftp.get(outputfile, os);
//        System.out.println(os.toString());

          sftp.get(outputfile, filepath+".out");

          session.close();
        }

        // Quit
        sftp.quit();
        ssh.disconnect();
      }
      else 
        JOptionPane.showMessageDialog(null, 
            "Problem logging in!\nCheck username and password.",
            "Authentication Problem",
            JOptionPane.ERROR_MESSAGE);


    } catch(IOException ioe){}
    finally
    {
      JOptionPane.showMessageDialog(null,
            "Finished \n" + program,
            "Process Finished",
            JOptionPane.INFORMATION_MESSAGE);
    }
  }

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
    return settings;
  }


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
