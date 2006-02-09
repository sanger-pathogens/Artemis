/* SshLogin.java
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
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;
import javax.swing.JPasswordField;
import javax.swing.SwingConstants;

import java.awt.Dimension;
import java.awt.GridLayout;
import java.io.IOException;
import java.util.Properties;

import java.util.logging.FileHandler;
import java.util.logging.Handler;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;

import com.sshtools.j2ssh.SshClient;
import com.sshtools.j2ssh.authentication.AuthenticationProtocolState;
import com.sshtools.j2ssh.authentication.PasswordAuthenticationClient;
import com.sshtools.j2ssh.transport.IgnoreHostKeyVerification;


/**
*
* Client to use ssh connection to server to run blast/fasta
* remotely. 
*
*/
public class SshLogin
{
  // defaults
  private String logfile  = null;

  // login variables
  private String hostname = null;
  private String user     = null;
  private int port        = -1;
  private static JPasswordField pfield = new JPasswordField(16);
  private static JTextField portfield  = new JTextField(16);
  private static JTextField hostfield  = new JTextField(16);
  private static JTextField ufield  = new JTextField(16);
  private static SshClient ssh;
  private static Properties settings;

  public SshLogin()
  {
    try
    {
      logfile = System.getProperty("logfile");
      // Setup a logfile
      if(logfile != null)
      {
        Handler fh = new FileHandler(logfile);
        fh.setFormatter(new SimpleFormatter());
        Logger.getLogger("com.sshtools").setUseParentHandlers(false);
        Logger.getLogger("com.sshtools").addHandler(fh);
        Logger.getLogger("com.sshtools").setLevel(Level.WARNING);
      }
      else
        Logger.getLogger("com.sshtools").setLevel(Level.OFF);
    }
    catch(IOException ioe){}

   if(settings == null)
      settings = setProperties();
  }

  public SshClient getSshClient()
  {
    return getSshClient(false);
  }

  public SshClient getSshClient(final boolean fail)
  {
    if(ssh == null || !ssh.isConnected())
    {
      try
      {
        ssh = login(fail);
      }
      catch(java.net.ConnectException ce)
      {
        ce.printStackTrace(); 
      }
      catch(IOException ioe){}
    }
    return ssh;
  }

  public String getUser()
  {
    return ufield.getText().trim();
  }

  public static String getHostname()
  {
    return hostfield.getText().trim();
  }

  public static String getPort()
  {
    return portfield.getText().trim();
  } 

  public static Properties getProperties()
  {
    return settings;
  }

  /**
  *
  * Log the user in.
  *
  */
  private SshClient login(final boolean fail)
            throws IOException
  {
    SshClient ssh = null;
    int result = AuthenticationProtocolState.FAILED;
    int count  = 0;

    while(result != AuthenticationProtocolState.COMPLETE)
    {
      if( !(count == 0 && pfield.getPassword().length > 0) ) 
      {
        if(!setLogin())
          return null;
      }

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
        ssh.connect(hostname, new IgnoreHostKeyVerification());
      else
        ssh.connect(hostname,port, new IgnoreHostKeyVerification());

      // Try the authentication
      result = ssh.authenticate(pwd);
      if(fail && result == AuthenticationProtocolState.FAILED)
        return null;

      count++;
    }
    return ssh;
  }

  public JPanel getLogin()
  {
    JPanel promptPanel = new JPanel(new GridLayout(4,2));

    if(hostname != null && hostfield.getText().equals(""))
      hostfield.setText(hostname);

    if(port >-1 && portfield.getText().equals(""))
      portfield.setText(Integer.toString(port));

    if(user != null && ufield.getText().equals(""))
      ufield.setText(user);

    JLabel hostlab = new JLabel(" Hostname:  ", SwingConstants.RIGHT);
    JLabel portlab = new JLabel("     Port:  ", SwingConstants.RIGHT);

    JLabel ulab = new JLabel(" Username:  ", SwingConstants.RIGHT);
    JLabel plab = new JLabel(" Password:  ", SwingConstants.RIGHT);

    //add labels etc
    promptPanel.add(hostlab);
    promptPanel.add(hostfield);

    promptPanel.add(portlab);
    promptPanel.add(portfield);

    promptPanel.add(ulab);
    promptPanel.add(ufield);

    promptPanel.add(plab);
    promptPanel.add(pfield);
    return promptPanel;
  }

  /**
  *
  * Get password field
  *
  */
  public JPasswordField getJPasswordField()
  {
    return pfield;
  }

  /**
  *
  * Set the login information.
  *
  */
  private boolean setLogin()
  {
    JPanel promptPanel = getLogin();

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
  private Properties setProperties()
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
    if(user == null)
      user = System.getProperty("user.name");

    return settings;
  }

}
