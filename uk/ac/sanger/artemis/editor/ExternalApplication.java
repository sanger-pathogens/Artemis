/*
 *
 * created: Wed Aug 3 2004
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2000  Genome Research Limited
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
 */

package uk.ac.sanger.artemis.editor;

import java.io.*;

/**
*
* Used to run an external command and get the stdout
* and stderr.
*
*/
public class ExternalApplication
{

  /** running process */
  private Process p;
  /** standard out */
  private StringBuffer stdout = new StringBuffer();
  /** standard error */
  private StringBuffer stderr = new StringBuffer();
  /** running directory */
  private File project;
  /** process status */
  private String status;
  private StdoutHandler stdouth;
  private StderrHandler stderrh;

  /**
  *
  * @param cmd            command to run
  * @param envp           environment
  * @param project        running directory
  *
  */
  public ExternalApplication(String[] cmd, 
                       String[] envp, File project)
  {
    this.project = project;
    status = "0";

    Runtime cmdRun = Runtime.getRuntime();
    try
    {
      p = cmdRun.exec(cmd,envp,project);

      // 2 threads to read in stdout & stderr buffers 
      // to prevent blocking
      stdouth = new StdoutHandler(this);
      stderrh = new StderrHandler(this);
      stdouth.start();
      stderrh.start();
    }
    catch(IOException ioe)
    {
      ioe.printStackTrace();
      System.out.println("ExternalApplication Error executing: "+
                          cmd);
      status = "1";
    }
  }


  /**
  *
  * Read in the process stderr.
  *
  */
  private void readProcessStderr()
  {

    BufferedInputStream stderrStream = null;
    BufferedReader stderrRead = null;
    try
    {
      String line;
      stderrStream =
         new BufferedInputStream(p.getErrorStream());
      stderrRead =
         new BufferedReader(new InputStreamReader(stderrStream));
      char c[] = new char[100];
      int nc = 0;

      while((nc = stderrRead.read(c,0,100)) != -1)
        stderr = stderr.append(new String(c,0,nc));
    }
    catch (IOException io)
    {
      System.err.println("ExternalApplication: Error in "+
                                "collecting standard out");
    }
    finally
    {
      try
      {
        if(stderrStream!=null)
          stderrStream.close();
      }
      catch(IOException ioe)
      {
        System.err.println("ExternalApplication: Error closing stream");
      }
      try
      {
        if(stderrRead!=null)
          stderrRead.close();
      }
      catch(IOException ioe)
      {
        System.err.println("ExternalApplication: Error closing reader");
      }
    }

    return;
  }

  /**
  *
  * Read in the process stdout.
  *
  */
  private void readProcessStdout()
  {
    
    BufferedInputStream stdoutStream = null;
    BufferedReader stdoutRead = null;
    try
    {
      String line;
      stdoutStream =
         new BufferedInputStream(p.getInputStream());
      stdoutRead =
         new BufferedReader(new InputStreamReader(stdoutStream));
 
      
      char c[] = new char[100];
      int nc = 0;
      String chunk;

      while((nc = stdoutRead.read(c,0,100)) != -1)
      {
        chunk  = new String(c,0,nc);
        stdout = stdout.append(chunk);
      }

    }
    catch (IOException io)
    {
      System.err.println("ExternalApplication: Error in "+ 
                                "collecting standard out");
    }
    finally
    {
      try
      {
        if(stdoutStream!=null)
          stdoutStream.close();
      }
      catch(IOException ioe)
      {
        System.err.println("ExternalApplication: Error closing stream");
      } 
      try
      {
        if(stdoutRead!=null)
          stdoutRead.close();
      }
      catch(IOException ioe)
      {
        System.err.println("ExternalApplication: Error closing reader");
      }
    }
 
    return;
  }

  /**
  *
  * Get the stdout for the process.
  * @return standard out.
  *
  */
  public String getProcessStdout()
  {
    try
    {
      // make sure we hang around for stdout
      while(stdouth.isAlive())
        Thread.currentThread().sleep(10);
    }
    catch(InterruptedException ie)
    {
      ie.printStackTrace();
    }
                                                                                
    return new String(stdout.toString().trim());
  }


  /**
  *
  * Get the stderr for the process.
  * @return standard error.
  *
  */
  public String getProcessStderr()
  {
    try
    {
      // make sure we hang around for stderr
      while(stderrh.isAlive())
        Thread.currentThread().sleep(10);
    }
    catch(InterruptedException ie)
    {
      ie.printStackTrace();
    }
                                                                                
    return new String(stderr.toString().trim());
  }

  /**
  *
  * Wait for the process to end
  *
  */
  public int waitFor()
  {
    try
    {
      return p.waitFor();
    }
    catch(InterruptedException ie)
    {
      ie.printStackTrace();
    }
    return -1;
  }

  /**
  *
  * @return process
  *
  */
  public Process getProcess()
  {
    return p;
  }

  /**
  *
  * @return status
  *
  */
  public String getStatus()
  {
    return status;
  }

  class StdoutHandler extends Thread
  {
    ExternalApplication rea;

    protected StdoutHandler(ExternalApplication rea)
    {
      this.rea = rea;
    }

    public void run()
    {
      rea.readProcessStdout();
    }
  }

  class StderrHandler extends Thread
  {
    ExternalApplication rea;

    protected StderrHandler(ExternalApplication rea)
    {
      this.rea = rea;
    }

    public void run()
    {
      rea.readProcessStderr();
    }
  }

}

