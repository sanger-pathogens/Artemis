/* RunSamTools
 *
 * created: 2009
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2009  Genome Research Limited
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
package uk.ac.sanger.artemis.components.alignment;

import java.io.*;
import java.util.List;

import htsjdk.samtools.SAMRecord;

/**
* Used to run an samtools process this reads stdout and 
* stderr in separate threads.
*/
public class RunSamTools
{

  /** running process */
  private Process p;
  /** standard out */
  private List<SAMRecord> readsInView;
  /** standard error */
  private StringBuffer stderr = new StringBuffer();
  private StringBuffer stdout = new StringBuffer();
  private String initialIOError = null;
  
  /** running directory */
  //private File project;
  /** process status */
  private String status;
  private StdoutHandler stdouth;
  private StderrHandler stderrh;

  /**
  * @param cmd              command to run
  * @param env              environment
  * @param project          running directory
  */
  public RunSamTools(String cmd[], 
                     String[] envp,
                     File project,
                     List<SAMRecord> readsInView)
  {
    //this.project = project;
    this.readsInView = readsInView;
    status = "0";

    Runtime run = Runtime.getRuntime();
    try
    {
      p = run.exec(cmd,envp,project);

      // 2 threads to read in stdout & stderr buffers 
      // to prevent blocking
      stdouth = new StdoutHandler(this);
      stderrh = new StderrHandler(this);
      stdouth.start();
      stderrh.start();
    }
    catch(IOException ioe)
    {
      System.err.println("Error executing: "+
                          cmd[0]);
      initialIOError = ioe.getMessage();
      status = "1";
    }
  }

  /**
  * Read in the process stderr.
  */
  private void readProcessStderr()
  {
    BufferedInputStream stderrStream = null;
    BufferedReader stderrRead = null;
    try
    {
      //String line;
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
      System.err.println("RunEmbossApplication2: Error in "+
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
        System.err.println("RunEmbossApplication2: Error closing stream");
      }
      try
      {
        if(stderrRead!=null)
          stderrRead.close();
      }
      catch(IOException ioe)
      {
        System.err.println("RunEmbossApplication2: Error closing reader");
      }
    }

    return;
  }

  /**
  * Read in the process stdout.
  */
  private void readProcessStdout()
  {
    InputStreamReader stdoutStream = null;
    BufferedReader stdoutRead = null;
    try
    {
      //String line;
      stdoutStream =
         new InputStreamReader(p.getInputStream());
      stdoutRead =
         new BufferedReader(stdoutStream);    
      
      String line;
      while((line = stdoutRead.readLine()) != null)
      {
        if(readsInView != null)
        {
          String fields[] = line.split("\t");

          SAMRecord pread = new SAMRecord(null);
          pread.setReadName(fields[0]);
          pread.setFlags(Integer.parseInt(fields[1]));
          pread.setReferenceName(fields[2]);
          pread.setAlignmentStart(Integer.parseInt(fields[3]));
          pread.setMappingQuality(Integer.parseInt(fields[4]));
          pread.setCigarString(fields[5]);
          pread.setMateReferenceName(fields[6]);
          pread.setMateAlignmentStart( Integer.parseInt(fields[7]));
          pread.setInferredInsertSize(Integer.parseInt(fields[8]));
          pread.setReadString(fields[9]);
          
          readsInView.add(pread);
        }
        else
          stdout.append(line+"\n");
      }
    }
    catch (IOException io)
    {
      System.err.println("RunEmbossApplication2: Error in "+ 
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
        System.err.println("RunEmbossApplication2: Error closing stream");
      } 
      try
      {
        if(stdoutRead!=null)
          stdoutRead.close();
      }
      catch(IOException ioe)
      {
        System.err.println("RunEmbossApplication2: Error closing reader");
      }
    }
 
    return;
  }

  /**
   * @return standard out
   */
   public void waitForStdout()
   {
     try
     {
       // make sure we hang around for stdout
       while(stdouth.isAlive())
         Thread.sleep(10);
     }
     catch(InterruptedException ie)
     {
       ie.printStackTrace();
     }
                                                                                 
     return;
   }
   
   /**
    * @return standard out
    */
    public String getProcessStdout()
    {
      try
      {
        // make sure we hang around for stdout
        while(stdouth.isAlive())
          Thread.sleep(10);
      }
      catch(InterruptedException ie)
      {
        ie.printStackTrace();
      }
                                                                                  
      return stdout.toString();
    }
   
  /**
  * @return stderr
  */
  public String getProcessStderr()
  {
    try
    {
      // make sure we hang around for stdout
      while(stderrh.isAlive())
        Thread.sleep(10);
    }
    catch(InterruptedException ie)
    {
      ie.printStackTrace();
    }
                                                                                
    return new String(stderr.toString().trim());
  }

  /**
  * Wait for the process to end
  */
  public void waitFor()
  {
	try
	{
	  int exitVal = p.waitFor();
	  
	  if(exitVal != 0)
	    System.out.println("Exit value:: "+exitVal);
	}
	catch(InterruptedException ie)
	{
	  ie.printStackTrace();
	}
  }
  
  /**
  * @return process
  */
  public Process getProcess()
  {
    return p;
  }

  /**
  * @return status
  */
  public String getStatus()
  {
    return status;
  }

  class StdoutHandler extends Thread
  {
    RunSamTools rea;

    protected StdoutHandler(RunSamTools rea)
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
    RunSamTools rea;

    protected StderrHandler(RunSamTools rea)
    {
      this.rea = rea;
    }

    public void run()
    {
      rea.readProcessStderr();
    }
  }

  public String getInitialIOError() 
  {
    return initialIOError;
  }

}

