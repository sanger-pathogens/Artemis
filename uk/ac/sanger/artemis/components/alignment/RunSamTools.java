
package uk.ac.sanger.artemis.components.alignment;

import java.io.*;
import javax.swing.JTextArea;

/**
* Used to run an samtools process this reads stdout and 
* stderr in separate threads.
*/
public class RunSamTools
{

  /** running process */
  private Process p;
  /** standard out */
  private StringBuffer stdout = new StringBuffer();
  /** standard error */
  private StringBuffer stderr = new StringBuffer();
  
  private String initialIOError = null;
  
  /** running directory */
  //private File project;
  /** process status */
  private String status;
  private StdoutHandler stdouth;
  private StderrHandler stderrh;
  private JTextArea textArea;


  /**
  * @param cmd              command to run
  * @param env              environment
  * @param project          running directory
  */
  public RunSamTools(String cmd[], 
                     String[] envp, File project)
  {
    //this.project = project;
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
    BufferedInputStream stdoutStream = null;
    BufferedReader stdoutRead = null;
    try
    {
      //String line;
      stdoutStream =
         new BufferedInputStream(p.getInputStream());
      stdoutRead =
         new BufferedReader(new InputStreamReader(stdoutStream));
    
      char c[] = new char[200];
      int nc = 0;
      String chunk;

      while((nc = stdoutRead.read(c,0,200)) != -1)
      {
        chunk  = new String(c,0,nc);
        stdout = stdout.append(chunk);
        if(textArea != null)
        {
          textArea.append(chunk);
          textArea.setCaretPosition(textArea.getDocument().getLength());
        }
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
  * This method can be called after the process has completed
  * to write the stdout to the project directory.
  */
//  private void writeStdout()
//  {
//    if(project != null)
//    {
//      PrintWriter out = null;
//      String fname = "";
//      try
//      {
//        fname = project.getCanonicalPath() + 
//                             "/application_stdout";
//        File so = new File(fname);
//
//        if(!so.exists())
//          so.createNewFile();
//
//        out = new PrintWriter(new FileWriter(fname,true));
//        out.println(stdout);
//      }
//      catch(IOException ioe)
//      {
//        System.err.println("RunEmbossApplication2: Error writing" +
//                            fname);
//      }
//      finally
//      {
//        if(out!=null)
//          out.close();
//      }
//    }
//  }

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

