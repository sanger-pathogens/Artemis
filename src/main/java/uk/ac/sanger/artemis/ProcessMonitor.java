/* ProcessMonitor.java
 *
 * created: Wed Aug  6 2003
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 2003  Genome Research Limited
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/ProcessMonitor.java,v 1.1 2004-06-09 09:45:01 tjc Exp $
 */

package uk.ac.sanger.artemis;

import java.io.*;

/**
 *  Objects of this class watch a Process object and send a
 *  ExternalProgramEvent when the process status changes (eg. the end of a
 *  process).
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: ProcessMonitor.java,v 1.1 2004-06-09 09:45:01 tjc Exp $
 **/

public class ProcessMonitor 
    extends SimpleExternalProgramMonitor
    implements ExternalProgramMonitor 
{

  /** The Process that was passed to the constructor. */
  final public Process process;

  /**
   *  Create a new ProcessMonitor object for the given Process.
   *  @param process The Process to monitor.
   *  @param name The name of the process that is being monitored.
   *  @param logger The log for errors, STDOUT and STDERR of the Process.
   **/
  public ProcessMonitor(final Process process, final String name,
                        final Logger logger) 
  {
    super(name, logger);
    this.process = process;
  }

  /**
   *  This code will wait for a Process to change status then log it and call
   *  sendEvent().
   **/
  public void run() 
  {
    try 
    {
      final Reader reader = new InputStreamReader (process.getErrorStream ());
      getLogger().log (reader);
    }
    catch(IOException e) 
    {
      getLogger().log("cannot read the errer stream from " +
                       getProgramName ());
    }

    try 
    {
      final Reader reader = new InputStreamReader (process.getInputStream ());
      getLogger().log (reader);
    } 
    catch (IOException e) 
    {
      getLogger().log("cannot read the ouput stream from " +
                       getProgramName ());
    }
    
    while(true) 
    {
      try
      {
        final int return_value = process.waitFor();

        final boolean core_dumped = (return_value & 0x80) != 0;

        getLogger ().log ("\n--------------------------------------" +
                          "---------------------\n\n");

        final String log_message;

        if(core_dumped) 
          log_message = getProgramName() + " process dumped core";
        else 
        {
          final int sig_number = return_value & 0x7f;

          if (sig_number > 0) 
            log_message =
              getProgramName() + " process received signal: " + sig_number;
          else 
          {
            final int exit_code = return_value >> 8;

            if(exit_code == 0) 
              log_message = getProgramName() + " process completed";
            else 
              log_message = getProgramName() +
                " process finished with exit code: " + exit_code;
          }
        }

        final ExternalProgramEvent event =
          new ExternalProgramEvent(ExternalProgramEvent.FINISHED,
                                    log_message, process);

        sendEvent(event);
        getLogger().log (log_message + "\n");

        return;
      } 
      catch (InterruptedException e) 
      {
        // go around the loop again
      }
    }
  }

}
