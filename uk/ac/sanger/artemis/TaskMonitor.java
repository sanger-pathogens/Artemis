/* TaskMonitor.java
 *
 * created: Thu Aug  7 2003
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/TaskMonitor.java,v 1.1 2004-06-09 09:45:11 tjc Exp $
 */

package uk.ac.sanger.artemis;

import uk.ac.sanger.jcon.job.*;

/**
 *  Objects of this class watch a Task object and then display a
 *  MessageFrame window when the process status changes (eg. the end of a
 *  process) and send a ExternalProgramEvent.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: TaskMonitor.java,v 1.1 2004-06-09 09:45:11 tjc Exp $
 **/

public class TaskMonitor
    extends SimpleExternalProgramMonitor
    implements ExternalProgramMonitor {
  /**
   *  Create a new ProcessMonitor object for the given Process.
   *  @param task The Task to monitor.
   *  @param name The name of the Task that is being monitored (eg. blastn).
   *  @param logger The log for errors, STDOUT and STDERR of the task.
   **/
  public TaskMonitor (final Task task, final String name,
                      final Logger logger) {
    super (name, logger);
    this.task = task;
    SLEEP_TIME =
      1000 *
      Options.getOptions ().getIntegerProperty ("jcon_task_check_time").intValue ();
  }

  /**
   *  The length of time (in seconds) to sleep for.
   **/
  private final int SLEEP_TIME;

  /**
   *  Wait for the status of a Task to change then log it and call sendEvent().
   **/
  public void run () {
    while (true) {
      try {
        final Status status = task.getStatus ();
        String message = null;

        if (status == null) {
          message = getProgramName () + " disappeared";
        } else {
          final int status_id = status.getId ();
          
          switch (status_id) {
          case Status.COMPLETED:
            message = getProgramName () + " completed";
            break;
          case Status.FAILED:
            message = getProgramName () + " failed";
            break;
          case Status.CANCELLED:
            message = getProgramName () + " cancelled";
            break;
          case Status.SKIPPED:
            message = getProgramName () + " skipped";
            break;
          default:
            // keep going
          }
        }

        if (message != null) {
          final ExternalProgramEvent event =
            new ExternalProgramEvent (ExternalProgramEvent.FINISHED,
                                      message, task);

          sendEvent (event);
          getLogger ().log (message + "\n");

          return;
        }

        Thread.sleep (SLEEP_TIME);
      } catch (IllegalStateException e) {
        return;
      } catch (InterruptedException _) {
        // ignore interrupt during a sleep()
      }
    }
  }

  /**
   *  The Task that was passed to the constructor.
   **/
  final public Task task;
}
