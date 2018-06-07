/* ProcessWatcher.java
 *
 * created: Mon Oct  4 1999
 *
 * This file is part of Artemis
 *
 * Copyright (C) 1999  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/ProcessWatcher.java,v 1.1 2004-06-09 09:47:13 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.Logger;

import java.io.*;

/**
 *  Objects of this class watch a Process object and then display a
 *  MessageFrame window when the process finishes.
 *
 *  @author Kim Rutherford
 *  @version $Id: ProcessWatcher.java,v 1.1 2004-06-09 09:47:13 tjc Exp $
 **/

public class ProcessWatcher
    implements Runnable  {
  /**
   *  Create a new ProcessWatcher object for the given Process.  When the
   *  process finishes a MessageFrame will alert the user and all
   *  ProcessWatcherListeners will be informed.
   *  @param process The Process to watch.
   *  @param name The name of the process to watch.
   **/
  public ProcessWatcher (final Process process, final String name) {
    this (process, name, true);
  }

  /**
   *  Create a new ProcessWatcher object for the given Process.  When the
   *  process finishes a MessageFrame will alert the user if and only if the
   *  alert_user argument is true and all ProcessWatcherListeners will be
   *  informed.
   *  @param process The Process to watch.
   *  @param name The name of the process to watch.
   *  @param alert_user The user will be informed when a process ends if and
   *    only if this is true.
   **/
  public ProcessWatcher (final Process process, final String name,
                         final boolean alert_user) {
    this.process = process;
    this.name = name;
    this.alert_user = alert_user;
  }

  /**
   *  This code will wait for the Process to finish then display the return
   *  value in a MessageFrame.
   **/
  public void run () {
    try {
      final Reader reader = new InputStreamReader (process.getErrorStream ());
      getLogger ().log (reader);
    } catch (IOException e) {
      final MessageFrame message_frame =
        new MessageFrame ("error while reading output of: " + name);

      message_frame.setVisible (true);
    }

    try {
      final Reader reader = new InputStreamReader (process.getInputStream ());
      getLogger ().log (reader);
    } catch (IOException e) {
      final MessageFrame message_frame =
        new MessageFrame ("error while reading output of: " + name);

      message_frame.setVisible (true);
    }

    while (true) {
      try {
        final int return_value = process.waitFor ();

        for (int i = 0 ; i < listeners.size () ; ++i) {
          final ProcessWatcherEvent event =
            new ProcessWatcherEvent (process, return_value);
          final ProcessWatcherListener listener =
            (ProcessWatcherListener) listeners.elementAt (i);
          listener.processFinished (event);
        }

        final boolean core_dumped = (return_value & 0x80) != 0;

        getLogger ().log ("\n--------------------------------------" +
                          "---------------------\n\n");

        final String log_message;

        if (core_dumped) {
          log_message = name + " process dumped core";
          new MessageFrame (log_message +
                            " - check the log window").setVisible (true);
        } else {
          final int sig_number = return_value & 0x7f;

          if (sig_number > 0) {
            log_message = name + " process received signal: " + sig_number;
            new MessageFrame (log_message +
                              " - check the log window").setVisible (true);
          } else {
            final int exit_code = return_value >> 8;

            if (exit_code == 0) {
              log_message = name + " process completed";
              if (alert_user) {
                final MessageFrame message_frame =
                  new MessageFrame (log_message);
                message_frame.setVisible (true);
              }
            } else {
              log_message =
                name + " process finished with exit code: " + exit_code;
              new MessageFrame (log_message +
                                " - check the log window").setVisible (true);
            }
          }
        }

        getLogger ().log (log_message + "\n");

        return;
      } catch (InterruptedException e) {
       // go around the loop again
      }
    }
  }

  /**
   *  Return the global Logger object.
   **/
  private static Logger getLogger () {
    return Splash.getLogger ();
  }

  /**
   *  Add the given object as a ProcessWatcherListener for this ProcessWatcher.
   **/
  public void addProcessWatcherListener (final ProcessWatcherListener l) {
    listeners.addElement (l);
  }

  /**
   *  Remove the given object as a ProcessWatcherListener for this
   *  ProcessWatcher.
   **/
  public void removeProcessWatcherListener (final ProcessWatcherListener l) {
    listeners.removeElement (l);
  }

  /**
   *  The Process reference that was passed to the constructor.
   **/
  private Process process;

  /**
   *  The program name that was passed to the constructor.
   **/
  private String name;

  /**
   *  The alert_user argument that was passed to the constructor.
   **/
  private boolean alert_user;

  /**
   *  The list of objects that are listening for ProcessWatcher events.
   **/
  private final java.util.Vector listeners = new java.util.Vector ();
}
