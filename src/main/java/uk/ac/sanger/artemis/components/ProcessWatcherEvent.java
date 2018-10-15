/* ProcessWatcherEvent.java
 *
 * created: Tue Feb 29 2000
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 2000  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/ProcessWatcherEvent.java,v 1.1 2004-06-09 09:47:14 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

/**
 *  This event is sent when the process that is watched by ProcessWatcher
 *  finishes.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: ProcessWatcherEvent.java,v 1.1 2004-06-09 09:47:14 tjc Exp $
 **/

public class ProcessWatcherEvent {
  /**
   *  Create a new ProcessWatcherEvent object.
   *  @param process The Process that has finished.
   *  @param exit_code The exit code of the process that has finished.
   **/
  public ProcessWatcherEvent (final Process process, final int exit_code) {
    this.process = process;
    this.exit_code = exit_code;
  }

  /**
   *  Return the process that was passed to the constructor.
   **/
  public Process getProcess () {
    return process;
  }

  /**
   *  Return the exit_code that was passed to the constructor.
   **/
  public int getExitCode () {
    return exit_code;
  }

  /**
   *  The process that was passed to the constructor.
   **/
  final Process process;

  /**
   *  The exit_code that was passed to the constructor.
   **/
  final int exit_code;
}
