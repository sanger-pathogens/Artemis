/* ExternalProgramEvent.java
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/ExternalProgramEvent.java,v 1.2 2007-11-09 10:36:30 tjc Exp $
 */

package uk.ac.sanger.artemis;

//import uk.ac.sanger.jcon.job.*;

/**
 *  An event that represents a change in the state of an ExternalProgram.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: ExternalProgramEvent.java,v 1.2 2007-11-09 10:36:30 tjc Exp $
 **/

public class ExternalProgramEvent {
  /**
   *  Create a new ExternalProgramEvent object of the given type.
   *  @param message A summary of the event.
   *  @param task The Task that caused this event.
   **/
  /*
  public ExternalProgramEvent (final int type, final String message,
                               final Task task) {
    this.type = type;
    this.message = message;
    this.task = task;
  }
  */
  
  /**
   *  Create a new ExternalProgramEvent object of the given type.
   *  @param message A summary of the event.
   *  @param process The Process that caused this event.
   **/
  public ExternalProgramEvent (final int type, final String message,
                               final Process process) {
    this.type = type;
    this.message = message;
    this.process = process;
  }

  /**
   *  The type of event sent when the ExternalProgram finishes.
   **/
  final public static int FINISHED = 1;

  /**
   *  The type of event sent when the ExternalProgram starts.
   **/
  final public static int STARTED = 1;

  /**
   *  Return the type of this event.
   **/
  public int getType () {
    return type;
  }

  /**
   *  Return the message that was passed to the constructor.
   **/
  public String getMessage () {
    return message;
  }

  /**
   *  Return the Task that was passed to the constructor or null if a Process
   *  was passed to the constructor.
   **/
  /*public Task getTask () {
    return task;
  }*/

  /**
   *  Return the Process that was passed to the constructor or null if a
   *  Task was passed to the constructor.
   **/
  public Process getProcess () {
    return process;
  }

  /**
   *  The type that was passed to the constructor.
   **/
  final private int type;

  /**
   *  The message that was passed to the constructor.
   **/
  final private String message;

  /**
   *  The Task that was passed to the constructor or null if a Process was
   *  passed to the constructor.
   **/
  //private Task task = null;

  /**
   *  The Process that was passed to the constructor or null if a Task was
   *  passed to the constructor.
   **/
  private Process process = null;
}
