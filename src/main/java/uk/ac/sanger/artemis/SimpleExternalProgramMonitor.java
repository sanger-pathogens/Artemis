/* SimpleExternalProgramMonitor.java
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/SimpleExternalProgramMonitor.java,v 1.1 2004-06-09 09:45:09 tjc Exp $
 */

package uk.ac.sanger.artemis;

import java.util.Vector;

/**
 *  This class contains methods common to all implementations of
 *  ExternalProgramMonitor
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: SimpleExternalProgramMonitor.java,v 1.1 2004-06-09 09:45:09 tjc Exp $
 **/

public class SimpleExternalProgramMonitor
    extends Thread
    implements ExternalProgramMonitor 
{

  /** Process name that was passed to the constructor. */
  final public String name;

  /** The Logger that was passed to the constructor. */
  final private Logger logger;

  /**
   *  Create a new SimpleExternalProgramMonitor.
   *  @param name The name of the external program that is being monitored.
   *  @param logger The log for errors, STDOUT and STDERR of the external
   *    program.
   **/
  public SimpleExternalProgramMonitor(final String name,
                                      final Logger logger)
  {
    this.name = name;
    this.logger = logger;
  }

  /**
   *  Add a listener for ExternalProgramEvents.
   **/
  public void addExternalProgramListener (ExternalProgramListener l) 
  {
    listeners.add (l);
  }

  /**
   *  Send the given event to all the listeners.
   **/
  protected void sendEvent(ExternalProgramEvent e) 
  {
    for(int i = 0 ; i < listeners.size () ; ++i) 
      ((ExternalProgramListener)listeners.elementAt(i)).statusChanged (e);
  }

  private Vector listeners = new Vector ();

  /**
   *  Return the Logger that was passed to the constructor.
   **/
  public Logger getLogger()
  {
    return logger;
  }

  /**
   *  Return the name of the external program we are monitoring.
   **/
   public String getProgramName()
   {
     return name;
   }

}
