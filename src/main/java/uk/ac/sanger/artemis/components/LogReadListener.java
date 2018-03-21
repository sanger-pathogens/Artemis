/* LogReadListener.java
 *
 * created: Fri Nov 28 2003
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/LogReadListener.java,v 1.2 2004-12-03 17:47:04 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.io.ReadListener;
import uk.ac.sanger.artemis.io.ReadEvent;

import javax.swing.*;

/**
 *  A class that implements ReadListener by logging all ReadEvents.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: LogReadListener.java,v 1.2 2004-12-03 17:47:04 tjc Exp $
 **/

public class LogReadListener implements ReadListener {
  /**
   *  Create a new DialogReadListener.
   *  @param source The source we are reading from - generally a file name or
   *    URL.  This file name is prepended to the log entries.
   **/
  public LogReadListener (final String source) {
    this.source = source;
  }
  
  /**
   *  Implementation of ReadListener.
   **/
  public void notify (final ReadEvent event) {

    if (source == null) {
      Splash.getLogger ().log (event.getMessage () + "\n");
    } else {
      Splash.getLogger ().log ("while reading from " + source + ": " +
                               event.getMessage () + "\n");
    }
    seen_message = true;
  }
  
  /**
   *  Return true if and only if notify() has been called at least once.
   **/
  public boolean seenMessage () {
    return seen_message;
  }
  
  /**
   *  The name that was passed to the constructor.
   **/
  private String source;

  private boolean seen_message = false;
}
