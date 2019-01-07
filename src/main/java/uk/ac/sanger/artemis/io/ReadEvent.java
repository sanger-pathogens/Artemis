/* EntryStreamEvent.java
 *
 * created: Fri Aug 15 2003
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/ReadEvent.java,v 1.1 2004-06-09 09:50:21 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

/**
 *  An Event that holds a message, a warning or an error generated while
 *  reading an Entry.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: ReadEvent.java,v 1.1 2004-06-09 09:50:21 tjc Exp $
 **/

public class ReadEvent extends java.util.EventObject {
  /**
   *  Create a new ReadEvent with the given source and message.
   **/
  public ReadEvent (final Object source, final String message) {
    super (source);
    this.message = message;
  }

  /**
   *  Return the message that was passed to the constructor.
   **/
  public String getMessage () {
    return message;
  }

  /**
   *  The message that was passed to the constructor.
   **/
  final String message;
}
