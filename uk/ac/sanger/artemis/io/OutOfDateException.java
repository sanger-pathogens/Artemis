/* OutOfDateException.java
 *
 * created: Mon May 24 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/OutOfDateException.java,v 1.1 2004-06-09 09:50:02 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

/**
 *  This exception is thrown when an attempt is made to change an entry or
 *  feature that has changed since the last call to Entry.getDateStamp ().
 *
 *  @author Kim Rutherford
 *  @version $Id: OutOfDateException.java,v 1.1 2004-06-09 09:50:02 tjc Exp $
 **/

public class OutOfDateException extends Exception {
  /**
   *  This constructor creates a OutOfDateException with no detail message.
   **/
  public OutOfDateException () {
    super ();
  }
}

