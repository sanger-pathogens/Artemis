/* EntryInformationException.java
 *
 * created: Fri Feb 11 2000
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/EntryInformationException.java,v 1.1 2004-06-09 09:49:17 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

/**
 *  This exception is thrown when a Key, Qualifier or Key/Qualifier
 *  combination is used in a way that is not allowed in the current
 *  EntryInformation environment.  This could happen if a required qualifier
 *  is removed.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: EntryInformationException.java,v 1.1 2004-06-09 09:49:17 tjc Exp $
 **/

public class EntryInformationException extends Exception {
  /**
   *  Create a new EntryInformationException object with the given String as
   *  the message.
   *  @param message the detail message
   **/
  public EntryInformationException (final String message) {
    super (message);
  }

}
