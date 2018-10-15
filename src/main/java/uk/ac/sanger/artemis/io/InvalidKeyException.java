/* InvalidKeyException.java
 *
 * created: Sun Jan  3 1999
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 1998,1999,2000  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/InvalidKeyException.java,v 1.1 2004-06-09 09:49:40 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

/**
 *  This exception is thrown when a String passed to the Key constructor is
 *  not a valid feature key.
 *
 *  @author Kim Rutherford
 *  @version $Id: InvalidKeyException.java,v 1.1 2004-06-09 09:49:40 tjc Exp $
 **/

public class InvalidKeyException extends EntryInformationException {
  /**
   *  Create a new InvalidKeyException object with the given String as
   *  the message.
   *  @param message the detail message
   **/
  public InvalidKeyException (String message, Key key) {
    super (message);
    this.key = key;
  }

  /**
   *  Return the Key that was passed to the constructor.
   **/
  public Key getKey () {
    return key;
  }

  /**
   *  The Key that was passed to the constructor.
   **/
  final private Key key;
}
