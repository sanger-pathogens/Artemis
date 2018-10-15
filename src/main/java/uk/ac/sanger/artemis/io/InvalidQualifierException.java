/* InvalidQualifierException.java
 *
 * created: Wed Jan  6 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/InvalidQualifierException.java,v 1.1 2004-06-09 09:49:41 tjc Exp $
 **/

package uk.ac.sanger.artemis.io;

/**
 *  This exception is thrown when a String passed as the name to the Qualifier
 *  constructor is not a valid feature qualifier name or the the value passed
 *  to the qualifier is invalid.
 *
 *  @author Kim Rutherford
 *  @version $Id: InvalidQualifierException.java,v 1.1 2004-06-09 09:49:41 tjc Exp $
 **/

public class InvalidQualifierException extends EntryInformationException {
  /**
   *  Create a new InvalidQualifierException object with the given String as
   *  the message.
   *  @param message the detail message
   *  @param qualifier The Qualifier that caused the exception.
   **/
  public InvalidQualifierException (final String message,
                                    final Qualifier qualifier) {
    super (message);

    this.qualifier = qualifier;
  }

  /**
   *  Return the Qualifier that was passed to the constructor.
   **/
  public Qualifier getQualifier () {
    return qualifier;
  }

  /**
   *  The Qualifier that was passed to the constructor.
   **/
  final private Qualifier qualifier;
}


