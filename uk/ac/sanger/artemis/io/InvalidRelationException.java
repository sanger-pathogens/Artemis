/* InvalidRelationException.java
 *
 * created: Fri Jan  1 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/InvalidRelationException.java,v 1.1 2004-06-09 09:49:42 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

/**
 *  This exception is thrown when an inappropriate Qualifier is added to or
 *  accessed in a Feature.  For example adding a /codon_start qualifier to a
 *  misc_feature Feature.
 *
 *  @author Kim Rutherford
 *  @version $Id: InvalidRelationException.java,v 1.1 2004-06-09 09:49:42 tjc Exp $
 **/

public class InvalidRelationException extends EntryInformationException {
  /**
   *  Create a new InvalidRelationException object with the given String as
   *  the message.
   *  @param message the detail message
   *  @param key The key that caused the exception.
   *  @param qualifier The Qualifier that caused the exception.
   **/
  public InvalidRelationException (String message, Key key,
                                   Qualifier qualifier) {
    super (message);
    this.key = key;
    this.qualifier = qualifier;
  }

  /**
   *  Create a new InvalidRelationException object with the given String as
   *  the message.  This constructor should be used only when the qualifier
   *  that caused the problem is not known.
   *  @param message the detail message
   *  @param key The key that caused the exception.
   **/
  public InvalidRelationException (String message, Key key) {
    super (message);
    this.key = key;
    this.qualifier = null;
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

  /**
   *  Return the Qualifier that was passed to the constructor (which could be
   *  null if the second qualifier was called).
   **/
  public Qualifier getQualifier () {
    return qualifier;
  }

  /**
   *  The Qualifier that was passed to the constructor.
   **/
  final private Qualifier qualifier;
}
