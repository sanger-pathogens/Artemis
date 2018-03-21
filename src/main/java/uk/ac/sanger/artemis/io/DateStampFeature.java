/* DateStampFeature.java
 *
 * created: Mon Jul 12 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/DateStampFeature.java,v 1.1 2004-06-09 09:49:02 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.*;

import java.util.Date;

/**
 *  This interface is a Feature with set () methods that take a datestamp as
 *  an argument.
 *
 *  @author Kim Rutherford
 *  @version $Id: DateStampFeature.java,v 1.1 2004-06-09 09:49:02 tjc Exp $
 **/

public interface DateStampFeature extends Feature {

  /**
   *  Set the value of this object.
   *  @param key The new feature key
   *  @param location The Location object for the new feature
   *  @param qualifiers The qualifiers for the new feature
   *  @exception InvalidRelationException Thrown if this Feature cannot contain
   *    any of the given Qualifier objects.
   *  @exception OutOfDate If the key has changed in the server since the time
   *    given by datestamp.
   *  @exception ReadOnlyException Thrown if this Feature cannot be changed.
   *  @exception OutOfRangeException Thrown if any part of the location is out
   *    of range for the sequence.
   **/
  void set (final Date datestamp,
            final Key key,
            final Location location,
            final QualifierVector qualifiers)
      throws InvalidRelationException, OutOfRangeException, ReadOnlyException,
             OutOfDateException;
  /**
   *  Set the key field of this object.
   *  @param key The new feature key
   *  @exception ReadOnlyException Thrown if this Feature cannot be changed.
   *  @exception OutOfDate If the key has changed in the server since the time
   *    given by datestamp.
   **/
  void setKey (final Date datestamp, final Key key)
      throws InvalidRelationException, ReadOnlyException,
             OutOfDateException;

    /**
   *  Set the location of this object.
   *  @param location The Location object for the new feature
   *  @exception ReadOnlyException Thrown if this Feature cannot be changed.
   *  @exception OutOfDate If the key has changed in the server since the time
   *    given by datestamp.
   *  @exception OutOfRangeException Thrown if any part of the location is out
   *    of range for this sequence.
   **/
  void setLocation (final Date datestamp,
                    final Location location)
      throws ReadOnlyException, OutOfDateException, OutOfRangeException;

  /**
   *  Set the qualifiers of this object discarding all the current qualifiers.
   *  @param qualifiers The qualifiers for the new feature
   *  @exception InvalidRelationException Throw if this Feature cannot contain
   *    any of the given Qualifier objects.
   *  @exception ReadOnlyException Thrown if this Feature cannot be changed.
   *  @exception OutOfDateException If the key has changed in the server since
   *    the time given by datestamp.
   **/
  void setQualifiers (final Date datestamp,
                      final QualifierVector qualifiers)
      throws InvalidRelationException, ReadOnlyException, OutOfDateException;

  /**
   *  Return the datestamp of this feature - that is, the time it was last
   *  changed.
   **/
  Date getDatestamp ();
}

