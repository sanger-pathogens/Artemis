/* Feature.java
 *
 * created: Mon Oct 12 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/Feature.java,v 1.2 2005-11-28 16:46:38 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.*;

/**
 *  Feature interface
 *
 *  @author Kim Rutherford
 *  @version $Id: Feature.java,v 1.2 2005-11-28 16:46:38 tjc Exp $
 *
 */

public interface Feature {
  /**
   *  Set the value of this object.
   *  @param key The new feature key
   *  @param location The Location object for the new feature
   *  @param qualifiers The qualifiers for the new feature
   *  @exception EntryInformationException Thrown if this Feature cannot contain
   *    the given Key, Qualifier or Key/Qualifier combination.
   *  @exception InvalidRelationException Thrown if this Feature cannot contain
   *    any of the given Qualifier objects.
   **/
  void set (final Key key,
                   final Location location,
                   final QualifierVector qualifiers)
      throws EntryInformationException, ReadOnlyException, OutOfRangeException;

  /**
   *  Set the key field of this object.
   *  @param key The new feature key
   *  @exception InvalidKeyException Thrown if the given Key is not a valid
   *    Key for this type of feature.
   *  @exception EntryInformationException Thrown if this Feature cannot contain
   *    the given Key, Qualifier or Key/Qualifier combination.
   *  @exception InvalidRelationException Thrown if this Feature cannot contain
   *    any of the given Qualifier objects.
   **/
  void setKey (final Key key)
      throws EntryInformationException, ReadOnlyException;

  /**
   *  Set the location of this object.
   *  @param location The Location object for the new feature
   **/
  void setLocation (final Location location)
      throws OutOfRangeException, ReadOnlyException;

  void setLocation (final Location location, Entry saved_entry)
      throws OutOfRangeException, ReadOnlyException;

  /**
   *  Set the qualifiers of this object discarding all the current qualifiers.
   *  @param qualifiers The qualifiers for the new feature
   *  @exception EntryInformationException Thrown if this Feature cannot contain
   *    the given Key, Qualifier or Key/Qualifier combination.
   *  @exception InvalidRelationException Throw if this Feature cannot contain
   *    any of the given Qualifier objects.
   **/
  void setQualifiers (final QualifierVector qualifiers)
      throws EntryInformationException, ReadOnlyException;

  /**
   *  Add the given Qualifier to this Feature.  If this Feature contains a
   *  Qualifier with the same name as the new Qualifier it will be replaced.
   *  @param qualifier The new qualifier to add.
   *  @exception EntryInformationException Thrown if this Feature cannot contain
   *    the given Key, Qualifier or Key/Qualifier combination.
   **/
  void setQualifier (final Qualifier qualifier)
      throws EntryInformationException, ReadOnlyException;

  /**
   *  Remove the Qualifier with the given name.  If there is no Qualifier with
   *  that name then return immediately.
   *  @exception EntryInformationException Thrown if this Feature cannot contain
   *    the given Key, Qualifier or Key/Qualifier combination.
   *  @param name The qualifier name to look for.
   **/
  void removeQualifierByName (final String name)
      throws EntryInformationException, ReadOnlyException;

  /**
   *  Return the key of this feature, as passed to the constructor.
   *  @exception EntryInformationException Thrown if this Feature cannot
   *    contain the given Key, Qualifier or Key/Qualifier combination.
   **/
  Key getKey ();

  /**
   *  Return the Location of this Feature, as passed to the constructor.
   **/
  Location getLocation ();

  /**
   *  Return a QualifierVector object containing the qualifiers for this
   *  feature.  This method does not return a copy of the qualifier vector so
   *  changing the vector will change the feature.
   **/
  QualifierVector getQualifiers ();

  /**
   *  Return a copy of the Qualifier in this Feature with the given name or
   *  null if there no such Qualifier.
   **/
  Qualifier getQualifierByName (final String name)
      throws InvalidRelationException;

  /**
   *  Return the first base of this feature.
   **/
  int getFirstBase ();

  /**
   *  Return the last base of this feature.
   **/
  int getLastBase ();

  /**
   *  Return the Entry object that contains this Feature.
   **/
  Entry getEntry ();

  /**
   *  Return the reference of a new copy of this Feature.  The new Feature
   *  will not be in any Entry, so it should be explicitly added with a call
   *  to Entry.add ().
   **/
  Feature copy ();

  /**
   *  Set the user data reference for this object.  The user data object is
   *  chosen by the user of this class.  This facility is provided to make it
   *  easier to wrap this class in more powerful classes.
   *  @param user_data The new value of the user data reference of this object.
   *    The existing user data reference (if any) is lost.
   **/
  void setUserData (final Object user_data);

  /**
   *  Return the user data reference for this object, as set by a call to
   *  setUserData ().
   **/
  Object getUserData ();

  /**
   *  Returns true if and only if this Feature can't be changed or can't be
   *  removed from it's entry.
   **/
  boolean isReadOnly ();
}
