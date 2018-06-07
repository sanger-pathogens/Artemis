/* ReadOnlyEmblStreamFeature.java
 *
 * created: Mon Aug 13 2001
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2001  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/ReadOnlyEmblStreamFeature.java,v 1.1 2004-06-09 09:50:24 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.*;

/**
 *  A read-only version of EmblStreamFeature - all set() methods throw
 *  ReadOnlyException.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: ReadOnlyEmblStreamFeature.java,v 1.1 2004-06-09 09:50:24 tjc Exp $
 **/

public class ReadOnlyEmblStreamFeature extends EmblStreamFeature {
  /**
   *  Create a new EmblStreamFeature object.
   *  @param key The new feature key
   *  @param location The Location object for the new feature
   *  @param qualifiers The qualifiers for the new feature
   *  @exception InvalidRelationException Thrown if this Feature cannot contain
   *    the given Qualifier.
   **/
  public ReadOnlyEmblStreamFeature (Key key,
                                    Location location,
                                    QualifierVector qualifiers)
      throws InvalidRelationException {
    super (key, location, qualifiers);

    finished_constructor = true;
  }

  /**
   *  Set the owning DocumentEntry of this Feature.
   *  @param entry The Entry that now owns this Feature.
   **/
  public void setDocumentEntry (final DocumentEntry entry)
      throws ReadOnlyException {
    if (finished_constructor) {
      throw new ReadOnlyException ();
    } else {
      super.setDocumentEntry (entry);
    }
  }

  /**
   *  Set the value of this object.
   *  @param key The new feature key
   *  @param location The Location object for the new feature
   *  @param qualifiers The qualifiers for the new feature
   *  @exception EntryInformationException Thrown if this Feature cannot
   *    contain the given Key, Qualifier or Key/Qualifier combination.
   **/
  public void set (final Key key,
                   final Location location,
                   final QualifierVector qualifiers)
      throws OutOfRangeException, EntryInformationException,
             ReadOnlyException {
    if (finished_constructor) {
      throw new ReadOnlyException ();
    } else {
      super.set (key, location, qualifiers);
    }
  }

  /**
   *  Set the key field of this object.
   *  @param key The new feature key
   *  @exception EntryInformationException Thrown if this Feature cannot
   *    contain the given Key, Qualifier or Key/Qualifier combination.
   **/
  public void setKey (final Key key)
      throws EntryInformationException, ReadOnlyException {
    if (finished_constructor) {
      throw new ReadOnlyException ();
    } else {
      super.setKey (key);
    }
  }

  /**
   *  Set the location of this object.
   *  @param location The Location object for the new feature
   **/
  public void setLocation (final Location location)
      throws ReadOnlyException, OutOfRangeException {
    if (finished_constructor) {
      throw new ReadOnlyException ();
    } else {
      super.setLocation (location);
    }
  }

  /**
   *  Set the qualifiers of this object.
   *  @param qualifiers The new qualifiers for the feature.  null means remove
   *    all qualifiers.
   *  @exception EntryInformationException Thrown if this Feature cannot
   *    contain the given Key, Qualifier or Key/Qualifier combination.
   **/
  public void setQualifiers (final QualifierVector qualifiers)
      throws EntryInformationException, ReadOnlyException {
    if (finished_constructor) {
      throw new ReadOnlyException ();
    } else {
      super.setQualifiers (qualifiers);
    }
  }

  /**
   *  Add the given Qualifier to this Feature.  If this Feature contains a
   *  Qualifier with the same name as the new Qualifier it will be replaced.
   *  @param qualifier The new qualifier to add.
   *  @exception EntryInformationException Thrown if this Feature cannot
   *    contain the given Key, Qualifier or Key/Qualifier combination.
   **/
  public void setQualifier (final Qualifier qualifier)
      throws EntryInformationException, ReadOnlyException {
    if (finished_constructor) {
      throw new ReadOnlyException ();
    } else {
      super.setQualifier (qualifier);
    }
  }

  /**
   *  Remove the Qualifier with the given name.  If there is no Qualifier with
   *  that name then return immediately.
   *  @exception EntryInformationException Thrown if a required qualifier is
   *    removed.
   *  @param name The qualifier name to look for.
   **/
  public void removeQualiferByName (final String name)
      throws EntryInformationException, ReadOnlyException {
    if (finished_constructor) {
      throw new ReadOnlyException ();
    } else {
      super.removeQualifierByName (name);
    }
  }

  /**
   *  Add the values from the given qualifier to the Qualifier object with the
   *  same name in this Feature or otherwise add a copy of the argument.
   *  @param qualifier This object contains name and values to add.
   *  @return The Qualifier that was changed or created.
   *  @exception EntryInformationException Thrown if this Feature cannot
   *    contain the given Key, Qualifier or Key/Qualifier combination.
   **/
  public Qualifier addQualifierValues (final Qualifier qualifier)
      throws EntryInformationException, ReadOnlyException {
    if (finished_constructor) {
      throw new ReadOnlyException ();
    } else {
      return super.addQualifierValues (qualifier);
    }
  }

  private boolean finished_constructor = false;
}
