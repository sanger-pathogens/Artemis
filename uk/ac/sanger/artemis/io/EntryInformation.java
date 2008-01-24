/* EntryInformation.java
 *
 * created: Tue Jan  5 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/EntryInformation.java,v 1.3 2008-01-24 14:15:30 tjc Exp $
 **/

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.StringVector;


/**
 *  This interface contains methods that give information about the possible
 *  keys and qualifiers of an entry.
 *
 *  @author Kim Rutherford
 *  @version $Id: EntryInformation.java,v 1.3 2008-01-24 14:15:30 tjc Exp $
 **/

public interface EntryInformation {
  /**
   *  Add a new feature key to the list of keys returned by getValidKeys ()
   *  and getSortedValidKeys ().
   **/
  void addKey (final Key key);

  /**
   *  Add a QualifierInfo object to this EntryInformation object.  This will
   *  change the return values of getSortedValidKeys() and
   *  getValidKeys() if the new QualifierInfo object refers to a new key.
   *  @exception QualifierInfoException Thrown if a QualifierInfo object with
   *    the same name, but conflicting types has already been added to this
   *    EntryInformation object.
   **/
  void addQualifierInfo (final QualifierInfo qualifier_info)
      throws QualifierInfoException;

  /**
   *  Return a vector containing the valid feature keys.  The returned
   *  Vector is a copy.
   *  @return null if and only if all keys are valid.
   **/
  KeyVector getValidKeys ();

  /**
   *  Return a alphanumerically sorted vector of the valid keys.
   *  @return null if and only if all keys are valid.
   **/
  KeyVector getSortedValidKeys ();

  /**
   *  Returns the keys that where added with addKey ()
   **/
  KeyVector getUserKeys ();

  /**
   *  Return the Key that should be used when a new (empty) feature is created.
   **/
  Key getDefaultKey ();

  /**
   *  Return a vector containing all valid qualifiers for a feature with
   *  the given key.  The returned StringVector is a copy.
   **/
  StringVector getValidQualifierNames (final Key key);

 /**
   *  Return a vector containing the required feature qualifiers for the given
   *  key.
   **/
  StringVector getRequiredQualifiers (final Key key);

  /**
   *  Return true if and only if the given String contains a valid
   *  qualifier name.
   **/
  boolean isValidQualifier (final String name);

  /**
   *  Return true if and only if given qualifier_name is a legal qualifier
   *  name for the given key.
   **/
  boolean isValidQualifier (final Key key,
                                   final String qualifier_name);

  /**
   *  Return true if and only if given qualifier_name is a legal qualifier
   *  name for the given key and must be present in a feature with that key.
   **/
  boolean isRequiredQualifier (final Key key,
                                      final String qualifier_name);

  /**
   *  Return true if and only if given Key is a valid.
   **/
  boolean isValidKey (final Key key);

  /**
   *  Return the QualifierInfo object for the given qualifier name, or null if
   *  there are no legal qualifiers with the given name.
   **/
  QualifierInfo getQualifierInfo (final String qualifier_name);

  /**
   *  Return a Vector contains all the QualifierInfo objects that this
   *  EntryInformation knows about.
   **/
  QualifierInfoHash getAllQualifierInfo ();

  /**
   *  Returns true if strict EMBL format will be used when writing.  If true
   *  qualifiers will always wrap at 80 columns.
   *  EMBL format locations have the complement (if any) inside the join.  eg:
   *  join(complement(1..100),complement(200..300)).
   *  The default format has the complement on the outside, which is much more
   *  readable.  eg: complement(join(1..100,200..300))
   **/
  boolean useEMBLFormat ();

  /**
   *  Set the use_embl_format flag (see useEMBLFormat ()).  true means the
   *  strict EMBL format should be used when writing.
   **/
  void setEMBLFormat (final boolean use_embl_format);

  /**
   *  Fix this EntryInformation so that the given exception won't happen
   *  again.
   **/
  void fixException (final EntryInformationException exception);
}
