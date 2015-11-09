/* SimpleEntryInformation.java
 *
 * created: Wed Feb  9 2000
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/SimpleEntryInformation.java,v 1.5 2008-06-26 09:38:54 tjc Exp $
 */


package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.*;


/**
 *  A SimpleEntryInformation is an EntryInformation object allows any key or
 *  qualifier.
 *
 *  @author Kim Rutherford
 *  @version $Id: SimpleEntryInformation.java,v 1.5 2008-06-26 09:38:54 tjc Exp $
 **/

public class SimpleEntryInformation
  implements EntryInformation {

  /**
   *  Create a new SimpleEntryInformation object.
   **/
  public SimpleEntryInformation () {

  }

  /**
   *  Create a new SimpleEntryInformation object which is a copy of the given
   *  object.
   **/
  public SimpleEntryInformation (final EntryInformation new_info) {
    if (new_info.getAllQualifierInfo () != null) {
      qualifier_info_hash = new_info.getAllQualifierInfo ();
    }
    if (new_info.getValidKeys () != null) {
      valid_keys = new_info.getValidKeys ();
    }
    if (new_info.getUserKeys () != null) {
      user_keys = new_info.getUserKeys ().copy ();
    }

    use_embl_format = new_info.useEMBLFormat ();
  }

  /**
   *  Add a new feature key to the list of keys returned by getValidKeys()
   *  and getSortedValidKeys().  Any key added with this method can have any
   *  qualifier.
   **/
  public void addKey (final Key key) {
    if (user_keys != null && user_keys.contains (key)) {
      return;
    }

    if (user_keys == null) {
      user_keys = new KeyVector ();
    }

    user_keys.add (key);
  }

  /**
   *  Add a QualifierInfo object to this EntryInformation object.  This will
   *  change the return values of getSortedValidKeys() and
   *  getValidKeys() if the new QualifierInfo object refers to a new key.
   *  @exception QualifierInfoException Thrown if a QualifierInfo object with
   *    the same name, but conflicting types has already been added to this
   *    EntryInformation object.
   **/
  public void addQualifierInfo (final QualifierInfo qualifier_info)
      throws QualifierInfoException {

    if (qualifier_info_hash == null) {
      qualifier_info_hash = new QualifierInfoHash ();
    }

    final QualifierInfo current_qualifier_info =
      getQualifierInfo (qualifier_info.getName ());

    if (current_qualifier_info != null) {
      if (qualifier_info.getType () != qualifier_info.UNKNOWN &&
          current_qualifier_info.getType () != qualifier_info.getType ()) {
        final String message =
          "qualifier " + qualifier_info.getName () + " used with " +
          "conflicting types";
        throw new QualifierInfoException (message);
      }

      if (qualifier_info.getValidKeys () == null ||
          qualifier_info.getRequiredKeys () == null) {

        qualifier_info_hash.put (qualifier_info);
      }
    } else {
      qualifier_info_hash.put (qualifier_info);
    }

    final KeyVector qualifier_valid_keys = qualifier_info.getValidKeys ();

    if (qualifier_valid_keys != null) {
      for (int i = 0 ; i < qualifier_valid_keys.size () ; ++i) {
        final Key this_key = (Key)qualifier_valid_keys.get(i);

        if (valid_keys != null && valid_keys.contains (this_key)) {
          continue;
        }

        if (valid_keys == null) {
          valid_keys = new KeyVector ();
        }

        valid_keys.add (this_key);
      }
    }
  }

  /**
   *  Return a vector containing the valid feature keys.  The returned
   *  Vector is a copy.
   *  @return null if and only if all keys are valid.
   **/
  public KeyVector getValidKeys () {
    if (valid_keys == null) {
      if (user_keys == null) {
        return null;
      } else {
        return user_keys.copy ();
      }
    } else {
      final KeyVector return_keys = valid_keys.copy ();

      if (user_keys != null) {
        for (int i = 0 ; i < user_keys.size () ; ++i) {
          if (!return_keys.contains ((Key)user_keys.get (i))) {
            return_keys.add ((Key)user_keys.get (i));
          }
        }
      }

      return return_keys;
    }
  }

  /**
   *  Return a alphanumerically sorted vector of the valid keys.
   *  @return null if and only if all keys are valid.
   **/
  public KeyVector getSortedValidKeys () {
    final KeyVector return_vector = getValidKeys ();

    if (return_vector == null) {
      return null;
    }

    return_vector.mysort ();
    return return_vector;
  }

  /**
   *  Return the Key that should be used when a new (empty) feature is created.
   **/
  public Key getDefaultKey () {
    Key misc_feature_key = new Key ("misc_feature");

    if (isValidKey (misc_feature_key))
      return misc_feature_key;

    misc_feature_key = new Key ("region");
    if (isValidKey (misc_feature_key))
      return misc_feature_key;

    return (Key)getValidKeys ().get (0);
  }

  /**
   *  Return a vector containing all valid feature qualifier names for
   *  features with the given key.  The returned StringVector is a copy.
   *  @return null if and only if any Qualifier is legal.
   **/
  public StringVector getValidQualifierNames (final Key key) {
    if (getQualifierInfoHash () == null) {
      return null;
    }

    final StringVector all_names = getQualifierInfoHash ().names ();

    final StringVector return_names = new StringVector ();

    for (int i = 0 ; i < all_names.size () ; ++i) {
      final QualifierInfo this_qualifier_info =
        getQualifierInfoHash ().get ((String)all_names.elementAt (i));

      if (this_qualifier_info.isValidFor (key)) {
        return_names.add (this_qualifier_info.getName ());
      }
    }

    if (return_names.size () == 0) {
      return null;
    } else {
      return_names.sort ();
      return return_names;
    }
  }

 /**
   *  Return a vector containing the names of the required feature qualifiers
   *  for the given key.
   *  @return null if and only if no Qualifier is necessary.
   **/
  public StringVector getRequiredQualifiers (final Key key) {
    if (getQualifierInfoHash () == null) {
      return null;
    }

    final StringVector all_names = getQualifierInfoHash ().names ();

    final StringVector return_names = new StringVector ();

    for (int i = 0 ; i < all_names.size () ; ++i) {
      final QualifierInfo this_qualifier_info =
        getQualifierInfoHash ().get ((String)all_names.elementAt (i));

      if (this_qualifier_info.isRequiredFor (key)) {
        return_names.add (this_qualifier_info.getName ());
      }
    }

    if (return_names.size () == 0) {
      return null;
    } else {
      return_names.sort ();
      return return_names;
    }
  }

  /**
   *  Return true if and only if the given String contains a valid
   *  qualifier name.
   **/
  public boolean isValidQualifier (final String name) {
    if (getQualifierInfoHash () == null) {
      return true;
    }

    if (getQualifierInfoHash ().get (name) == null) {
      return false;
    } else {
      return true;
    }
  }

  /**
   *  Return true if and only if given qualifier_name is a legal qualifier
   *  name for the given key.
   **/
  public boolean isValidQualifier (final Key key,
                                   final String qualifier_name) {
    if (getUserKeys () != null && getUserKeys ().contains (key)) {
      // any qualifier is valid for a user key
      return true;
    }

    if (getQualifierInfoHash () == null) {
      return true;
    }

    final QualifierInfo qualifier_info =
      getQualifierInfoHash ().get (qualifier_name);

    if (qualifier_info == null) {
      return false;
    } else {
      return qualifier_info.isValidFor (key);
    }
  }

  /**
   *  Return true if and only if given qualifier_name is a legal qualifier
   *  name for the given key and must be present in a feature with that key.
   **/
  public boolean isRequiredQualifier (final Key key,
                                      final String qualifier_name) {
    if (getUserKeys () != null && getUserKeys ().contains (key)) {
      // there are no required qualifiers for a user key
      return false;
    }

    if (getQualifierInfoHash () == null) {
      return false;
    }

    final QualifierInfo qualifier_info =
      getQualifierInfoHash ().get (qualifier_name);

    if (qualifier_info == null) {
      return false;
    } else {
      return qualifier_info.isRequiredFor (key);
    }
  }

  /**
   *  Return true if and only if the given key is a valid in this
   *  EntryInformation.
   **/
  public boolean isValidKey (final Key key) {
    if (valid_keys == null) {
      return true;
    } else {
      if (valid_keys.contains (key)) {
        return true;
      } else {
        if (getUserKeys () != null && getUserKeys ().contains (key)) {
          return true;
        } else {
          return false;
        }
      }
    }
  }

  /**
   *  Return the QualifierInfo object for the given qualifier name, or null if
   *  there are no qualifiers with the given name associated with this
   *  EntryInformation object.
   **/
  public QualifierInfo getQualifierInfo (final String qualifier_name) {
    if (getQualifierInfoHash () == null) {
      return null;
    } else {
      return getQualifierInfoHash ().get (qualifier_name);
    }
  }

  /**
   *  Return a default EntryInformation object.  It will allow any key or
   *  qualifier.
   **/
  public static EntryInformation getDefaultEntryInformation () {
    return new SimpleEntryInformation ();
  }

  /**
   *  Returns true if strict EMBL format will be used when writing.  If true
   *  qualifiers will always wrap at 80 columns.
   *  EMBL format locations have the complement (if any) inside the join.  eg:
   *  join(complement(1..100),complement(200..300)).
   *  The default format has the complement on the outside, which is much more
   *  readable.  eg: complement(join(1..100,200..300))
   **/
  public boolean useEMBLFormat () {
    return use_embl_format;
  }

  /**
   *  Set the use_embl_format flag (see useEMBLFormat ()).  true means the
   *  strict EMBL format should be used when writing.
   **/
  public void setEMBLFormat (final boolean use_embl_format) {
    this.use_embl_format = use_embl_format;
  }

  /**
   *  Return a Vector contains all the QualifierInfo objects that this
   *  EntryInformation knows about.
   *  @return null if and only if there are no QualifierInfo objects in this
   *    object.
   **/
  public QualifierInfoHash getAllQualifierInfo () {
    if (getQualifierInfoHash () == null) {
      return null;
    } else {
      return getQualifierInfoHash ().copy ();
    }
  }

  /**
   *  Fix this EntryInformation so that the given exception won't happen
   *  again.
   **/
  public void fixException (final EntryInformationException exception) {
    if (exception instanceof InvalidKeyException) {
      final InvalidKeyException key_exception =
        (InvalidKeyException) exception;

      addKey (key_exception.getKey ());
    }

    if (exception instanceof InvalidRelationException) {
      final InvalidRelationException relation_exception =
        (InvalidRelationException) exception;

      if (relation_exception.getQualifier () == null) {
        final Key key_to_add = relation_exception.getKey ();

        addKey (key_to_add);
      } else {
        final Qualifier qualifier = relation_exception.getQualifier ();

        // change the QualifierInfo object in the EntryInformation so that it
        // can handle this key/qualifier combination

        final QualifierInfo qualifier_info =
          getQualifierInfo (qualifier.getName ());

        final QualifierInfo new_qualifier_info;

        if (qualifier_info == null) {
          new_qualifier_info =
            new QualifierInfo (qualifier.getName (),
                               QualifierInfo.QUOTED_TEXT, null, null, false);
        } else {
          new_qualifier_info =
            new QualifierInfo (qualifier_info.getName (),
                               qualifier_info.getType (), null, null, false);
        }

        try {
          addQualifierInfo (new_qualifier_info);
        } catch (QualifierInfoException e) {
          throw new Error ("internal error - unexpected exception: " + e);
        }
      }
    }

    if (exception instanceof InvalidQualifierException) {
      final String exception_qualifier_name =
        ((InvalidQualifierException) exception).getQualifier ().getName ();

      // add it again, but make it less restrictive
      final QualifierInfo new_info =
        new QualifierInfo (exception_qualifier_name, QualifierInfo.UNKNOWN,
                           null, null, false);

      try {
        addQualifierInfo (new_info);
      } catch (QualifierInfoException e) {
        throw new Error ("internal error - unexpected exception: " + e);
      }
    }
  }

  /**
   *  Returns the keys that where added with addKey ().
   *  @return null if and only if there are no user keys yet.
   **/
  public KeyVector getUserKeys () {
    return user_keys;
  }

  /**
   *  Return qualifier_info_hash.
   **/
  private QualifierInfoHash getQualifierInfoHash () {
    return qualifier_info_hash;
  }

  /**
   *  Returns valid_keys.
   **/
  private KeyVector getKeys () {
    return valid_keys;
  }

  /**
   *  This is set by the constructor and addQualifierInfo().  null means
   *  any qualifier is valid
   **/
  private QualifierInfoHash qualifier_info_hash = null;

  /**
   *  This is set by the constructor and addQualifierInfo().  null means any
   *  key is valid
   **/
  private KeyVector valid_keys = null;

  /**
   *  This is set by the constructor and by addUserKey().  null means any
   **/
  private KeyVector user_keys = null;

  /**
   *  true if the EMBL location format should be used when writing.  (See
   *  useEMBLLocationFormat () ).
   **/
  private boolean use_embl_format = false;
}
