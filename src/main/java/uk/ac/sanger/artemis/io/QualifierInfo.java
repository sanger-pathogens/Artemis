/* QualifierInfo.java
 *
 * created: Sun Feb 21 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/QualifierInfo.java,v 1.1 2004-06-09 09:50:07 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

/**
 *  Each object of this class contains the information about one qualifier
 *  type.
 *
 *  @author Kim Rutherford
 *  @version $Id: QualifierInfo.java,v 1.1 2004-06-09 09:50:07 tjc Exp $
 **/

public class QualifierInfo {
  /**
   *  Create a new QualifierInfo object.
   *  @param name The name of the qualifier.
   *  @param type The type of the qualifier (one of QUOTED_TEXT, TEXT,
   *    NO_VALUE or OPTIONAL).
   *  @param keys The possible keys this qualifier can be associated
   *    with. (null means the qualifier is valid for all features).
   *  @param required_keys The keys of the features in which this qualifier is
   *    compulsary. (null means the qualifier is not required by any feature).
   *  @param once_only If true this qualifier may occur only once in a
   *    feature (ie. the qualifier can only have one value).
   **/
  public QualifierInfo (final String name, final int type,
                        final KeyVector keys, final KeyVector required_keys,
                        final boolean once_only) {
    this.name = name;
    this.type = type;
    this.once_only = once_only;

    if (keys != null) {
      this.keys = keys.copy ();
    }

    if (required_keys != null) {
      this.required_keys = required_keys.copy ();
    }
  }
  
  /**
   *  Used if the value of the qualifier should a text string quoted with
   *  double quotes: "some text".
   **/
  public final static int QUOTED_TEXT          = 1;

  /**
   *  Used if the value of the qualifier should be an unquoted text string
   *  (containing no spaces).
   **/
  public final static int TEXT                 = 2;

  /**
   *  Used if this qualifier does not take a value (for example /partial)
   **/
  public final static int NO_VALUE             = 3;

  /**
   *  Used if the value of this qualifier should be quoted but is optional.
   **/
  public final static int OPTIONAL_QUOTED_TEXT = 4;

  /**
   *  Used if the value of this qualifier can have any type.
   **/
  public final static int ANY                  = 5;

  /**
   *  Used for other types.
   **/
  public final static int UNKNOWN             = 0;

  /**
   *  Returns the type ID of this QualifierInfo object.
   **/
  public int getType () {
    return type;
  }

  /**
   *  Return the name of the QualifierInfo object ("label", "note", etc.)
   **/
  public String getName () {
    return name;
  }

  /**
   *  Return if and only if a qualifier of this type is allowed in a feature
   *  with the given key.
   **/
  public boolean isValidFor (final Key key) {
    if (keys == null) {
      return true;
    } else {
      return keys.contains (key);
    }
  }

  /**
   *  Return if and only if a qualifier of this type is required in a feature
   *  with the given key.
   **/
  public boolean isRequiredFor (final Key key) {
    if (required_keys == null) {
      return false;
    } else {
      return required_keys.contains (key);
    }
  }

  /**
   *  Return a list of the Keys for which a qualifier of this type is
   *  valid. (null means the qualifier is valid for all features)
   **/
  public KeyVector getValidKeys () {
    return keys;
  }
  
  /**
   *  Return a list of the Keys for which a qualifier of this type is
   *  required.  (null means the qualifier is not required by any feature)
   **/
  public KeyVector getRequiredKeys () {
    return required_keys;
  }
  
  /**
   *  Return true if and only if this qualifier can appear only once in a
   *  feature.
   **/
  public boolean isOnceOnly () {
    return once_only;
  }

  /**
   *  A String array containing the possible type strings.  This is used by
   *  getTypeID () to check that a valid type string was passed to the
   *  constructor.
   **/
  private static final String [] possible_type_strings = {
    "\"text\"", "text", "\"list\"", "\"opt\"", "(::)", "EC", "feature", "list",
    "location", "modbase", "none", "number", "real", "ref", "item"
  };

  /**
   *  Return the type ID for the given type_string.  "text" returns
   *  QualifierInfo.QUOTED_TEXT, text returns TEXT, none returns NO_VALUE and
   *  "opt" returns OPTIONAL_QUOTED_TEXT.
   **/
  public static int getQualifierTypeID (final String type_string) {
    for (int i = 0 ; i < possible_type_strings.length ; ++i) {
      final String this_type_string = possible_type_strings[i];

      if (this_type_string.equals (type_string)) {

        // found a match so return the correct ID

        if (type_string.equals ("\"text\"") ||
            type_string.equals ("\"list\"")) {
          return QualifierInfo.QUOTED_TEXT;
        } else {
          if (type_string.equals ("text")) {
            return QualifierInfo.TEXT;
          } else {
            if (type_string.equals ("none")) {
              return QualifierInfo.NO_VALUE;
            } else {
              if (type_string.equals ("\"opt\"")) {
                return QualifierInfo.OPTIONAL_QUOTED_TEXT;
              } else {
                // default is text which covers all the non quoted types
                return QualifierInfo.TEXT;
              }
            }
          }
        }
      }
    }

    throw new Error ("unknown type string: " + type_string);
  }

  /**
   *  The name string that was passed to the constructor.
   **/
  private String name;

  /**
   *  The type ID of this object (as passed to the constructor).
   **/
  private int type;

  /**
   *  The valid Keys for this object (as passed to the constructor).
   **/
  private KeyVector keys = null;

  /**
   *  The required Keys for this object (as passed to the constructor).
   **/
  private KeyVector required_keys = null;

  /**
   *  (See constructor)
   **/
  private boolean once_only;
}


