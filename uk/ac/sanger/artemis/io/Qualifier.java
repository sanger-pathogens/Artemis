/* Qualifier.java
 *
 * created: Tue Oct 13 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/Qualifier.java,v 1.2 2005-10-11 14:20:31 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.*;

/**
 *  Each object of this class represents a group of embl qualifiers
 *  (name/value pairs) with the same name, so in effect the object contains
 *  one name and many values.
 *
 *  @author Kim Rutherford
 *  @version $Id: Qualifier.java,v 1.2 2005-10-11 14:20:31 tjc Exp $
 * */

public class Qualifier {
  /**
   *  Create a new Qualifier object.  This object consists of a name and
   *  possibly some values.  In the raw embl file we have /name=value.
   *  @param name The name of this qualifier (ie. the text immediately after
   *    the / in the qualifier)
   * @param values The values of this qualifier (ie the text immediately after
   *    the = in the qualifiers).  This argument should be null if the
   *    qualifier has no values.  The values String objects should not include
   *    the quote characters (),"" or [].  For example if the original
   *    qualifier was /citation=[3] then the value String should be: 3.  This
   *    argument is copied by the constructor.
   **/
  public Qualifier (final String name, final StringVector values) {
    initialise (name, values);
  }

  /**
   *  Create a new Qualifier object.  This object consists of a name and
   *  possibly some values.  In the raw embl file we have /name=value.
   *  @param name The name of this qualifier (ie. the text immediately after
   *    the / in the qualifier)
   *  @param value The value of this qualifier (see the other constructor for
   *    details).  Other values may be added later.
   **/
  public Qualifier (final String name, final String value) {
    initialise (name, new StringVector (value));
  }

  /**
   *  Create a new Qualifier object.  This object will consist of just a name.
   *  In the raw embl file we have /name.
   *  @param name The name of this qualifier (ie. the text immediately after
   *    the / in the qualifier)
   **/
  public Qualifier (final String name) {
    initialise (name, (StringVector) null);
  }

  /**
   *  Called by the constructor to initialise the object.
   *  @param name The name of this qualifier (ie. the text immediately after
   *    the / in the qualifier)
   *  @param values The values of this qualifier (see the other constructor for
   *    details).  Other values may be added later.
   **/
  private void initialise (final String name, final StringVector values) {
    this.name = name;
    if (values == null) {
      this.values = null;
    } else {
      if (values.size () == 0) {
        throw new Error ("internal error - zero length values vector");
      }

      this.values = values.copy ();
//    this.values = values;
    }
  }

  /**
   *  Return true if and only if all the values of this Qualifier are null.
   **/
  private boolean allValuesNull () {
    if (values == null) {
      return true;
    } else {
      for (int i = 0 ; i < values.size () ; ++i) {
        if (values.elementAt (i) != null) {
          return false;
        }
      }

      return true;
    }
  }

  /**
   *  Return true if and only if all the values of this Qualifier are not null.
   **/
  private boolean allValuesNotNull () {
    if (values == null) {
      return false;
    } else {
      for (int i = 0 ; i < values.size () ; ++i) {
        if (values.elementAt (i) == null) {
          return false;
        }
      }

      return true;
    }
  }
  /**
   *  Return the name of this qualifier.
   **/
  public String getName () {
    return name;
  }

  /**
   *  Return the values of this qualifier.  The StringVector that is returned
   *  is a copy (use addValues (), removeValue () etc. to change the
   *  Qualifier), but the String objects inside the vector are not copied so
   *  the same String references will occur in all copies.
   **/
  public StringVector getValues () {
    if (values == null) {
      return null;
    } else {
      return values.copy ();
    }
  }

  /**
   *  Add the given values to this object.
   **/
  public void addValues (final StringVector new_values) {
    if (values == null) {
      values = new StringVector ();
      values.add ((String)null);
    }
    if (new_values == null) {
      // if new_values is null then we have a qualifier with no value (like
      // /pseudo).  we add one null value for each occurrence of the
      // qualifier.
      values.add ((String)null);
    } else {
      for (int i = 0 ; i < new_values.size () ; ++i) {
        values.add (new_values.elementAt (i));
      }
    }
  }

  /**
   *  Add the given value to this object.
   **/
  public void addValue (final String new_value) {
    if (values == null) {
      values = new StringVector ();
    }
    values.add (new_value);
  }

  /**
   *  Remove the given value from this object.  The reference of the argument
   *  String must have previously been passed to addValue () or addValues ()
   *  for a removal to occur.
   **/
  public void removeValue (final String value) {
    values.remove (value);

    if (values.size () == 0) {
      values = null;
    }
  }

  /**
   *  Return the reference of a new copy of this Qualifier.
   **/
  public Qualifier copy () {
    return new Qualifier (getName (), getValues ());
  }

  /**
   *  Return true if object has the same name and values as this Qualifier.
   *  The values can be in any order.
   **/
  public boolean equals (final Object object) {
    if (object instanceof Qualifier) {
      final Qualifier comp_qualifier = (Qualifier) object;

      if (comp_qualifier.getName ().equals (getName ())) {
        if (comp_qualifier.values == null && values == null) {
          return true;
        } else {
          if (comp_qualifier.values == null || values == null) {
            return false;
          } else {
            if (comp_qualifier.values.size () != values.size ()) {
              return false;
            }

            final StringVector this_string_vector = values.copy ();
            final StringVector comp_string_vector =
              comp_qualifier.values.copy ();

            this_string_vector.sort ();
            comp_string_vector.sort ();

            for (int i = 0 ; i < values.size () ; ++i) {
              final String this_string = (String)this_string_vector.elementAt (i);

              if (this_string != null &&
                  !this_string.equals (comp_string_vector.elementAt (i))) {
                return false;
              }
            }
            return true;
          }
        }
      } else {
        return false;
      }
    } else {
      return false;
    }
  }

  /**
   *  The name that was passed to the constructor.
   **/
  private String name;

  /**
   *  The values that were passed to the constructor.
   **/
  private StringVector values;
}


