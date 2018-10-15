/* Key.java
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/Key.java,v 1.3 2006-09-01 10:05:04 tjc Exp $
 */

package uk.ac.sanger.artemis.io;


/**
 *  Each object in this class represents a feature key.  A key is a String
 *  with a limited number of possible values.
 *
 *  @author Kim Rutherford
 *  @version $Id: Key.java,v 1.3 2006-09-01 10:05:04 tjc Exp $
 **/

public class Key 
{
  /** A convenience copy of the CDS Key. */
  public static final Key CDS = new Key ("CDS", true);

  /** The String that was passed to the constructor. */
  private String key_string;

  /**
   *  Create a new Key object from the given String.
   *  @param key_string The text of the new Key.
   **/
  public Key (final String key_string) 
  {
    this.key_string = key_string;
  }

  /**
   *  Create a new Key object from the given String without checking to see if
   *  the key_string is valid.
   *  @param key_string The text of the new Key.
   *  @param dummy This is present only to distinguish this constructor from
   *    the previous one.
   **/
  Key (final String key_string, final boolean dummy) 
  {
    this.key_string = key_string;
  }

  /**
   *  Return the String reference that was passed to the constructor.
   **/
  public String getKeyString () 
  {
    return key_string;
  }

  /**
   *  Return a String representation of this Key.  This currently does this
   *  same as getKeyString ().
   **/
  public String toString () 
  {
    return getKeyString ();
  }

  /**
   *  Compares this Key to the given Object.  Returns true if and only if the
   *  test_object is a String or Key with the same value as this Key.
   **/
  public boolean equals (final Object test_object) 
  {
    if(test_object instanceof String) 
      return key_string.equals(test_object);
    else if(test_object instanceof Key) 
      return key_string.equals (((Key) test_object).getKeyString ());
    else 
      return false;
  }

  /**
   * Returns a hash code for this Object.
   *
   * @return  a hash code value for this object.
   */
  public int hashCode() 
  {
    return getKeyString().hashCode();
  }

  /**
   *  Return the length of the String that was passed to the constructor.
   **/
  public int length () 
  {
    return getKeyString ().length ();
  }

  /**
   *  Compares two strings lexicographically. The comparison is based on the
   *  Unicode value of each character in the strings.
   *  @param another_key The Key to be compared.
   *  @return The value 0 if the argument Key is equal to this Key; a
   *    value less than 0 if this Key is lexicographically less than
   *    the Key argument; and a value greater than 0 if this Key is
   *    lexicographically greater than the Key argument.
   **/
  public int compareTo (final Key another_key) 
  {
    return key_string.compareTo (another_key.key_string);
  }

}

