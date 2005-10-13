/* StringVector.java
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/util/StringVector.java,v 1.5 2005-10-13 12:06:12 tjc Exp $
 */

package uk.ac.sanger.artemis.util;

import java.util.Vector;
import java.util.Collections;
import java.util.Comparator;

/**
 *  This class implements a Vector of String objects.
 *
 *  @author Kim Rutherford
 *  @version $Id: StringVector.java,v 1.5 2005-10-13 12:06:12 tjc Exp $
 **/

public class StringVector extends Vector
{
  /**
   *  Create a new vector of String objects.
   **/
  public StringVector()
  {
    super();
  }

  /**
   *  Create a new vector which contains the given Strings.
   **/
  public StringVector(final String[] new_strings) 
  {
    super(new_strings.length);
    int len = new_strings.length;
    for(int i = 0; i < len; ++i) 
      add(new_strings[i]);
  }

  public StringVector(final String new_string)
  {
    super();
    add(new_string);
  }

  /**
   *  Create a new vector which contains the given Strings.
   **/
  public StringVector(final StringVector new_strings) 
  {
    super(new_strings);
  }

  /**
   *  Call add() on each of the String objects in the given StringVector.
   **/
  public void add(final StringVector new_strings) 
  {
    for (int i = 0; i < new_strings.size(); ++i)
      add (new_strings.elementAt(i));
  }

  /**
   *  Sorts the elements of the vector using quicksort from the collections
   *  package.
   */
  public void sort()
  {
    final Comparator comparator = new Comparator()
    {
      public int compare(Object fst, Object snd) 
      {
        if(fst == null) 
        {
          if(snd == null)
            return 0;
          else
            return -1;
        } 
        else 
        {
          if(snd == null)
            return 1;
        }
        return ((String)fst).compareTo((String) snd);
      }
    };

    Collections.sort(this, comparator);
  }

  /**
   *  Return a new copy of this object.
   **/
  public StringVector copy() 
  {
    return new StringVector(this);
  }

  /**
   *  Return a StringVector containing the values of the given String after
   *  splitting using the given characters.  If the argument String is zero
   *  length or it consists only of the characters used to split, the return
   *  vector will be zero length.
   *  @param keep_zero_char_tokens If true then zero width tokens will be
   *    returned.  eg. when spliting on tabs if this parameter is true then
   *    splitting this "\t\tfoo" will return "" and "foo".  If this flag is
   *    false then the split_characters will be treated as a block (and "foo"
   *    would be returned in the example.
   **/
  public static StringVector getStrings(final String argument,
                                        String split_characters,
                                        final boolean keep_zero_char_tokens) 
  {
    final StringVector return_vector = new StringVector();

    String value;
    String last_value = null;

    int ind1 = 0;
    int ind2;
    int argLen  = argument.length();

    while(ind1 < argLen)
    {
      ind2 = argument.indexOf(split_characters,ind1);
      if(ind2 == ind1)
      {
        ind1++;
        continue;
      }

      if(ind2 < 0)
        ind2 = argLen;
 
      value = argument.substring(ind1,ind2);
      ind1 = ind2+1;

      if(value.length() == 1 &&
         split_characters.indexOf(value.charAt(0)) != -1) 
      {
        // ignore the split characters

        if(keep_zero_char_tokens &&
           (last_value == null ||
            last_value != null && last_value.length () == 1 &&
            split_characters.indexOf (last_value) != -1)) 
        {
          // we need to add a space because of two split_characters in a row
          return_vector.add("");
        }
      } 
      else
        return_vector.add(value);

      last_value = value;
    }

    return return_vector;
  }

  /**
   *  Return a StringVector containing the values of the given String after
   *  splitting using the given characters.  If the argument String is zero
   *  length or it consists only of the characters used to split, the return
   *  vector will be zero length.
   **/
  public static StringVector getStrings(final String argument,
                                        final String split_characters) 
  {
    return getStrings(argument, split_characters, false);
  }

  /**
   *  Return a StringVector containing the values of the given String after
   *  splitting on whitespace.  The return object contains one String for each
   *  sequence of non-whitespace characters in the argument.  If the argument
   *  String is zero length or it consists only of whitespace, the return
   *  vector will be zero length.
   **/
  public static StringVector getStrings(final String argument) 
  {
    return getStrings(argument, " ", false);
  }

  public static void main(String args[])
  {
    String argument = "a a g g g c a c g t c g c a t c g a c t c";
    long startTime = System.currentTimeMillis();

    for(int i=0; i<10000000; i++)
      getStrings(argument, " ", true);

    long endTime = System.currentTimeMillis();

    System.out.println("TIME TAKEN "+  Long.toString(endTime-startTime));
  }

}
