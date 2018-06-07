/* LocationLexer.java
 *
 * created: Tue Oct  6 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/LocationLexer.java,v 1.1 2004-06-09 09:49:50 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

/**
 *  The LocationLexer class provides methods for breaking a EMBL feature
 *  location string into tokens.  The complete list of possible tokens is
 *  given below.
 *
 *
 *  @author Kim Rutherford
 *  @version $Id: LocationLexer.java,v 1.1 2004-06-09 09:49:50 tjc Exp $
 *
 */

public class LocationLexer {
  /**
   *  Create a new LocationLexer object that can be used to tokenise
   *  location_string.
   */
  public LocationLexer (String location_string) {
    this.location_string = location_string;
  }

  /**
   *  Return a TokenEnumeration containing all the tokens in this string.
   */
  public TokenEnumeration getTokens () {
    return new TokenEnumeration (location_string);
  }

  /**
   *  This is a helper class for LocationLexer - see LocationLexer.getTokens ()
   *
   */
  public class TokenEnumeration {
    public TokenEnumeration (String location_string) {
      this.location_string = location_string;
    }

    /**
     *  Return the next token but don't remove it from the enumeration.
     *  @return The next token or null if there are no more.
     */
    public Object peekElement () {
      if (peeked_object == null) {
        peeked_object = removeNextToken ();
      }

      return peeked_object;
    }

    /**
     *  Return the next token and remove it from the enumeration
     */
    public Object nextElement () {
      Object o = removeNextToken ();
      return o;
    }

    /**
     *  Try to "eat" the first token in the enumeration.  If the next token is
     *  the same as the argument String then the next token will be removed
     *  from the enum and it will return true, otherwise the enum will not
     *  change and it wil return false.
     */
    public boolean eatToken (String token) {
      if (peekElement () instanceof String &&
          ((String)peekElement ()).equals (token)) {
        nextElement ();
        return true;
      } else {
        return false;
      }
    }

    /**
     *  Try to "eat" the first token in the enumeration.  If the next token is
     *  the same as the argument Character then the next token will be removed
     *  from the enum and it will return true, otherwise the enum will not
     *  change and it wil return false.
     */
    public boolean eatToken (final char token) {
      if (peekElement () instanceof Character &&
          ((Character)peekElement ()).charValue () == token) {
        nextElement ();
        return true;
      } else {
        return false;
      }
    }

    /**
     *  Return a String contains all the remaining tokens concatenated
     *  together.
     **/
    public String toString () {
      // the number we pick for the initial StringBuffer size is not critical,
      // but should cover most possibilities
      final StringBuffer spare_tokens_string = new StringBuffer (100);

      Object next_token = nextElement ();
      
      while (next_token != null) {
        spare_tokens_string.append (next_token.toString ());
        next_token = nextElement ();
      }

      return spare_tokens_string.toString ();
    }

    /**
     *  Return the next token and logically remove it from the start of the
     *  enumeration.  Returns null when there are no more tokens.
     *  removeNextToken () will return peeked_object (and set it to null)
     *  rather than removing a token iff peeked_object is not null.
     */
    private Object removeNextToken () {
      if (peeked_object == null) {

        // loop until we get to the end of location_string or until we return
        // a token.
        while (true) {
          if (next_char_index == location_string.length ()) {
            // all tokens have been read
            return null;
          }

          final char current_char = location_string.charAt (next_char_index);

          switch (current_char) {

          case ' ': case '\t':
            // go around the loop again
            next_char_index++;
            continue;

          case '(': case ')': case ',': case '^': case ':':
            // handle single character tokens (except ".")
            next_char_index++;
            return new Character (current_char);

          case '>':
            next_char_index++;
            if (next_char_index < location_string.length () &&
                Character.isDigit (location_string.charAt (next_char_index))) {
              return new UpperInteger (removeInteger ());
            } else {
              return new Character ('>');
            }

          case '<':
            next_char_index++;
            if (next_char_index < location_string.length () &&
              Character.isDigit (location_string.charAt (next_char_index))) {
              return new LowerInteger (removeInteger ());
            } else {
              return new Character ('<');
            }

          case '0': case '1': case '2': case '3': case '4':
          case '5': case '6': case '7': case '8': case '9':
            return removeInteger ();

          case '.':
            if (next_char_index + 1 == location_string.length ()) {
              // special case we have a "." at the end of the location string -
              // the parser will catch this are report the error
              next_char_index++;
              return new Character ('.');
            } else {
              if (location_string.charAt (next_char_index + 1) == '.') {
                next_char_index += 2;
                return "..";
              }
              else {
                next_char_index ++;
                return new Character ('.');
              }
            }

          default:
            {
              // everything else is a label (eg AF009694), functional name
              // (ie complement, join or order) or garbage
              final String label = removeLabel ();

              if (label.equals ("")) {
                // couldn't read a label so just return the current character
                // and let the parser sort it out

                next_char_index++;
                return new String ("" + current_char);
              } else {
                return label;
              }
            }
          }
        }
      } else {
        final Object tmp_object = peeked_object;
        peeked_object = null;
        return tmp_object;
      }
    }


    /**
     *  Reads an integer from the current position (next_char_index) in
     *  location_string and increments next_char_index to point to the next
     *  non digit in location_string.
     */
    private Integer removeInteger () {
      String integer_string = "";

      char current_char = location_string.charAt (next_char_index);

      while (Character.isDigit (current_char)) {
        integer_string += current_char;

        next_char_index++;

        if (next_char_index >= location_string.length ()) {
          break;
        }

        current_char = location_string.charAt (next_char_index);
      }

      return new Integer (integer_string);
    }


    /**
     *  Remove a string or letters, digits and colons from location_string and
     *  adjust next_char_index appropriately.  Returns an empty String if the
     *  next character is not alphanumeric.
     */
    private String removeLabel () {
      String return_string = "";

      char current_char = location_string.charAt (next_char_index);

      if (!Character.isLetter (current_char)) {
        // first character must be a letter
        return "";
      }

      while (Character.isLetterOrDigit (current_char) ||
             current_char == '.' || current_char == '_' ||
             current_char == '*' || current_char == '\'' ||
             current_char == '-') {
        return_string += current_char;

        next_char_index++;

        if (next_char_index >= location_string.length ()) {
          break;
        }

        current_char = location_string.charAt (next_char_index);
      }
      
      return return_string;
    }

    /**
     *  Contains string passed to the constructor.
     */
    private String location_string;


    /**
     *  A pointer into location_string indicating the next character we should
     *  read.
     */
    private int next_char_index = 0;

    /**
     *  If peekElement () has been called then the token we read from the
     *  remaining_string is stored here. (see removeNextToken ()).
     */
    private Object peeked_object = null;
  }


  /**
   *  This contains the String that was passed to the constructor
   */
  private String location_string;
}


