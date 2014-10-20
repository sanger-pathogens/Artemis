/* StreamQualifier.java
 *
 * created: Wed Dec 30 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/StreamQualifier.java,v 1.3 2008-11-07 17:54:26 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.*;

import java.io.BufferedReader;
import java.io.IOException;

/**
 *  This class contains routines for reading and writing Qualifiers.
 *
 *  @author Kim Rutherford
 *  @version $Id: StreamQualifier.java,v 1.3 2008-11-07 17:54:26 tjc Exp $
 **/

public // XXX

class StreamQualifier {
  /**
   *  Create a new Qualifier object by unquoting the value part of the
   *  Qualifier and the calling the Qualifier constructor.  This object
   *  consists of a name and a value.  In the raw embl file we have
   *  /name=value.
   *  @param name The name of this qualifier (ie. the text immediately after
   *    the / in the qualifier)
   *  @param value The value of this qualifier (ie the text immediately after
   *    the = in the qualifier).  This argument may be null if the qualifier
   *    has no value.  Unlike the Qualifier constructor the value String
   *    should include the quote characters (),"" or [] if the original
   *    qualifier contains them.  For example if the original qualifier was
   *    /citation=[3] then the value String should be: [3].
   *  @exception QualifierParseException Thrown if the value String is
   *    incorrectly quoted.
   **/
  public static Qualifier
    makeStreamQualifier (final String name,
                         final String value,
                         final EntryInformation entry_information)
      throws QualifierParseException {

    if (!entry_information.isValidQualifier (name)) {
      // use this qualifier value to decide how qualifiers with this name
      // should be quoted
      final QualifierInfo new_qualifier_info;
      if (value.startsWith ("\"")) {
        new_qualifier_info =
          new QualifierInfo (name, QualifierInfo.QUOTED_TEXT,
                             null, null, false);
      } else {
        new_qualifier_info =
          new QualifierInfo (name, QualifierInfo.TEXT, null, null, false);
      }

      try {
        entry_information.addQualifierInfo (new_qualifier_info);
      } catch (QualifierInfoException e) {
        // this shouldn't happen because we have just checked that there is no
        // qualifier with this name
        throw new Error ("internal error - unexpected exception: " + e);
      }
    }

    return new Qualifier (name, unquote (value));
  }

  /**
   *  Return a String version of the given Qualifier.
   *  @param qualifier_info Used to determine how to quote the qualifiers
   **/
  public static String toString (final QualifierInfo qualifier_info,
                                 final Qualifier qualifier) {
    final StringVector values = qualifier.getValues ();

    if (values == null) {
      return '/' + qualifier.getName ();
    } else {
      // the number we pick for the initial StringBuffer size is not critical,
      // but should cover most possibilities
      final StringBuffer buffer = new StringBuffer (50);

      for (int i = 0 ; i < values.size () ; ++i) {
        buffer.append ('/');
        buffer.append (qualifier.getName ());
        if (values.elementAt (i) != null) {
          buffer.append ('=');
          buffer.append (quotedValue (qualifier_info,
                                      qualifier.getName (),
                                      (String)values.elementAt (i).replaceAll("(^|[^\"])\"([^\"]|$)","$1\"\"$2")));
        }
      }

      return buffer.toString ();
    }
  }

  /**
   *  Return a StringVector containing one String (of the form /name=value)
   *  for each of the values of the given qualifier.
   *  @param qualifier_info Used to determine how to quote the qualifiers
   **/
  public static StringVector
    toStringVector (final QualifierInfo qualifier_info,
                    final Qualifier qualifier) {
    final StringVector values = qualifier.getValues ();

    final StringVector return_vector = new StringVector ();

    if (values == null) {
      return_vector.add ('/' + qualifier.getName ());
    } else {
      for (int i = 0 ; i < values.size () ; ++i) {
        // the number we pick for the initial StringBuffer size is not
        // critical
        final StringBuffer buffer = new StringBuffer (50);

        buffer.append ('/');
        buffer.append (qualifier.getName ());
        if (values.elementAt (i) != null) {
          buffer.append ('=');
          buffer.append (quotedValue (qualifier_info,
                                      qualifier.getName (),
                                      (String)values.elementAt (i).replaceAll("(^|[^\"])\"([^\"]|$)","$1\"\"$2")));
        }
        return_vector.add (buffer.toString ());
      }
    }

    return return_vector;
  }

  /**
   *  This is used by readFromStream () as temporary storage.  It is a class
   *  member rather than a local variable so that we don't need to allocate a
   *  object for each call.  The number we pick for the initial StringBuffer
   *  size is not critical, but should cover most possibilities to prevent
   *  reallocation.
   **/
  final private static StringBuffer read_name_string_buffer =
    new StringBuffer (20);

  /**
   *  Read a qualifier name from a stream.
   *  @param buffered_reader the stream to read from
   *  @return the qualifier name if successful, otherwise null
   */
  static String readName (final BufferedReader buffered_reader)
      throws QualifierParseException, IOException {

    int current_char;

    while ((current_char = buffered_reader.read ()) != -1 &&
           0 != current_char    // Kaffe 1.00 returns 0 at end of string
           ) {
      if (' ' == current_char ||
          '\n' == current_char ||
          '\r' == current_char ||
          '\t' == current_char) {
        // read a whitespace character so go back to the top of the loop
        continue;
      } else {
        if ('/' == current_char) {
          // we have found the start of the qualifier name
          break;
        } else {
          // if the character isn't a / or space then something is wrong
          throw new QualifierParseException ("failed to read a qualifier " +
                                             "name from this string: " +
                                             (char)current_char +
                                             buffered_reader.readLine ());
        }
      }
    }

    if (-1 == current_char ||
        0 == current_char       // Kaffe 1.00 returns 0 at end of string
        ) {
      // end of file
      return null;
    }

    buffered_reader.mark (1);

    read_name_string_buffer.setLength (0);

    while ((current_char = buffered_reader.read ()) != -1) {
      if (Character.isLetter ((char) current_char) ||
          Character.isDigit ((char) current_char) ||
          '_' == current_char ||
          '+' == current_char) {

        read_name_string_buffer.append ((char) current_char);

        // save the new position and go around the loop again
        buffered_reader.mark (1);
        continue;
      } else {
        // we have read one character too many
        buffered_reader.reset ();
        break;
      }
    }

    final String return_string = read_name_string_buffer.toString ();

    if (return_string.length () == 0) {
      throw new QualifierParseException ("zero length qualifier name read " +
                                         "from this string: " +
                                         buffered_reader.readLine ());
    } else {
      return return_string;
    }
  }

  /**
   *  This is used by readFromStream () as temporary storage.  It is a class
   *  member rather than a local variable so that we don't need to allocate a
   *  object for each call.  The number we pick for the initial array
   *  size is not critical, but should cover most possibilities to prevent
   *  reallocation.
   **/
  private static char [] read_value_buffer = new char [5000];

  /**
   *  The index into read_value_buffer - used by readValue () to keep track of
   *  where to put the next character.
   **/
  private static int buffer_index = 0;


  /**
   *  Append the given char to read_value_buffer (at the position
   *  buffer_index), reallocating the buffer if necessary
   **/
  private static void appendToValueBuffer (final char new_char) {
    if (buffer_index >= read_value_buffer.length) {
      // reallocate as the buffer is full

      final char [] temp_buffer = new char [read_value_buffer.length*2];

      System.arraycopy (read_value_buffer, 0,
                        temp_buffer, 0,
                        read_value_buffer.length);
      read_value_buffer = temp_buffer;
    }

    read_value_buffer [buffer_index++] = (char) new_char;
  }

  /**
   *  Read a qualifier value from a stream.
   *  @param buffered_reader the stream to read from
   *  @return the qualifier value if successful, otherwise null
   *  @exception QualifierParseException Thrown if the format of the
   *    value String is not appropriate for a Qualifier with the given name or
   *    if the qualifier can't be read.
   *    Each qualifier has a specific format for the value part which depends
   *    on the name, for example the value part of /codon_start qualifier must
   *    be a number: 1, 2 or 3.
   */
  static synchronized String readValue (final BufferedReader buffered_reader)
      throws QualifierParseException, IOException {

    buffer_index = 0;

    buffered_reader.mark (1);

    int current_char = buffered_reader.read ();

    if (-1 == current_char) {
      return "";
    }

    // this is the character the marks the end of the value string.  the
    // default value of 0 means a '/' should end the string.
    char final_char = 0;

    // this will be set to ", [ or ( if the value starts with one of those
    // characters
    char start_char = 0;

    // this is is used to balance the round or square brackets.  it is
    // incremented each time an open bracket is seen (after the first one) and
    // decremented each time a close bracket is seen (except for the last).
    //
    int bracket_count = 0;

    if ('"' == current_char) {
      final_char = '"';
    }
    if ('[' == current_char) {
      final_char = ']';
      start_char = '[';
      ++bracket_count;
    }
    if ('(' == current_char) {
      final_char = ')';
      start_char = '(';
      ++bracket_count;
    }

    if (0 == final_char) {
      // the character we read isn't one of the delimiter characters so put it
      // back
      buffered_reader.reset ();
    } else {
      // append the char now so that loop doesn't stop immediately in the '"'
      // case
      appendToValueBuffer ((char) current_char);
    }

    buffered_reader.mark (1);

    while ((current_char = buffered_reader.read ()) != -1) {

      // change newlines and other control characters to spaces
      if (Character.isISOControl ((char)current_char) &&
          current_char != '\t') {
        current_char = ' ';
      }
      
      if (current_char != '"') {
        if (current_char == start_char) {
          ++bracket_count;
        } else {
          if (current_char == final_char) {
            --bracket_count;
          }
        }
      }

      if (current_char == final_char && bracket_count == 0) {

        if (current_char == '"') {
          // check for two quotes in a row

          // since the current character is a quote we know we can change the
          // mark
        
          buffered_reader.mark (1);

          final int next_char = buffered_reader.read ();

          if (next_char == '"') {
            // we have hit a quoted quote
            appendToValueBuffer ('"');
            appendToValueBuffer ('"');
            continue;
          } else {
            // end of line or next qualifier
            
            if (next_char != -1) {
              buffered_reader.reset ();
            }

            appendToValueBuffer ('"');
            break;
          }
        } else {
          // end of value
          appendToValueBuffer ((char) current_char);
          break;
        }
      } else {

        if (0 == final_char && '/' == current_char) {
          // in this case '/' marks the end of the value. we need to push back
          // the '/' so that reading the next qualifier will work
          buffered_reader.reset ();
          break;
        } else {
          appendToValueBuffer ((char) current_char);

          // save the new position and go around the loop again
          buffered_reader.mark (1);
          continue;
        }
      }
    }

    if (bracket_count > 0) {
      throw new QualifierParseException ("hit the end of line while looking " +
                                         "for a \"" + final_char + "\"");

    }

    // move buffer_index back past any whitespace
    while (buffer_index > 0 &&
           Character.isWhitespace (read_value_buffer[buffer_index-1])) {
      --buffer_index;
    }

    return new String (read_value_buffer, 0, buffer_index);
  }

  /**
   *  Return the value part of a Qualifier correctly quoted for insertion into
   *  a embl entry.
   *  @param qualifier_info The type of the qualifier that we will quote.
   *  @param name The name part of qualifier.  The quote characters to check
   *    for depend on this name.
   *  @param value Quote this value.
   **/
  private static String quotedValue (final QualifierInfo qualifier_info,
                                     final String name, final String value) {
    if (qualifier_info != null &&
        (qualifier_info.getType () == QualifierInfo.QUOTED_TEXT ||
         qualifier_info.getType () == QualifierInfo.OPTIONAL_QUOTED_TEXT)) {
      return '"' + value + '"';
    } else {
      if (value.indexOf ('/') != -1) {
        // quote it anyway
        return '"' + value + '"';
      } else {
        return value;
      }
    }
  }

  /**
   *  Return the value part of a qualifier with any quote characters removed.
   *  @param name The name part of qualifier.  The quote characters to check
   *    for depend on this name.
   *  @param value This is the value String to strip the quote characters
   *    from.  This may be null if this qualifier has no value part (for
   *    example /partial).
   *  @return The unquoted version of the qualifier value or null if the value
   *    passed to unquote() is null.
   *  @exception QualifierParseException Thrown if the value String is
   *    incorrectly quoted for a qualifier with the given name.  For example
   *    there is a quote at one end of the value and not the other.
   **/
  private static String unquote (final String value)
      throws QualifierParseException {
    if (value.length () >= 2) {
      final char first_char = value.charAt (0);
      final char last_char = value.charAt (value.length () - 1);

      if (first_char == '"' && last_char == '"') {
        return value.substring (1, value.length () - 1);
      }
      if (first_char != '"' && last_char != '"') {
        return value;
      }

      throw new QualifierParseException ("unbalanced quotes: " + value);
    } else {
      return value;
    }
  }
}
