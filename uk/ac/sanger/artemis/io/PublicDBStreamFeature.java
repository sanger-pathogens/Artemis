/* PublicDBStreamFeature.java
 *
 * created: Tue Sep 14 1999
 *
 * This file is part of Artemis
 *
 * Copyright (C) 1999  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/PublicDBStreamFeature.java,v 1.4 2007-04-12 10:26:40 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.*;
import java.io.*;

/**
 *  This is an implementation of Feature that can read and write itself to a
 *  EMBL or GENBANK stream.
 *
 *  @author Kim Rutherford
 *  @version $Id: PublicDBStreamFeature.java,v 1.4 2007-04-12 10:26:40 tjc Exp $
 **/

abstract public class PublicDBStreamFeature
    extends SimpleDocumentFeature
    implements DocumentFeature, StreamFeature, ComparableFeature {
  /**
   *  Create a new PublicDBStreamFeature object.  The feature should be added
   *  to an Entry (with Entry.add ()).
   *  @param key The new feature key
   *  @param location The Location object for the new feature
   *  @param qualifiers The qualifiers for the new feature
   **/
  public PublicDBStreamFeature (final Key key,
                                final Location location,
                                final QualifierVector qualifiers) {
    super (null);
    try {
      setKey (key);
      setLocation (location);
      setQualifiers (qualifiers);
    } catch (EntryInformationException e) {
      // this should never happen because the feature will not be in an Entry
      throw new Error ("internal error - unexpected exception: " + e);
    } catch (ReadOnlyException e) {
      // this should never happen because the feature will not be in an Entry
      throw new Error ("internal error - unexpected exception: " + e);
    } catch (OutOfRangeException e) {
      // this should never happen because the feature will not be in an Entry
      throw new Error ("internal error - unexpected exception: " + e);
    }
  }

  /**
   *  Create a new PublicDBStreamFeature with the same key, location and
   *  qualifiers as the given feature.  The feature should be added to an
   *  Entry (with Entry.add ()).
   *  @param feature The feature to copy.
   **/
  public PublicDBStreamFeature (final Feature feature) {
    super (null);
    try {
      setKey (feature.getKey ());
      setLocation (feature.getLocation ());
      setQualifiers (feature.getQualifiers ());
    } catch (EntryInformationException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    } catch (ReadOnlyException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    } catch (OutOfRangeException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    }
  }

  /**
   *  Write this Feature to the given stream.
   *  @param writer The stream to write to.
   *  @exception IOException thrown if there is an io problem while writing
   *    the Feature.
   **/
  public synchronized void writeToStream (final Writer writer)
      throws IOException {
    writeKey (writer);
    writeLocation (writer);
    writeQualifiers (writer);
  }

  /**
   *  Read and return a PublicDBStreamFeature from a stream.  A feature must
   *  be the next thing in the stream.
   *  @param stream the Feature is read from this stream
   *  @param feature_type this flag indicates whether to read the feature as
   *    an EMBL feature (flag == LineGroup.EMBL_FEATURE_TABLE) or as a GENBANK
   *    feature (flag == LineGroup.GENBANK_FEATURE_TABLE).
   *  @exception IOException thrown if there is a problem reading the Feature -
   *    most likely ReadFormatException.
   *  @return null if in_stream is at the end of file when the method is
   *  called
   **/
  protected static PublicDBStreamFeature
    readFromStream (LinePushBackReader stream,
                    int feature_type)
      throws IOException {

    final PublicDBStreamFeature new_feature;

    if (feature_type == LineGroup.EMBL_FEATURE) {
      new_feature = new EmblStreamFeature ();
    } else {
      new_feature = new GenbankStreamFeature ();
    }

    try {
      new_feature.setFromStream (stream);
    } catch (EOFException e) {
      return null;
    }

    return new_feature;
  }

  /**
   *  This is used by readFromStream () as temporary storage.  It is a class
   *  member rather than a local variable so that we don't need to allocate a
   *  object for each call.  The number we pick for the initial StringBuffer
   *  size is not critical, but should cover most possibilities
   **/
  final static private StringBuffer qualifier_string_buffer =
    new StringBuffer (1500);

  /**
   *  Read the details of a feature from an EMBL stream into the current
   *  object.  (Called only by readFromStream ()).
   *  @param in_stream the Feature is read from this stream
   *  @exception IOException thrown if there is a problem reading the Feature -
   *    most likely ReadFormatException.
   **/
  private void setFromStream (final LinePushBackReader in_stream)
      throws IOException {

    final String first_line = in_stream.readLine ();

    if (first_line == null) {
      // we tried to read a Feature starting at the end of file
      throw new EOFException ("while reading a feature");
    }

    final int line_type = getLineType (first_line);

    if (this instanceof EmblStreamFeature &&
        line_type != LineGroup.EMBL_FEATURE ||
        this instanceof GenbankStreamFeature &&
        line_type != LineGroup.GENBANK_FEATURE) {
      // this line is not the first line of a feature

      in_stream.pushBack (first_line);

      throw new EOFException ("end of feature table");
    }

    if (first_line.length () < 15) {
      throw new ReadFormatException ("line too short",
                                     in_stream.getLineNumber ());
    }

    final String key_string =
      getKeyStringFromLine (first_line, in_stream.getLineNumber ());

    if (key_string == null) {
      throw new ReadFormatException ("expected the first line of a " +
                                     "feature",
                                     in_stream.getLineNumber ());
    }

    // set to true when we see the first qualifier.
    // we need this so we know when to stop adding to location_string and
    // start adding to qualifier_string.
    boolean location_string_finished = false;

    String location_string = getRestOfFeatureLine (first_line);

    qualifier_string_buffer.setLength (0);

    String line;

    // the line of the input where the qualifiers start - used for error
    // reporting
    //int qualifier_start_line = -1;

    // the line of the input where this feature starts - used for error
    // reporting
    final int feature_start_line = in_stream.getLineNumber ();

    // loop until there are no more lines in the file or we hit the start of
    // the next feature or the start of the next line group
    while (true) {
      line = in_stream.readLine ();

      if (line == null) {
        // we have reached the end of file - break out of loop
        break;
      }

      final int current_line_type = getLineType (line);

      if (this instanceof EmblStreamFeature &&
          current_line_type == LineGroup.EMBL_FEATURE ||
          this instanceof GenbankStreamFeature &&
          current_line_type == LineGroup.GENBANK_FEATURE) {

        if (getKeyStringFromLine (line, in_stream.getLineNumber ()) == null) {

          // read the text in this line after the key
          final String rest_of_line = getRestOfFeatureLine (line);

          if (rest_of_line == null) {
            throw new ReadFormatException ("line too short while reading " +
                                           "feature",
                                           in_stream.getLineNumber ());
          }

          if (location_string_finished) {

            final int qualifier_string_length =
              qualifier_string_buffer.length ();

            final char last_char =
              qualifier_string_buffer.charAt (qualifier_string_length - 1);

            if (last_char != '"') {
              // we put space in the string to handle those cases when a
              // qualifier wraps over the end of a line, but "" needs to be
              // kept as a single token
              qualifier_string_buffer.append (" ");
            }

            qualifier_string_buffer.append (rest_of_line);
          } else {
            // the last line(s) we read was part of location string

            if (rest_of_line.startsWith ("/")) {

              // we have now seen the first qualifier line so we can use the
              // location_string to create a Location object

              location_string_finished = true;

              qualifier_string_buffer.append (rest_of_line);

              //qualifier_start_line = in_stream.getLineNumber ();

            } else {
              location_string = location_string + rest_of_line;
            }
          }
          // continue with loop
        } else {

          // line has a feature key then it must be the start of a new
          // feature
          in_stream.pushBack (line);
          break;

        }
      } else {
        // this line is not a feature line so return now
        in_stream.pushBack (line);
        break;
      }
    }


    final Key key = new Key (key_string);

    final Location location;

    try {
      location = new Location (location_string);
    } catch (LocationParseException exception) {
      // re-throw the exception with the line number added

      final String new_error_string = exception.getMessage ();

      // subtract 1 because the error was on the previous line
      throw new ReadFormatException (new_error_string,
                                     feature_start_line);
    }


    final QualifierVector qualifiers;

    final String qualifier_string = qualifier_string_buffer.toString ();

    try {
      qualifiers = getQualifiersFromString (qualifier_string,
                                            getEntryInformation ());
    } catch (QualifierParseException exception) {
      // re-throw the exception with the line number added
      final String new_error_string = exception.getMessage ();

      // subtract 1 because the error was on the previous line
      throw new ReadFormatException (new_error_string,
                                     feature_start_line);

    }

    try {
      set (key, location, qualifiers);
    } catch (EntryInformationException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    } catch (OutOfRangeException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    }

    setDirtyFlag ();
  }


  /**
   *  Get the text in a line from a feature minus the FT tag and the key
   *  @param line the entry line to process
   */
  private static String getRestOfFeatureLine (String line) {

    if (line.length () < LOCATION_START_COLUMN) {
      return null;
    } else {
      return line.substring (LOCATION_START_COLUMN).trim ();
    }
  }

  private final static int KEY_FIELD_WIDTH = 16;

  /**
   *  Return the key from a embl entry line.
   *  @param line_string the text of the entry line to process
   *  @param line_number the line number of line_string in the input stream
   *  @return null if this isn't the first line of a feature, otherwise the
   *    key of this feature
   */
  private static String getKeyStringFromLine (final String line_string,
                                              final int line_number)
      throws ReadFormatException {
    // the first line of a feature starts with "FT   " - remove that first
    final String rest_of_line = getRestOfLine (line_string);

    if (rest_of_line == null ||
        rest_of_line.startsWith (" ")) {
      // this isn't the first line of a feature
      return null;
    } else {
      if (rest_of_line.length () < KEY_FIELD_WIDTH) {
        return null;
      } else {
        if (rest_of_line.charAt (KEY_FIELD_WIDTH - 1) != ' ') {
          throw new ReadFormatException ("column " +
                                         LOCATION_START_COLUMN + " must "
                                         + "be empty", line_number);
        }

        final String key_field = rest_of_line.substring (0,15);

        return key_field.trim ();
      }
    }
  }

  /**
   *  Read some embl feature qualifiers from a stream into a QualifierVector
   *  object.  The stream should contain qualifiers in this form:
   *  <PRE>  /name1=value1/name2="value2"/name3=[value3]  </PRE>
   *  @param in_stream the qualifiers are read from this stream
   *  @exception IOException thrown if there is a problem reading the
   *    qualifiers, such as end of file.
   *  @exception QualifierParseException Thrown if the format of the value
   *    String is not appropriate for a Qualifier with the given name.  Each
   *    qualifier has a specific format for the value part which depends on
   *    the name, for example the value part of /codon_start qualifier must be
   *    a number: 1, 2 or 3.
   *  @return A Vector containing one Qualifier object for each name/value
   *    pair read from the stream.
   **/
  public static QualifierVector
    readQualifiers (final Reader in_stream,
                    final EntryInformation entry_information)
      throws QualifierParseException, IOException {

    QualifierVector return_vector = new QualifierVector ();

    BufferedReader buffered_reader = new BufferedReader (in_stream);

    String name;
    String value;

    // loop until end of file
    while (true) {

      name = StreamQualifier.readName (buffered_reader);

      if (name == null) {
        // end of file/stream
        break;
      }

      // save one character in case the next char is not a '='
      buffered_reader.mark (1);

      final int next_char = buffered_reader.read ();

      if (next_char == -1) {
        value = null;
      } else {
        if (next_char == '=') {
          value = StreamQualifier.readValue (buffered_reader);
        } else {
          // this qualifier doesn't have a value
          value = null;
          buffered_reader.reset ();
        }
      }

      final Qualifier new_qualifier;

      if (value == null) {
        new_qualifier = new Qualifier (name);
      } else {
        new_qualifier =
          StreamQualifier.makeStreamQualifier (name, value,
                                               entry_information);
      }

      return_vector.addQualifierValues (new_qualifier);
    }

    return return_vector;
  }

  /**
   *  Return a QualifierVector containing the qualifiers from a String.
   *  @param qual_string contains the qualifiers to parse
   */
  private static QualifierVector
    getQualifiersFromString (final String qual_string,
                             final EntryInformation entry_information)
      throws QualifierParseException {

    final StringReader string_reader = new StringReader (qual_string);

    final QualifierVector qualifiers;

    try {
      qualifiers = readQualifiers (string_reader, entry_information);
    } catch (IOException exception) {
      throw (new QualifierParseException (exception.getMessage ()));
    }

    string_reader.close ();

    return qualifiers;
  }

  /**
   *  The column of the output where we should start writting the key.
   **/
  private final static int KEY_START_COLUMN = 5;

  /**
   *  The column of the output where we should start writting the location.
   **/
  private final static int LOCATION_START_COLUMN = 21;

  /**
   *  We should wrap at this column - no characters will be put in or after
   *  this column
   **/
  private final static int LAST_COLUMN = 81;

  /**
   *  Write the key of this Feature to the given stream.
   **/
  private void writeKey (Writer writer)
      throws IOException {
    if (this instanceof EmblStreamFeature) {
      writer.write ("FT   " + getKey ());
    } else {
      writer.write ("     " + getKey ());
    }

    final StringBuffer spaces = new StringBuffer (20);

    // add spaces so that the location starts at coloumn 21
    for (int i = 0 ;
         i < LOCATION_START_COLUMN - KEY_START_COLUMN - getKey ().length () ;
         ++i) {
      spaces.append (' ');
    }

    writer.write (spaces.toString ());
  }


  /**
   *  Write the location of this feature to a stream.  It is written in the
   *  usual EMBL format.  Line that are more than 79 characters wide are
   *  wrapped.  The wrapped lines start with "FT                   ", the
   *  first line doesn't.
   *  @param writer The stream to write to.
   *  @exception IOException thrown if there is an io problem while writing
   *    the Feature.
   **/
  private void writeLocation (final Writer writer)
      throws IOException {

    if (getEntryInformation ().useEMBLFormat ()) {
      String location_string = getLocation ().toString ();

      // write the first line without the "FT ..." bit as this will already
      // have been written by writeKey ()

      while (true) {
        final int wrap_position =
          getWrapPosition (location_string,
                           ',',
                           LAST_COLUMN - LOCATION_START_COLUMN - 1);
        if (wrap_position == -1) {
          // write out the remainer of the location_string
          writer.write (location_string);
          writer.write ("\n");
          break;
        } else {
          // write the start of the string (the bit before wrap_position) and
          // save the rest for the next time around the loop

          writer.write (location_string.substring (0, wrap_position + 1));
          if (this instanceof EmblStreamFeature) {
            writer.write ("\nFT                   ");
          } else {
            writer.write ("\n                     ");
          }

          location_string = location_string.substring (wrap_position + 1);
        }
      }
    } else {
      String location_string = getLocation ().toStringShort ();

      // write it out on one line
      writer.write (location_string);
      writer.write ("\n");
    }
  }

  /**
   *  The value to pass as the line_length argument of getWrapPosition ().
   **/
  final static private int QUALIFIER_WRAP_LENGTH =
    LAST_COLUMN - LOCATION_START_COLUMN - 1;

  /**
   *  Write the qualifiers of this feature to a stream.  The qualifiers are
   *  written in EMBL format: ie. <p>
   *  FT                   /codon_start=1
   *  <p> etc.
   *  @param writer The stream to write to.
   *  @exception IOException thrown if there is an io problem while writing
   *    the Feature.
   **/
  private void writeQualifiers (final Writer writer)
      throws IOException {
    for (int i = 0 ; i < getQualifiers ().size () ; ++i) {
      final Qualifier current_qualifier = (Qualifier)getQualifiers ().elementAt (i);

      //final String qualifier_name = current_qualifier.getName ();

      final EntryInformation entry_information = getEntryInformation ();
      
      final QualifierInfo qualifier_info =
        entry_information.getQualifierInfo (current_qualifier.getName ());

      // this will contain one String for each /name=value pair
      final StringVector qualifier_strings =
        StreamQualifier.toStringVector (qualifier_info, current_qualifier);

      for (int value_index = 0 ;
           value_index < qualifier_strings.size () ;
           ++value_index) {

        String qualifier_string = (String)qualifier_strings.elementAt (value_index);

        while (true) {
          int wrap_position = getWrapPosition (qualifier_string, ' ',
                                               QUALIFIER_WRAP_LENGTH);

          if (entry_information.useEMBLFormat () &&
              (wrap_position == -1 || wrap_position > QUALIFIER_WRAP_LENGTH)) {
            // EMBL entries cannot go over 80 columns
            if (qualifier_string.length () > QUALIFIER_WRAP_LENGTH) {
              wrap_position = QUALIFIER_WRAP_LENGTH;
            }
          }
 
          final String this_string;

          if (wrap_position == -1) {
            this_string = qualifier_string; 
          } else {
            // write the start of the string (the bit before wrap_position -
            // ignoring the space that we wrapped on)
            this_string = qualifier_string.substring (0, wrap_position);
          }

          if (this instanceof EmblStreamFeature) {
            writer.write ("FT                   ");
            writer.write (this_string);
            writer.write ("\n");
          } else {
            writer.write ("                     ");
            writer.write (this_string);
            writer.write ("\n");
          }

          if (wrap_position == -1) {
            break;
          }

          // don't write the spaces that start the next line
          while (wrap_position < qualifier_string.length () &&
                 qualifier_string.charAt (wrap_position) == ' ') {
            wrap_position += 1;
          }

          //save the rest for the next time around the loop
          qualifier_string = qualifier_string.substring (wrap_position);

          if (qualifier_string.length () == 0) {
            break;
          }
        }
      }
    }
  }

  /**
   *  Find a suitable character position to wrap at.
   *  (Helper method for writeLocation () and writeQualifiers ().)
   *  @param string_to_wrap This is the string to search for the wrap
   *    character.
   *  @param wrap_character This is the character after which we should wrap
   *    the line.
   *  @param line_length The maximum line length before wrapping.
   *  @return The character position in the string_to_wrap of wrap_character
   *    or -1 if there is no need to wrap. 
   **/
  private int getWrapPosition (String string_to_wrap,
                               char wrap_character,
                               int line_length) {
    //final int return_index;

    if (string_to_wrap.length () < line_length) {
      return -1;                // no need to wrap
    }

    final int wrap_char_index =
      string_to_wrap.substring (0, line_length).lastIndexOf (wrap_character);

    if (wrap_char_index == -1) {
      // there is no wrap_character within the line_length characters - so
      // just return the first occurrence in the string (or -1).
      return string_to_wrap.indexOf (wrap_character);
    } else {
      return wrap_char_index;
    }
  }
}
