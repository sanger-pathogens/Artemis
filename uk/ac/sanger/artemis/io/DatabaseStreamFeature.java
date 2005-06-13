/* DatbaseStreamFeature.java
 *
 * created: Mar 2005
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2005  Genome Research Limited
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
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.*;
import java.io.*;

/**
 *  This is an implementation of Feature that can read and write itself to a
 *  CHADO stream.
 *
 *  @version $Id: DatabaseStreamFeature.java,v 1.2 2005-06-13 18:36:49 tjc Exp $
 **/

public class DatabaseStreamFeature
    extends GFFStreamFeature
    implements DocumentFeature, StreamFeature, ComparableFeature
{
  /**
   *  Create a new DatabaseStreamFeature object.  The feature should be added
   *  to an Entry (with Entry.add ()).
   *  @param key The new feature key
   *  @param location The Location object for the new feature
   *  @param qualifiers The qualifiers for the new feature
   **/
  public DatabaseStreamFeature(final Key key,
                               final Location location,
                               final QualifierVector qualifiers)
  {
    super (null);
    try 
    {
      setKey (key);
      setLocation (location);
      setQualifiers (qualifiers);
    } 
    catch (EntryInformationException e) 
    {
      // this should never happen because the feature will not be in an Entry
      throw new Error ("internal error - unexpected exception: " + e);
    } 
    catch (ReadOnlyException e) 
    {
      // this should never happen because the feature will not be in an Entry
      throw new Error ("internal error - unexpected exception: " + e);
    } 
    catch (OutOfRangeException e) 
    {
      // this should never happen because the feature will not be in an Entry
      throw new Error ("internal error - unexpected exception: " + e);
    }
  }

  /**
   *  Create a new DatabaseStreamFeature with the same key, location and
   *  qualifiers as the given feature.  The feature should be added to an
   *  Entry (with Entry.add ()).
   *  @param feature The feature to copy.
   **/
  public DatabaseStreamFeature (final Feature feature) 
  {
    super (null);
    try 
    {
      setKey (feature.getKey ());
      setLocation (feature.getLocation ());
      setQualifiers (feature.getQualifiers ());
    } 
    catch (EntryInformationException e) 
    {
      throw new Error ("internal error - unexpected exception: " + e);
    } 
    catch (ReadOnlyException e) 
    {
      throw new Error ("internal error - unexpected exception: " + e);
    } 
    catch(OutOfRangeException e) 
    {
      throw new Error ("internal error - unexpected exception: " + e);
    }
  }


  public DatabaseStreamFeature(String line)
  {
    super(null);
    try
    {
      setFromString(line);
    }
    catch(IOException ioe)
    {
      throw new Error ("internal error - unexpected exception: " + ioe);
    }
  }

  /**
   *  Return the reference of a new copy of this Feature.
   **/
  public Feature copy() 
  {
    final Feature return_value = new DatabaseStreamFeature(this);
    return return_value;
  }

  /**
   *  Write this Feature to the given stream.
   *  @param writer The stream to write to.
   *  @exception IOException thrown if there is an io problem while writing
   *    the Feature.
   **/
  public synchronized void writeToStream (final Writer writer)
      throws IOException 
  {
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
//protected static DatabaseStreamFeature
//  readFromStream(LinePushBackReader stream)
//    throws IOException
//{
//  String line = stream.readLine ();
//  return readFromStream(line);
//}

//protected static DatabaseStreamFeature
//  readFromStream(String line)
//    throws IOException
//{
//  if(line == null)
//    return null;

//  return new DatabaseStreamFeature(line);
//}

  /**
   *  This is used by readFromStream() as temporary storage.  It is a class
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
  private void setFromStream(final LinePushBackReader in_stream)
      throws IOException
  {
    final String first_line = in_stream.readLine();
    setFromString(first_line);
  }

  /**
   *  Read the details of a feature from an EMBL stream into the current
   *  object.  (Called only by readFromStream ()).
   *  @param in_stream the Feature is read from this stream
   *  @exception IOException thrown if there is a problem reading the Feature -
   *    most likely ReadFormatException.
   **/
  private void setFromString(String first_line)
      throws IOException 
  {
    if(first_line == null) 
      throw new EOFException ("while reading a feature");

    final String key_string = getKeyStringFromLine(first_line);

    if(key_string == null)
      throw new ReadFormatException("expected the first line of a " +
                                    "feature");

    final Key key = new Key(key_string);
    final Location location;

    int ind1 = first_line.indexOf(" ")+1;
    int ind2 = first_line.indexOf(";");

    try 
    {
      location = new Location(first_line.substring(ind1,ind2));
    } 
    catch(LocationParseException exception) 
    {
      // re-throw the exception with the line number added

      final String new_error_string = exception.getMessage();

      // subtract 1 because the error was on the previous line
      throw new ReadFormatException(new_error_string);
    }


    final QualifierVector qualifiers;
    final String qualifier_string = first_line.substring(ind2+1);

    try 
    {
      qualifiers = getQualifiersFromString(qualifier_string,
                                           getEntryInformation());
    } 
    catch (QualifierParseException exception) 
    {
      // re-throw the exception with the line number added
      final String new_error_string = exception.getMessage ();

      // subtract 1 because the error was on the previous line
      throw new ReadFormatException(new_error_string);
    }

    try 
    {
      set(key, location, qualifiers);
    }
    catch(EntryInformationException e) 
    {
      throw new Error("internal error - unexpected exception: " + e);
    } 
    catch (OutOfRangeException e)
    {
      throw new Error("internal error - unexpected exception: " + e);
    }

    setDirtyFlag();
  }


  /**
   *  Return the key from a embl entry line.
   *  @param line_string the text of the entry line to process
   *  @return null if this isn't the first line of a feature, otherwise the
   *    key of this feature
   */
  private static String getKeyStringFromLine(final String line_string)
      throws ReadFormatException 
  {
    if(line_string == null ||
       line_string.startsWith (" ")) 
      return null;
    else
    {
      int ind1 = line_string.indexOf("=")+1;
      int ind2 = line_string.indexOf(" ");
      final String key_field = line_string.substring(ind1,ind2);
      return key_field.trim();
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
    readQualifiers(final Reader in_stream,
                   final EntryInformation entry_information)
      throws QualifierParseException, IOException 
  {

    QualifierVector return_vector = new QualifierVector ();

    BufferedReader buffered_reader = new BufferedReader (in_stream);

    String name;
    String value;

    // loop until end of file
    while( (name = StreamQualifier.readName (buffered_reader)) != null)
    {
      // save one character in case the next char is not a '='
      buffered_reader.mark(1);

      final int next_char = buffered_reader.read();

      if(next_char == -1) 
        value = null;
      else
      {
        if (next_char == '=') 
          value = StreamQualifier.readValue(buffered_reader);
        else 
        {
          // this qualifier doesn't have a value
          value = null;
          buffered_reader.reset();
        }
      }

      final Qualifier new_qualifier;

      if(value == null)
        new_qualifier = new Qualifier(name);
      else
      {
        new_qualifier =
          StreamQualifier.makeStreamQualifier(name, value,
                                              entry_information);
      }

      return_vector.addQualifierValues(new_qualifier);
    }

    return return_vector;
  }

  /**
   *  Return a QualifierVector containing the qualifiers from a String.
   *  @param qual_string contains the qualifiers to parse
   */
  private static QualifierVector
    getQualifiersFromString(final String qual_string,
                            final EntryInformation entry_information)
      throws QualifierParseException 
  {
    final StringReader string_reader = new StringReader(qual_string);
    final QualifierVector qualifiers;

    try 
    {
      qualifiers = readQualifiers (string_reader, entry_information);
    } 
    catch (IOException exception) 
    {
      throw (new QualifierParseException (exception.getMessage ()));
    }

    string_reader.close ();

    return qualifiers;
  }

  /**
   *  Write the key of this Feature to the given stream.
   **/
  private void writeKey (Writer writer)
      throws IOException 
  {
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
      throws IOException 
  {
  }

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
      throws IOException 
  {
  }

}
