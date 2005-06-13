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
 *  @version $Id: DatabaseStreamFeature.java,v 1.3 2005-06-13 18:51:26 tjc Exp $
 **/

public class DatabaseStreamFeature
    extends GFFStreamFeature
    implements DocumentFeature, StreamFeature, ComparableFeature
{

  /**
   *  Create a new DatabaseStreamFeature with the same key, location and
   *  qualifiers as the given feature.  The feature should be added to an
   *  Entry (with Entry.add ()).
   *  @param feature The feature to copy.
   **/
  public DatabaseStreamFeature (final Feature feature) 
  {
    super(feature);
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
   *  This is used by readFromStream() as temporary storage.  It is a class
   *  member rather than a local variable so that we don't need to allocate a
   *  object for each call.  The number we pick for the initial StringBuffer
   *  size is not critical, but should cover most possibilities
   **/
  final static private StringBuffer qualifier_string_buffer =
    new StringBuffer (1500);

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
