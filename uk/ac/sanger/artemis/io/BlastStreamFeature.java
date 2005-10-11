/* BlastStreamFeature.java
 *
 * created: Wed May  8 2002
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/BlastStreamFeature.java,v 1.2 2005-10-11 14:20:31 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.*;

import java.io.*;

/**
 *  A StreamFeature that thinks it is a Blast feature.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: BlastStreamFeature.java,v 1.2 2005-10-11 14:20:31 tjc Exp $
 **/

public class BlastStreamFeature
    extends SimpleDocumentFeature
    implements DocumentFeature, StreamFeature, ComparableFeature {
  /**
   *  Create a new BlastStreamFeature object.  The feature should be added
   *  to an Entry (with Entry.add ()).
   *  @param key The new feature key
   *  @param location The Location object for the new feature
   *  @param qualifiers The qualifiers for the new feature
   **/
  public BlastStreamFeature (final Key key,
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
   *  Create a new BlastStreamFeature with the same key, location and
   *  qualifiers as the given feature.  The feature should be added to an
   *  Entry (with Entry.add ()).
   *  @param feature The feature to copy.
   **/
  public BlastStreamFeature (final Feature feature) {
    super (null);

    if (feature instanceof BlastStreamFeature) {
      blast_line = ((BlastStreamFeature)feature).blast_line;
    }

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
   *  Return the reference of a new copy of this Feature.
   **/
  public Feature copy () {
    final Feature return_value = new BlastStreamFeature (this);

    return return_value;
  }

  /**
   *  Create a new BlastStreamFeature from the given line.  The String
   *  should be in gene finder format.
   **/
  private BlastStreamFeature (final String line)
      throws ReadFormatException {
    super (null);

    final StringVector line_bits = StringVector.getStrings (line, "\t");

    if (line_bits.size () < 12) {
      throw new ReadFormatException ("invalid Blast line (not enough " +
                                     "fields): " + line);
    }

    try {
      int query_start = Integer.valueOf ((String)line_bits.elementAt (6)).intValue ();
      int query_end = Integer.valueOf ((String)line_bits.elementAt (7)).intValue ();

      final String percent_id = (String)line_bits.elementAt (2);

      final String query_id = (String)line_bits.elementAt (0);
      final String subject_id = (String)line_bits.elementAt (1);
      final String subject_start_string = (String)line_bits.elementAt (8);
      final String subject_end_string = (String)line_bits.elementAt (9);

      final String score = (String)line_bits.elementAt (11);
      final String e_value = (String)line_bits.elementAt (10);

      final Qualifier blast_score_qualifier =
        new Qualifier ("blast_score", score);
      // score qualifier must be 1-100
      final Qualifier score_qualifier =
        new Qualifier ("score", percent_id);
      final Qualifier percent_id_qualifier =
        new Qualifier ("percent_id", percent_id);
      final Qualifier query_id_qualifier =
        new Qualifier ("query_id", query_id);
      final Qualifier subject_start_qualifier =
        new Qualifier ("subject_start", subject_start_string);
      final Qualifier subject_end_qualifier =
        new Qualifier ("subject_end", subject_end_string);
      final Qualifier subject_id_qualifier =
        new Qualifier ("subject_id", subject_id);

      setQualifier (blast_score_qualifier);
      setQualifier (score_qualifier);
      setQualifier (percent_id_qualifier);
      setQualifier (query_id_qualifier);
      setQualifier (subject_start_qualifier);
      setQualifier (subject_end_qualifier);
      setQualifier (subject_id_qualifier);

      int subject_start = Integer.valueOf (subject_start_string).intValue ();
      int subject_end = Integer.valueOf (subject_end_string).intValue ();

      final Key key = new Key ("BLASTCDS");

      setKey (key);

      final StringVector note_values = new StringVector ();

      note_values.add ("hit to " + subject_id + " " + subject_start +
                       ".." + subject_end + "  score: " + score +
                        "  percent id: " + percent_id + "  e-value: " +
                       e_value);

      final Qualifier note_qualifier = new Qualifier ("note", note_values);

      setQualifier (note_qualifier);


      boolean complement_flag;

      if (subject_end < subject_start) {
        complement_flag = true;
      } else {
        complement_flag = false;
      }

      if (query_start > query_end) {
        final int tmp = query_end;
        query_end = query_start;
        query_start = tmp;

        complement_flag = !complement_flag;
      }

      final RangeVector ranges =
        new RangeVector (new Range (query_start, query_end));

      setLocation (new Location (ranges, complement_flag));
    } catch (ReadOnlyException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    } catch (EntryInformationException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    } catch (OutOfRangeException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    }

    this.blast_line = line;
  }

  /**
   *  Read and return a BlastStreamFeature from a stream.  A feature must
   *  be the next thing in the stream.
   *  @param stream the Feature is read from this stream
   *  @exception IOException thrown if there is a problem reading the Feature -
   *    most likely ReadFormatException.
   *  @exception InvalidRelationException Thrown if this Feature cannot contain
   *    the given Qualifier.
   *  @return null if in_stream is at the end of file when the method is
   *    called
   **/
  protected static BlastStreamFeature
    readFromStream (LinePushBackReader stream)
      throws IOException, InvalidRelationException {

    String line = stream.readLine ();

    if (line == null) {
      return null;
    }

    try {
      final BlastStreamFeature new_feature =
        new BlastStreamFeature (line);

      return new_feature;
    } catch (ReadFormatException exception) {
      // re-throw the exception with the line number added

      final String new_error_string = exception.getMessage ();

      throw new ReadFormatException (new_error_string,
                                     stream.getLineNumber ());
    }
  }

  /**
   *  Read the details of a feature from an EMBL stream into the current
   *  object.
   *  @param entry_information The EntryInformation object of the Entry that
   *    will contain the Feature.
   *  @param in_stream the Feature is read from this stream
   *  @exception IOException thrown if there is a problem reading the Feature -
   *    most likely ReadFormatException if the stream does not contain Blast
   *    feature.
   **/
  public void setFromStream (final EntryInformation entry_information,
                             final LinePushBackReader in_stream)
      throws IOException, InvalidRelationException, ReadOnlyException {
    throw new ReadOnlyException ();
  }

  /**
   *  Write this Feature to the given stream.
   *  @param writer The stream to write to.
   *  @exception IOException thrown if there is an io problem while writing
   *    the Feature.
   **/
  public void writeToStream (final Writer writer)
      throws IOException {

    // for now Blast features are read-only so just write what we read
    writer.write (blast_line + "\n");
  }

  /**
   *  The DocumentEntry object that contains this Feature as passed to the
   *  constructor.
   **/
  private DocumentEntry entry;

  /**
   *  This is the line of blastall -m 8 input that was read to get this
   *  BlastStreamFeature.
   **/
  private String blast_line = null;
}
