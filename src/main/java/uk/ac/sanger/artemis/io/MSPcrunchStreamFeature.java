/* MSPcrunchStreamFeature.java
 *
 * created: Sat Apr 15 2000
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/MSPcrunchStreamFeature.java,v 1.7 2009-04-20 14:41:35 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.*;

import java.io.*;

/**
 *  A StreamFeature that thinks it is a MSPcrunch feature.
 *
 *  @author Kim Rutherford
 *  @version $Id: MSPcrunchStreamFeature.java,v 1.7 2009-04-20 14:41:35 tjc Exp $
 **/

public class MSPcrunchStreamFeature
    extends SimpleDocumentFeature
    implements DocumentFeature, StreamFeature, ComparableFeature {
  
  private static Key KEY_X = new Key ("CRUNCH_X");
  private static Key KEY_D = new Key ("CRUNCH_D");

  /**
   *  Create a new MSPcrunchStreamFeature object.  The feature should be added
   *  to an Entry (with Entry.add ()).
   *  @param key The new feature key
   *  @param location The Location object for the new feature
   *  @param qualifiers The qualifiers for the new feature
   **/
  public MSPcrunchStreamFeature (final Key key,
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
   *  Create a new MSPcrunchStreamFeature with the same key, location and
   *  qualifiers as the given feature.  The feature should be added to an
   *  Entry (with Entry.add ()).
   *  @param feature The feature to copy.
   **/
  public MSPcrunchStreamFeature (final Feature feature) {
    super (null);

    if (feature instanceof MSPcrunchStreamFeature) {
      mspcrunch_line = ((MSPcrunchStreamFeature)feature).mspcrunch_line;
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
    final Feature return_value = new MSPcrunchStreamFeature (this);

    return return_value;
  }

  /**
   *  Create a new MSPcrunchStreamFeature from the given line.  The String
   *  should be in gene finder format.
   **/
  private MSPcrunchStreamFeature (final String line)
      throws ReadFormatException {
    super (null);

//  final StringVector line_bits = StringVector.getStrings (line, " ");

    int index = 0;
    int lastIndex = 0;
    int count = 1;

    while( (index = line.indexOf(" ", index)) > -1)
    {
      count++;
      index++;
    }

    if(count > 8)
      count = 9;
    else if(count < 7) 
      throw new ReadFormatException ("invalid MSPcrunch line (not enough " +
                                     "fields): " + line);

    final String line_bits[] = new String[count];
    count = 0;
    index = 0;

    while( (index = line.indexOf(" ", index)) > -1 &&
           count < line_bits.length-1)
    {
      line_bits[count] = line.substring(lastIndex, index).trim();
      if(line_bits[count].equals(""))
        count--;

      count++;
      index++;
      lastIndex = index;  
    }
    if(lastIndex < line.length())
    {
      line_bits[count] = line.substring(lastIndex);
      count++;
    }

    try
    {
      int query_start = Integer.valueOf (line_bits[2]).intValue ();
      int query_end = Integer.valueOf (line_bits[3]).intValue ();

      final boolean crunch_x;
      final boolean complement_flag;
      
      if (line_bits[1].startsWith("(+")) 
      {
        crunch_x = true;
        complement_flag = false;
      } 
      else 
      {
        if (line_bits[1].startsWith("(-")) 
        {
          crunch_x = true;
          complement_flag = true;
        } 
        else 
        {
          if (line_bits[1].charAt (0) == '.' ||
              Character.isDigit (line_bits[1].charAt (0))) 
          {
            crunch_x = false;
            if (query_start > query_end) 
              complement_flag = true;
            else 
              complement_flag = false;
            
          } 
          
          else 
          {
            final String message =
              "invalid MSPcrunch line - column 2 should be a " +
              "number, (-1) or (+1): " + line;
            throw new ReadFormatException (message);
          }
        }
      }

      
      if(query_start > query_end) 
      {
        final int tmp = query_end;
        query_end = query_start;
        query_start = tmp;
      }

      if(crunch_x) 
        line_bits[1] = null;

      final Qualifier query_id_qualifier;

      final String subject_start;
      final String subject_end;
      final String subject_id;
      final String description;

      if(crunch_x) 
      {
        query_id_qualifier = new Qualifier ("query_id", "unknown");
        subject_start = line_bits[4];
        subject_end = line_bits[5];
        subject_id = line_bits[6];
        final StringBuffer desc_buffer = new StringBuffer ();

        for (int i = 7 ; i < line_bits.length ; ++i) 
        {
          desc_buffer.append (line_bits[i]);
          if(i < line_bits.length - 1)
            desc_buffer.append (" ");
        }
        description = desc_buffer.toString ();
        
        setKey(KEY_X);
      } 
      else 
      {
        query_id_qualifier = new Qualifier ("query_id", line_bits[4]);
        subject_start = line_bits[5];
        subject_end = line_bits[6];
        subject_id  = line_bits[7];
        description = line_bits[8];
        
        setKey(KEY_D);
      }
      
      final Qualifier blast_score_qualifier =
        new Qualifier ("blast_score", line_bits[0]);

      final Qualifier subject_start_qualifier =
        new Qualifier ("subject_start", subject_start);

      final Qualifier subject_end_qualifier =
        new Qualifier ("subject_end", subject_end);

      final Qualifier subject_id_qualifier =
        new Qualifier ("subject_id", subject_id);

      setQualifier (blast_score_qualifier);

      if(line_bits[1] != null) 
      {
        // score qualifier must be 1-100
        final Qualifier score_qualifier =
          new Qualifier ("score", line_bits[1]);
        
        final Qualifier percent_id_qualifier =
          new Qualifier ("percent_id", line_bits[1]);
        
        setQualifier (score_qualifier);
        setQualifier (percent_id_qualifier);
      }

      setQualifier (query_id_qualifier);
      setQualifier (subject_start_qualifier);
      setQualifier (subject_end_qualifier);
      setQualifier (subject_id_qualifier);

//    final StringVector note_values = new StringVector ();
//    note_values.add (description);
//
// This is more memory intensive than just using the description:
//    note_values.add ("hit to " + subject_id + " " + subject_start +
//                     ".." + subject_end + "  score: " + line_bits[0] +
//                     (line_bits[1] == null ?
//                      "" :
//                      "  percent id: " + line_bits[1]) +
//                     "  " + description);

      final Qualifier note_qualifier = new Qualifier ("note", description);
      setQualifier (note_qualifier);

      final RangeVector ranges =
        new RangeVector (new Range (query_start, query_end));

      setLocation (new Location (ranges, complement_flag));
    } catch (ReadOnlyException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    } catch (EntryInformationException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    } catch (OutOfRangeException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    } catch (LocationParseException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    }

    this.mspcrunch_line = line;
  }

  /**
   *  Read and return a MSPcrunchStreamFeature from a stream.  A feature must
   *  be the next thing in the stream.
   *  @param stream the Feature is read from this stream
   *  @exception IOException thrown if there is a problem reading the Feature -
   *    most likely ReadFormatException.
   *  @exception InvalidRelationException Thrown if this Feature cannot contain
   *    the given Qualifier.
   *  @return null if in_stream is at the end of file when the method is
   *    called
   **/
  protected static MSPcrunchStreamFeature
    readFromStream (LinePushBackReader stream)
      throws IOException, InvalidRelationException {

    String line = stream.readLine ();

    if (line == null) {
      return null;
    }

    try {
      final MSPcrunchStreamFeature new_feature =
        new MSPcrunchStreamFeature (line);

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
   *    most likely ReadFormatException if the stream does not contain MSPcrunch
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

    // for now MSPcrunch features are read-only so just write what we read
    writer.write (mspcrunch_line + "\n");

    /*

    final RangeVector ranges = getLocation ().getRanges ();

    for (int i = 0 ; i < ranges.size () ; ++i) {
      final Range this_range        = ranges.elementAt (i);
      Qualifier   source            = getQualifierByName ("mspcrunch_source");
      Qualifier   MSPcrunch_feature = getQualifierByName ("mspcrunch_feature");
      Qualifier   score             = getQualifierByName ("score");
      Qualifier   MSPcrunch_group   = getQualifierByName ("mspcrunch_group");
      Qualifier   note              = getQualifierByName ("note");

      if (source == null) {
        source = new Qualifier ("source", "");
      }

      if (MSPcrunch_feature == null) {
        MSPcrunch_feature = new Qualifier ("MSPcrunch_feature", "");
      }

      if (score == null) {
        score = new Qualifier ("score", "");
      }

      if (MSPcrunch_group == null || MSPcrunch_group.getValues () == null ||
          MSPcrunch_group.getValues ().elementAt (0).equals ("")) {
        final Qualifier gene = getQualifierByName ("gene");

        if (gene == null) {
          MSPcrunch_group = new Qualifier ("MSPcrunch_group", "");
        } else {
          MSPcrunch_group = gene;
        }
      }

      if (note == null) {
        note = new Qualifier ("note", "");
      }

      String frame = ".";

      final Qualifier codon_start = getQualifierByName ("codon_start");

      if (codon_start != null && i == 0) {
        frame = codon_start.getValues ().elementAt (0);
      }

      writer.write (getKey () + "\t" +
                    source.getValues ().elementAt (0) + "\t" +
                    MSPcrunch_feature.getValues ().elementAt (0) + "\t" +
                    this_range.getStart () + "\t" +
                    this_range.getEnd () + "\t" +
                    score.getValues () .elementAt (0)+ "\t" +
                    (getLocation ().isComplement () ? "\t-\t" : "\t+\t") +
                    frame + "\t" +
                    MSPcrunch_group.getValues ().elementAt (0) + "\t" +
                    note.getValues ().elementAt (0) + "\n");
    }
    */

  }

  /**
   *  The DocumentEntry object that contains this Feature as passed to the
   *  constructor.
   **/
  private DocumentEntry entry;

  /**
   *  This is the line of MSPcrunch -d input that was read to get this
   *  MSPcrunchStreamFeature.
   **/
  private String mspcrunch_line = null;
}
