/* EmblStreamSequence.java
 *
 * created: Mon Jun 14 1999
 *
 * This file is part of Artemis
 * 
 * Copyright(C) 1999  Genome Research Limited
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or(at your option) any later version.
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/EmblStreamSequence.java,v 1.1 2004-06-09 09:49:13 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.LinePushBackReader;

import java.io.IOException;
import java.io.Writer;

/**
 *  This is a subclass of StreamSequence containing EMBL sequence.
 *
 *  @author Kim Rutherford
 *  @version $Id: EmblStreamSequence.java,v 1.1 2004-06-09 09:49:13 tjc Exp $
 **/

public class EmblStreamSequence extends StreamSequence 
{

  /** The header line(s) of this sequence(set by setHeader()). */
  private String header_line = null;

  /**
   *  Create a new EmblStreamSequence object from a stream.  The next line to
   *  read from the stream should be the header line of the sequence.
   *  @param in_stream The stream to read from.  When the constructor returns
   *    the stream will at the next line after the sequence.
   **/
  public EmblStreamSequence(final LinePushBackReader in_stream)
      throws IOException 
  {
    readHeader(in_stream);
    readSequence(in_stream);
  }

  /**
   *  Make a new EmblStreamSequence containing the same sequence as the given
   *  Sequence.
   **/
  public EmblStreamSequence(final Sequence sequence) 
  {
    setFromString(sequence.toString());
  }

  /**
   *  Make a new EmblStreamSequence containing the same sequence as the given
   *  String. 
   *  @param sequence_string The String containing the sequence for the new
   *    EmblStreamSequence.
   **/
  public EmblStreamSequence(final String sequence_string)
  {
    setFromString(sequence_string);
  }

  /**
   *  Return a new StreamSequence object that is a copy of this one.
   **/
  public StreamSequence copy() 
  {
    return new EmblStreamSequence(this);
  }

  /**
   *  Return the sequence type(EMBL_FORMAT for this class).
   **/
  public int getFormatType() 
  {
    return StreamSequenceFactory.EMBL_FORMAT;
  }

  /**
   *  Read the header for this sequence(if any).
   **/
  protected void readHeader(final LinePushBackReader in_stream)
      throws IOException 
  {
    final String seq_header_line = in_stream.readLine();

    setHeader(seq_header_line);
  }

  /**
   *  Return the header line(s) of this Sequence or null if there is no
   *  header.  The returned String will not end in a newline character.
   **/
  public String getHeader() 
  {
    return header_line;
  }

  /**
   *  Set the header line(s) of this Sequence.  The argument should not end in
   *  a newline character.
   **/
  public void setHeader(final String sequence_string) 
  {
    header_line = sequence_string;
  }

  /**
   *  This method will read EMBL sequence from the stream in this object
   *  removing whitespace as it goes.
   **/
  protected void readSequence(final LinePushBackReader in_stream)
      throws IOException 
  {
    final StringBuffer this_line_sequence_buffer = new StringBuffer(81);
    final String seq_header_line = getHeader();
    final int header_base_count  = getHeaderBaseCount(seq_header_line);

    final int buffer_capacity;

    if(header_base_count > 50000) 
    {
      // add 100 for good luck
      buffer_capacity = header_base_count + 100;
    }
    else 
    {
      // try to pick an initial that covers a large number of cases
      buffer_capacity = 50000;
    }

    // we buffer up the sequence bases then assign them to sequence once all
    // bases are read
    StringBuffer sequence_buffer = new StringBuffer(buffer_capacity);

    String line;

    while((line = in_stream.readLine()) != null) 
    {
      if(line.startsWith("//") || !line.startsWith("     ")) 
      {
        // end of Sequence -
        in_stream.pushBack(line);
        break;
      }

      if(line.length() < 72) 
        throw new ReadFormatException("line too short while reading " +
                                       "embl sequence data",
                                       in_stream.getLineNumber());

      this_line_sequence_buffer.setLength(0);

      // if the input file is a well formatted embl entry, then the bases
      // should be in 6 columns of 10 bases each (with a single space
      // separating them). The first column starts at index 5 of the input
      // line.

      final char [] this_line_char_array = line.toCharArray();

      for(int i = 0 ; i < 6 ; ++i) 
      {
        final int first_base = 5 + i * 11;

        this_line_sequence_buffer.append(this_line_char_array,
                                         first_base,
                                         10);
      }

      // add the sequence from this line to the sequence we already have.

      String this_line_sequence =
        this_line_sequence_buffer.toString();

      if(this_line_char_array[59] == ' ') 
      {
        // call trim() because this line(the last line) ends in spaces.
        sequence_buffer.append(this_line_sequence.trim().toLowerCase());
      } 
      else 
        sequence_buffer.append(this_line_sequence.toLowerCase());
    }

    setFromString(sequence_buffer.toString());
  }

  /**
   *  Write this Sequence to the given stream.
   *  @param writer The stream to write to.
   **/
  public synchronized void writeToStream(final Writer writer)
      throws IOException 
  {

    final String sequence = toString();

    // first count A,C,G,T and other bases

    final int SEQUENCE_LINE_BASE_COUNT = 60;

    writer.write("SQ   Sequence " +
                  sequence.length() + " BP; " +
                  getACount() + " A; " +
                  getCCount() + " C; " +
                  getGCount() + " G; " +
                  getTCount() + " T; " +
                  getOtherCount() + " other;\n");

    int line_length_so_far = 0;

    for(int i = 0 ; i < sequence.length() ; i += SEQUENCE_LINE_BASE_COUNT) 
    {
      // get the bases in chunks of at most 60
      final int this_line_length;

      if(sequence.length() - i < SEQUENCE_LINE_BASE_COUNT) 
        this_line_length = sequence.length() - i;
      else 
        this_line_length = SEQUENCE_LINE_BASE_COUNT;

      writer.write("    ");
      line_length_so_far += 4;

      final int BLOCK_LENGTH = 10;
      final int BLOCK_COUNT  = 6;

      for(int j = 0 ; j < this_line_length ; j += BLOCK_LENGTH)
      {
        final int this_block_length;
        writer.write(' ');

        if(this_line_length - j < BLOCK_LENGTH) 
          this_block_length = this_line_length - j;
        else 
          this_block_length = BLOCK_LENGTH;

        writer.write(sequence.substring(i + j, i + j + this_block_length));

        line_length_so_far += this_block_length + 1;
      }

      // the base counter to write at the end of each line
      final int base_count = i + this_line_length;
      final String string_base_count = String.valueOf(base_count);

      // now pad the line with spaces
      final int count_width = string_base_count.length();
      final int LINE_WIDTH = 80;

      for(int char_index = 0 ;
           char_index < LINE_WIDTH - count_width - line_length_so_far ;
           ++char_index) 
        writer.write(' ');

      writer.write(string_base_count);
      line_length_so_far = 0;
      writer.write("\n");
    }
  }

  /**
   *  Take a EMBL sequence header line and extract and return the base count.
   *  If the line cannot be parsed as a header line then -1 is returned.
   **/
  private int getHeaderBaseCount(final String line) 
  {
    if(line.startsWith("SQ   Sequence ")) 
    {
      final String temp_line = line.substring(14);
      final int space_index = temp_line.indexOf(' ');

      if(space_index == -1) 
        return -1;

      final String count_string = temp_line.substring(0, space_index);

      try 
      {
        return Integer.parseInt(count_string);
      } 
      catch(NumberFormatException e) 
      {
        return -1;
      }
    }
    else 
      return -1;
  }

}

