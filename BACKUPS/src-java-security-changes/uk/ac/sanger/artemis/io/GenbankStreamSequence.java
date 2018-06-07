/* GenbankStreamSequence.java
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/GenbankStreamSequence.java,v 1.6 2008-03-25 13:49:30 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.LinePushBackReader;

import java.io.IOException;
import java.io.Writer;

/**
 *  This is a subclass of StreamSequence containing GENBANK format sequence.
 *
 *  @author Kim Rutherford
 *  @version $Id: GenbankStreamSequence.java,v 1.6 2008-03-25 13:49:30 tjc Exp $
 **/

public class GenbankStreamSequence extends StreamSequence
{

  /** The header line(s) of this sequence(set by setHeader()). */
  private String header_line = null;

  /**
   *  Create a new GenbankStreamSequence object from a stream.  The next line
   *  to read from the stream should be the header line of the sequence.
   *  @param in_stream The stream to read from.  When the constructor returns
   *    the stream will at the next line after the sequence.
   **/
  public GenbankStreamSequence(final LinePushBackReader in_stream)
      throws IOException 
  {
    readHeader(in_stream);
    readSequence(in_stream);
  }

  /**
   *  Make a new GenbankStreamSequence containing the same sequence as the
   *  given Sequence.
   **/
  public GenbankStreamSequence(final Sequence sequence) 
  {
    setFromChar(((StreamSequence)sequence).getCharSequence());
  }

  /**
   *  Make a new GenbankStreamSequence containing the same sequence as the
   *  given String.
   *  @param sequence_string The String containing the sequence for the new
   *    GenbankStreamSequence.
   **/
  public GenbankStreamSequence(final String sequence_string) 
  {
    setFromChar(sequence_string.toCharArray());
  }

  /**
   *  Return a new StreamSequence object that is a copy of this one.
   **/
  public StreamSequence copy() 
  {
    return new GenbankStreamSequence(this);
  }

  /**
   *  Return the sequence type(GENBANK_FORMAT for this class).
   **/
  public int getFormatType() 
  {
    return StreamSequenceFactory.GENBANK_FORMAT;
  }

  /**
   *  Read the header for this sequence(if any).
   **/
  protected void readHeader(final LinePushBackReader in_stream)
      throws IOException 
  {
    String sequence_header = in_stream.readLine();

    // GenBank format has changed - BASE COUNT no longer appears
    if(sequence_header.startsWith("BASE COUNT")) 
    {
      final String next_line = in_stream.readLine();

      if(next_line.startsWith("ORIGIN")) 
        sequence_header = sequence_header + next_line;
      else 
        throw new ReadFormatException("Genbank sequence data should have " +
                                       "\"ORIGIN\" on the second line");
    }
    else 
    {
      if(!sequence_header.startsWith("ORIGIN"))
        throw new ReadFormatException("Genbank sequence data should begin " +
                                       "with \"BASE COUNT\" or \"ORIGIN\"");
    }

    setHeader(sequence_header);
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
   *  This method will read raw sequence from the stream in this object
   *  removing whitespace as it goes.  No checks are made on the format of the
   *  sequence, apart from checking that the stream contains only letters.
   **/
  protected void readSequence(final LinePushBackReader in_stream)
      throws IOException 
  {

    final int buffer_capacity = 50000;

    // we buffer up the sequence bases then assign them to sequence once all
    // bases are read
    String line;
    setSequencePackingCapacity(buffer_capacity);
    final StringBuffer this_line_sequence_buffer = new StringBuffer(100);
    while((line = in_stream.readLine()) != null) 
    {
      if(line.equals("//")) 
      {
        // end of the sequence
        in_stream.pushBack(line);
        break;
      }
      
      // set to true when we see a base - numbers and spaces are allowed before
      // the first base is seen in the line
      boolean seen_base = false;
      this_line_sequence_buffer.setLength(0);

      for(int i = 0 ; i < line.length() ; ++i) 
      {
        final char this_char = line.charAt(i);
  
        if(!seen_base)  
        {
          if(Character.isDigit(this_char) ||
             Character.isSpaceChar(this_char)) 
          {
              continue;
            // ok - ignore it
          }
        }

        if(Character.isLetter(this_char) ||
           this_char == '.' ||
           this_char == '-' ||
           this_char == '*') 
        {
          seen_base = true;
          this_line_sequence_buffer.append(this_char);
        } 
        else
        {
          if(Character.isSpaceChar(this_char)) {
            // just ignore it
          } 
          else 
            throw new ReadFormatException("GENBANK sequence file contains " +
                                           "a character that is not a " +
                                           "letter: " + this_char,
                                           in_stream.getLineNumber());
        }
      }
      
      if(this_line_sequence_buffer.length() > 0)
        appendChar(this_line_sequence_buffer.toString().toLowerCase().toCharArray());
    }
    setCounts();
  }

  /**
   *  Write this Sequence to the given stream.
   *  @param writer The stream to write to.
   **/
  public synchronized void writeToStream(final Writer writer)
      throws IOException 
  {
    final StringBuffer line_buffer = new StringBuffer(90);

    // first count A,C,G,T and other bases

    final int SEQUENCE_LINE_BASE_COUNT = 60;

    line_buffer.setLength(0);
    line_buffer.ensureCapacity(90);

    line_buffer.append("BASE COUNT  ");

    final int MAX_WIDTH = 7;

    appendAndPad(line_buffer, MAX_WIDTH, String.valueOf(getACount()));
    line_buffer.append(" a");

    appendAndPad(line_buffer, MAX_WIDTH, String.valueOf(getCCount()));
    line_buffer.append(" c");

    appendAndPad(line_buffer, MAX_WIDTH, String.valueOf(getGCount()));
    line_buffer.append(" g");

    appendAndPad(line_buffer, MAX_WIDTH, String.valueOf(getTCount()));
    line_buffer.append(" t");
  
    writer.write(line_buffer + "\nORIGIN\n");

    for(int i = 0 ; i < length() ; i += SEQUENCE_LINE_BASE_COUNT) 
    {
      // get the bases in chunks of at most 60
      final int this_line_length;

      if(length() - i < SEQUENCE_LINE_BASE_COUNT)
        this_line_length = length() - i;
      else 
        this_line_length = SEQUENCE_LINE_BASE_COUNT;

      line_buffer.setLength(0);
      line_buffer.ensureCapacity(90);

      // the base counter to write at the end of each line
      final int base_count = i + 1;

      final String string_base_count = String.valueOf(base_count);

      // now pad the line with spaces
      final int LINE_WIDTH = 80;

      final int MAX_COUNT_WIDTH = 9;

      appendAndPad(line_buffer, MAX_COUNT_WIDTH, string_base_count);

      final int BLOCK_LENGTH = 10;
      final int BLOCK_COUNT  = 6;

      for(int j = 0 ; j < this_line_length ; j += BLOCK_LENGTH) 
      {
        final int this_block_length;

        line_buffer.append(' ');

        if(this_line_length - j < BLOCK_LENGTH) 
          this_block_length = this_line_length - j;
        else 
          this_block_length = BLOCK_LENGTH;

        line_buffer.append(getCharSubSequence(i + j + 1, i + j + this_block_length));

//sequence.substring(i + j,
//                                            i + j + this_block_length));
      }
   
      line_buffer.append("\n");
      writer.write(line_buffer.toString());
    }
  }

  /**
   *  Append string to line_buffer, but pad it (in front) with spaces so that
   *  it takes max_width characters.
   **/
  private void appendAndPad(final StringBuffer line_buffer,
                             final int max_width,
                             final String string) 
  {
    final int string_width = string.length();

    for(int char_index = 0; char_index < max_width - string_width;
        ++char_index) 
      line_buffer.append(' ');

    line_buffer.append(string);
  }

}

