/* RawStreamSequence.java
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/RawStreamSequence.java,v 1.2 2005-07-08 15:11:12 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.LinePushBackReader;

import java.io.IOException;
import java.io.Writer;
import java.util.Vector;

/**
 *  This is a subclass of StreamSequence containing raw sequence.
 *
 *  @author Kim Rutherford
 *  @version $Id: RawStreamSequence.java,v 1.2 2005-07-08 15:11:12 tjc Exp $
 **/

public class RawStreamSequence extends StreamSequence 
{

  /**
   *  The _character_ positions of the header lines relative to the sequence.
   *  This is in the same order as fasta_header_strings.
   **/
  private Vector fasta_header_positions = null;

  /**
   *  This Vector holds the fasta headers (if any) that where seen while
   *  reading the sequence.  This is in the same order as
   *  fasta_header_positions.
   **/
  private Vector fasta_header_strings = null;

  /**
   *  Create a new RawStreamSequence object from a stream containing raw
   *  bases.
   *  @param in_stream The stream to read from.  When the constructor returns
   *    the stream will at the next line after the sequence.
   **/
  public RawStreamSequence(final LinePushBackReader in_stream)
      throws IOException 
  {
    readSequence(in_stream);
  }

  /**
   *  Make a new RawStreamSequence containing the same sequence as the given
   *  Sequence.
   **/
  public RawStreamSequence(final Sequence sequence) 
  {
    this(sequence.toString());

    if(sequence instanceof RawStreamSequence) 
    {
      final RawStreamSequence raw_stream_sequence =
                                      (RawStreamSequence) sequence;
      this.fasta_header_positions =
        (Vector)raw_stream_sequence.fasta_header_positions.clone();
      this.fasta_header_strings =
        (Vector)raw_stream_sequence.fasta_header_strings.clone();
    }
  }

  /**
   *  Make a new RawStreamSequence containing the same sequence as the given
   *  String
   **/
  public RawStreamSequence(final String sequence_string) 
  {
    setFromChar(sequence_string.toCharArray());
  }

  /**
   *  Return a new StreamSequence object that is a copy of this one.
   **/
  public StreamSequence copy() 
  {
    return new RawStreamSequence(this);
  }

  /**
   *  Return the sequence type (RAW_FORMAT for this class).
   **/
  public int getFormatType() 
  {
    return StreamSequenceFactory.RAW_FORMAT;
  }

  /**
   *  Read the header for this sequence(if any).
   **/
  protected void readHeader(final LinePushBackReader in_stream)
      throws IOException 
  {
    // no header for raw sequence
  }

  /**
   *  This method will read raw sequence from the stream removing whitespace
   *  as it goes.  No checks are made on the format of the sequence, apart
   *  from checking that the stream contains only letters.
   **/
  protected void readSequence(final LinePushBackReader in_stream)
      throws IOException 
  {

    // initialise these here because this method is called before class
    // variables are initialised
    fasta_header_positions = new Vector();
    fasta_header_strings   = new Vector();

    final int buffer_capacity = 5000;

    // we buffer up the sequence bases then assign them to sequence once all
    // bases are read
    setSequencePackingCapacity(buffer_capacity);

    String line;
    int nbase = 0;
    while((line = in_stream.readLine()) !=null) 
    {
      if(line.startsWith(">")) 
      {
        final String header_string = line.substring(1);
        final Integer header_position = new Integer(nbase);

        fasta_header_strings.addElement(header_string);
        fasta_header_positions.addElement(header_position);

        // ignore header lines
        continue;
      }

      line = line.trim().toLowerCase();
//    System.out.println(line);
      appendChar(line.toCharArray());
      nbase += line.length();
    }

  }

  /**
   *  Write this Sequence to the given stream.
   *  @param writer The stream to write to.
   **/
  public void writeToStream(final Writer writer)
      throws IOException 
  {
    final StringBuffer line_buffer = new StringBuffer(90);
    final String sequence = toString();

    final int SEQUENCE_LINE_BASE_COUNT = 60;

    int header_counter = 0;

    final int [] header_positions = getFastaHeaderPositions();
    final String [] header_strings = getFastaHeaderStrings();

    int i = 0;

    while(i < sequence.length()) 
    {
      if(header_counter < header_positions.length) 
      {
        if(i == header_positions[header_counter]) 
        {
          // dump a header
          writer.write(">" + header_strings[header_counter] + "\n");
          ++header_counter;
        }
      }

      // get the bases in chunks of at most 60

      int this_line_length;

      if(sequence.length() - i < SEQUENCE_LINE_BASE_COUNT) 
        this_line_length = sequence.length() - i;
      else 
        this_line_length = SEQUENCE_LINE_BASE_COUNT;

      if(header_counter < header_positions.length) 
      {
        if(i + this_line_length > header_positions[header_counter]) 
          this_line_length = header_positions[header_counter] - i;
      }

      line_buffer.setLength(0);
      line_buffer.append(sequence.substring(i, i + this_line_length));
      line_buffer.append("\n");

      writer.write(line_buffer.toString());

      if((i / SEQUENCE_LINE_BASE_COUNT) % 100 == 0) 
        Thread.yield();

      if(header_counter < header_positions.length) 
      {
        if(header_positions[header_counter] < i + SEQUENCE_LINE_BASE_COUNT) 
        {
          i = header_positions[header_counter];
          continue;
        }
      }

      i += SEQUENCE_LINE_BASE_COUNT;
    }
  }

  /**
   *  Return an array containing the character positions of the fasta header
   *  lines.  The positions are returned in the same order as the strings from
   *  getFastaHeaderStrings().
   **/
  int [] getFastaHeaderPositions() 
  {
    if(fasta_header_positions != null && fasta_header_positions.size() > 0) 
    {
      final int [] return_value = new int [fasta_header_positions.size()];

      for(int i = 0 ; i < return_value.length ; ++i) 
        return_value[i] =
                ((Integer)fasta_header_positions.elementAt(i)).intValue();
      
      return return_value;
    } 
    else 
      return new int [0];
  }

  /**
   *  Return an array containing the fasta headers for the input sequence.
   *  The positions are returned in the same order as the strings from
   *  getFastaHeaderPositions().
   **/
  String [] getFastaHeaderStrings() 
  {
    if(fasta_header_strings != null && fasta_header_strings.size() > 0) 
    {
      final String [] return_value = new String [fasta_header_strings.size()];

      for(int i = 0 ; i < return_value.length ; ++i) 
        return_value[i] =(String)fasta_header_strings.elementAt(i);

      return return_value;
    } 
    else 
      return new String [0];
  }

}
