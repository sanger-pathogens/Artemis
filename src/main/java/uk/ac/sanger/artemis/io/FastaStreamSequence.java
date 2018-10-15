/* FastaStreamSequence.java
 *
 * created: Mon Jun 14 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/FastaStreamSequence.java,v 1.1 2004-06-09 09:49:19 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.LinePushBackReader;

import java.io.IOException;
import java.io.Writer;

/**
 *  This is a subclass of StreamSequence containing Fasta sequence.
 *
 *  @author Kim Rutherford
 *  @version $Id: FastaStreamSequence.java,v 1.1 2004-06-09 09:49:19 tjc Exp $
 **/

public class FastaStreamSequence extends RawStreamSequence {
  /**
   *  Create a new FastaStreamSequence object from a stream.  The next line to
   *  read from the stream should be the header line of the sequence.
   *  @param in_stream The stream to read from.  When the constructor returns
   *    the stream will at the next line after the sequence.
   **/
  public FastaStreamSequence (final LinePushBackReader in_stream)
      throws IOException {
    super (in_stream);
  }

  /**
   *  Make a new FastaStreamSequence containing the same sequence as the given
   *  Sequence.  The header will be ">artemis_sequence" unless the Sequence is
   *  a FastaStreamSequence.
   **/
  public FastaStreamSequence (final Sequence sequence) {
    super (sequence);
  }

  /**
   *  Make a new FastaStreamSequence containing the same sequence as the given
   *  String.
   *  @param sequence_string The String containing the sequence for the new
   *    FastaStreamSequence.
   *  @param header The header to use for the new object.
   **/
  public FastaStreamSequence (final String sequence_string,
                              final String header) {
    super (sequence_string);
    this.header = header;
  }

  /**
   *  Make a new FastaStreamSequence containing the same sequence as the given
   *  String.  The header will be ">"
   *  @param sequence_string The String containing the sequence for the new
   *    FastaStreamSequence.
   **/
  public FastaStreamSequence (final String sequence_string) {
    this (sequence_string, "");
  }

  /**
   *  Return a new StreamSequence object that is a copy of this one.
   **/
  public StreamSequence copy () {
    return new FastaStreamSequence (this);
  }

  /**
   *  Return the sequence type (EMBL_FORMAT for this class).
   **/
  public int getFormatType () {
    return StreamSequenceFactory.FASTA_FORMAT;
  }

  /**
   *  Write this Sequence to the given stream.  The output will be in FASTA
   *  format.
   *  @param writer The stream to write to.
   **/
  public void writeToStream (final Writer writer)
      throws IOException {
    if (getHeader () != null) {
      writer.write (">" + getHeader () + "\n");
    }

    super.writeToStream (writer);
  }

  /**
   *  Return the header that was passed to the FastaStreamSequence (final
   *  String sequence_string, final String header) constructor.
   **/
  private String getHeader () {
    return header;
  }

  /**
   *  The header to use in writeToStream() if the FastaStreamSequence was
   *  created with the (final String sequence_string, final String header)
   *  constructor.
   **/
  private String header = null;
}
