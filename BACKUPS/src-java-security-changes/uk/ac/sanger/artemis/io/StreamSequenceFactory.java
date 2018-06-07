/* StreamSequenceFactory.java
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/StreamSequenceFactory.java,v 1.1 2004-06-09 09:50:38 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.LinePushBackReader;

import java.io.IOException;

/**
 *  This class contains the method makeStreamSequence (), which reads a
 *  StreamSequence object from a LinePushBackReader object. 
 *
 *  @author Kim Rutherford
 *  @version $Id: StreamSequenceFactory.java,v 1.1 2004-06-09 09:50:38 tjc Exp $
 **/

abstract public class StreamSequenceFactory 
{

  /** The tag used for sequence that is in EMBL format. */
  final public static int EMBL_FORMAT = 1;

  /** The tag used for sequence that is in GENBANK format. */
  final public static int GENBANK_FORMAT = 2;

  /** The tag use for sequence that is in raw format. */
  final public static int RAW_FORMAT = 3;

  /** The tag use for sequence that is in FASTA or similar format. */
  final public static int FASTA_FORMAT = 4;
  
  final public static int INDEXED_FASTA_FORMAT = 5;

  /** 
   *  Read a StreamSequence object from a LinePushBackReader object.
   **/
  public static StreamSequence makeStreamSequence(final LinePushBackReader
                                                   in_stream, Entry entry)
      throws IOException 
  {
    final int sequence_type = getSequenceType(in_stream);

    switch (sequence_type) 
    {
      case EMBL_FORMAT:
        return new EmblStreamSequence(in_stream);
      case FASTA_FORMAT:
      {
        if(IndexFastaStream.isIndexed(entry))
        {
          IndexFastaStream ifs = new IndexFastaStream(entry);
          if(ifs.useIndex())
            return ifs;
        }
        return new FastaStreamSequence(in_stream);
      }
      case GENBANK_FORMAT:
        return new GenbankStreamSequence(in_stream);
      case RAW_FORMAT:
      default:
        return new RawStreamSequence(in_stream);
    }
  }

  /**
   *  Return the sequence type that can be read from the given
   *  LinePushBackReader.  ie EMBL_FORMAT, FASTA_FORMAT, etc.
   **/
  public static int getSequenceType(LinePushBackReader in_stream)
      throws IOException 
  {
    final String seq_header_line = in_stream.readLine ();

    in_stream.pushBack (seq_header_line);

    if(seq_header_line.startsWith ("SQ   ")) 
      return EMBL_FORMAT;
    else
    {
      if(seq_header_line.startsWith (">")) 
        return FASTA_FORMAT;
      else 
      {
        if (seq_header_line.startsWith ("BASE COUNT") ||
            seq_header_line.startsWith ("ORIGIN")) 
          return GENBANK_FORMAT;
        else 
          return RAW_FORMAT;
      }
    } 
  }

}

