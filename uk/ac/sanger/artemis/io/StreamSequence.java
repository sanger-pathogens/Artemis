/* StreamSequence.java (formally ReaderSequence.java)
 *
 * created: Wed Dec 30 1998
 *
 * This file is part of Artemis
 *
 * Copyright (C) 1998,1999,2000  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/StreamSequence.java,v 1.3 2004-12-23 15:33:46 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import java.io.IOException;
import java.io.Writer;

/**
 *  This is an implementation of Sequence that can read and write itself to a
 *  stream.
 *
 *  @author Kim Rutherford
 *  @version $Id: StreamSequence.java,v 1.3 2004-12-23 15:33:46 tjc Exp $
 **/

public abstract class StreamSequence
    extends LineGroup implements Sequence {
  /**
   *  Create a new StreamSequence object that contains no sequence.
   **/
  protected StreamSequence () {
    setFromString ("");
  }
  
  /**
   *  Return a new StreamSequence object that is a copy of this one.
   **/
  abstract public StreamSequence copy ();

  /**
   *  Return the sequence type (one of EMBL_FORMAT, RAW_FORMAT, FASTA_FORMAT,
   *  etc.)
   **/
  abstract public int getFormatType ();

  /**
   *  Return the type of this LineGroup object (SEQUENCE, FEATURE_TABLE, etc.)
   *  Returns SEQUENCE for objects of this class.
   **/
  public int getType () {
    return SEQUENCE;
  }

  /**
   *  Return the contents of this Sequence as a String.  This currently does
   *  not return a copy so any changes will change the Sequence object.
   **/
  public String toString () {
    return sequence;
  }

  /**
   *  Return a the given range of bases as a String.  Returns an empty
   *  sequence if the end position is less than the start position.
   *  @param start The start base of the range.
   *  @param end The end base of the range.
   **/
  public String getSubSequence (int start, int end) {
    if (end < start) {
      // empty range
      return "";
    } else {
      if (start == 1 && end == length ()) {
        return sequence;
      }

      return sequence.substring (start - 1, end);
    }
  }

  public char[] getCharSubSequence (int start, int end) 
  {
    char[] dst = new char[end-start+1];
    sequence.getChars(start-1, end, dst, 0);
    return dst;
  }


  /**
   *  Set this sequence to hold the bases in the given String.
   **/
  public void setFromString (final String new_sequence) {
    sequence = new_sequence;
    setCounts ();
  }

  /**
   *  Write this Sequence to the given stream.
   *  @param writer The stream to write to.
   **/
  public abstract void writeToStream (final Writer writer)
      throws IOException;

  /**
   *  Returns the length of the sequence in bases.
   **/
  public int length () {
    return sequence.length ();
  }

  /**
   *  Return the count of c bases in the whole sequence.
   **/
  public int getCCount () {
    return c_count;
  }

  /**
   *  Return the count of g bases in the whole sequence.
   **/
  public int getGCount () {
    return g_count;
  }

  /**
   *  Return the count of a bases in the whole sequence.
   **/
  public int getACount () {
    return a_count;
  }

  /**
   *  Return the count of t bases in the whole sequence.
   **/
  public int getTCount () {
    return t_count;
  }

  /**
   *  Return the count of non-g,c,t,a bases in the whole sequence.
   **/
  public int getOtherCount () {
    return
      length () - (getCCount () + getACount () + getTCount () + getGCount ());
  }

  /**
   *  Set the a_count, c_count, t_count and g_count variables.
   **/
  private void setCounts () {
    a_count = c_count = t_count = g_count = 0;

    final int sequence_length = sequence.length ();

    final char [] sequence_chars = new char [sequence_length];

    sequence.getChars (0, sequence_length, sequence_chars, 0);

    for (int i = 0 ; i < sequence_length ; ++i) {
      switch (sequence_chars[i]) {
      case 'a':
        ++a_count;
        break;
      case 'c':
        ++c_count;
        break;
      case 'g':
        ++g_count;
        break;
      case 't':
        ++t_count;
        break;
      default:
        break;
      }
    }
  }

  /**
   *  Contains the sequence data for this object.  It will contain the bases
   *  of the sequence with no spaces after the Feature constructor finishes.
   **/
  private String sequence;

  /**
   *  Count of the a bases in the sequence.
   **/
  private int a_count = 0;

  /**
   *  Count of the c bases in the sequence.
   **/
  private int c_count = 0;

  /**
   *  Count of the g bases in the sequence.
   **/
  private int g_count = 0;

  /**
   *  Count of the t bases in the sequence.
   **/
  private int t_count = 0;
}
