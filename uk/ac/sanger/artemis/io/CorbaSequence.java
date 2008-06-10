/* CorbaSequence.java
 *
 * created: Mon Feb  8 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/CorbaSequence.java,v 1.4 2008-06-10 15:33:13 tjc Exp $
 **/

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.ReadOnlyException;

import java.io.IOException;

/**
 *  This is a subclass of Sequence that can read itself from a Corba embl
 *  object.
 *
 *  @author Kim Rutherford
 *  @version $Id: CorbaSequence.java,v 1.4 2008-06-10 15:33:13 tjc Exp $
 **/

public class CorbaSequence implements Sequence {
  /**
   *  Create a new CorbaSequence object from the given handle.
   *  @param data This is the corba object that we will read from.
   **/
  public CorbaSequence (final nsdb.EmblSeq corba_handle) {
    this.corba_handle = corba_handle;

    sequence  = corba_handle.getSeq();
  }

  /**
   *  Return this Sequence as a String.
   **/
  public String toString () {
    return sequence;
  }

  /**
   *  Return a the given range of bases as a String.
   *  @param start The start base of the range.
   *  @param end The end base of the range.
   **/
  public String getSubSequence (int start, int end) {
    if (start == 1 && end == length ()) {
      return sequence;
    }

    if (end < start) {
      // empty range
      return "";
    } else {
      return sequence.substring (start - 1, end);
    }
  }

  public char[] getCharSubSequence (int start, int end)
  {
    char[] dst = new char[end-start+1];
    StringBuffer buff = new StringBuffer(sequence);
    buff.getChars(start-1, end, dst, 0);
    return dst;
  }

  public char charAt(int i)
  {
    return sequence.charAt(i);
  }

  /**
   *  Set this sequence to hold the bases in the given String - throws
   *  ReadOnlyException for CorbaSequence objects.
   *  @exception ReadOnlyException If this Sequence cannot be changed.
   **/
  public void setFromChar(final char[] new_sequence)
      throws ReadOnlyException {
    throw new ReadOnlyException ();
  }

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
    return corba_handle.getCountC();
  }

  /**
   *  Return the count of g bases in the whole sequence.
   **/
  public int getGCount () {
    return corba_handle.getCountG();
  }

  /**
   *  Return the count of a bases in the whole sequence.
   **/
  public int getACount () {
    return corba_handle.getCountA();
  }

  /**
   *  Return the count of t bases in the whole sequence.
   **/
  public int getTCount () {
    return corba_handle.getCountT();
  }

  /**
   *  Return the count of non-g,c,t,a bases in the whole sequence.
   **/
  public int getOtherCount () {
    return
      length () - (getCCount () + getACount () + getTCount () + getGCount ());
  }

  public void clear() {}

  /**
   *  The sequence that was read from the Corba object by the constructor.
   **/
  private String sequence;
  
  /**
   *  The corba object that was passed to the constructor.
   **/
  private nsdb.EmblSeq corba_handle;
}


