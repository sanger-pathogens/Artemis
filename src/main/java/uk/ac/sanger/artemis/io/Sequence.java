/* Sequence.java
 *
 * created: Mon Oct 12 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/Sequence.java,v 1.5 2008-12-11 16:54:23 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.ReadOnlyException;
import org.biojava.bio.symbol.IllegalSymbolException;

/**
 *  Sequence interface
 *
 *  @author Kim Rutherford
 *  @version $Id: Sequence.java,v 1.5 2008-12-11 16:54:23 tjc Exp $
 *
 */

public interface Sequence 
{
  /**
   *  Return a the given range of bases as a String.
   *  @param start The start base of the range.
   *  @param end The end base of the range.
   **/
  String getSubSequence(int start, int end);
  char[] getCharSubSequence(int start, int end);

  char charAt(int i);

  /**
   *  Set this sequence to hold the bases in the given String.
   *  @exception ReadOnlyException If this Sequence cannot be changed.
   **/
  void setFromChar(final char sequence[])
      throws ReadOnlyException, IllegalSymbolException;

  /**
   *  Returns the length of the sequence in bases.
   **/
  int length();

  /**
   *  Return the count of c bases in the whole sequence.
   **/
  int getCCount();

  /**
   *  Return the count of g bases in the whole sequence.
   **/
  int getGCount();

  /**
   *  Return the count of a bases in the whole sequence.
   **/
  int getACount();

  /**
   *  Return the count of t bases in the whole sequence.
   **/
  int getTCount();

  /**
   *  Return the count of non-g,c,t,a bases in the whole sequence.
   **/
  int getOtherCount();
  
  //void clear();
}
