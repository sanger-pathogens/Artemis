/* ComparisonData.java
 *
 * created: Mon Jul 12 1999
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 1999,2000  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/ComparisonData.java,v 1.2 2004-12-14 10:41:42 tjc Exp $
 */

package uk.ac.sanger.artemis;

import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.io.Range;

import uk.ac.sanger.artemis.sequence.*;

/**
 *  Objects that implement this interface provide data (AlignMatch objects)
 *  for the alignment of two sequences.
 *
 *  @author Kim Rutherford
 *  @version $Id: ComparisonData.java,v 1.2 2004-12-14 10:41:42 tjc Exp $
 **/

public interface ComparisonData 
{
  /**
   *  Return an array containing all the AlignMatch objects for this
   *  comparison.
   **/
  public AlignMatch[] getMatches();

  /**
   *  Return all the AlignMatch objects in this comparison which overlap
   *  first_seq_range on the first sequence or second_seq_range on the second
   *  sequence.
   **/
//public AlignMatch[] getMatchesInRange(final Range first_seq_range,
//                                      final Range second_seq_range);

  /**
   *  If this object contain valid matches for a comparison between
   *  first_sequence and second_sequence return null (first_sequence is the
   *  subject of the comparison second_sequence is the query).  If the
   *  comparison would be valid if the data for the ends of the matches were
   *  swapped, then return a copy of this object with all the matches flipped.
   *  (For now, valid means that none of the matches goes over the end of the
   *  sequence.)
   *  @exception OutOfRangeException Thrown if the data in this object is not
   *    valid for either orientation.
   **/
  public ComparisonData flipMatchesIfNeeded(final Bases first_sequence,
                                            final Bases second_sequence)
      throws OutOfRangeException;

  /**
   *  Return the maximum score of all the AlignMatch objects in this object.
   **/
  public int getMaximumScore();

  /**
   *  Return the minimum score of all the AlignMatch objects in this object.
   **/
  public int getMinimumScore();
}
