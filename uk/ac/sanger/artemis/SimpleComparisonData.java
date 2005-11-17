/* SimpleComparisonData.java
 *
 * created: Wed May 17 2000
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/SimpleComparisonData.java,v 1.3 2005-11-17 16:50:50 tjc Exp $
 */

package uk.ac.sanger.artemis;

import uk.ac.sanger.artemis.sequence.*;

import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.util.LinePushBackReader;

import java.io.*;
import java.util.StringTokenizer;
import java.util.Vector;
import java.util.Hashtable;

/**
 *  This class contains methods that are common to all ComparisonData
 *  objects.  In particular it has methods for managing AlignMatch objects.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: SimpleComparisonData.java,v 1.3 2005-11-17 16:50:50 tjc Exp $
 **/

abstract class SimpleComparisonData implements ComparisonData 
{
  /** array of matches created by the constructor */
  private AlignMatch [] matches;

  /** array is used as a buffer */
  private AlignMatch [] match_buffer;

  /** Set by the constructor and returned by getMaximumScore() */
  private int max_score = -1;

  /** Set by the constructor and returned by getMinimumScore() */
  private int min_score = 999999999;

  /** Set by setMatches() to be the highest base we see in the subject */
  private int subject_sequence_max_base = -1;

  /** Set by setMatches() to be the highest base we see in the query   */
  private int query_sequence_max_base = -1;

  /**
   *  Create a new SimpleComparisonData by reading from the given
   *  LinePushBackReader.
   **/
  public SimpleComparisonData(final LinePushBackReader stream)
      throws IOException 
  {
    final Vector align_match_vector = new Vector();

    String line;

    while( (line = stream.readLine()) != null )
    {
      if(line.trim().length() == 0) 
        continue;

      final AlignMatch new_match = makeMatchFromString(line);

      // not a blank line or a comment
      if(new_match != null) 
        align_match_vector.addElement(new_match);
    }

    final AlignMatch[] matches = new AlignMatch[align_match_vector.size()];

    for(int i = 0; i < matches.length; ++i) 
      matches[i] = (AlignMatch)align_match_vector.elementAt(i);

    setMatches(matches);
  }

  /**
   *  Create a new, empty instance of SimpleComparisonData.
   **/
  protected SimpleComparisonData() 
  {
  }

  /**
   *  Return an array containing all the AlignMatch objects for this
   *  comparison.
   **/
  public AlignMatch[] getMatches() 
  {
    return matches;
  }

  
  /**
   *  If this object contains only valid matches for a comparison between
   *  subject_sequence and query_sequence return null (subject_sequence is the
   *  subject of the comparison query_sequence is the query).  If the
   *  comparison would be valid if the data for the ends of the matches were
   *  swapped, then return a copy of this object with all the matches flipped.
   *  (For now, valid means that none of the matches goes over the end of the
   *  sequence.)
   *  @exception OutOfRangeException Thrown if the data in this object is not
   *    valid for either orientation.
   **/
  public ComparisonData flipMatchesIfNeeded(final Bases subject_sequence,
                                            final Bases query_sequence)
      throws OutOfRangeException 
  {
    final AlignMatch forward_error_match =
      checkMatches(subject_sequence, query_sequence);

    if(forward_error_match == null) 
      return null;
    else 
    {
      final AlignMatch reverse_error_match =
        checkMatches(query_sequence, subject_sequence);

      if(reverse_error_match == null)
      {
        final SimpleComparisonData new_comparison_data =
                                   getNewSimpleComparisonData();

        int length = matches.length;
        final AlignMatch[] new_matches = new AlignMatch[length];

        for(int i = 0; i < length; ++i)
        {
          final AlignMatch this_match = matches[i];

          final AlignMatch new_match =
            new AlignMatch(this_match.getQuerySequenceRange(),
                           this_match.getSubjectSequenceRange(),
                           this_match.isRevMatch(),
                           this_match.getScore(),
                           this_match.getPercentID());

          new_matches[i] = new_match;
        }

        new_comparison_data.setMatches(new_matches);

        return new_comparison_data;
      } 
      else
      {
        final String message;

        if(forward_error_match.getSubjectSequenceStart() >
           subject_sequence.getLength()) 
          message = "match goes off end of subject sequence: " +
                    forward_error_match.getSubjectSequenceStart();
        else 
        {
          if(forward_error_match.getSubjectSequenceEnd() >
             subject_sequence.getLength()) 
            message = "match goes off end of subject sequence: " +
                      forward_error_match.getSubjectSequenceEnd();
          else
          {
            if(forward_error_match.getQuerySequenceStart() >
               query_sequence.getLength())
              message = "match goes off end of query sequence: " +
                         forward_error_match.getQuerySequenceStart();
            else 
            {
              if(forward_error_match.getQuerySequenceEnd() >
                 query_sequence.getLength())
                message = "match goes off end of query sequence: " +
                           forward_error_match.getQuerySequenceEnd();
              else 
                throw new Error("internal error - unreachable code");
            }
          }
        }
        throw new OutOfRangeException(message);
      }
    }
  }

  /**
   *  Returns a new, empty instance of this type of object;
   **/
  abstract protected SimpleComparisonData getNewSimpleComparisonData();

  /**
   *  Make an AlignMatch object from the given String.  The String must be in
   *  a format appropriate for this object.
   **/
  abstract protected AlignMatch makeMatchFromString(final String line)
      throws IOException;

  /**
   *  Return null if and only if this object contains only valid matches for a
   *  comparison between subject_sequence and query_sequence.  The first
   *  invalid AlignMatch is returned otherwise.
   **/
  private AlignMatch checkMatches(final Bases subject_sequence,
                                  final Bases query_sequence) 
  {
    int length = matches.length;
    for(int i = 0; i < length; ++i) 
    {
      final AlignMatch match = matches[i];

      if(match.getSubjectSequenceStart() > subject_sequence.getLength() ||
         match.getSubjectSequenceEnd() > subject_sequence.getLength()) 
        return match;

      if(match.getQuerySequenceStart() > query_sequence.getLength() ||
         match.getQuerySequenceEnd() > query_sequence.getLength())
        return match;
    }

    return null;
  }


  /**
   *  Set the array of AlignMatch objects.
   **/
  protected void setMatches(final AlignMatch[] matches) 
  {
    this.matches = matches;

    int length = matches.length;
    match_buffer = new AlignMatch[length];

    for(int i = 0 ; i < length ; ++i) 
    {
      final AlignMatch this_match = matches[i];

      final int score = this_match.getScore();

      final int this_match_subject_sequence_end =
        this_match.getSubjectSequenceEnd();

      final int this_match_query_sequence_end =
        this_match.getQuerySequenceEnd();

      if(this_match_subject_sequence_end > subject_sequence_max_base) 
        subject_sequence_max_base = this_match_subject_sequence_end;

      if(this_match_query_sequence_end > query_sequence_max_base) 
        query_sequence_max_base = this_match_query_sequence_end;
    }
  }

  /**
   *  Make and return a new AlignMatch.
   **/
  static protected AlignMatch makeAlignMatch(int subject_sequence_start,
                                             int subject_sequence_end,
                                             int query_sequence_start,
                                             int query_sequence_end,
                                             final int score,
                                             final int percent_id) 
  {
    try 
    {
      // true if and only if the query hits the reverse complement of the
      // subject
      boolean rev_match = false;

      if(subject_sequence_end < subject_sequence_start) 
      {
        final int tmp = subject_sequence_start;
        subject_sequence_start = subject_sequence_end;
        subject_sequence_end = tmp;
        rev_match = !rev_match;
      }

      if(query_sequence_end < query_sequence_start) 
      {
        final int tmp = query_sequence_start;
        query_sequence_start = query_sequence_end;
        query_sequence_end = tmp;
        rev_match = !rev_match;
      }

      return new AlignMatch(new Range(subject_sequence_start,
                                      subject_sequence_end),
                            new Range(query_sequence_start,
                                      query_sequence_end),
                            rev_match, score, percent_id);
    }
    catch(OutOfRangeException e) 
    {
      throw new Error("internal error - unexpected exception: " + e);
    }
  }

  /**
   *  Set the values of min_score and max_score.
   **/
  private void setMinMaxScore()
  {
    int length = matches.length;
    for(int i = 0; i < length; ++i) 
    {
      final AlignMatch this_match = matches[i];

      final int score = this_match.getScore();

      if(score > -1) 
      {
        if(score > max_score) 
          max_score = score;

        if(score < min_score)
          min_score = score;
      }
    }
  }

  /**
   *  Return the maximum score of all the AlignMatch objects in this object.
   **/
  public int getMaximumScore() 
  {
    if(max_score == -1)
      setMinMaxScore();

    return max_score;
  }

  /**
   *  Return the minimum score of all the AlignMatch objects in this object.
   **/
  public int getMinimumScore() 
  {
    if(max_score == -1) 
      setMinMaxScore();

    return min_score;
  }

}
