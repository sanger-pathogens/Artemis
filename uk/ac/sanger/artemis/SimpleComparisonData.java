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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/SimpleComparisonData.java,v 1.1 2004-06-09 09:45:07 tjc Exp $
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
 *  @version $Id: SimpleComparisonData.java,v 1.1 2004-06-09 09:45:07 tjc Exp $
 **/

abstract class SimpleComparisonData implements ComparisonData {
  /**
   *  Create a new SimpleComparisonData by reading from the given
   *  LinePushBackReader.
   **/
  public SimpleComparisonData (final LinePushBackReader stream)
      throws IOException {
    final Vector align_match_vector = new Vector ();

    while (true) {
      final String line = stream.readLine ();

      if (line == null) {
        break;
      }

      if (line.trim ().length () == 0) {
        continue;
      }

      final AlignMatch new_match = makeMatchFromString (line);

      if (new_match == null) {
        // hit a blank line or a comment - loop again
      } else {
        align_match_vector.addElement (new_match);
      }
    }

    final AlignMatch [] matches = new AlignMatch [align_match_vector.size ()];

    for (int i = 0 ; i < matches.length ; ++i) {
      matches[i] = (AlignMatch) align_match_vector.elementAt (i);
    }

    setMatches (matches);
  }

  /**
   *  Create a new, empty instance of SimpleComparisonData.
   **/
  protected SimpleComparisonData () {

  }

  /**
   *  Return an array containing all the AlignMatch objects for this
   *  comparison.
   **/
  public AlignMatch [] getMatches () {
    return matches;
  }

  /**
   *  Return all the AlignMatch objects in this comparison which overlap
   *  subject_seq_range on the subject sequence or query_seq_range on the query
   *  sequence.
   **/
  public AlignMatch [] getMatchesInRange (final Range subject_seq_range,
                                             final Range query_seq_range) {

    // a count of how many objects we have put into match_buffer so far.
    int match_buffer_count = 0;

    for (int i = 0 ; i < spare_buckets.size () ; ++i) {
      final AlignMatch this_match = (AlignMatch) spare_buckets.elementAt (i);

      if (matchInRange (this_match, subject_seq_range, query_seq_range)) {
        match_buffer[match_buffer_count] = this_match;
        ++match_buffer_count;
      }
    }

    // used to make sure we don't return any duplicates
    final Hashtable table = new Hashtable (100);

    for (int bucket_index = subject_seq_range.getStart () / BUCKET_SIZE ;
         bucket_index < subject_seq_range.getEnd () / BUCKET_SIZE ;
         ++bucket_index) {
      for (int i = 0 ;
           i < subject_sequence_buckets[bucket_index].size () ;
           ++i) {
        final AlignMatch this_match =
          (AlignMatch) subject_sequence_buckets[bucket_index].elementAt (i);

        if (this_match.getSubjectSequenceRange ().overlaps (subject_seq_range)) {
          match_buffer[match_buffer_count] = this_match;
          ++match_buffer_count;
          table.put (this_match, this_match);
        }
      }
    }

    for (int bucket_index = query_seq_range.getStart () / BUCKET_SIZE ;
         bucket_index < query_seq_range.getEnd () / BUCKET_SIZE ;
         ++bucket_index) {
      for (int i = 0 ;
           i < query_sequence_buckets[bucket_index].size () ;
           ++i) {
        final AlignMatch this_match =
          (AlignMatch) query_sequence_buckets[bucket_index].elementAt (i);

        if (table.containsKey (this_match)) {
          continue;
        }

        if (this_match.getQuerySequenceRange ().overlaps (query_seq_range)) {
          match_buffer[match_buffer_count] = this_match;
          ++match_buffer_count;
        }
      }
    }

    final AlignMatch [] return_matches = new AlignMatch [match_buffer_count];

    System.arraycopy (match_buffer, 0,
                      return_matches, 0,
                      return_matches.length);

    return return_matches;
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
  public ComparisonData flipMatchesIfNeeded (final Bases subject_sequence,
                                             final Bases query_sequence)
      throws OutOfRangeException {
    final AlignMatch forward_error_match =
      checkMatches (subject_sequence, query_sequence);

    if (forward_error_match == null) {
      return null;
    } else {
      final AlignMatch reverse_error_match =
        checkMatches (query_sequence, subject_sequence);

      if (reverse_error_match == null) {
        final SimpleComparisonData new_comparison_data =
          getNewSimpleComparisonData ();

        final AlignMatch [] new_matches = new AlignMatch [matches.length];

        for (int i = 0 ; i < matches.length ; ++i) {
          final AlignMatch this_match = matches[i];

          final AlignMatch new_match =
            new AlignMatch (this_match.getQuerySequenceRange (),
                            this_match.getSubjectSequenceRange (),
                            this_match.isRevMatch (),
                            this_match.getScore (),
                            this_match.getPercentID ());

          new_matches [i] = new_match;
        }

        new_comparison_data.setMatches (new_matches);

        return new_comparison_data;
      } else {
        final String message;

        if (forward_error_match.getSubjectSequenceStart () >
            subject_sequence.getLength ()) {
          message = "match goes off end of subject sequence: " +
            forward_error_match.getSubjectSequenceStart ();
        } else {
          if (forward_error_match.getSubjectSequenceEnd () >
              subject_sequence.getLength ()) {
            message = "match goes off end of subject sequence: " +
              forward_error_match.getSubjectSequenceEnd ();
          } else {
            if (forward_error_match.getQuerySequenceStart () >
                query_sequence.getLength ()) {
              message = "match goes off end of query sequence: " +
                forward_error_match.getQuerySequenceStart ();
            } else {
              if (forward_error_match.getQuerySequenceEnd () >
                  query_sequence.getLength ()) {
                message = "match goes off end of query sequence: " +
                  forward_error_match.getQuerySequenceEnd ();
              } else {
                throw new Error ("internal error - unreachable code");
              }
            }
          }
        }

        throw new OutOfRangeException (message);
      }
    }
  }

  /**
   *  Returns a new, empty instance of this type of object;
   **/
  abstract protected SimpleComparisonData getNewSimpleComparisonData ();

  /**
   *  Make an AlignMatch object from the given String.  The String must be in
   *  a format appropriate for this object.
   **/
  abstract protected AlignMatch makeMatchFromString (final String line)
      throws IOException;

  /**
   *  Return null if and only if this object contains only valid matches for a
   *  comparison between subject_sequence and query_sequence.  The first
   *  invalid AlignMatch is returned otherwise.
   **/
  private AlignMatch checkMatches (final Bases subject_sequence,
                                final Bases query_sequence) {
    for (int i = 0 ; i < matches.length ; ++i) {
      final AlignMatch match = matches[i];

      if (match.getSubjectSequenceStart () > subject_sequence.getLength () ||
          match.getSubjectSequenceEnd () > subject_sequence.getLength ()) {
        return match;
      }

      if (match.getQuerySequenceStart () > query_sequence.getLength () ||
          match.getQuerySequenceEnd () > query_sequence.getLength ()) {
        return match;
      }
    }

    return null;
  }

  /**
   *  Return true if and only if the given AlignMatch object overlaps
   *  subject_seq_range on the subject sequence or query_seq_range on the
   *  query sequence.
   **/
  private boolean matchInRange (final AlignMatch match,
                                final Range subject_seq_range,
                                final Range query_seq_range) {
    if (match.getSubjectSequenceRange ().overlaps (subject_seq_range) ||
        match.getQuerySequenceRange ().overlaps (query_seq_range)) {
      return true;
    } else {
      return false;
    }
  }

  /**
   *  Set the array of AlignMatch objects.
   **/
  protected void setMatches (final AlignMatch [] matches) {
    this.matches = matches;

    match_buffer = new AlignMatch [matches.length];

    for (int i = 0 ; i < matches.length ; ++i) {
      final AlignMatch this_match = matches[i];

      final int score = this_match.getScore ();

      final int this_match_subject_sequence_end =
        this_match.getSubjectSequenceEnd ();

      final int this_match_query_sequence_end =
        this_match.getQuerySequenceEnd ();

      if (this_match_subject_sequence_end > subject_sequence_max_base) {
        subject_sequence_max_base = this_match_subject_sequence_end;
      }

      if (this_match_query_sequence_end > query_sequence_max_base) {
        query_sequence_max_base = this_match_query_sequence_end;
      }
    }
  }

  /**
   *  The number of base per bucket.
   **/
  final private int BUCKET_SIZE = 1000;

  /**
   *  Create subject_sequence_buckets, query_sequence_buckets and
   *  spare_buckets.
   **/
  private void makeBuckets () {
    subject_sequence_buckets =
      new Vector [subject_sequence_max_base / BUCKET_SIZE + 1];
    query_sequence_buckets =
      new Vector [query_sequence_max_base / BUCKET_SIZE + 1];

    for (int i = 0 ; i < matches.length ; ++i) {
      final AlignMatch match = matches[i];

      if (match.getSubjectSequenceRange ().getCount () > BUCKET_SIZE ||
          match.getQuerySequenceRange ().getCount () > BUCKET_SIZE) {
        spare_buckets.addElement (match);
      } else {
        final int match_subject_sequence_start =
          match.getSubjectSequenceStart ();

        final int match_query_sequence_start =
          match.getQuerySequenceStart ();

        final int subject_buckets_index =
          match_subject_sequence_start / BUCKET_SIZE;
        subject_sequence_buckets[subject_buckets_index].addElement (match);

        final int query_buckets_index =
          match_query_sequence_start / BUCKET_SIZE;
        query_sequence_buckets[query_buckets_index].addElement (match);
      }
    }
  }

  /**
   *  Make and return a new AlignMatch.
   **/
  static protected AlignMatch makeAlignMatch (int subject_sequence_start,
                                              int subject_sequence_end,
                                              int query_sequence_start,
                                              int query_sequence_end,
                                              final int score,
                                              final int percent_id) {
    try {
      // true if and only if the query hits the reverse complement of the
      // subject
      boolean rev_match = false;

      if (subject_sequence_end < subject_sequence_start) {
        final int tmp = subject_sequence_start;
        subject_sequence_start = subject_sequence_end;
        subject_sequence_end = tmp;
        rev_match = !rev_match;
      }

      if (query_sequence_end < query_sequence_start) {
        final int tmp = query_sequence_start;
        query_sequence_start = query_sequence_end;
        query_sequence_end = tmp;
        rev_match = !rev_match;
      }

      return new AlignMatch (new Range (subject_sequence_start,
                                        subject_sequence_end),
                             new Range (query_sequence_start,
                                        query_sequence_end),
                             rev_match, score, percent_id);
    } catch (OutOfRangeException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    }
  }

  /**
   *  Set the values of min_score and max_score.
   **/
  private void setMinMaxScore () {
    for (int i = 0 ; i < matches.length ; ++i) {
      final AlignMatch this_match = matches[i];

      final int score = this_match.getScore ();

      if (score > -1) {
        if (score > max_score) {
          max_score = score;
        }

        if (score < min_score) {
          min_score = score;
        }
      }
    }
  }

  /**
   *  Return the maximum score of all the AlignMatch objects in this object.
   **/
  public int getMaximumScore () {
    if (max_score == -1) {
      setMinMaxScore ();
    }

    return max_score;
  }

  /**
   *  Return the minimum score of all the AlignMatch objects in this object.
   **/
  public int getMinimumScore () {
    if (max_score == -1) {
      setMinMaxScore ();
    }

    return min_score;
  }

  /**
   *  This is the array of matches created by the constructor.
   **/
  private AlignMatch [] matches;

  /**
   *  This is the array is used as a buffer by getMatchesInRange ().
   **/
  private AlignMatch [] match_buffer;

  /**
   *  Set by the constructor and returned by getMaximumScore ().
   **/
  private int max_score = -1;

  /**
   *  Set by the constructor and returned by getMinimumScore ().
   **/
  private int min_score = 999999999;

  /**
   *  Set by setMatches () to be the highest base we see in the subject
   *  sequence.
   **/
  private int subject_sequence_max_base = -1;

  /**
   *  Set by setMatches () to be the highest base we see in the query
   *  sequence.
   **/
  private int query_sequence_max_base = -1;

  /**
   *  This array contains a Vector for each BUCKET_SIZE bases in the subject
   *  sequence.  All AlignMatch objects where the match start in the subject
   *  sequence (ie AlignMatch.getSubjectSequenceStart ()) is >= 1 and <=
   *  BUCKET_SIZE will be put in the subject bucket.  If >= BUCKET_SIZE + 1
   *  and <= BUCKET_SIZE * 2 it will be in the query bucket, etc.
   **/
  private Vector [] subject_sequence_buckets = null;

  /**
   *  This array contains a Vector for each BUCKET_SIZE bases in the query
   *  sequence.
   **/
  private Vector [] query_sequence_buckets = null;

  /**
   *  This Vector contains the AlignMatch objects where the match is bigger
   *  than BUCKET_SIZE in either of the sequences.
   **/
  private Vector spare_buckets = new Vector ();
}
