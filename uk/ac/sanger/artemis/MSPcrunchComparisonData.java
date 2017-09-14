/* MSPcrunchComparisonData.java
 *
 * created: Tue Apr  4 2000
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/MSPcrunchComparisonData.java,v 1.2 2005-11-21 15:49:02 tjc Exp $
 */

package uk.ac.sanger.artemis;

import uk.ac.sanger.artemis.util.LinePushBackReader;

import java.io.IOException;
import java.io.Writer;
import java.util.StringTokenizer;

/**
 *  This class implements the ComparisonData interface for MSPcrunch -d
 *  output.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: MSPcrunchComparisonData.java,v 1.2 2005-11-21 15:49:02 tjc Exp $
 **/

public class MSPcrunchComparisonData extends SimpleComparisonData
    implements ComparisonData {
  /**
   *  Create a new MSPcrunchComparisonData by reading from the given
   *  LinePushBackReader.
   **/
  public MSPcrunchComparisonData (final LinePushBackReader stream)
      throws IOException {
    super (stream);
  }

  /**
   *  Create a new, empty instance of MSPcrunchComparisonData.
   **/
  protected MSPcrunchComparisonData () {
  }

  /**
   *  Returns a new, empty instance of this type of object;
   **/
  protected SimpleComparisonData getNewSimpleComparisonData () {
    return new MSPcrunchComparisonData ();
  }

  /**
   *  Make an AlignMatch object from the given String.
   **/
  private static AlignMatch makeMatchFromStringStatic (final String line)
      throws IOException {
    final StringTokenizer tokenizer = new StringTokenizer (line, " ");

    if (tokenizer.countTokens () < 8) {
      final String message = "while reading MSPcrunch data: " +
        "not enough fields in this line: " + line;
      throw new ComparisonDataParseException (message);
    }

    final String score_token = tokenizer.nextToken ();
    final String percent_ident_token = tokenizer.nextToken ();
    final String q_start_token = tokenizer.nextToken ();
    final String q_end_token = tokenizer.nextToken ();

    tokenizer.nextToken ();

    final String s_start_token = tokenizer.nextToken ();
    final String s_end_token = tokenizer.nextToken ();

    try {
      int score;
      try {
        score = Integer.valueOf (score_token).intValue ();
      } catch (NumberFormatException e) {
        score = Float.valueOf (score_token).intValue (); // blast+
      }
      final int percent_ident =
        (int)(Float.valueOf (percent_ident_token).floatValue ());
      final int q_start = Integer.valueOf (q_start_token).intValue ();
      final int q_end   = Integer.valueOf (q_end_token).intValue ();
      final int s_start = Integer.valueOf (s_start_token).intValue ();
      final int s_end   = Integer.valueOf (s_end_token).intValue ();

      final AlignMatch new_match =
        makeAlignMatch (s_start, s_end, q_start, q_end, score, percent_ident);

      return new_match;
    } catch (NumberFormatException e) {
      throw new IOException ("while reading MSPcrunch data: " +
                             "failed to parse a number from this string: " +
                             e.getMessage ());
    }
  }

  public static void writeMatchFromAlignMatch(final AlignMatch match,
                                              final String query, final String subject,
                                              final Writer writer)
     throws IOException 
  {
    writer.write( match.getScore() + " " +
                  match.getPercentID() + " " +
                  match.getQuerySequenceStart() + " " +
                  match.getQuerySequenceEnd() + " " +
                  query + " " +
                  match.getSubjectSequenceStart() + " " +
                  match.getSubjectSequenceEnd() + " " +
                  subject + "\n" );
  }

  /**
   *  Make an AlignMatch object from the given String.  The String must be in
   *  a format appropriate for this object.
   **/
  protected AlignMatch makeMatchFromString (final String line)
      throws IOException {
    return makeMatchFromStringStatic (line);
  }

  /**
   *  Returns true if and only if the given line is in the correct format for
   *  this type of ComparisonData.  This should be as strict as possible.
   **/
  public static boolean formatCorrect (final String line) {
    try {
      makeMatchFromStringStatic (line);
    } catch (IOException e) {
      return false;
    }
    return true;
  }
}
