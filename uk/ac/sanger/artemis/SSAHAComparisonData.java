/* SSAHAComparisonData.java
 *
 * created: Mon Jul  9 2001
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 2001  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/SSAHAComparisonData.java,v 1.1 2004-06-09 09:45:02 tjc Exp $
 */

package uk.ac.sanger.artemis;

import uk.ac.sanger.artemis.util.LinePushBackReader;

import java.io.*;
import java.util.StringTokenizer;
import java.util.Vector;

/**
 *  This class implements the ComparisonData interface for ssaha -pf
 *  output.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: SSAHAComparisonData.java,v 1.1 2004-06-09 09:45:02 tjc Exp $
 **/

public class SSAHAComparisonData extends SimpleComparisonData {
  /**
   *  Create a new SSAHAComparisonData by reading from the given
   *  LinePushBackReader.
   **/
  public SSAHAComparisonData (final LinePushBackReader stream)
      throws IOException {
    super (stream);
  }

  /**
   *  Create a new, empty instance of SSAHAComparisonData.
   **/
  public SSAHAComparisonData () {
    
  }

  /**
   *  Returns a new, empty instance of this type of object;
   **/
  protected SimpleComparisonData getNewSimpleComparisonData () {
    return new SSAHAComparisonData ();
  }

  /**
   *  Make an AlignMatch object from the given String.
   **/
  public static AlignMatch makeMatchFromStringStatic (final String line)
      throws IOException {

    final StringTokenizer tokenizer = new StringTokenizer (line, "\t");

    if (tokenizer.countTokens () != 9) {
      final String message = "while reading SSAHA data: " +
        "not enough fields in this line: " + line;
      throw new ComparisonDataParseException (message);
    }

    final String direction_token = tokenizer.nextToken ();
    // throw away the query name
    tokenizer.nextToken ();
    final String q_start_token = tokenizer.nextToken ();
    final String q_end_token = tokenizer.nextToken ();
    // throw away the subject name
    tokenizer.nextToken ();
    final String s_start_token = tokenizer.nextToken ();
    final String s_end_token = tokenizer.nextToken ();

    final String score_token = tokenizer.nextToken ();
    final String percent_ident_token = tokenizer.nextToken ();

    try {
      final int score   = Integer.valueOf (score_token).intValue ();
      final int percent_ident =
        (int)(Float.valueOf (percent_ident_token).floatValue ());
      final int q_start = Integer.valueOf (q_start_token).intValue ();
      final int q_end   = Integer.valueOf (q_end_token).intValue ();
      final int s_start = Integer.valueOf (s_start_token).intValue ();
      final int s_end   = Integer.valueOf (s_end_token).intValue ();

      if (direction_token.equals ("F")) {
        return makeAlignMatch (s_start, s_end, q_start, q_end, score,
                               percent_ident);
      } else {
        return makeAlignMatch (s_start, s_end, q_end, q_start, score,
                               percent_ident);
      }
    } catch (NumberFormatException e) {
      throw new IOException ("while reading SSAHA data: " +
                             "failed to parse a number from this string: " +
                             e.getMessage ());
    }
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
