/* MegaBlastComparisonData.java
 *
 * created: Tue Apr 30 2002
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/MegaBlastComparisonData.java,v 1.1 2004-06-09 09:44:57 tjc Exp $
 */

package uk.ac.sanger.artemis;

import uk.ac.sanger.artemis.util.LinePushBackReader;

import java.io.*;
import java.util.StringTokenizer;

/**
 *  This class implements the ComparisonData interface for MegaBlast output.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: MegaBlastComparisonData.java,v 1.1 2004-06-09 09:44:57 tjc Exp $
 **/

public class MegaBlastComparisonData extends SimpleComparisonData
    implements ComparisonData {
  /**
   *  Create a new MegaBlastComparisonData by reading from the given
   *  LinePushBackReader.
   **/
  public MegaBlastComparisonData (final LinePushBackReader stream)
      throws IOException {
    super (stream);
  }

  /**
   *  Create a new, empty instance of MegaBlastComparisonData.
   **/
  public MegaBlastComparisonData () {
    
  }

  /**
   *  Returns a new, empty instance of this type of object;
   **/
  protected SimpleComparisonData getNewSimpleComparisonData () {
    return new MegaBlastComparisonData ();
  }


  /**
   *  Make an AlignMatch object from the given String.
   **/
  private static AlignMatch makeMatchFromStringStatic (final String line)
      throws IOException {
    final StringTokenizer tokenizer = new StringTokenizer (line, " ");

    if (tokenizer.countTokens () != 6) {
      final String message = "while reading MegaBlast data: " +
        "not enough fields in this line: " + line;
      throw new ComparisonDataParseException (message);
    }

    // first token looks like: 'all'=='-selected' - just check that we have
    // something that looks sensible
    final String file_names = tokenizer.nextToken ();

    if (!file_names.startsWith ("'") || !file_names.endsWith ("'") ||
        file_names.indexOf ("==") == -1) {
      final String message = "while reading MegaBlast data: " +
        "first field (" + file_names + ") is badly formatted in this line: " +
        line;
      throw new ComparisonDataParseException (message);
    }

    String s_start_token = tokenizer.nextToken ();
    final String q_start_token = tokenizer.nextToken ();
    final String s_end_token = tokenizer.nextToken ();
    String q_end_token = tokenizer.nextToken ();

    if (s_start_token.startsWith ("(")) {
      s_start_token = s_start_token.substring (1);
    } else {
      final String message = "while reading MegaBlast data: " +
        "second field " +
        s_start_token + " " + q_start_token + " " +
        s_end_token + " " + q_end_token +
        " is badly formatted in this line: " + line;
      throw new ComparisonDataParseException (message);
    }

    if (q_end_token.endsWith (")")) {
      q_end_token = q_end_token.substring (0, q_end_token.length () - 1);
    } else {
      final String message = "while reading MegaBlast data: " +
        "second field (" +
        s_start_token + " " + q_start_token + " " +
        s_end_token + " " + q_end_token +
        " is badly formatted in this line: " + line;
      throw new ComparisonDataParseException (message);
    }

    // dunno what it is - throw it away
    tokenizer.nextToken ();

    try {
      final int q_start = Integer.valueOf (q_start_token).intValue ();
      final int q_end   = Integer.valueOf (q_end_token).intValue ();
      final int s_start = Integer.valueOf (s_start_token).intValue ();
      final int s_end   = Integer.valueOf (s_end_token).intValue ();

      return makeAlignMatch (s_start, s_end, q_start, q_end, 100,
                             100);
    } catch (NumberFormatException e) {
      throw new IOException ("while reading blast -m 8 data: " +
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
