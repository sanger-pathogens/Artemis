/* SequenceViewer.java
 *
 * created: Sat Dec 19 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/SequenceViewer.java,v 1.1 2004-06-09 09:47:47 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import java.io.StringWriter;
import java.io.PrintWriter;

/**
 *  This component provides a viewer for dna or amino acid sequences.  The
 *  units are numbered automatically, given a view like this:
 *  <pre>
 *  ATGATGATGATGATATGCATGATCG
 *           10        20
 *  </pre>
 *  @author Kim Rutherford
 *  @version $Id: SequenceViewer.java,v 1.1 2004-06-09 09:47:47 tjc Exp $
 **/

public class SequenceViewer extends FileViewer {
  /**
   *  Create a new SequenceViewer component.
   *  @param title The name to attach to the new JFrame.
   *  @param include_numbers If true then the sequence will be numbered
   *    (every second line of the display will be numbers rather than
   *    sequence).
   **/
  public SequenceViewer (final String title,
                         final boolean include_numbers) {
    super (title);

    this.include_numbers = include_numbers;
  }

  /**
   *  Set the sequence to display.
   *  @param comment_line A comment to put on the first line.  This may be
   *    null in which case there is no comment line.  This String should not
   *    be "\n" terminated.
   *  @param sequence The sequence to display.  This may be null in which case
   *    there is no comment line.  This String should not be "\n" terminated.
   **/
  public void setSequence (final String comment_line, final String sequence) {
    this.comment_line = comment_line;
    this.sequence = sequence;

    setView ();
  }

  /**
   *  Set the sequence to display with no comment line.
   *  @param sequence The sequence to display.  This may be null in which case
   *    there is no comment line.  This String should not be "\n" terminated.
   **/
  public void setSequence (final String sequence) {
    this.comment_line = null;
    this.sequence = sequence;

    setView ();
  }

  /**
   *  Write the sequence and comment_line in the view.
   **/
  private void setView () {
    final StringWriter string_writer = new StringWriter ();

    final PrintWriter writer = new PrintWriter (string_writer);

    if (comment_line != null) {
      writer.println (comment_line);
    }

    String rest_of_sequence = sequence;

    final int UNITS_PER_LINE = 60;

    int line_count = 0;

    while (rest_of_sequence != null) {
      final String this_line;

      if (rest_of_sequence.length () < UNITS_PER_LINE) {
        this_line = rest_of_sequence;
        rest_of_sequence = null;
      } else {
        this_line = rest_of_sequence.substring (0, UNITS_PER_LINE);
        rest_of_sequence = rest_of_sequence.substring (UNITS_PER_LINE);
      }


      writer.println (this_line);

      if (include_numbers) {
        final int COUNT_INTERVAL = 10;

        for (int i = 1 ;
             i <= this_line.length () / COUNT_INTERVAL;
             ++i) {
          final int this_unit_count =
            i * COUNT_INTERVAL + line_count * UNITS_PER_LINE;

          final int number_of_spaces =
            COUNT_INTERVAL - String.valueOf (this_unit_count).length ();

          for (int space_index = 0 ;
               space_index < number_of_spaces ;
               ++space_index) {
            writer.print (' ');
          }

          writer.print (this_unit_count);
        }

        writer.println ();
      }

      ++line_count;
    }

    writer.flush ();

    setText (string_writer.toString ());
  }

  /**
   *  The sequence to display.
   **/
  private String sequence = "";

  /**
   *  A comment to put on the first line.  If null don't display any comment.
   **/
  private String comment_line = null;

  /**
   *  If true then the amino acids will be numbered (every second line of the
   *  display will be numbers rather than sequence).
   **/
  private final boolean include_numbers;
}


