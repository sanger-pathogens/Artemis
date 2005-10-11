/* UserDataAlgorithm.java
 *
 * created: Wed May 10 2000
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/plot/UserDataAlgorithm.java,v 1.3 2005-10-11 14:20:31 tjc Exp $
 */

package uk.ac.sanger.artemis.plot;

import uk.ac.sanger.artemis.sequence.*;

import uk.ac.sanger.artemis.util.*;
import uk.ac.sanger.artemis.io.ReadFormatException;

import java.io.*;

/**
 *  Objects of this class have one useful method - getValues (), which takes a
 *  range of bases and returns a single floating point number.  The number is
 *  calculated by averaging the values from a data file.  The Strand to use is
 *  set in the constructor.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: UserDataAlgorithm.java,v 1.3 2005-10-11 14:20:31 tjc Exp $
 **/

public class UserDataAlgorithm extends BaseAlgorithm {
  /**
   *  Create a new UserDataAlgorithm object.
   *  @param strand The strand to do the calculation on.
   *  @param document The Document to read the data from.
   **/
  public UserDataAlgorithm (final Strand strand, final Document document)
      throws IOException {
    super (strand, "User algorithm from " + document.getName (), "user");

    final Reader document_reader = document.getReader ();

    LinePushBackReader pushback_reader =
      new LinePushBackReader (document_reader);

    final String first_line = pushback_reader.readLine ();

    final StringVector tokens = StringVector.getStrings (first_line, " ");

    if (tokens.size () < 1) {
      throw new ReadFormatException ("unknown file type");
    }

    this.number_of_values = tokens.size ();

    pushback_reader.pushBack (first_line);

    data = new float [strand.getSequenceLength ()][tokens.size ()];

    readData (pushback_reader);
  }

  /**
   *  Read all from buffered_reader into data.
   **/
  private void readData (final LinePushBackReader pushback_reader)
      throws IOException {
    String line = null;

    int count = 0;

    while ((line = pushback_reader.readLine ()) != null) {

      if (count >= getStrand ().getSequenceLength ()) {
        throw new ReadFormatException ("too many values in input file");
      }

      final StringVector tokens = StringVector.getStrings (line, " ");

      final int read_base;

      if (tokens.size () == data[0].length) {
        for (int i = 0 ; i < tokens.size () ; ++i) {
          try {
            final float value =
              Float.valueOf ((String)tokens.elementAt (i)).floatValue ();

            if (value > data_max) {
              data_max = value;
            }

            if (value < data_min) {
              data_min = value;
            }

            data[count][i] = value;

            average_value += value;
          } catch (NumberFormatException e) {
            throw new ReadFormatException ("cannot understand this number: " +
                                           tokens.elementAt (i) + " - " +
                                           e.getMessage ());
          }
        }
      } else {
        throw new ReadFormatException ("line has the wrong number of fields");
      }

      ++count;
    }

    average_value /= data[0].length * getStrand ().getSequenceLength ();
  }

  /**
   *  Return the value of the function between a pair of bases.
   *  @param start The start base (included in the range).
   *  @param end The end base (included in the range).
   *  @param values The one return value for this algorithm is returned in
   *    this array.
   **/
  public void getValues (int start, int end, final float [] values) {
    final int value_count = getValueCount ();

    for (int i = 0 ; i < value_count ; ++i) {
      values [i] = 0;
      for (int base = start ; base <= end ; ++base) {
        values [i] += data[base - 1][i] / (end - start + 1);
      }
    }
  }

  /**
   *  Return the number of values a call to getValues () will return - one
   *  in this case.
   **/
  public int getValueCount () {
    return number_of_values;
  }

  /**
   *  Return the default or optimal window size.
   *  @return null is returned if this algorithm doesn't have optimal window
   *    size.
   **/
  public Integer getDefaultWindowSize () {
    return new Integer (3);
  }

  /**
   *  Return the default maximum window size for this algorithm.
   *  @return null is returned if this algorithm doesn't have maximum window
   *    size.
   **/
  public Integer getDefaultMaxWindowSize () {
    return new Integer (100);
  }

  /**
   *  Return the default minimum window size for this algorithm.
   *  @return null is returned if this algorithm doesn't have minimum window
   *    size.
   **/
  public Integer getDefaultMinWindowSize () {
    return new Integer (1);
  }

  /**
   *  Return the default or optimal step size.
   *  @return null is returned if this algorithm doesn't have optimal step
   *    size.
   **/
  public Integer getDefaultStepSize (int window_size) {
    if (window_size > 10) {
      return new Integer (window_size / 10);
    } else {
      return null;
    }
  }

  /**
   *  Return the maximum value of this algorithm.
   **/
  protected Float getMaximumInternal () {
    return new Float (data_max);
  }

  /**
   *  Return the minimum value of this algorithm.
   **/
  protected Float getMinimumInternal () {
    return new Float (data_min);
  }

  /**
   *  Return the average value of function over the whole strand.
   **/
  public Float getAverage () {
    return new Float (average_value);
  }

  /**
   *  The data that was read by the constructor.
   **/
  private float data[][] = null;

  /**
   *  The maximum value in the data array.
   **/
  private float data_max = -9999999;

  /**
   *  The minimum value in the data array.
   **/
  private float data_min = 9999999;

  /**
   *  The average calculated by readData ().
   **/
  private float average_value = 0;

  /**
   *  The value returned by getValueCount ().
   **/
  private int number_of_values;
}
