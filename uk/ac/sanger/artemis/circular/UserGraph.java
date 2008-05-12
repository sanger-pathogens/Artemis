/*
 * Copyright (C) 2008  Genome Research Limited
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
 *  @author: Tim Carver
 */


package uk.ac.sanger.artemis.circular;

import java.io.IOException;
import java.io.Reader;

import uk.ac.sanger.artemis.io.ReadFormatException;
import uk.ac.sanger.artemis.util.Document;
import uk.ac.sanger.artemis.util.LinePushBackReader;
import uk.ac.sanger.artemis.util.StringVector;

public class UserGraph extends Graph
{
  private static final long serialVersionUID = 1L;

  /**
   *  The data that was read by the constructor.
   **/
  private float data[] = null;

  /**
   *  The maximum value in the data array.
   **/
  private float data_max = Float.MIN_VALUE;

  /**
   *  The minimum value in the data array.
   **/
  private float data_min = Float.MAX_VALUE;

  /**
   *  The average calculated by readData ().
   **/
  private float average_value = 0;

  
  public UserGraph(DNADraw currentDna, final Document document)
         throws IOException
  {
    super(currentDna);
    
    final Reader document_reader = document.getReader();
    LinePushBackReader pushback_reader =
      new LinePushBackReader(document_reader);

    final String first_line = pushback_reader.readLine ();
    final StringVector tokens = StringVector.getStrings (first_line, " ");

    if (tokens.size () < 1)
      throw new ReadFormatException ("unknown file type");

    pushback_reader.pushBack (first_line);
    data = new float [getBases().getLength()];
    readData (pushback_reader);
  }
 
 /**
  * Return the value between the given pair of bases.
  */
  protected float calculateValue(int start, int end)
  {
    float value = 0.f;

    for (int base = start ; base <= end ; ++base) 
      value += data[base - 1] / (end - start + 1);
      
    return value;
  }
  
  /**
   *  Read all from buffered_reader into data.
   **/
  private void readData (final LinePushBackReader pushback_reader)
      throws IOException
  {
    String line = null;
    int count = 0;

    while ((line = pushback_reader.readLine ()) != null) 
    {
      if(count >= getBases().getLength()) 
        throw new ReadFormatException ("too many values in input file");

      final StringVector tokens = StringVector.getStrings (line, " ");

      if(tokens.size () == 1) 
      {
        try 
        {
          final float value =
              Float.valueOf ((String)tokens.elementAt (0)).floatValue ();

          if (value > data_max)
            data_max = value;

          if (value < data_min)
            data_min = value;

          data[count] = value;
          average_value += value;
        } 
        catch (NumberFormatException e) 
        {
          throw new ReadFormatException ("cannot understand this number: " +
                                         tokens.elementAt (0) + " - " +
                                         e.getMessage ());
        }
      } 
      else 
      {
        throw new ReadFormatException ("line has the wrong number of fields");
      }
      ++count;
    }

    average_value /= getBases().getLength();
  }
}