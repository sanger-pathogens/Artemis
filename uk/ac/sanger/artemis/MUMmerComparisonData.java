/* MUMmerComparisonData.java
 *
 * created: Thu Jul 15 1999
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
 */

package uk.ac.sanger.artemis;

import uk.ac.sanger.artemis.util.LinePushBackReader;

import java.io.*;
import java.util.StringTokenizer;
import java.util.Vector;

/**
 *  This class implements the ComparisonData interface for MUMmer output.
 *
 *  @author Kim Rutherford
 *  @version $Id: MUMmerComparisonData.java,v 1.1 2006-11-07 17:29:30 tjc Exp $
 **/

public class MUMmerComparisonData extends SimpleComparisonData
    implements ComparisonData 
{
  /**
   *  Create a new MUMmerComparisonData object by reading from the given
   *  LinePushBackReader.
   **/
  public MUMmerComparisonData (final LinePushBackReader stream)
      throws IOException 
  {
    super (stream);
  }
  
  public MUMmerComparisonData ()
  {
    
  }
  
  /**
   *  Make an AlignMatch object from the given String.
   **/
  protected AlignMatch makeMatchFromString(final String line)
      throws IOException 
  {
    if(line == null || line.startsWith(">"))
      return null;

    final StringTokenizer tokenizer = new StringTokenizer(line, " ");

    if(tokenizer.countTokens() < 3)
    {
      final String message = "unable to understand this line: " + line;
      throw new ComparisonDataParseException(message);
    }

    final String first_token = tokenizer.nextToken();
    final String second_token = tokenizer.nextToken();
    final String third_token = tokenizer.nextToken();

    String forth_token = "+";

    if(tokenizer.countTokens() > 0)
      forth_token = tokenizer.nextToken();

    try
    {
      final int first_number = Integer.valueOf(first_token).intValue();
      final int second_number = Integer.valueOf(second_token).intValue();
      final int match_length = Integer.valueOf(third_token).intValue();

      //        System.err.println ("token: " + forth_token + "  " + line);

      int score = -1;

      if((forth_token.equals("+") || forth_token.equals("-"))
          && tokenizer.countTokens() > 0)
      {
        score = Integer.valueOf(tokenizer.nextToken()).intValue();

        if(score < -1)
          score = -1;
      }

      if(forth_token.equals("-"))
      {
        return makeAlignMatch(first_number, first_number + match_length - 1,
            second_number + match_length - 1, second_number, score, -1);
      }
      else
      {
        return makeAlignMatch(first_number, first_number + match_length - 1,
            second_number, second_number + match_length - 1, score, -1);
      }
    }
    catch(NumberFormatException e)
    {
      return null;
/*      throw new IOException("while reading MUMmer file: failed to "
          + "parse a number from this string: " + e.getMessage());*/
    }
  }


  protected SimpleComparisonData getNewSimpleComparisonData()
  {
    return new MUMmerComparisonData();
  }



}
