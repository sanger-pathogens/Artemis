/* MarkerRangeRequesterEvent.java
 *
 * created: Mon Jul 10 2000
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/MarkerRangeRequesterEvent.java,v 1.2 2004-10-29 09:36:24 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.sequence.*;
import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.io.LocationLexer;
import uk.ac.sanger.artemis.io.LocationLexer.*;
import uk.ac.sanger.artemis.util.*;

import javax.swing.*;

/**
 *  This event is sent when the user presses OK or Cancel in a
 *  MarkerRangeRequester component.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: MarkerRangeRequesterEvent.java,v 1.2 2004-10-29 09:36:24 tjc Exp $
 **/

public class MarkerRangeRequesterEvent extends TextRequesterEvent {
  /**
   *  Create a new MarkerRangeRequesterEvent object.
   *  @param source The MarkerRangeRequester that generated the event.
   *  @param requester_text The contents of the TextField in the
   *    MarkerRangeRequester.
   *  @param type The type of event.
   **/
  public MarkerRangeRequesterEvent(final MarkerRangeRequester source,
                                   final String requester_text,
                                   final int type)
  {
    super(source, requester_text, type);
  }


  /**
   *  Parse a Range and return the start and end positions in an array that is
   *  two elements long.
   *  @return null if and only if there was a parse error.  If the range is on
   *    the forward strand the first element will be less than or equal to the
   *    second, otherwise the first will be greater than the second.
   **/
  private int[] getRangeInternal() 
  {
    if(getRequesterText().length() == 0) 
      return null;

    final LocationLexer lexer = new LocationLexer(getRequesterText());

    final TokenEnumeration enumTk = lexer.getTokens();

    boolean complement_flag = false;

    if(enumTk.peekElement() instanceof String &&
        ((String)enumTk.peekElement()).equals("complement")) 
    {
      complement_flag = true;
      enumTk.nextElement();
    }

    enumTk.eatToken('(');

    if(enumTk.peekElement() instanceof Integer) 
    {
      int first_base = ((Integer)enumTk.nextElement()).intValue();

      if(enumTk.peekElement() instanceof Integer ||
         (enumTk.eatToken("..") || enumTk.eatToken('.') ||
          enumTk.eatToken("-")) &&
          enumTk.peekElement() instanceof Integer)
      {
        int last_base = ((Integer)enumTk.nextElement()).intValue();

        enumTk.eatToken(')');

        if(enumTk.peekElement() == null) 
        {
          if(complement_flag) 
          {
            final int temp = first_base;
            first_base = last_base;
            last_base = temp;
          }

          return new int[]
          {
            first_base, last_base
          };
        }
        else
        {
          new MessageDialog((JFrame) getSource(),
                            "garbage at the end of the range: " + enumTk);
          return null;
        }
      }
    }

    // if we get to here then there was a parse error
    new MessageDialog((JFrame) getSource(),
                      "error in range: the range should have this " +
                      "form: 100..200 - error at: " + enumTk);

    return null;
  }

  /**
   *  Parse the contents of the String that was passed to the constructor as a
   *  Range and return it.
   *  @param bases The Bases object to use to find the Strand that is passed
   *    to the MarkerRange constructor.
   *  @return The Range or null if the Range can't be parsed.
   **/
  public MarkerRange getMarkerRange(final Bases bases) 
  {
    try 
    {
      final MarkerRange marker_range;

      final int [] return_values = getRangeInternal();

      if(return_values == null) 
        return null;

      final int first_base = return_values[0];
      final int last_base  = return_values[1];

      if(first_base <= last_base) 
      {
        // forward strand
        final Strand strand;
        strand = bases.getForwardStrand();

        marker_range =
          strand.makeMarkerRangeFromPositions(first_base,
                                              last_base);
      }
      else
      {
        // reverse strand
        final Strand strand = bases.getReverseStrand();
        final int raw_first = bases.getComplementPosition(first_base);
        final int raw_last  = bases.getComplementPosition(last_base);

        marker_range =
          strand.makeMarkerRangeFromPositions(raw_last,
                                              raw_first);
      }

      return marker_range;
    } 
    catch(OutOfRangeException e)
    {
      new MessageDialog((JFrame)getSource(),
                         "the bases are out of range for this " +
                         "sequence");
      return null;
    }
  }

  /**
   *  Parse and return a raw Range object from the text.
   *  @return The Range or null if the range could not be parsed.
   **/
  public Range getRawRange()
  {
    final int [] return_values = getRangeInternal();

    if(return_values == null) 
      return null;

    final int first_base = return_values[0];
    final int last_base  = return_values[1];

    try 
    {
      if(first_base < last_base) 
        return new Range(first_base, last_base);
      else 
        return new Range(last_base, first_base);
    } 
    catch(OutOfRangeException e) 
    {
      throw new Error("internal error - unexpected exception: " + e);
    }
  }
}

