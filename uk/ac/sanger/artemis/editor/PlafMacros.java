/*
 *
 * created: Wed Aug 3 2004
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2000  Genome Research Limited
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or(at your option) any later version.
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

package uk.ac.sanger.artemis.editor;

import javax.swing.*;
import java.util.*;
import java.io.*;
import java.lang.*;
import java.awt.*;

/**
*
* Macros class
*
*/
public class PlafMacros implements SwingConstants 
{
  // don't make these final, since the value is 
  // different on each platform
  private static String LINE_SEPARATOR = 
                  System.getProperty("line.separator");
  private static int LINE_SEPARATOR_LEN = 
                  LINE_SEPARATOR.length();

  /**
  *
  * Break text into separate lines
  * @param text		text to break up
  * @return		multiple lines in string array
  *
  */
  public static String[] breakupLines(String text) 
  {
    int len = text.length();
    if (len == 0)
      return new String[] {""}; 
    else 
    {
      Vector data = new Vector(10);
      int start=0;
      int i=0;
      while (i<len) {
        if (text.startsWith(LINE_SEPARATOR,i)) 
        {
          data.addElement(text.substring(start,i));
          start=i+LINE_SEPARATOR_LEN;
          i=start;
        } 
        else if (text.charAt(i)=='\n')
        {
          data.addElement(text.substring(start,i));
          start=i+1;
          i=start;
        } 
        else { i++; }
      }
      if (start != len)
        data.addElement(text.substring(start));
      int numlines = data.size();
      String lines[] = new String[numlines];
      data.copyInto(lines);
      return lines;
    }
  }

  /**
  *
  * Get the line separator string
  * @return 	line separator
  *
  */
  public static String getLineSeparator()
  {
    return LINE_SEPARATOR;
  }
}

