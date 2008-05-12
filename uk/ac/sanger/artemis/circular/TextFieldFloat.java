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

import java.awt.Toolkit;
import java.text.NumberFormat;
import java.text.ParseException;
import java.util.Locale;
import javax.swing.text.AttributeSet;
import javax.swing.text.BadLocationException;
import javax.swing.text.Document;
import javax.swing.text.PlainDocument;

/**
*
* JTextfield for float fields in the EMBOSS form
*
*/
public class TextFieldFloat extends TextFieldSink 
{

  private Toolkit toolkit;
  private NumberFormat decimalFormatter;

  public TextFieldFloat() 
  {
    super();
    toolkit = Toolkit.getDefaultToolkit();
    decimalFormatter = NumberFormat.getNumberInstance(Locale.UK);
  }

  public double getValue() 
  {
    double retVal = 0;
    try 
    {
      retVal = decimalFormatter.parse(getText()).doubleValue();
    } 
    catch (ParseException e) 
    {
      // This should never happen because insertString allows
      // only properly formatted data to get in the field.
      toolkit.beep();
//    System.err.println("TextFieldFloat getValue: " + retVal);
    }
    return retVal;
  }

  public void setValue(double value) 
  {
    setText(decimalFormatter.format(value));
  }

  protected Document createDefaultModel() 
  {
    return new DecimalNumberDocument();
  }

  protected class DecimalNumberDocument extends PlainDocument 
  {
    public void insertString(int offs, 
                             String str,
                             AttributeSet a) 
              throws BadLocationException 
    {
      char[] source = str.toCharArray();
      char[] result = new char[source.length];
      int j = 0;

      for (int i = 0; i < result.length; i++) 
      {
        if (Character.isDigit(source[i]) ||
            source[i] == '.' || source[i] == '-')
          result[j++] = source[i];
        else 
        {
          if(source[i] != ',')
            toolkit.beep();
//        System.err.println("insertString: " + source[i]);
        }
      }
      super.insertString(offs, new String(result, 0, j), a);
    }
  }
}

