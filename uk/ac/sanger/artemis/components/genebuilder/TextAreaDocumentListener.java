/* TextAreaDocumentListener
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2009  Genome Research Limited
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
package uk.ac.sanger.artemis.components.genebuilder;

import java.awt.Dimension;
import java.awt.FontMetrics;
import java.util.StringTokenizer;

import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;

import uk.ac.sanger.artemis.components.QualifierTextArea;


public class TextAreaDocumentListener implements DocumentListener
{
  private QualifierTextArea qta;

  public TextAreaDocumentListener(QualifierTextArea qta)
  {
    this.qta = qta;
    setQualifierTextAreaSize();
  }

  public void insertUpdate(DocumentEvent e)
  {
    updateSize(e);
  }

  public void removeUpdate(DocumentEvent e)
  {
    updateSize(e);
  }

  // Plain text components do not fire these events
  public void changedUpdate(DocumentEvent e)
  {
  }

  private void updateSize(DocumentEvent e)
  {
    setQualifierTextAreaSize();
  }

/*  private int getLineCount(Document doc)
  {
    // get last visible character's offset
    int lineEndOffset = 0;
    int line = 0;
    // go to end of each line until last character of last line is reached
    try
    {
      while (lineEndOffset < qta.getDocument().getEndPosition().getOffset() && line < 50)
      {
        lineEndOffset = Utilities.getRowEnd(qta, lineEndOffset+1);
        line++;
      }
    }
    catch (BadLocationException e){}

    if(line == 50)
    {
      Element root = doc.getDefaultRootElement();
      line = root.getElementCount();
    }
    return line;   
  }*/
  
  /**
   * Calculate the number of lines, taking into account line wrapping at
   * word boundaries (whitespace).
   * @param fm
   * @param text
   * @param width
   * @return
   */
  private int getNumberOfLines(final int width)
  {
    String text = qta.getText();
    String lines[] = text.split("\n");
    int lineCount = lines.length;
    for(int i=0; i<lines.length; i++)
      lineCount += getWrappedLineCount(lines[i], width);
    
    if(lineCount < 1)
      lineCount = 1;
    return lineCount;
  }

  /**
   * For a given line count how many times it is wrapped.
   * @param text
   * @param width
   * @return
   */
  private int getWrappedLineCount(final String text, final int width)
  {
    String delim = " \t\n\f\r";
    StringTokenizer tok = new StringTokenizer(text, delim, true);
    FontMetrics fm = qta.getFontMetrics(qta.getFont());

    int lineOffset = 0;
    int lineNumber = 0;
    while(tok.hasMoreTokens())
    {
      int thisWordLength = fm.stringWidth(tok.nextToken());
      lineOffset+=thisWordLength;
      if(lineOffset>width)
      {
        lineNumber++;
        lineOffset = thisWordLength;
      }
    }
    return lineNumber;
  }

  /**
   * Set the size from the number of lines.
   * @param doc
   */
  private void setQualifierTextAreaSize()
  {
    int lines = getNumberOfLines(qta.getPreferredSize().width);
    int lineHeight = qta.getFontMetrics(qta.getFont()).getHeight();

    qta.setPreferredSize(new Dimension(qta.getPreferredSize().width, 
                                       lineHeight * lines));
  }
}