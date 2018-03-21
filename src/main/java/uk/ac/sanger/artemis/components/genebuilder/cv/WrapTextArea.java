/* WrapTextArea.java
 *
 * created: 2009
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2009  Genome Research Limited
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
 **/
package uk.ac.sanger.artemis.components.genebuilder.cv;

import java.awt.Dimension;
import java.awt.FontMetrics;
import java.util.StringTokenizer;

import javax.swing.JTextArea;

  /** 
   * Text component used for product & controlled_curation qualifier value.
   */
  class WrapTextArea extends JTextArea
  {
    private static final long serialVersionUID = 1L;
    private int labelWidth;
    
    public WrapTextArea(final String text, 
                        final Dimension go_dimension,
                        int width)
    {
      super(text);
      setOpaque(false);
      setEditable(false);
      setLineWrap(true);
      setWrapStyleWord(true);
      FontMetrics fm  = getFontMetrics(getFont());
      
      
      if(go_dimension != null)
        labelWidth = go_dimension.width;
      else
        labelWidth= fm.stringWidth("GO:0001234 [F] ");
      
      width = labelWidth+width;
      int rows = getNumberOfLines(fm, text, width);
      final Dimension d = new Dimension(width, (int) (getRowHeight()*rows) );
      setPreferredSize(d);
      setMaximumSize(d);
    }
    
    /**
     * Calculate the number of lines, taking into account line wrapping at
     * word boundaries (whitespace).
     * @param fm
     * @param text
     * @param width
     * @return
     */
    private int getNumberOfLines(FontMetrics fm, final String text, final int width)
    {
      String delim = " \t\n\f\r";
      StringTokenizer tok = new StringTokenizer(text, delim, true);
      //final String words[] = text.split("\\s");
      int lineOffset = 0;
      int lineNumber = 1;
      int w2 = (fm.stringWidth(" ")+1)/2;
     
      while(tok.hasMoreTokens())
      {
        int thisWordLength = fm.stringWidth(tok.nextToken());
        lineOffset+=thisWordLength;
        if(lineOffset>=width-w2)
        {
          lineNumber++;
          lineOffset = thisWordLength;
        }
      }

      /*int stringWidth = fm.stringWidth(text);
      int rows = Math.round((stringWidth/width)+.5f);*/
      return lineNumber;
    }
    
    protected int getLabelWidth()
    {
      return labelWidth;
    }
  }