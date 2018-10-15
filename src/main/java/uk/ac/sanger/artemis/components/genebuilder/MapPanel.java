/* MapPanel.java
 *
 * created: 2008
 *
 * This file is part of Artemis
 *
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
 **/

package uk.ac.sanger.artemis.components.genebuilder;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.FontMetrics;
import java.awt.GradientPaint;
import java.awt.Graphics2D;
import java.awt.geom.RoundRectangle2D;

import javax.swing.JPanel;

import uk.ac.sanger.artemis.Selection;
import uk.ac.sanger.artemis.io.ChadoCanonicalGene;


public class MapPanel extends JPanel
{
  private static final long serialVersionUID = 1L;
  protected static int border = 15;
  protected ChadoCanonicalGene chado_gene;
  protected Selection selection;
  
  /**
   * Draw rectangular box for a feature.
   * @param g2d
   * @param start   start of feature
   * @param end     end of feature
   * @param ypos    y position
   * @param colour  feature colour
   * @param size    parameter to control the height of the feature
   */
  protected static void drawFeature(Graphics2D g2d, int start, int end, 
                           int ypos, Color colour, float size,
                           boolean selected, float selected_size,
                           int fontHeight)
  {
    drawFeature(g2d, start, end, ypos, colour, Color.BLACK, size,
        selected, selected_size, fontHeight);
  }
  
  /**
   * Draw rectangular box for a feature.
   * @param g2d
   * @param start   start of feature
   * @param end     end of feature
   * @param ypos    y position
   * @param colour  feature colour
   * @param size    parameter to control the height of the feature
   */
  protected static void drawFeature(Graphics2D g2d, int start, int end, 
                           int ypos, Color colour, Color borderColor, float size,
                           boolean selected, float selected_size,
                           int fontHeight)
  {
    RoundRectangle2D e = new RoundRectangle2D.Float(start, ypos, 
        end-start,
        fontHeight*size, 0, ypos);

    if(colour == null)
      colour = Color.BLACK;
    
    GradientPaint gp = new GradientPaint(start, ypos, 
        colour,
        start, ypos+( (fontHeight/2) * size ), 
        Color.white, true);
    g2d.setPaint(gp); 
    g2d.fill(e);
    
    if(selected)
      g2d.setStroke(new BasicStroke(selected_size));
    else
      g2d.setStroke(new BasicStroke(1.f));
    
    // draw boundary
    g2d.setColor(borderColor);
    g2d.draw(e);
  }
  
  protected int getFontHeight()
  {
    final FontMetrics fm = this.getFontMetrics(getFont());
    return fm.getHeight();  
  }
  
  /**
   * Macro for getting the size of the transcipt and
   * exon image.
   * @return
   */
  protected int getTranscriptSize()
  {
    return (2 * border) + (getFontHeight() * 3);  
  }
  
  protected int getViewerBorder()
  {
    return border; 
  }
}