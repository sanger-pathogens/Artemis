/* SegmentBorder.java
 *
 * created: Fri Nov 19 2004
 *
 * This file is part of Artemis
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

package uk.ac.sanger.artemis.components;

import java.awt.*;

/**
 *
 **/

public class SegmentBorder
{

  private boolean highlight_feature;
  private boolean highlight_segment;
  private boolean draw_arrow;
  private int x;
  private int y;
  private int width;
  private int feature_direction;
  public static Color HIGHLIGHT_BORDER_COLOUR = new Color(140,25,25);
 
  /**
  * Information stored in this object is used to draw the feature
  * segment borders.
  * @param highlight_feature  true if the feature is highlighted
  * @param highlight_segment  true if the segment is highlighted
  * @param draw_arrow         true if an arrow is to be drawn on this segment
  * @param x                  top right hand x position of the segment
  * @param y                  top right hand y position of the segment
  * @param width              segment width
  * @param feature_direction  feature direction 
  *
  */
  public SegmentBorder(final boolean highlight_feature, 
                       final boolean highlight_segment, 
                       final boolean draw_arrow, int x, int y, int width,
                       int feature_direction)
  {
    this.highlight_feature = highlight_feature;
    this.highlight_segment = highlight_segment;
    this.draw_arrow = draw_arrow;

    this.x = x;
    this.y = y;
    this.width  = width;

    this.feature_direction = feature_direction;
  }


  protected void drawSegmentBorder(Graphics g, int height, int arrowWidth)
  {
    final Graphics2D g2d = (Graphics2D)g;
    if(highlight_feature)  // highlight selected features
    {
      // selected - highlight by drawing a thicker line
      final BasicStroke stroke = (BasicStroke)g2d.getStroke();

      if(highlight_segment)
      {
        g2d.setColor(HIGHLIGHT_BORDER_COLOUR);
        g2d.setStroke(new BasicStroke(4.f));
      }
      else
        g2d.setStroke(new BasicStroke(3.f));

      g2d.drawRect(x, y, width, height);
      g2d.setColor(Color.black);
      g2d.setStroke(stroke);
    }
    else
      g.drawRect(x, y, width, height);

    // draw the arrow point
    if(draw_arrow)
    {
      int xpos = x;
      int arrow_tip_x = x + feature_direction * arrowWidth;
      if(feature_direction ==1)
      {
        xpos += width;
        arrow_tip_x += width;
      }

      final int arrow_tip_y = y + (height/2);

      g.drawLine(xpos, y, arrow_tip_x, arrow_tip_y);
      g.drawLine(arrow_tip_x, arrow_tip_y, xpos, y+height);
    }
  }

}
