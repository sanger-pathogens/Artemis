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
  private int height;
  private int feature_direction;
 
  /**
  *
  * Information stored in this object is used to draw the feature
  * segment borders.
  *
  * @param highlight_feature  true if the feature is highlighted
  * @param highlight_segment  true if the segment is highlighted
  * @param draw_arrow         true if an arrow is to be drawn on this segment
  * @param x                  top right hand x position of the segment
  * @param y                  top right hand y position of the segment
  * @param width              segment width
  * @param height             segment height
  * @param feature_direction  feature direction 
  *
  */
  public SegmentBorder(boolean highlight_feature, boolean highlight_segment, 
                       boolean draw_arrow, int x, int y, int width, int height,
                       int feature_direction)
  {
    this.highlight_feature = highlight_feature;
    this.highlight_segment = highlight_segment;
    this.draw_arrow = draw_arrow;

    this.x = x;
    this.y = y;
    this.width  = width;
    this.height = height;

    this.feature_direction = feature_direction;
  }

 
  protected int getXPoint()
  {
    return x;
  }

  protected int getYPoint()
  {
    return y;
  }

  protected int getWidth()
  {
    return width;
  }

  protected int getHeight()
  {
    return height;
  }

  protected int getFeatureDirection()
  {
    return feature_direction;
  }

  protected boolean isFeatureHighlight()
  {
    return highlight_feature;
  }

  protected boolean isSegmentHighlight()
  {
    return highlight_segment;
  }

  protected boolean showArrow()
  {
    return draw_arrow;
  }

}
