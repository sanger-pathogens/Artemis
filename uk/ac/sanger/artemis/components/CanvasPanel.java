/* CanvasPanel.java
 *
 * created: Sat Jun 17 2000
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/CanvasPanel.java,v 1.1 2004-06-09 09:46:07 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.Options;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

/**
 *  This is a JPanel that contains a JPanel containing a JComponent.  Both Panels
 *  have BorderLayout.  The JComponent is added at "Center".
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: CanvasPanel.java,v 1.1 2004-06-09 09:46:07 tjc Exp $
 **/

abstract public class CanvasPanel extends JPanel 
{

  /**
   *  Off screen image used for double buffering when drawing the canvas.
   **/
  private Image offscreen;

  /** Contains the canvas. */
  private JPanel mid_panel = null;

  /** The drawing area for this component. */
  private JComponent canvas = null;

  /** The height of the font used in this component. */
  private int font_ascent;

  /** maximum height of the font used in this component. */
  private int font_max_ascent;

  /** descent of the font used in this component. */
  private int font_descent;

  /** The(maximum) width of the font used in this component. */
  private int font_width;

  /** base line of the font used in this component. */
  private int font_base_line;

  /**
   *  Create a new JPanel(mid_panel) and a JComponent.
   **/
  public CanvasPanel() 
  {
    setLayout(new BorderLayout());
    setFontInfo();
    createCanvas();
  }

  /**
   *  Call repaint() on the canvas object.
   **/
  protected void repaintCanvas() 
  {
    getCanvas().repaint();
  }

  /**
   *  Create a JPanel(mid_panel) and a JComponent object.
   **/
  private void createCanvas() 
  {
    mid_panel = new JPanel();
    mid_panel.setLayout(new BorderLayout());

    canvas = new JComponent() 
    {
      /**
       *  Set the offscreen buffer to null as part of invalidation.
       **/
      public void invalidate() 
      {
        super.invalidate();
        offscreen = null;
      }

      /**
       *  Override update to *not* erase the background before painting
       */
      public void update(final Graphics g) 
      {
        paint(g);
      }

      /**
       *  Paint the canvas.
       */
      public void paint(final Graphics g) 
      {
        final int canvas_width = canvas.getSize().width;
        final int canvas_height = canvas.getSize().height;

        // there is no point drawing into a zero width canvas
        if(canvas_height <= 0 || canvas_width <= 0) 
          return;

        if(!canvas.isVisible()) 
          return;

        if(offscreen == null) 
          offscreen = canvas.createImage(canvas_width, canvas_height);

        Graphics og = offscreen.getGraphics();
        og.setClip(0, 0, canvas_width, canvas_height);

        paintCanvas(og);
        g.drawImage(offscreen, 0, 0, null);
        g.dispose();
      }
    };

    mid_panel.add(canvas, "Center");
    add(mid_panel, "Center");
  }

  /**
   *  Return the JComponent that was created by createCanvas().
   **/
  protected JComponent getCanvas() 
  {
    return canvas;
  }

  /**
   *  Returns the sub-JPanel that contains the JComponent.
   **/
  protected JPanel getMidPanel() 
  {
    return mid_panel;
  }

  /**
   *  Called by canvas.paint() when
   **/
  abstract protected void paintCanvas(final Graphics graphics);

  /**
   *  Set font_width and font_ascent from the default font.
   **/
  private void setFontInfo() 
  {
    FontMetrics fm = getFontMetrics(getFont());

    // find the width of a wide character
    font_width = fm.charWidth('M');
    font_ascent = fm.getAscent();
    font_max_ascent = fm.getMaxAscent();
    font_descent = fm.getDescent();
  }

  /**
   *  Return the width of the canvas.
   **/
  public int getCanvasWidth() 
  {
    return getCanvas().getSize().width;
  }

  /**
   *  Return the height of the canvas.
   **/
  public int getCanvasHeight() 
  {
    return getCanvas().getSize().height;
  }

  /**
   *  Return the width of our font, as calculated by setFontInfo().
   **/
  public int getFontWidth() 
  {
    return font_width;
  }

  /**
   *  Return the ascent(height above the baseline) of our font, as calculated
   *  by setFontInfo().
   **/
  public int getFontAscent() 
  {
    return font_ascent;
  }

  /**
   *  Return the max ascent(height above the baseline) of our font, as
   *  calculated by setFontInfo().
   **/
  public int getFontMaxAscent() 
  {
    return font_ascent;
  }

  /**
   *  Return the descent of our font, as calculated by setFontInfo().
   **/
  public int getFontDescent() 
  {
    return font_descent;
  }

  /**
   *  The max ascent + descent of the default font.
   **/
  public int getFontHeight() 
  {
    return getFontMaxAscent() + getFontDescent();
  }

}
