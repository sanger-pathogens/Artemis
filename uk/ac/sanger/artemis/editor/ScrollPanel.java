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
import java.awt.LayoutManager;
import java.awt.Rectangle;
import java.awt.Dimension;
import java.awt.Color;

/**
*
* Extends JPanel to implement Scrollable to speed scroll pane
* scrolling
*
*/
public class ScrollPanel extends JPanel implements Scrollable
{

  public ScrollPanel()
  {
    super();
  }

  public ScrollPanel(LayoutManager l)
  {
    super(l);
  }

  public Dimension getPreferredScrollableViewportSize()
  {
    return getPreferredSize();
  }

  public boolean getScrollableTracksViewportHeight()
  {
    return false;
  }

  public boolean getScrollableTracksViewportWidth()
  {
    return false;
  }

  public int getScrollableBlockIncrement(Rectangle r,
                    int orientation, int direction)
  {
    return 60;
  }

  public int getScrollableUnitIncrement(Rectangle r,
                    int orientation, int direction)
  {
    return 10;
  }


}

