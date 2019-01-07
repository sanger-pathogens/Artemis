/*
 *
 * created: Wed Sep 7 2004
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
import java.awt.*;
import java.awt.event.*;

/**
*
* Extend JButton to show mouse over colour change.
*
*/
public class MouseOverButton extends JButton       
{
  /** */
  private static final long serialVersionUID = 1L;
  private boolean over = false;
  private HitInfo hit;

  public MouseOverButton()
  {
    super();
  }

  public MouseOverButton(String s)
  {
    super(s);
  }

  public MouseOverButton(HitInfo hit)
  {
    super(hit.getID());
    this.hit = hit;
  }

  public void paintComponent(Graphics g)
  {
    super.paintComponent(g);

    if(!getText().equals(""))
      return;

    Graphics2D g2 = (Graphics2D)g;

    if(over)
      g2.setColor(new Color(100,100,200));
    else
      g2.setColor(Color.blue);
    g2.setStroke(new BasicStroke(1.5f));
    g2.drawLine(1,3,11,3);
    g2.drawLine(1,7,11,7);

    setSize(12,12);
  }

  public String getToolTipText()
  {
    if(hit == null)
      return null;

    return hit.getOrganism();
  }

  protected void processMouseEvent(MouseEvent evt)
  {
    super.processMouseEvent(evt);
    switch (evt.getID())
    {
      case MouseEvent.MOUSE_ENTERED:
        over = true; 
        setForeground(new Color(100,100,200));
	repaint();
	break;
      case MouseEvent.MOUSE_EXITED:
        over = false;
        setForeground(Color.blue);
        repaint();
        break;
    }
  }
}

