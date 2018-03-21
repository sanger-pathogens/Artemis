/* ZoomScrollBar.java
 *
 * created: Fri Oct  9 1998
 *
 * This file is part of Artemis
 *
 * Copyright(C) 1998,1999,2000,2001,2002  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/ZoomScrollBar.java,v 1.2 2008-11-27 10:48:57 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import java.awt.Component;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionListener;

import javax.swing.JScrollBar;

class ZoomScrollBar extends JScrollBar
{
  private static final long serialVersionUID = 1L;
  private FeatureDisplay display;
  private MouseMotionListener mouseMotionListener = null;
  
  public ZoomScrollBar(final FeatureDisplay display)
  {
    super(JScrollBar.VERTICAL);
    this.display = display;
    
    //  try to arrange for the scrollbar to have a maximum value big enough
    // that the whole sequence can be visible at once
    final int MAX_FACTOR =
      (int)Math.round(Math.log(display.getSequenceLength()/20) /  Math.log(3));
    setValues(display.getScaleFactor(), 1, 0, MAX_FACTOR);
    setBlockIncrement(1);
    setUnitIncrement(1);
    addAdjustmentListener(new AdjustmentListener() 
    {
      public void adjustmentValueChanged(AdjustmentEvent e) 
      {
        display.setScaleFactor(e.getValue());
      }
    });

    if(display.getScaleFactor() >= MAX_FACTOR) 
    {
      display.setScaleFactor(MAX_FACTOR - 1);
      setValue(display.getScaleFactor());
    }
  }
  
  /**
   * Add MouseMotionListener to auto hide scroll
   */
  protected void addMouseMotionListenerToFeatureDisplay()
  {
    if(mouseMotionListener == null)
      mouseMotionListener = new MouseMotionListener()
      {
        public void mouseDragged(MouseEvent e){}

        public void mouseMoved(MouseEvent e)
        {
          int thisWidth = WIDTH;
          if(thisWidth < 5)
            thisWidth = 15;
          if(e.getX() > (display.getSize().width - thisWidth))
          {
            if(!containsComponent())
              display.add(ZoomScrollBar.this, "East");
          }
          else
          {
            if(containsComponent())
              display.remove(ZoomScrollBar.this);
          }
          display.repaint();
          display.revalidate();
        }
      };
    
    display.addMouseMotionListener(mouseMotionListener);
  }
  
  /**
   * Remove MouseMotionListener
   */
  protected void removeMouseMotionListenerFromFeatureDisplay()
  {
    if(mouseMotionListener != null)
      display.removeMouseMotionListener(mouseMotionListener); 
  }
  
  /**
   * Check to see if this component is contained by the display
   * (FeatureDisplay) component.
   * @return
   */
  private boolean containsComponent()
  {
    Component[] c = display.getComponents();
    for(int i=0; i<c.length; i++)
    {
      if(c[i].equals(this))
        return true;
    }
    
    return false;
  }
}