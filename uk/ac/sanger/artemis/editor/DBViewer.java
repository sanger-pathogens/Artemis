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

import java.awt.*;
import javax.swing.*;
import java.util.Vector;
import java.util.Enumeration;

public class DBViewer extends JPanel
{

  Vector hitInfoCollection = null;

  public DBViewer(FastaTextPane fastaPane)
  {
    super();
    hitInfoCollection = fastaPane.getHitCollection();
  }

  /**
  *
  * Override paintComponent
  * @param g    graphics
  *
  */
  public void paintComponent(Graphics g)
  {
// let UI delegate paint first (incl. background filling)
    super.paintComponent(g);
    
    Graphics2D g2 = (Graphics2D)g;
    int bound = 10;
    int width = getWidth()-(2*bound);
    g2.setColor(Color.black);
    g2.setStroke(new BasicStroke(3.f));
    g2.drawLine(bound,bound,bound+width,bound);

    int ydisp = bound+5;
    Enumeration enumHits = hitInfoCollection.elements();
    while(enumHits.hasMoreElements())
    {
      ydisp += 5;
      HitInfo hit = (HitInfo)enumHits.nextElement();
 
      g2.setColor(Color.red);
      g2.setStroke(new BasicStroke(1.f));
 
      float hit_unit = (float)width/(float)hit.getQueryLength();     
      int start = (int)(bound+(hit_unit*hit.getQueryStart()));
      int end   = (int)(bound+(hit_unit*hit.getQueryEnd()));
      g2.drawLine(start,bound+ydisp,end,bound+ydisp);

      System.out.println(hit.getID()+" "+hit.getQueryStart()+" -> "+hit.getQueryEnd()+
                         " start "+start+" end "+end+" width "+width);
//                       ", E = "+hit.getEValue()+"  Length = "+hit.getQueryLength());
    }

  }

}

