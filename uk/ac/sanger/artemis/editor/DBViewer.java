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
import java.awt.event.*;
import javax.swing.*;
import java.util.Vector;
import java.util.Enumeration;

public class DBViewer extends ScrollPanel
{

  /** collection of hits */
  private Vector hitInfoCollection = null;
  /** query length */
  private int qlen;
  /** viewer boundary */
  private int bound = 10;
  /** y displacement for each hit */
  private int ydisp = 5;
  /** maximum hit score */
  private float max_score = 0;
  /** colour the hits by score or by evalue */
  private boolean colourByScore = false;
  /** popup menu */
  private JPopupMenu popup;
  /** number height */
  private int hgtNumber;

  public DBViewer(FastaTextPane fastaPane)
  {
    super();
    hitInfoCollection = fastaPane.getHitCollection();
    qlen = fastaPane.getQueryLength();
    
    Dimension d = new Dimension(500,
                               (hitInfoCollection.size()*ydisp)+(6*bound));
    setPreferredSize(d);
    setToolTipText("");   //enable tooltip display

    Enumeration enumHits = hitInfoCollection.elements();
    while(enumHits.hasMoreElements())
    {
      HitInfo hit = (HitInfo)enumHits.nextElement();
      float score = Float.parseFloat(hit.getScore());
      if(score > max_score)
        max_score = score;
    }

    // Popup menu
    addMouseListener(new PopupListener());
    popup = new JPopupMenu();

    JRadioButtonMenuItem colourScore = new JRadioButtonMenuItem("Colour by Score");
    colourScore.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        setColourByScore(true);
      }
    });
    popup.add(colourScore);

    JRadioButtonMenuItem colourEval  = new JRadioButtonMenuItem("Colour by E-value");
    colourEval.setSelected(true);
    colourEval.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        setColourByScore(false);
      }
    });
    popup.add(colourEval);

    ButtonGroup butt = new ButtonGroup();
    butt.add(colourScore);
    butt.add(colourEval);
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
 
    FontMetrics metrics = g.getFontMetrics();
    hgtNumber = metrics.getAscent();

//  Dimension d = new Dimension(500,
//                             (hitInfoCollection.size()*ydisp)+(4*bound)+hgtNumber);
//  setPreferredSize(d);

    Graphics2D g2 = (Graphics2D)g;
    int width = (int)getPreferredSize().getWidth()-(2*bound);
    g2.setColor(Color.black);
    g2.setStroke(new BasicStroke(3.f));
    g2.drawLine(bound,bound+hgtNumber,bound+width,bound+hgtNumber);

    g2.setStroke(new BasicStroke(1.f));
    g2.drawLine(bound,bound+hgtNumber,bound,bound+hgtNumber-6);
    g2.drawString("0",bound,hgtNumber+3);
    
// draw hits
    int ypos = bound+ydisp+hgtNumber;
    Enumeration enumHits = hitInfoCollection.elements();
    while(enumHits.hasMoreElements())
    {
      ypos += ydisp;
      HitInfo hit = (HitInfo)enumHits.nextElement();
 
      if(colourByScore)
      {
        float score = Float.parseFloat(hit.getScore());
        if(score > max_score/2)
          g2.setColor(Color.red);
        else if(score > max_score/4)
          g2.setColor(Color.blue);
        else
          g2.setColor(Color.cyan);
      }
      else
      {
        double evalue = Double.parseDouble(hit.getEValue());
        if(evalue < 0.005)
          g2.setColor(Color.red);
        else if(evalue < 0.)
          g2.setColor(Color.blue);
        else
          g2.setColor(Color.cyan);
      }

      g2.setStroke(new BasicStroke(1.f));
 
      float hit_unit = (float)width/(float)qlen;     
      int start = (int)(bound+(hit_unit*hit.getQueryStart()));
      int end   = (int)(bound+(hit_unit*hit.getQueryEnd()));
      g2.drawLine(start,bound+ypos,end,bound+ypos);

//    System.out.println(hit.getID()+" "+hit.getQueryStart()+" -> "+hit.getQueryEnd()+
//                       " start "+start+" end "+end+" width "+width);
//                       ", E = "+hit.getEValue()+"  Length = "+hit.getQueryLength());
    }
  
  }

  protected void setColourByScore(boolean colourByScore)
  {
    this.colourByScore = colourByScore;
    repaint();
  }
  

  /**
  *
  * Popup menu listener
  *
  */
  class PopupListener extends MouseAdapter
  {
    public void mousePressed(MouseEvent e)
    {
      maybeShowPopup(e);
    }

    public void mouseReleased(MouseEvent e)
    {
      maybeShowPopup(e);
    }

    private void maybeShowPopup(MouseEvent e)
    {
      if(e.isPopupTrigger())
        popup.show(e.getComponent(),
                e.getX(), e.getY());
    }
  }


  /**
  *
  * Determine the tool tip to display
  * @param e    mouse event
  * @return     tool tip
  *
  */
  public String getToolTipText(MouseEvent e)
  {
    Point loc = e.getPoint();
    int seqPos = (int)((loc.y-bound-ydisp-bound-ydisp-hgtNumber)/ydisp);

    if(seqPos >= 0 && seqPos<=hitInfoCollection.size())
    {
      HitInfo hit = (HitInfo)hitInfoCollection.get(seqPos);
      return hit.getID();
    }
    return null;
  }

}

