/* BaseAlgorithm.java
 *
 * created: Wed Dec 16 1998
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 1998,1999,2000  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/plot/BaseAlgorithm.java,v 1.9 2009-06-26 15:52:48 tjc Exp $
 */

package uk.ac.sanger.artemis.plot;

import uk.ac.sanger.artemis.sequence.*;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.FontMetrics;
import java.awt.Color;
import java.awt.BasicStroke;

/**
 *  The BaseAlgorithm class is the base class for algorithms that work
 *  directly on bases.  A BaseAlgorithm has a name and is specific to one
 *  Strand of DNA, meaning the algorithm can't change strand part way along.
 *
 *  @author Kim Rutherford
 *  @version $Id: BaseAlgorithm.java,v 1.9 2009-06-26 15:52:48 tjc Exp $
 **/

public abstract class BaseAlgorithm extends Algorithm 
{
  private Strand strand;

  /**
   *  Create a new BaseAlgorithm object.
   *  @param strand The strand to do the calculation on.
   *  @param algorithm_name A String used to identify this algorithm to the
   *    user.
   *  @param algorithm_short_name A String used to identify this algorithm
   *    internally.  See the Algorithm constructor for more details.
   **/
  public BaseAlgorithm (final Strand strand, final String algorithm_name,
                        final String algorithm_short_name) 
  {
    super (algorithm_name, algorithm_short_name);
    this.bases = strand.getBases();
    this.strand = strand;

    if (strand.isForwardStrand ()) {
      forward_flag = true;
    } else {
      forward_flag = false;
    }
  }

  /**
   *  Return the Bases object of the Strand that was passed to the
   *  constructor.
   **/
  public Bases getBases () {
    return bases;
  }

  /**
   *  Returns the strand we will do the calculation on.
   **/
  public Strand getStrand () {
    if (forward_flag ^ rev_comp_display) {
      return getBases ().getForwardStrand ();
    } else {
      return getBases ().getReverseStrand ();
    }
  }
  
  /**
   *  If rev_comp_display is true all calculations will be performs on the
   *  opposite Strand to the strand that was passed to the constructor.
   **/
  public void setRevCompDisplay (final boolean rev_comp_display) {
    this.rev_comp_display = rev_comp_display;
  }

  /**
   *  Returns true if the FeatureDisplay is reverse complemented.  All
   *  calculations should be performed on the opposite Strand to the strand
   *  that was passed to the constructor.
   **/
  public boolean isRevCompDisplay () {
    return rev_comp_display;
  }


  /**
  *  Draw in a legend
  */
  public void drawLegend(Graphics g, int font_height,
                         int font_width, LineAttributes[] lines,
                         int numPlots)
  {
    Graphics2D g2d = (Graphics2D)g;

    FontMetrics fm = g2d.getFontMetrics();
    int lineHgt    = 3 * font_height/4; 

    if( (strand.isForwardStrand() && !isRevCompDisplay()) ||
        (!strand.isForwardStrand() && isRevCompDisplay()))
    {
      int width = 0;
      for(int i=0; i<numPlots; i++)
      {
        g2d.setColor(Color.black);
        
        if(lines[i].getLabel() != null) // user defined label
        {
          g2d.drawString(lines[i].getLabel(),width,font_height);
          width += lines[i].getLabelWidth(fm);
        }
        else
        {
          g2d.drawString(Integer.toString(i+1),width,font_height);
          width += 5*font_width;
        }
        
        BasicStroke stroke = (BasicStroke)g2d.getStroke();
        g2d.setStroke(new BasicStroke(3.f));
        g2d.setColor(lines[i].getLineColour());
        g2d.drawLine(width - (font_width*1), lineHgt, width - (font_width*3), lineHgt);
        g2d.setStroke(stroke);
      }
    }
    else
    {
      g2d.setColor(Color.black);
      g2d.drawString("4",0,font_height);
      g2d.drawString("5",font_width*5,font_height);
      g2d.drawString("6",font_width*10,font_height);

      BasicStroke stroke = (BasicStroke)g2d.getStroke();
      g2d.setStroke(new BasicStroke(3.f));
      int frame = strand.getSequenceLength() % 3;
      
      //System.out.println("FRAME "+frame+"  length="+strand.getSequenceLength());
      Color col4 = null;
      Color col5 = null;
      Color col6 = null;
       
      switch(frame)
      {
         case 0:
           col4 = lines[1].getLineColour();
           col5 = lines[2].getLineColour();
           col6 = lines[0].getLineColour();
           break;
         case 1:
           col4 = lines[2].getLineColour();
           col5 = lines[0].getLineColour();
           col6 = lines[1].getLineColour();
           break;
         case 2:
           col4 = lines[0].getLineColour();
           col5 = lines[1].getLineColour();
           col6 = lines[2].getLineColour();
           break;
      }
    
      g2d.setColor(col4);
      g2d.drawLine(font_width*2, lineHgt, font_width*4, lineHgt);
  
      g2d.setColor(col5);
      g2d.drawLine(font_width*7, lineHgt, font_width*9, lineHgt);

      g2d.setColor(col6);
      g2d.drawLine(font_width*12, lineHgt, font_width*14, lineHgt);
      g2d.setStroke(stroke);
    }
  }

  /**
   *  Return the value of the function between a pair of bases.
   *  @param start The start base (included in the range).
   *  @param end The end base (included in the range).
   *  @param values The results are returned in this array, hence it should be
   *    allocated at the size given by getValueCount ().
   **/
  public abstract void getValues (int start, int end, final float [] values);

  /**
   *  Return the number of values a call to getValues () will return.
   **/
  public abstract int getValueCount ();
  
  /**
   *  The Bases we will do the calculation on.
   **/
  private Bases bases;

  /**
   *  If rev_comp_display is true all calculations will be performed on the
   *  opposite Strand to the strand that was passed to the constructor.
   **/
  private boolean rev_comp_display = false;

  /**
   *  true if and only if the calculations should be done on the forward
   *  Strand.
   **/
  private boolean forward_flag;
}
