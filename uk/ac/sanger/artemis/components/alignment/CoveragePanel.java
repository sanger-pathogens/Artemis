/* CoveragePanel
 *
 * created: 2009
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2009  Genome Research Limited
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

package uk.ac.sanger.artemis.components.alignment;

import java.awt.Color;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.geom.GeneralPath;
import java.util.List;

import javax.swing.JPanel;

import net.sf.samtools.SAMRecord;

  public class CoveragePanel extends JPanel
  {
    private static final long serialVersionUID = 1L;

    private int start;
    private int end;
    private float pixPerBase;
    private Color lightBlue = new Color(30,144,255);
    private BamView jamView;
    
    public CoveragePanel(BamView jamView)
    {
      super();
      setBackground(Color.white);
      this.jamView = jamView;
    }
    
    /**
     * Override
     */
    protected void paintComponent(Graphics g)
    {
      super.paintComponent(g);
      Graphics2D g2 = (Graphics2D)g;
      
      List<SAMRecord> readsInView = jamView.getReadsInView();
      if(readsInView == null)
        return;
      
      int windowSize = 10;
      int nBins = Math.round((end-start+1.f)/windowSize);
      int coverage[] = new int[nBins];
      for(int i=0; i<coverage.length; i++)
        coverage[i] = 0;
      
      int max = 0;
      for(int i=0; i<readsInView.size(); i++)
      {
        SAMRecord thisRead = readsInView.get(i);
        int offset = jamView.getSequenceOffset(thisRead.getReferenceName());
        int length = thisRead.getReadLength();
        
        for(int j=0; j<length;j++)
        {
          int bin = 
            (int)(((thisRead.getAlignmentStart()-start) + j + offset) / windowSize);

          if(bin < 0 || bin > coverage.length-1)
            continue;   
          coverage[bin]+=1;
          if(coverage[bin] > max)
            max = coverage[bin];
        }
      }
      
      g2.setColor(lightBlue);
      GeneralPath shape = new GeneralPath();
      shape.moveTo(0,getHeight());
      for(int i=0; i<coverage.length; i++)
      {
        float xpos = ((i*(windowSize)) - windowSize/2.f)*pixPerBase;
        shape.lineTo(xpos,
            getHeight() - (((float)coverage[i]/(float)max)*getHeight()));
      }

      shape.lineTo(getWidth(),getHeight());
      g2.fill(shape);

      String maxStr = Float.toString(max/windowSize);
      FontMetrics fm = getFontMetrics(getFont());
      g2.setColor(Color.black);
      
      int xpos = getWidth() - fm.stringWidth(maxStr) - 
                 jamView.getJspView().getVerticalScrollBar().getWidth();
      g2.drawString(maxStr, xpos, fm.getHeight());
    }
    
    protected void setStartAndEnd(int start, int end)
    {
      this.start = start;
      this.end = end;
    }

    protected void setPixPerBase(float pixPerBase)
    {
      this.pixPerBase = pixPerBase;
    }
  }