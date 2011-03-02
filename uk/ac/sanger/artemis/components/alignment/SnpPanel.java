/* SnpPanel
 *
 * created: 2010
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2010  Genome Research Limited
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

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.geom.GeneralPath;
import java.util.List;

import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JTextField;

import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.sequence.Bases;
import uk.ac.sanger.artemis.util.OutOfRangeException;

import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMRecord;

  public class SnpPanel extends AbstractGraphPanel
  {
    private static final long serialVersionUID = 1L;
    private BamView bamView;

    private Bases bases;
    private float minBaseQualityFilter = 0;
       
    public SnpPanel(final BamView bamView, Bases bases)
    {
      super();
      setBackground(Color.white);
      this.bamView = bamView;
      this.bases = bases;
      initPopupMenu(this, popup);
      
      JMenuItem configure = new JMenuItem("Filter by Base Quality...");
      configure.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          // filter by base quality
          JTextField filterField = new JTextField(Float.toString(minBaseQualityFilter));

          int status = JOptionPane.showConfirmDialog(SnpPanel.this, 
              filterField, "Base Quality Filter", 
              JOptionPane.OK_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE);
          if(status == JOptionPane.CANCEL_OPTION)
            return;
          try
          {
            minBaseQualityFilter = Float.parseFloat(filterField.getText());
          }
          catch(NumberFormatException nfe)
          {
            JOptionPane.showMessageDialog(SnpPanel.this, nfe.getMessage(), 
                "Number Format", JOptionPane.WARNING_MESSAGE);
          }
        }
      });
      popup.add(configure);
    }
    
    /**
     * Override
     */
    protected void paintComponent(Graphics g)
    {
      super.paintComponent(g);
      Graphics2D g2 = (Graphics2D)g;
      
      if(bases == null)
        return;
      
      List<SAMRecord> readsInView = bamView.getReadsInView();
      if(readsInView == null)
        return;
      
      int windowSize;
      
      if(autoWinSize)
      {
        windowSize = (bamView.getBasesInView()/300);
        userWinSize = windowSize;
      }
      else
        windowSize = userWinSize;
      
      if(windowSize < 1)
        windowSize = 1;

      int nBins = Math.round((end-start+1.f)/windowSize);
      int max = drawPlot(g2, nBins, windowSize);
      
      String maxStr = Float.toString(max/windowSize);
      FontMetrics fm = getFontMetrics(getFont());
      g2.setColor(Color.black);
      
      int xpos = getWidth() - fm.stringWidth(maxStr) - 
      bamView.getJspView().getVerticalScrollBar().getWidth();
      g2.drawString(maxStr, xpos, fm.getHeight());
    }
    
    
    private int drawPlot(Graphics2D g2, int nBins, int windowSize)
    {
      //lines = CoveragePanel.getLineAttributes(bamView.bamList.size());
      List<SAMRecord> readsInView = bamView.getReadsInView();
      
      int snpCount[] = new int[nBins];
      for(int i=0; i<snpCount.length; i++)
        snpCount[i] = 0;
      
      int max = 0;
      for(int i=0; i<readsInView.size(); i++)
      {
        SAMRecord thisRead = readsInView.get(i);
        max = calculateSNPs(thisRead, windowSize, nBins, snpCount, max);
      }

      g2.setColor(Color.red);
      g2.setStroke(new BasicStroke(1.f));
      
      if(windowSize == 1)
      {
        GeneralPath shape = new GeneralPath();
        shape.moveTo(0,getHeight());
        for(int i=0; i<snpCount.length; i++)
        {
          float xpos1 = ((i*windowSize) )*pixPerBase;
          float xpos2 = ((i*windowSize) + windowSize)*pixPerBase;
          
          shape.lineTo(xpos1,getHeight());
          shape.lineTo(xpos1,
              getHeight() - (((float)snpCount[i]/(float)max)*getHeight()));
          shape.lineTo(xpos2,
              getHeight() - (((float)snpCount[i]/(float)max)*getHeight()));
          
          shape.lineTo(xpos2,getHeight());
        }

        shape.lineTo(getWidth(),getHeight());
        g2.fill(shape);
      }
      else
      {
        for(int i=1; i<snpCount.length; i++)
        {
          int x0 = (int) (((i*windowSize) - windowSize/2.f)*pixPerBase);
          int y0 = (int) (getHeight() - (((float)snpCount[i-1]/(float)max)*getHeight()));
          int x1 = (int) ((((i+1)*windowSize) - windowSize/2.f)*pixPerBase);
          int y1 = (int) (getHeight() - (((float)snpCount[i]/(float)max)*getHeight()));
        
          g2.drawLine(x0, y0, x1, y1);
        }
      }

      return max;
    }
    
    /**
     * Display the SNPs for the given read.
     * @param g2
     * @param thisRead
     * @param pixPerBase
     * @param ypos
     */
    private int calculateSNPs(SAMRecord thisRead,
                               int windowSize, 
                               int nBins,
                               int[] snpCount,
                               int max)
    {
      int thisStart = thisRead.getAlignmentStart();
      int thisEnd   = thisRead.getAlignmentEnd();
      int offset    = bamView.getSequenceOffset(thisRead.getReferenceName());
      // use alignment blocks of the contiguous alignment of
      // subsets of read bases to a reference sequence
      List<AlignmentBlock> blocks = thisRead.getAlignmentBlocks();
      byte[] phredQuality = thisRead.getBaseQualities();
      try
      {
        char[] refSeq = bases.getSubSequenceC(
            new Range(thisStart+offset, thisEnd+offset), Bases.FORWARD);
        byte[] readSeq = thisRead.getReadBases();

        offset = offset - bamView.getBaseAtStartOfView();
        for(int i=0; i<blocks.size(); i++)
        {
          AlignmentBlock block = blocks.get(i);
          for(int j=0; j<block.getLength(); j++)
          {
            int readPos = block.getReadStart()-1+j;
            int refPos  = block.getReferenceStart()+j;

            if (Character.toUpperCase(refSeq[refPos-thisStart]) != readSeq[readPos])
            {
              if(phredQuality[readPos] < minBaseQualityFilter)
                continue;
              int bin = (int)((refPos+offset) / windowSize);

              if(bin < 0 || bin > nBins-1)
                continue;

              snpCount[bin]+=1;
              if(snpCount[bin] > max)
                max = snpCount[bin];
            }
          }
        }
      }
      catch (OutOfRangeException e)
      {
        e.printStackTrace();
      }
      return max;
    }
   
  }