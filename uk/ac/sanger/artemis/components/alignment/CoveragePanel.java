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

import java.awt.AlphaComposite;
import java.awt.Color;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.geom.GeneralPath;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.List;

import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;

import net.sf.samtools.SAMRecord;

  public class CoveragePanel extends JPanel
  {
    private static final long serialVersionUID = 1L;

    private int start;
    private int end;
    private float pixPerBase;
    private BamView jamView;
    private JPopupMenu popup;
    private LineAttributes lines[];
    
    public CoveragePanel(final BamView jamView)
    {
      super();
      setBackground(Color.white);
      this.jamView = jamView;
      
      popup = new JPopupMenu();
      JMenuItem configure = new JMenuItem("Configure...");
      configure.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          lines =
            LineAttributes.configurePlots(jamView.bamList, lines, CoveragePanel.this);
        }
      });
      popup.add(configure);
      
      addMouseListener(new PopupListener());
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
      
      int windowSize = (jamView.getBasesInView()/200);
      if(windowSize < 5)
        windowSize = 5;

      int nBins = Math.round((end-start+1.f)/windowSize);
      int max = drawPlot(g2, nBins, windowSize);
      
      String maxStr = Float.toString(max/windowSize);
      FontMetrics fm = getFontMetrics(getFont());
      g2.setColor(Color.black);
      
      int xpos = getWidth() - fm.stringWidth(maxStr) - 
                 jamView.getJspView().getVerticalScrollBar().getWidth();
      g2.drawString(maxStr, xpos, fm.getHeight());
    }
    
    
    private int drawPlot(Graphics2D g2, int nBins, int windowSize)
    {
      List<SAMRecord> readsInView = jamView.getReadsInView();
      List<String> bamList = jamView.bamList;
      final Hashtable<String, Integer[]> plots = new Hashtable<String, Integer[]>();
      
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

          if(bin < 0 || bin > nBins-1)
            continue;

          String fileName;
          if(bamList.size() > 1)
            fileName = bamList.get((Integer) thisRead.getAttribute("FL"));
          else
            fileName = bamList.get(0);
          Integer coverage[] = plots.get(fileName);
          
          if(coverage == null)
          {
            coverage = new Integer[nBins];
            for(int k=0; k<coverage.length; k++)
              coverage[k] = 0;
            plots.put(fileName, coverage);
          }         
 
          coverage[bin]+=1;
          if(coverage[bin] > max)
            max = coverage[bin];
        }
      }

      lines = getLineAttributes();
      Enumeration<String> plotEum = plots.keys();
      while(plotEum.hasMoreElements())
      {
        String fileName = (String) plotEum.nextElement();
        Integer[] thisPlot = plots.get(fileName);
        
        g2.setColor(lines[bamList.indexOf(fileName)].getLineColour());
        
        if(lines[bamList.indexOf(fileName)].getPlotType() == LineAttributes.PLOT_TYPES[0])
        {
          g2.setStroke(lines[bamList.indexOf(fileName)].getStroke());
          for(int i=1; i<thisPlot.length; i++)
          {
            int x0 = (int) ((((i-1)*(windowSize)) - windowSize/2.f)*pixPerBase);
            int y0 = (int) (getHeight() - (((float)thisPlot[i-1]/(float)max)*getHeight()));
            int x1 = (int) (((i*(windowSize)) - windowSize/2.f)*pixPerBase);
            int y1 = (int) (getHeight() - (((float)thisPlot[i]/(float)max)*getHeight()));
            
            g2.drawLine(x0, y0, x1, y1);
          }
        }
        else // filled plots
        {
          g2.setComposite(makeComposite(0.75f));

          GeneralPath shape = new GeneralPath();
          shape.moveTo(0,getHeight());
          for(int i=0; i<thisPlot.length; i++)
          {
            float xpos = ((i*(windowSize)) - windowSize/2.f)*pixPerBase;
            shape.lineTo(xpos,
                getHeight() - (((float)thisPlot[i]/(float)max)*getHeight()));
          }

          shape.lineTo(getWidth(),getHeight());
          g2.fill(shape);
        }
      }
      return max;
    }
    
    private AlphaComposite makeComposite(float alpha)
    {
      int type = AlphaComposite.SRC_OVER;
      return(AlphaComposite.getInstance(type, alpha));
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
    
    
    private LineAttributes[] getLineAttributes()
    {
      if(lines == null)
        lines = LineAttributes.init(jamView.bamList.size());
      return lines;
    }
  
  
  /**
   * Popup menu listener
   */
   class PopupListener extends MouseAdapter
   {
     JMenuItem gotoMateMenuItem;
     JMenuItem showDetails;
     
     public void mouseClicked(MouseEvent e)
     {
     }
     
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
       {
         popup.show(e.getComponent(),
                 e.getX(), e.getY());
       }
     }
   }
  }