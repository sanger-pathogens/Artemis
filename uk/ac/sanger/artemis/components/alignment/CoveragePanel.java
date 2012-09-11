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
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.geom.GeneralPath;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.List;

import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;

import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMRecord;

  public class CoveragePanel extends AbstractGraphPanel
  {
    private static final long serialVersionUID = 1L;
    private static LineAttributes lines[];
    private boolean includeCombined = false;
    private Hashtable<String, int[][]> plots;
    private int combinedCoverage[][];

    private static boolean redraw = false;
    private boolean setMaxBases = false;
    
    private boolean plotByStrand = false;

    protected CoveragePanel(final BamView bamView)
    {
      this();
      this.bamView = bamView;
      createMenus(popup);
      addMouseListener(new PopupListener());
    }
    
    protected CoveragePanel()
    {
      super();
      setMaxBases = true;
    }
    
    protected void createMenus(JComponent menu)
    {
      final JMenuItem configure = new JMenuItem("Configure Line(s)...");
      configure.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          int size = bamView.bamList.size();
          if(includeCombined)
            size++;
          lines =
            LineAttributes.configurePlots(bamView.bamList, 
                getLineAttributes(size), CoveragePanel.this);
          bamView.refreshColourOfBamMenu();
        }
      });
      menu.add(configure);
      
      final JMenuItem optMenu = new JMenuItem("Options...");
      optMenu.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          defineOpts();
          bamView.repaint();
        }
      });
      menu.add(optMenu);
    }
    
    /**
     * Override
     */
    protected void paintComponent(Graphics g)
    {
      super.paintComponent(g);
      Graphics2D g2 = (Graphics2D)g;
      if(plots == null)
        return;

      drawSelectionRange(g2, pixPerBase, start, end, getHeight(), Color.PINK);
      drawPlot(g2);
      drawMax(g2);
    }
    
    protected void init(BamView bamView, float pixPerBase, int start, int end)
    {
      this.bamView = bamView;
      setPixPerBase(pixPerBase);
      setStartAndEnd(start, end);
      init();
    }
    
    private void init()
    {
      if(autoWinSize)
      {
        windowSize = (bamView.getBasesInView()/300);
        userWinSize = windowSize;
      }
      else
        windowSize = userWinSize;
      
      if(windowSize < 1)
        windowSize = 1;
      nBins = Math.round((end-start+1.f)/windowSize);

      plots = new Hashtable<String, int[][]>();
      combinedCoverage = null;
      if(includeCombined)
      {
        combinedCoverage = new int[nBins][2];
        for(int k=0; k<combinedCoverage.length; k++)
          for(int l=0; l<2; l++)
            combinedCoverage[k][l] = 0;
        plots.put("-1", combinedCoverage);
      }
      max = 0;
    }

    private void drawPlot(Graphics2D g2)
    {
      max = 0;
      Enumeration<String> plotEum = plots.keys();
      while(plotEum.hasMoreElements())
      {
        String fileName = (String) plotEum.nextElement();
        int[][] thisPlot = plots.get(fileName);
        for(int i=1; i<thisPlot.length; i++)
        {
          if(plotByStrand)
          {
            for(int j=0; j<2; j++)
              if(max < thisPlot[i][j])
                max = thisPlot[i][j];
          }
          else if(max < thisPlot[i][0])
            max = thisPlot[i][0];
        }
      }
      draw(g2, getWidth(), getHeight());
    }
    
    protected void addRecord(SAMRecord thisRead, int offset, String fileName)
    {
      int coverage[][] = plots.get(fileName);
      if(coverage == null)
      {
        coverage = new int[nBins][2];
        for(int k=0; k<nBins; k++)
          for(int l=0; l<2; l++)
            coverage[k][l] = 0; 
        plots.put(fileName, coverage);
      }

      final int col;
      if(plotByStrand && thisRead.getReadNegativeStrandFlag())
        col = 1;
      else
        col = 0;
      final List<AlignmentBlock> blocks = thisRead.getAlignmentBlocks();
      for(int j=0; j<blocks.size(); j++)
      {
        AlignmentBlock block = blocks.get(j);
        int refStart = block.getReferenceStart();
        for(int k=0; k<block.getLength(); k++)
        {
          int pos = refStart + k + offset;
          int bin = pos/windowSize;
          if(bin < 0 || bin > nBins-1)
            continue;

          coverage[bin][col]+=1;
          if(coverage[bin][col] > max)
            max = coverage[bin][col];

          if(includeCombined)
          {
            combinedCoverage[bin][col]+=1;
            if(combinedCoverage[bin][col] > max)
              max = combinedCoverage[bin][col];
          }
        } 
      }
    }
    
    protected void draw(Graphics2D g2, int wid, int hgt)
    {
      int size = bamView.bamList.size();
      if(includeCombined)
      {
        lines = getLineAttributes(size+1);
        lines[size].setLineColour(Color.black);
      }
      else
        lines = getLineAttributes(size);

      int hgt2 = hgt/2;

      Enumeration<String> plotEum = plots.keys();
      while(plotEum.hasMoreElements())
      {
        String fileName = (String) plotEum.nextElement();
        int[][] thisPlot = plots.get(fileName);

        int index;
        if(fileName.equals("-1"))
          index = lines.length-1;
        else
          index = bamView.bamList.indexOf(fileName);
        
        g2.setColor(lines[index].getLineColour());
        if(lines[index].getPlotType() == LineAttributes.PLOT_TYPES[0])
        {
          g2.setStroke(lines[index].getStroke());
          for(int i=1; i<thisPlot.length; i++)
          {
            int x0 = (int) ((((i-1)*(windowSize)) - windowSize/2.f)*pixPerBase);
            int x1 = (int) (((i*(windowSize)) - windowSize/2.f)*pixPerBase);
            int y0, y1;
            if(plotByStrand)
            {
              for(int col=0; col<2; col++)
              {
                final int factor;
                if(col == 0)
                  factor = 1;   // fwd strand
                else
                  factor = -1;  // reverse strand
                
                y0 = (int) (hgt2 - (factor)*(((float)thisPlot[i-1][col]/(float)max)*hgt2));
                y1 = (int) (hgt2 - (factor)*(((float)thisPlot[i][col]/(float)max)*hgt2));
                
                g2.drawLine(x0, y0, x1, y1);
              }
            }
            else
            {
              y0 = (int) (hgt - (((float)(thisPlot[i-1][0])/(float)max)*hgt));
              y1 = (int) (hgt - (((float)(thisPlot[i][0])/(float)max)*hgt));
              g2.drawLine(x0, y0, x1, y1);
            }
          }
        }
        else // filled plots
        {
          g2.setComposite(makeComposite(0.75f));

          if(plotByStrand)
          {
            final GeneralPath shapeFwd = new GeneralPath();
            shapeFwd.moveTo(0,hgt2);
            final  GeneralPath shapeBwd = new GeneralPath();
            shapeBwd.moveTo(0,hgt2);
          
            for(int i=0; i<thisPlot.length; i++)
            {
              float xpos = ((i*(windowSize)) - windowSize/2.f)*pixPerBase;
              for(int col=0; col<2; col++)
              {
                if(col == 0)
                  shapeFwd.lineTo(xpos,
                    hgt2 - (((float)thisPlot[i][col]/(float)max)*hgt2));
                else
                  shapeBwd.lineTo(xpos,
                    hgt2 + (((float)thisPlot[i][col]/(float)max)*hgt2));
              }
            }

            shapeBwd.lineTo(wid,hgt2);
            shapeFwd.lineTo(wid,hgt2);
            g2.fill(shapeBwd);
            g2.fill(shapeFwd);
          }
          else
          {
            final GeneralPath shape = new GeneralPath();
            shape.moveTo(0,hgt);
            for(int i=0; i<thisPlot.length; i++)
            {
              float xpos = ((i*(windowSize)) - windowSize/2.f)*pixPerBase;
              shape.lineTo(xpos,
                  hgt - (((float)thisPlot[i][0]/(float)max)*hgt));
            }
            shape.lineTo(wid,hgt);
            g2.fill(shape);
          }
        }
      }

      if(plotByStrand)
      {
        g2.setColor(Color.GRAY);
        g2.drawLine(0, hgt2, wid+1, hgt2);
      }
    }
    
    private AlphaComposite makeComposite(float alpha)
    {
      int type = AlphaComposite.SRC_OVER;
      return(AlphaComposite.getInstance(type, alpha));
    }

    protected static LineAttributes[] getLineAttributes(int nsize)
    {
      if(lines == null)
        lines = LineAttributes.init(nsize);
      else if(lines.length < nsize)
      {
        LineAttributes tmpLines[] = LineAttributes.init(nsize);
        for(int i=0;i<lines.length;i++)
          tmpLines[i] = lines[i];
        lines = tmpLines;
      }
      return lines;
    }

    /**
     * @return the redraw
     */
    protected static boolean isRedraw()
    {
      if(redraw)
      {
        redraw = false;
        return true;
      }
      return redraw;
    }
    
    /**
     * @param plotByStrand the plotByStrand to set
     */
    protected void setPlotByStrand(boolean plotByStrand)
    {
      redraw = true;
      this.plotByStrand = plotByStrand;
    }
    
    private void defineOpts()
    {
      final JPanel opts = new JPanel(new GridBagLayout());
      final GridBagConstraints c = new GridBagConstraints();

      final JTextField newBaseMax = new JTextField(Integer.toString(bamView.getMaxBases()), 10);
      c.gridy = 0;
      if(setMaxBases)
      {
        final JLabel labMax1 = new JLabel("Zoom level before switching");
        final JLabel labMax2 = new JLabel("to coverage view (in bases):");
        c.anchor = GridBagConstraints.WEST;
        opts.add(labMax1, c);
        c.gridy = c.gridy+1;
        opts.add(labMax2, c);
        opts.add(newBaseMax, c);
      }
      
      final JTextField newWinSize = new JTextField(Integer.toString(userWinSize), 10);
      final JLabel lab = new JLabel("Window size:");
      lab.setEnabled(!autoWinSize);
      newWinSize.setEnabled(!autoWinSize);
      
      c.gridy = c.gridy+1;
      c.anchor = GridBagConstraints.EAST;
      opts.add(lab, c);
      opts.add(newWinSize, c);

      final JCheckBox autoSet = new JCheckBox("Automatically set window size", autoWinSize);
      autoSet.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          lab.setEnabled(!autoSet.isSelected());
          newWinSize.setEnabled(!autoSet.isSelected());
        }
      });
      c.anchor = GridBagConstraints.WEST;
      c.gridy = c.gridy+1;
      c.gridwidth = GridBagConstraints.REMAINDER;
      opts.add(autoSet, c);

      final JCheckBox showCombined = new JCheckBox("Show combined plot", includeCombined);
      if(bamView.bamList.size() == 1)
        showCombined.setEnabled(false);
      c.gridy = c.gridy+1;
      opts.add(showCombined, c);
      
      final JCheckBox byStrand = new JCheckBox("Plot by strand", plotByStrand);
      c.gridy = c.gridy+1;
      opts.add(byStrand, c);

      String window_options[] = { "OK", "Cancel" };
      int select = JOptionPane.showOptionDialog(null, opts, "Coverage Options",
          JOptionPane.DEFAULT_OPTION, JOptionPane.QUESTION_MESSAGE, null,
          window_options, window_options[0]);

      if(select == 1)
        return;
      
      redraw = true;
      autoWinSize = autoSet.isSelected();
      includeCombined = showCombined.isSelected();
      plotByStrand = byStrand.isSelected();
      
      try
      {
        userWinSize = Integer.parseInt(newWinSize.getText().trim());
        if(setMaxBases)
          bamView.setMaxBases(Integer.parseInt(newBaseMax.getText().trim()));
      }
      catch (NumberFormatException nfe)
      {
        return;
      }
    }

  }