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
    private Hashtable<String, int[]> plots;
    private int combinedCoverage[];

    private int nBins;
    private static boolean redraw = false;
    private boolean setMaxBases = false;
    
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
      JMenuItem configure = new JMenuItem("Configure Line(s)...");
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
        }
      });
      menu.add(configure);
      
      JMenuItem optMenu = new JMenuItem("Options...");
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

      plots = new Hashtable<String, int[]>();
      combinedCoverage = null;
      if(includeCombined)
      {
        combinedCoverage = new int[nBins];
        for(int k=0; k<combinedCoverage.length; k++)
          combinedCoverage[k] = 0;
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
        int[] thisPlot = plots.get(fileName);
        for(int i=1; i<thisPlot.length; i++)
          if(max < thisPlot[i])
            max = thisPlot[i];
      }
      draw(g2, getWidth(), getHeight());
    }
    
    protected void addRecord(SAMRecord thisRead, int offset, String fileName)
    {
      int coverage[] = plots.get(fileName);
      if(coverage == null)
      {
        coverage = new int[nBins];
        for(int k=0; k<coverage.length; k++)
          coverage[k] = 0;
        plots.put(fileName, coverage);
      }

      List<AlignmentBlock> blocks = thisRead.getAlignmentBlocks();
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

          coverage[bin]+=1;
          if(coverage[bin] > max)
            max = coverage[bin];

          if(includeCombined)
          {
            combinedCoverage[bin]+=1;
            if(combinedCoverage[bin] > max)
              max = combinedCoverage[bin];
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

      Enumeration<String> plotEum = plots.keys();
      while(plotEum.hasMoreElements())
      {
        String fileName = (String) plotEum.nextElement();
        int[] thisPlot = plots.get(fileName);

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
            int y0 = (int) (hgt - (((float)thisPlot[i-1]/(float)max)*hgt));
            int x1 = (int) (((i*(windowSize)) - windowSize/2.f)*pixPerBase);
            int y1 = (int) (hgt - (((float)thisPlot[i]/(float)max)*hgt));

            g2.drawLine(x0, y0, x1, y1);
          }
        }
        else // filled plots
        {
          g2.setComposite(makeComposite(0.75f));

          GeneralPath shape = new GeneralPath();
          shape.moveTo(0,hgt);
          for(int i=0; i<thisPlot.length; i++)
          {
            float xpos = ((i*(windowSize)) - windowSize/2.f)*pixPerBase;
            shape.lineTo(xpos,
                hgt - (((float)thisPlot[i]/(float)max)*hgt));
          }

          shape.lineTo(wid,hgt);
          g2.fill(shape);
        }
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

      final JCheckBox showCombined = new JCheckBox("Show Combined Plot", includeCombined);
      if(bamView.bamList.size() == 1)
        showCombined.setEnabled(false);
      c.gridy = c.gridy+1;
      opts.add(showCombined, c);
          
      String window_options[] = { "OK", "Cancel" };
      int select = JOptionPane.showOptionDialog(null, opts, "Coverage Options",
          JOptionPane.DEFAULT_OPTION, JOptionPane.QUESTION_MESSAGE, null,
          window_options, window_options[0]);

      if(select == 1)
        return;
      
      redraw = true;
      autoWinSize = autoSet.isSelected();
      includeCombined = showCombined.isSelected();
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