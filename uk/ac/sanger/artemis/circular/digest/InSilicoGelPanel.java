/*
 * Copyright (C) 2009  Genome Research Limited
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
 */
package uk.ac.sanger.artemis.circular.digest;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.List;
import java.util.Vector;

import javax.swing.ButtonGroup;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JRadioButtonMenuItem;

import uk.ac.sanger.artemis.components.FileViewer;

public class InSilicoGelPanel extends JPanel implements ActionListener
{
  private static final long serialVersionUID = 1L;
  private int marginHeight = 50;
  private int marginWidth = 55;
  private int panelHeight;
  private List<FragmentBand> genomeFragments = new Vector<FragmentBand>();
  private static int MAX_FRAGMENT_LENGTH = 0;
  private static int MIN_FRAGMENT_LENGTH = Integer.MAX_VALUE;
  private boolean drawLog = false;
  private JPopupMenu popup;
  private File restrictOutput;

  /**
   * @param genomeLength
   * @param cutSites
   * @param panelHeight
   * @param restrictOutput
   */
  public InSilicoGelPanel(final int genomeLength, final List<CutSite> cutSites,
      final int panelHeight, final File restrictOutput,
      final String name)
  {
    this.panelHeight = panelHeight;
    this.restrictOutput = restrictOutput;
    setToolTipText(name);

    Integer len;
    int firstSiteEnd = 0;
    int lastSite = 0;
    for (int i = 0; i < cutSites.size(); i++)
    {
      CutSite cutSite = cutSites.get(i);
      if (i == 0)
        firstSiteEnd = cutSite.getFivePrime() + 1;
      else
      {
        len = cutSite.getFivePrime() - lastSite;

        FragmentBand band = new FragmentBand();
        band.genomeFragmentLength = len;
        band.bandCutSite = cutSite;
        genomeFragments.add(band);
        if (len > MAX_FRAGMENT_LENGTH)
          MAX_FRAGMENT_LENGTH = len;
        if (len < MIN_FRAGMENT_LENGTH)
          MIN_FRAGMENT_LENGTH = len;
      }
      lastSite = cutSite.getFivePrime();
    }

    len = genomeLength - lastSite + firstSiteEnd;
    FragmentBand band = new FragmentBand();
    band.genomeFragmentLength = len;
    
    if(cutSites.size() > 0)
    {
      band.bandCutSite = cutSites.get(0);
      genomeFragments.add(band);
    }
    if (len > MAX_FRAGMENT_LENGTH)
      MAX_FRAGMENT_LENGTH = len;
    // System.out.println(len.toString());

    init();
  }
  
  public InSilicoGelPanel(final List<FragmentBand> genomeFragments,
      final int panelHeight, final File restrictOutput,
      final String name)
  {
    this.panelHeight = panelHeight;
    this.restrictOutput = restrictOutput;
    this.genomeFragments = genomeFragments;
   
    setToolTipText(name);
    for(int i=0; i<genomeFragments.size(); i++)
    {
      int len = genomeFragments.get(i).genomeFragmentLength;
      if (len > MAX_FRAGMENT_LENGTH)
        MAX_FRAGMENT_LENGTH = len;
      if (len < MIN_FRAGMENT_LENGTH)
        MIN_FRAGMENT_LENGTH = len;
    }
    init();
  }
  
  /**
   * Initialise the panel
   */
  private void init()
  {
    setBackground(Color.white);
    setPreferredSize(new Dimension(150, panelHeight));
    
    MouseListener popupListener = new PopupListener();
    addMouseListener(popupListener);

    popup = new JPopupMenu();
    JRadioButtonMenuItem linearScale = new JRadioButtonMenuItem("Linear scale");
    popup.add(linearScale);
    linearScale.addActionListener(this);
    JRadioButtonMenuItem logScale = new JRadioButtonMenuItem("Log scale");
    popup.add(logScale);
    logScale.addActionListener(this);
    ButtonGroup group = new ButtonGroup();
    group.add(linearScale);
    group.add(logScale);
    linearScale.setSelected(true);

    JMenuItem showCutSites = new JMenuItem("Show cut site details");
    popup.add(showCutSites);
    showCutSites.addActionListener(this);
  }

  /**
   * Override
   */
  public void paintComponent(Graphics g)
  {
    super.paintComponent(g);

    Graphics2D g2D = (Graphics2D) g;
    g2D.draw3DRect(marginWidth, (marginHeight / 2), marginWidth, panelHeight
        - (marginHeight), true);

    int gelHeight = (panelHeight - (2 * marginHeight));
    g2D.setColor(Color.blue);

    BasicStroke stroke = new BasicStroke(1.f);
    BasicStroke stroke2 = new BasicStroke(2.f);
    g2D.setStroke(stroke);

    for (int i = 0; i < genomeFragments.size(); i++)
    {
      final int y = getYPosition(genomeFragments.get(i), gelHeight);

      if ((genomeFragments.get(i).bandCutSite == null) || 
          !genomeFragments.get(i).bandCutSite.isHighlighted())
      {
        g2D.setStroke(stroke);
        g2D.setColor(Color.blue);
      }
      else
      {
        g2D.setColor(Color.black);
        g2D.drawString(genomeFragments.get(i).bandCutSite.getEnzymeName(),
            marginWidth + marginWidth + 2, y);
        g2D.setStroke(stroke2);
        g2D.setColor(Color.yellow);
      }
      g2D.drawLine(marginWidth, y, marginWidth + marginWidth, y);
    }

    drawScale(g2D, stroke, gelHeight);
  }
  
  private int getYPosition(final FragmentBand band, int gelHeight)
  {
    int fragmentLength = band.genomeFragmentLength;
    final int y;
    if (isDrawLog())
      y = getLogValue(fragmentLength, marginHeight, gelHeight);
    else
      y = gelHeight
          + marginHeight
          - (int) (((float) (gelHeight) / (float) (MAX_FRAGMENT_LENGTH - MIN_FRAGMENT_LENGTH)) * fragmentLength);
    return y;
  }

  /**
   * Draw the fragment length scale
   * 
   * @param g2D
   * @param stroke
   * @param gelHeight
   */
  private void drawScale(Graphics2D g2D, BasicStroke stroke, int gelHeight)
  {
    g2D.setColor(Color.black);
    g2D.setStroke(stroke);
    NumberFormat formatter = new DecimalFormat("#0.0");

    int nscale = 8;
    if (isDrawLog())
      nscale = 5;

    float range = (MAX_FRAGMENT_LENGTH - MIN_FRAGMENT_LENGTH) / (float) nscale;
    for (int i = 0; i < nscale + 1; i++)
    {
      float length = MIN_FRAGMENT_LENGTH + (range * i);

      int y;
      if (isDrawLog())
        y = getLogValue((int) length, marginHeight, gelHeight);
      else
        y = gelHeight
            + marginHeight
            - (int) (((float) (gelHeight) / (float) (MAX_FRAGMENT_LENGTH - MIN_FRAGMENT_LENGTH)) * length);

      g2D.drawLine(marginWidth, y, marginWidth - 10, y);

      g2D.drawString(formatter.format(length / 1000) + "kb", 0, y);
    }
  }

  private static final double LOG10SCALE = 1.d / Math.log(10);

  protected FragmentBand getBandAtLocation(Point loc)
  {
    int gelHeight = (panelHeight - (2 * marginHeight));
    for (int i = 0; i < genomeFragments.size(); i++)
    {
      int y = getYPosition(genomeFragments.get(i), gelHeight);

      if(loc.y == y)
        return genomeFragments.get(i);
    }
    return null;
  }
  
  /**
   * Get base 10 commons log
   * 
   * @param val
   * @return
   */
  private static double log10(double val)
  {
    return Math.log(val) * LOG10SCALE;
  }

  private int getLogValue(int val, int gelStart, int gelRange)
  {
    double log_low = log10(MIN_FRAGMENT_LENGTH);
    double log_high = log10(MAX_FRAGMENT_LENGTH);
    double log_val = log10(val);

    double log_unit = ((double) gelRange) / (log_high - log_low);
    return (int) ((double) (gelRange + gelStart) - ((log_val - log_low) * log_unit));
  }

  private boolean isDrawLog()
  {
    return drawLog;
  }

  private void setDrawLog(boolean drawLog)
  {
    this.drawLog = drawLog;
  }

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
      if (e.isPopupTrigger())
        popup.show(e.getComponent(), e.getX(), e.getY());
    }
  }

  public void actionPerformed(ActionEvent e)
  {
    if (e.getSource() instanceof JRadioButtonMenuItem)
    {
      JRadioButtonMenuItem radioButton = (JRadioButtonMenuItem) e.getSource();
      if (radioButton.isSelected())
      {
        setDrawLog(radioButton.getText().startsWith("Log"));
        repaint();
      }
    }
    else
    {
      final FileViewer viewer = new FileViewer(restrictOutput.getName(), true,
          false, false);
      BufferedReader br;
      try
      {
        br = new BufferedReader(new FileReader(restrictOutput));

        StringBuffer buff = new StringBuffer();
        String line;
        while ((line = br.readLine()) != null)
          buff.append(line + "\n");
        viewer.getTextPane().setText(buff.toString());
      }
      catch (Exception e2)
      {
        e2.printStackTrace();
      }
    }
  }
}