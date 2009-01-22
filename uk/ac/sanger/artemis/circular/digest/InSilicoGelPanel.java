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
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.util.List;
import java.util.Vector;

import javax.swing.ButtonGroup;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JRadioButtonMenuItem;

public class InSilicoGelPanel extends JPanel
                              implements ActionListener
{
	private static final long serialVersionUID = 1L;
	private int marginHeight = 50;
	private int marginWidth  = 60;
	private int panelHeight;
	private List<FragmentBand> genomeFragments = new Vector<FragmentBand>();
	private int MAX_FRAGMENT_LENGTH = 0;
	private int MIN_FRAGMENT_LENGTH = Integer.MAX_VALUE;
	private boolean drawLog = false;
	private JPopupMenu popup;


	public InSilicoGelPanel(final int genomeLength,
			                    final List<CutSite> cutSites,
			                    final int panelHeight)
	{
		this.panelHeight  = panelHeight;

		setBackground(Color.white);
		setPreferredSize(new Dimension(160,panelHeight));
		
		Integer len;
		int firstSiteEnd = 0;
		int lastSite = 0;
		for(int i=0; i<cutSites.size(); i++)
		{
			CutSite cutSite = cutSites.get(i);
			if(i == 0)
				firstSiteEnd = cutSite.getFivePrime()+1;
			else
			{
				len = cutSite.getFivePrime()-lastSite;
				System.out.println(i+" **** "+len.toString());
				FragmentBand band = new FragmentBand();
				band.genomeFragmentLength = len;
				band.bandCutSite = cutSite;
				genomeFragments.add(band);
				if(len > MAX_FRAGMENT_LENGTH)
					MAX_FRAGMENT_LENGTH = len;
				if(len < MIN_FRAGMENT_LENGTH)
					MIN_FRAGMENT_LENGTH = len;
			}
			lastSite = cutSite.getFivePrime();
		}

		len = genomeLength-lastSite+firstSiteEnd;
		FragmentBand band = new FragmentBand();
		band.genomeFragmentLength = len;
		band.bandCutSite = cutSites.get(0);
		genomeFragments.add(band);
		if(len > MAX_FRAGMENT_LENGTH)
			MAX_FRAGMENT_LENGTH = len;
		System.out.println(len.toString());
		
		MouseListener popupListener = new PopupListener();
		addMouseListener(popupListener);
		
		popup = new JPopupMenu();
		JRadioButtonMenuItem  linearScale = new JRadioButtonMenuItem ("Linear scale");
    popup.add(linearScale);
    linearScale.addActionListener(this);
    JRadioButtonMenuItem  logScale = new JRadioButtonMenuItem ("Log scale");
    popup.add(logScale);
    logScale.addActionListener(this);
    ButtonGroup group = new ButtonGroup();
    group.add(linearScale);
    group.add(logScale);
    linearScale.setSelected(true);
	}
	
	/**
	 * Override 
	 */
	public void paintComponent(Graphics g)
	{
		super.paintComponent(g);
		
		Graphics2D g2D = (Graphics2D)g;
		g2D.draw3DRect(marginWidth, (marginHeight/2), marginWidth, panelHeight-(marginHeight), true);
		
		int gelHeight = (panelHeight-(2*marginHeight));
		g2D.setColor(Color.blue);
	
		BasicStroke stroke  = new BasicStroke(1.f);
		BasicStroke stroke2 = new BasicStroke(2.f);
		g2D.setStroke(stroke);
		
		for(int i=0; i<genomeFragments.size(); i++)
		{
			int fragmentLength = genomeFragments.get(i).genomeFragmentLength;

			final int y;
			if(isDrawLog())
			  y = getLogValue(fragmentLength, marginHeight, gelHeight);
			else
				y = gelHeight+marginHeight-(int)( ((float)(gelHeight)/
						(float)(MAX_FRAGMENT_LENGTH-MIN_FRAGMENT_LENGTH)) * fragmentLength );
		
			if(!genomeFragments.get(i).bandCutSite.isHighlighted())
			{
				g2D.setStroke(stroke);
				g2D.setColor(Color.blue);
			}
			else
			{
				g2D.setStroke(stroke2);
				g2D.setColor(Color.yellow);
			}
			g2D.drawLine(marginWidth, y, marginWidth+marginWidth, y);
		}
	}
	
  private static final double LOG10SCALE = 1.d/Math.log(10);
  
  /**
   * Get base 10 commons log
   * @param val
   * @return
   */
  private static double log10(double val) { return Math.log(val) * LOG10SCALE; }
	
  private int getLogValue(int val, int gelStart, int gelRange) 
  {
      double log_low  = log10(MIN_FRAGMENT_LENGTH);
      double log_high = log10(MAX_FRAGMENT_LENGTH);
      double log_val  = log10(val);

      double log_unit = ((double)gelRange) / (log_high - log_low);
      return (int)((double)(gelRange + gelStart) - ((log_val - log_low) * log_unit));
  }
 
  private boolean isDrawLog()
	{
		return drawLog;
	}

  private void setDrawLog(boolean drawLog)
	{
		this.drawLog = drawLog;
	}
  
  class FragmentBand
  {
  	CutSite bandCutSite;
  	int genomeFragmentLength;
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
		JRadioButtonMenuItem radioButton = (JRadioButtonMenuItem)e.getSource();
		if(radioButton.isSelected())
		{
			setDrawLog(radioButton.getText().startsWith("Log"));
			repaint();
		}
	}
}