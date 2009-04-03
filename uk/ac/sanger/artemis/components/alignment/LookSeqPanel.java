/* LookSeqFrame.java
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

import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.net.MalformedURLException;
import java.net.URL;

import javax.swing.ButtonGroup;
import javax.swing.ImageIcon;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JFrame;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JSeparator;

import uk.ac.sanger.artemis.components.DisplayAdjustmentEvent;
import uk.ac.sanger.artemis.components.DisplayAdjustmentListener;
import uk.ac.sanger.artemis.components.FeatureDisplay;
import uk.ac.sanger.artemis.components.SwingWorker;

public class LookSeqPanel extends JPanel
                          implements DisplayAdjustmentListener
{
  private static final long serialVersionUID = 1L;
  /** image icon */
  private ImageIcon ii;
  private String urlStr;
  private String queryStr;
  private FeatureDisplay feature_display;
  private int lastStart = -1;
  private int lastEnd   = -1;
  final JPopupMenu popup  = new JPopupMenu("Plot Options");

  public LookSeqPanel()
  {
    super();
    setUpPopupMenu();
  }
  
  public LookSeqPanel(URL url)
  {
    super();
    ii = new ImageIcon(url);
    setUpPopupMenu();
  }
  
  /**
   * 
   * @param urlStr
   * @param query
   */
  public LookSeqPanel(final String urlStr, final String queryStr)
  {
    super();
    setImage(urlStr, queryStr);
    setUpPopupMenu();
  }
  
  private void setUpPopupMenu()
  {
    JCheckBoxMenuItem optionCoverage = new JCheckBoxMenuItem("Show coverage", false);
    optionCoverage.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        queryStr = queryStr.replaceFirst("view=\\w+", "view=coverage");
        setImage(urlStr, queryStr);
        revalidate();
      }
    });
    popup.add(optionCoverage);
    
    JCheckBoxMenuItem optionIndel = new JCheckBoxMenuItem("Show paired reads", true);
    
    optionIndel.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        queryStr = queryStr.replaceFirst("view=\\w+", "view=indel");
        setImage(urlStr, queryStr);
        revalidate();
      }
    });
    popup.add(optionIndel);
    
    ButtonGroup group = new ButtonGroup();
    group.add(optionCoverage);
    group.add(optionIndel);
    
    popup.add(new JSeparator());
    final JCheckBoxMenuItem perfectIndel = 
      new JCheckBoxMenuItem("Show perfect paired matches", true);
    perfectIndel.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        setDisplayOption(perfectIndel.isSelected(), "perfect");
      }
    });
    popup.add(perfectIndel);
    
    final JCheckBoxMenuItem snps = 
      new JCheckBoxMenuItem("Show paired reads with SNPs", true);
    snps.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        setDisplayOption(snps.isSelected(), "snps");
      }
    });
    popup.add(snps);
    
    final JCheckBoxMenuItem inversions = 
      new JCheckBoxMenuItem("Show inversions", true);
    inversions.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        setDisplayOption(inversions.isSelected(), "inversions");
      }
    });
    popup.add(inversions);
    
    final JCheckBoxMenuItem single = 
      new JCheckBoxMenuItem("Show single reads", true);
    single.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        setDisplayOption(single.isSelected(), "single");
      }
    });
    popup.add(single);
    
    final JCheckBoxMenuItem inversions_ext = 
      new JCheckBoxMenuItem("Show inversions_ext", true);
    inversions_ext.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        setDisplayOption(inversions_ext.isSelected(), "inversions_ext");
      }
    });
    popup.add(inversions_ext);
    
    final JCheckBoxMenuItem pairlinks = 
      new JCheckBoxMenuItem("Show link pairs", true);
    pairlinks.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        setDisplayOption(inversions_ext.isSelected(), "pairlinks");
      }
    });
    popup.add(pairlinks);
    
    final JCheckBoxMenuItem potsnps = 
      new JCheckBoxMenuItem("Show known SNPs", true);
    potsnps.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        setDisplayOption(inversions_ext.isSelected(), "potsnps");
      }
    });
    popup.add(potsnps);
    
    final JCheckBoxMenuItem uniqueness = 
      new JCheckBoxMenuItem("Show uniqueness", true);
    uniqueness.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        setDisplayOption(inversions_ext.isSelected(), "uniqueness");
      }
    });
    popup.add(uniqueness);
    
    addMouseListener(mouseListener);
  }

  public void setImage(final String urlStr, final String queryStr)
  {
    this.urlStr = urlStr;
    this.queryStr = queryStr;
    
    
    try
    {
      ii = new ImageIcon(new URL(urlStr+queryStr));
      setPreferredSize(new Dimension(ii.getIconWidth(),
          ii.getIconHeight()));
    }
    catch (MalformedURLException e)
    {
      e.printStackTrace();
    }
  }
  
  /**
  * Override the paintComponent
  */
  public void paintComponent(Graphics g)
  {
    super.paintComponent(g);
    if(ii != null)
      ii.paintIcon(this,g,0,0);
  }
  
  
  public void displayAdjustmentValueChanged(DisplayAdjustmentEvent event)
  {   
    if(lastStart == event.getStart() && lastEnd == event.getEnd())
      return;
    
    lastStart = event.getStart();
    lastEnd   = event.getEnd();
    queryStr = queryStr.replaceFirst("from=\\d+", "from="+Integer.toString(event.getStart()));
    queryStr = queryStr.replaceFirst("to=\\d+", "to="+Integer.toString(event.getEnd()));    

    if(feature_display != null)
    {
      int width = feature_display.getSize().width;
      queryStr = queryStr.replaceFirst("width=\\d+", "width="+Integer.toString(width));
    }
    //System.out.println("displayAdjustmentValueChanged() "+
    //    event.getStart()+".."+event.getEnd());
    displayImage();
  }
  
  private void displayImage()
  {
    SwingWorker worker = new SwingWorker()
    {
      public Object construct()
      {
        setImage(urlStr, queryStr);
        //System.out.println(queryStr);
        repaint();
        return null;
      }
    };
    worker.start();
  }
  
  private void setDisplayOption(boolean selected, String option)
  {
    if(selected)
      queryStr = queryStr.replaceFirst("display=", "display=\\|"+option);
    else
      queryStr = queryStr.replaceAll("\\|"+option, "");
    
    setImage(urlStr, queryStr);
    repaint();
    revalidate();
  }
  
  final MouseListener mouseListener = new MouseAdapter()
  {
    /**
     *  Listen for mouse press events.
     **/
    public void mousePressed(MouseEvent event)
    {
      if(!event.isPopupTrigger())
        return;
      popup.show(LookSeqPanel.this, event.getX(), event.getY());
    }
  };
  
  public static void main(String args[])
  {
    String urlStr = "http://www.sanger.ac.uk/cgi-bin/teams/team112/lookseq/get_data.pl?";
    String queryStr = "from=157682&to=479328&chr=MAL1&output=image&width=1024&lane=sample_2a&view=indel&display=|perfect|snps|inversions|pairlinks|potsnps|uniqueness|&debug=0";
    LookSeqPanel lookseq = new LookSeqPanel(urlStr, queryStr);
    
    JFrame f = new JFrame();
    f.getContentPane().add(lookseq);
    f.pack();
    f.setVisible(true);
  }

  public void setFeatureDisplay(FeatureDisplay feature_display)
  {
    this.feature_display = feature_display;
  }
}

