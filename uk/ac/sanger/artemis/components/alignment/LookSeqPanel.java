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
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.GridLayout;
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
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JSeparator;
import javax.swing.JTextField;

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
  private int drawCrossHairAt = -1;
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
    JCheckBoxMenuItem optionCoverage = new JCheckBoxMenuItem("Coverage", false);
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
    
    JCheckBoxMenuItem optionIndel = new JCheckBoxMenuItem("Paired reads", true);
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
    
    JCheckBoxMenuItem optionPileup = new JCheckBoxMenuItem("Pileup", false);
    optionPileup.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        queryStr = queryStr.replaceFirst("view=\\w+", "view=pileup");
        setImage(urlStr, queryStr);
        revalidate();
      }
    });
    popup.add(optionPileup);
    
    ButtonGroup group = new ButtonGroup();
    group.add(optionCoverage);
    group.add(optionIndel);
    group.add(optionPileup);
    
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
    
    if(drawCrossHairAt > -1)
      drawCrossHair(g);
  }
  
  private void drawCrossHair(Graphics g)
  {
    int xpos = (int) ((((float)drawCrossHairAt/(lastEnd-lastStart)))*getWidth());
    g.drawLine(xpos, 0, xpos, getHeight());
    g.drawString(Integer.toString(drawCrossHairAt+lastStart), xpos+1, 10);
  }
  
  public void displayAdjustmentValueChanged(final DisplayAdjustmentEvent event)
  {   
    SwingWorker worker = new SwingWorker()
    {
      public Object construct()
      {
        if(lastStart == event.getStart() && lastEnd == event.getEnd())
          return null;
        try
        {
          Thread.currentThread().sleep(500);
        }
        catch (InterruptedException e)
        {
          e.printStackTrace();
        }
        
        if(event.getStart() != ((FeatureDisplay)event.getSource()).getForwardBaseAtLeftEdge())
          return null;
      
        lastStart = event.getStart();
        lastEnd   = event.getEnd();
        queryStr = queryStr.replaceFirst("from=\\d+", "from="+Integer.toString(event.getStart()));
        queryStr = queryStr.replaceFirst("to=\\d+", "to="+Integer.toString(event.getEnd()));    

        if(feature_display != null)
        {
          int width = feature_display.getSize().width;
          queryStr = queryStr.replaceFirst("width=\\d+", "width="+Integer.toString(width));
        }
        
        setImage(urlStr, queryStr);
        drawCrossHairAt = -1;
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
  
  private int getBasePositionAt(int position)
  {
    return (int) ((((float)position)/getWidth())*(lastEnd-lastStart));
  }
  
  final MouseListener mouseListener = new MouseAdapter()
  {
    /**
     *  Listen for mouse press events.
     **/
    public void mousePressed(MouseEvent event)
    {
      if(event.isPopupTrigger())
        popup.show(LookSeqPanel.this, event.getX(), event.getY());
      else if(event.getClickCount() == 1 && 
              event.getButton() == MouseEvent.BUTTON1)
      {
        drawCrossHairAt = getBasePositionAt(event.getX());
        repaint();
      }
      else if(event.getClickCount() == 2)
      {
        drawCrossHairAt = -1;
        repaint();
      }
    }
  };

  public void setFeatureDisplay(FeatureDisplay feature_display)
  {
    this.feature_display = feature_display;
  }
  
  public void showOptions()
  {
    final JPanel optionsPanel = new JPanel(new GridBagLayout());
    GridBagConstraints c = new GridBagConstraints();
    final JTextField urlStrField = new JTextField(urlStr);
    c.gridy = 0;
    c.gridx = 0;
    c.anchor = GridBagConstraints.EAST;
    optionsPanel.add(new JLabel("Server:"), c);
    c.gridx = 1;
    c.anchor = GridBagConstraints.WEST;
    optionsPanel.add(urlStrField, c);
    
    String sampleStr = getText(queryStr,"chr=");
    final JTextField sampleField = new JTextField(sampleStr,40);
    c.gridy = 1;
    c.gridx = 0;
    c.anchor = GridBagConstraints.EAST;
    optionsPanel.add(new JLabel("Sample:"), c);
    c.gridx = 1;
    c.anchor = GridBagConstraints.WEST;
    optionsPanel.add(sampleField, c);
    
    String laneStr = getText(queryStr,"lane=");
    final JTextField laneField = new JTextField(laneStr,40);
    c.gridy = 2;
    c.gridx = 0;
    c.anchor = GridBagConstraints.EAST;
    optionsPanel.add(new JLabel("Lane:"), c);
    c.gridx = 1;
    c.anchor = GridBagConstraints.WEST;
    optionsPanel.add(laneField, c);
    
    String window_options[] = { "Display" };
    int select = JOptionPane.showOptionDialog(null, 
        optionsPanel,
        "Lookseq Options", JOptionPane.DEFAULT_OPTION,
        JOptionPane.QUESTION_MESSAGE, null, window_options,
        window_options[0]);
    
    urlStr = urlStrField.getText().trim();
    queryStr = queryStr.replaceFirst(
        "chr=\\w+", "chr="+sampleField.getText().trim());
    
    queryStr = queryStr.replaceFirst(
        "lane=\\w+", "lane="+laneField.getText().trim());
  }
  
  private String getText(String str, String subStr)
  {
    int beginIndex = str.indexOf(subStr);
    if(beginIndex == -1)
      return null;
    
    int endIndex = str.indexOf("&", beginIndex);
    if(endIndex == -1)
      endIndex = str.length();
    return str.substring(beginIndex+subStr.length(), endIndex);
  }
  
  public String getQueryStr()
  {
    return queryStr;
  }
  
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
}

