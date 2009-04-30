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

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.color.ColorSpace;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.image.BufferedImage;
import java.awt.image.ColorConvertOp;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;

import javax.imageio.ImageIO;
import javax.swing.ButtonGroup;
import javax.swing.ImageIcon;
import javax.swing.JCheckBox;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JSeparator;
import javax.swing.JTextField;

import uk.ac.sanger.artemis.circular.TextFieldInt;
import uk.ac.sanger.artemis.components.DisplayAdjustmentEvent;
import uk.ac.sanger.artemis.components.DisplayAdjustmentListener;
import uk.ac.sanger.artemis.components.FeatureDisplay;
import uk.ac.sanger.artemis.components.SwingWorker;

public class LookSeqPanel extends JPanel
                          implements DisplayAdjustmentListener
{
  private static final long serialVersionUID = 1L;
  /** image icon */
  private BufferedImage ii;
  /** URL lookseq service */
  private String urlStr;
  /** the query options (appended to lookseq service url) */
  private String queryStr;

  private FeatureDisplay feature_display;
  private int lastStart = -1;
  private int lastEnd   = -1;
  private int drawCrossHairAt = -1;
  final JPopupMenu popup  = new JPopupMenu("Plot Options");
  private boolean isGreyedOut = false;
  private static org.apache.log4j.Logger logger4j = 
    org.apache.log4j.Logger.getLogger(LookSeqPanel.class);

  public LookSeqPanel()
  {
    super();
    setUpPopupMenu();
  }
  
  public LookSeqPanel(URL url)
  {
    super();
    try
    {
      ii = ImageIO.read(url);
    }
    catch (IOException e)
    {
      e.printStackTrace();
    }
    setUpPopupMenu();
  }
  
  public LookSeqPanel(final String urlStr, final String queryStr)
  {
    super();
    setImage(urlStr, queryStr);
    setUpPopupMenu();
  }
  
  /**
   * Build the popup menu.
   */
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
      new JCheckBoxMenuItem("Show single reads", false);
    single.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        setDisplayOption(single.isSelected(), "single");
      }
    });
    popup.add(single);
    
    final JCheckBoxMenuItem inversions_ext = 
      new JCheckBoxMenuItem("Show inversions_ext", false);
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
        setDisplayOption(pairlinks.isSelected(), "pairlinks");
      }
    });
    popup.add(pairlinks);
    
    final JCheckBoxMenuItem potsnps = 
      new JCheckBoxMenuItem("Show known SNPs", true);
    potsnps.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        setDisplayOption(potsnps.isSelected(), "potsnps");
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
    popup.add(new JSeparator());
    
    final JMenuItem setIndelMaxSize = new JMenuItem("Maximum InDel Size...");
    setIndelMaxSize.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        setMaxIndelSize();
      }
    });
    popup.add(setIndelMaxSize);
    addMouseListener(mouseListener);
  }
  
  public void setUrl(final String urlStr, final String queryStr)
  {
    this.urlStr = urlStr;
    this.queryStr = queryStr;
  }

  public void setImage(final String urlStr, final String queryStr)
  {
    this.urlStr = urlStr;
    this.queryStr = queryStr;
    
    try
    {
      ii  = ImageIO.read(new URL(urlStr+queryStr));

      setPreferredSize(new Dimension(ii.getWidth(), ii.getHeight()));
      logger4j.debug("LookSeq URL    :: "+urlStr+queryStr);
      logger4j.debug("Proxy Settings :: "+System.getProperty("http.proxyHost")+":"+
                                          System.getProperty("http.proxyPort"));
    }
    catch (MalformedURLException e)
    {
      e.printStackTrace();
    }
    catch (IOException e)
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
    {
      if(isGreyedOut)
      {
        ColorSpace csGray = ColorSpace.getInstance(ColorSpace.CS_GRAY);
        new ImageIcon(new ColorConvertOp(csGray, null).filter(ii, null)).paintIcon(this,g,0,0); 
      }
      else
        new ImageIcon(ii).paintIcon(this,g,0,0);
    }
    
    if(drawCrossHairAt > -1)
      drawCrossHair(g);
  }
  
  private void drawCrossHair(Graphics g)
  {
    int xpos = (int) ((((float)drawCrossHairAt/(lastEnd-lastStart)))*getWidth());
    g.drawLine(xpos, 0, xpos, getHeight());
    g.drawString(Integer.toString(drawCrossHairAt+lastStart), xpos+1, 10);
  }
  
  /**
   * Set the maximum size.
   */
  private void setMaxIndelSize()
  {
    JPanel panel = new JPanel(new GridBagLayout());
    GridBagConstraints c = new GridBagConstraints();
   
    String maxdist = getText(queryStr,"maxdist=");
    c.gridx = 0;
    c.gridy = 0;
    panel.add(new JLabel("Set maximum indel size"), c);
    
    TextFieldInt maxIndel = new TextFieldInt();
    maxIndel.setColumns(10);
    c.gridx = 1;
    c.anchor = GridBagConstraints.WEST;
    panel.add(maxIndel, c);
    
    JCheckBox defaultValue = new JCheckBox("use default"); 
    c.gridy = 1;
    panel.add(defaultValue, c);
    
    if(maxdist == null)
      defaultValue.setSelected(true);
    else
    {
      defaultValue.setSelected(false);
      maxIndel.setValue(Integer.parseInt(maxdist));
    }
    
    String window_options[] = { "Set", "Cancel" };
    int select = JOptionPane.showOptionDialog(null, 
        panel, "Lookseq Options", JOptionPane.DEFAULT_OPTION,
        JOptionPane.QUESTION_MESSAGE, null, window_options,
        window_options[0]);

    if(select == 1)
      return;
    if(defaultValue.isSelected())
      queryStr = queryStr.replaceFirst(
          "&maxdist=\\w+", "");
    else
    {
      if(queryStr.indexOf("maxdist=") > -1)
        queryStr = queryStr.replaceFirst(
            "&maxdist=\\w+", "&maxdist="+maxIndel.getText().trim());
      else
        queryStr = queryStr.concat("&maxdist="+maxIndel.getText().trim());
    }
    
    setImage(urlStr, queryStr);
    repaint();
    revalidate();
  }
  
  /**
   * Called when the display is changed.
   */
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
        {
          if(!isGreyedOut)
          {
            isGreyedOut = true;
            repaint();
          }
          return null;
        }
      
        isGreyedOut = false;
        lastStart = event.getStart();
        lastEnd   = event.getEnd();
        setDisplay(lastStart, lastEnd, event);
        return null;
      }
    };
    worker.start();
  }
  
  /**
   * Set the start and end base positions to display.
   * @param start
   * @param end
   * @param event
   */
  public void setDisplay(int start,
                         int end,
                         DisplayAdjustmentEvent event)
  {
    queryStr = queryStr.replaceFirst("from=\\d+", "from="+Integer.toString(start));
    queryStr = queryStr.replaceFirst("to=\\d+", "to="+Integer.toString(end));    

    if(feature_display != null)
    {
      int width = feature_display.getSize().width;
      
      int possible_last_base =
        (int)(feature_display.getForwardBaseAtLeftEdge() + 
              feature_display.getMaxVisibleBases());
      
      // scale the image when scrolled passed the sequence length
      if(possible_last_base > end)
        width = (int)(width * (float)((float)(end-start)/
                                      (float)(possible_last_base-start)));
      
      queryStr = queryStr.replaceFirst("width=\\d+", "width="+Integer.toString(width));
    }
    
    setImage(urlStr, queryStr);
    drawCrossHairAt = -1;
    repaint();
    if(event == null || 
       event.getType() == DisplayAdjustmentEvent.SCALE_ADJUST_EVENT)
      revalidate();
  }
  
  /**
   * Utility for changing display options.
   * @param selected
   * @param option
   */
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
  
  /**
   * Get the base position from the pixel position on the panel
   * @param position
   * @return
   */
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
  
  /**
   * Show Lookseq options
   */
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
    
    
    
    String proxyHost = "";
    if(System.getProperty("http.proxyHost") != null)
      proxyHost = System.getProperty("http.proxyHost");
    
    final JTextField proxyHostField = new JTextField(proxyHost,40);
    c.gridy = 3;
    c.gridx = 0;
    c.anchor = GridBagConstraints.EAST;
    optionsPanel.add(new JLabel("Proxy Host:"), c);
    c.gridx = 1;
    c.anchor = GridBagConstraints.WEST;
    optionsPanel.add(proxyHostField, c);
    
    
    String proxyPort = "";
    if(System.getProperty("http.proxyPort") != null)
      proxyPort = System.getProperty("http.proxyPort");
    
    final JTextField proxyPortField = new JTextField(proxyPort,40);
    c.gridy = 4;
    c.gridx = 0;
    c.anchor = GridBagConstraints.EAST;
    optionsPanel.add(new JLabel("Proxy Port:"), c);
    c.gridx = 1;
    c.anchor = GridBagConstraints.WEST;
    optionsPanel.add(proxyPortField, c);
    
    
    String window_options[] = { "Display" };
    JOptionPane.showOptionDialog(null, 
        optionsPanel,
        "Lookseq Options", JOptionPane.DEFAULT_OPTION,
        JOptionPane.QUESTION_MESSAGE, null, window_options,
        window_options[0]);
    
    if(!proxyHostField.getText().trim().equals(""))
    {
      System.getProperties().put("http.proxyHost", proxyHostField.getText().trim());
      System.getProperties().put("http.proxyPort", proxyPortField.getText().trim());
      System.getProperties().put("proxySet","true");
    }
    
    urlStr = urlStrField.getText().trim();
    queryStr = queryStr.replaceFirst(
        "chr=[^&]+", "chr="+sampleField.getText().trim());
    
    queryStr = queryStr.replaceFirst(
        "lane=[^&]+", "lane="+laneField.getText().trim());
  }
  
  /**
   * Return the value of a sub-string in a string
   * @param str
   * @param subStr
   * @return
   */
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

