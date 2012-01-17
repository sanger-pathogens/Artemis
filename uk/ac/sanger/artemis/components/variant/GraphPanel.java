/*
 * created: 2011
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2011  Genome Research Limited
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
package uk.ac.sanger.artemis.components.variant;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.geom.GeneralPath;
import java.io.IOException;
import java.text.DecimalFormat;

import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JTextField;

import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.components.variant.BCFReader.BCFReaderIterator;

public class GraphPanel extends JPanel
{
  private static final long serialVersionUID = 1L; 
  protected boolean autoWinSize = true;
  protected int userWinSize = 1;
  protected JPopupMenu popup = new JPopupMenu();
  private VCFview vcfView;
  
  private static String TYPES[] = { "SNP", "DP" };
  private int type = 0;  // 0 - SNP; 1 - DP
  
  public GraphPanel(final VCFview vcfView)
  {
    super();
    this.vcfView = vcfView;
    setBackground(Color.white);
    initPopupMenu(this);
  }
  
  /**
   * Override
   */
  protected void paintComponent(Graphics g)
  {
    super.paintComponent(g);
    Graphics2D g2 = (Graphics2D)g;  
    
    float pixPerBase = vcfView.getPixPerBaseByWidth();
    int sbeg = vcfView.getBaseAtStartOfView();
    int send = sbeg+vcfView.getBasesInView();
    
    int windowSize;
    if(autoWinSize)
    {
      windowSize = (vcfView.getBasesInView()/300);
      userWinSize = windowSize;
    }
    else
      windowSize = userWinSize;
    if(windowSize < 1)
      windowSize = 1;

    int nBins = Math.round((send-sbeg+1.f)/windowSize);
    int snpCount[] = new int[nBins];
    int max = 0;
    for(int i=0; i<snpCount.length; i++)
      snpCount[i] = 0;
    
    AbstractVCFReader readers[] = vcfView.getVcfReaders();
    FeatureVector features = vcfView.getCDSFeaturesInRange(sbeg, send);
    for(int i=0; i<readers.length; i++)
    {
      max = countAll(i, readers[i], sbeg, send, windowSize, nBins, features, snpCount, max);
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

    DecimalFormat df = new DecimalFormat("0.00##");
    String maxStr = TYPES[type]+" : "+df.format((float)max/(float)windowSize);
    FontMetrics fm = getFontMetrics(getFont());
    g2.setColor(Color.black);
    
    int xpos = getWidth() - fm.stringWidth(maxStr) - 10;
    g2.drawString(maxStr, xpos, fm.getHeight());
  }
  
  private int countAll(int vcfIndex, AbstractVCFReader reader, int sbeg, int send, int windowSize, 
      int nBins, FeatureVector features, int snpCount[], int max)
  {
    if(vcfView.isConcatenate())
    {
      String[] contigs = reader.getSeqNames();
      for(int j=0; j<contigs.length; j++)
      {
        int offset = vcfView.getSequenceOffset(contigs[j]);
        int nextOffset;
        if(j<contigs.length-1)
          nextOffset = vcfView.getSequenceOffset(contigs[j+1]);
        else
          nextOffset = vcfView.seqLength;
        
        if( (offset >= sbeg && offset < send) ||
            (offset < sbeg && sbeg < nextOffset) )
        {
          int thisStart = sbeg - offset;
          if(thisStart < 1)
            thisStart = 1;
          int thisEnd   = send - offset;
          
          max = count(vcfIndex, reader, contigs[j], thisStart, thisEnd, windowSize, nBins, features, snpCount, max);
        }
      }
      return max;
    }
    else
      return count(vcfIndex, reader, vcfView.getChr(), sbeg, send, windowSize, nBins, features, snpCount, max);
  }
  
  private int count(int vcfIndex, AbstractVCFReader reader, String chr, int sbeg, int send, int windowSize, 
      int nBins, FeatureVector features, int snpCount[], int max)
  {
    if(reader instanceof BCFReader)
    {
      try
      {
        BCFReader bcfReader = (BCFReader)reader;
        BCFReaderIterator it = bcfReader.query(chr, sbeg, send);

        VCFRecord record;
        while((record = it.next()) != null)
        {
          max = calc(vcfIndex, record, features, reader, windowSize, nBins, snpCount, max);
        }
      }
      catch (IOException e)
      {
        e.printStackTrace();
      }
    }
    else
    {
      TabixReader.Iterator iter = 
        ((TabixReader)reader).query(chr+":"+sbeg+"-"+send); // get the iterator
      if (iter == null)
        return max;
      try
      {
        String s;
        while ((s = iter.next()) != null)
        {
          VCFRecord record = VCFRecord.parse(s, reader.getNumberOfSamples());
          max = calc(vcfIndex, record, features, reader, windowSize, nBins, snpCount, max);
        }
      }
      catch (IOException e)
      {
        e.printStackTrace();
      }
    }
    return max;
  }
  
  private int calc(int vcfIndex, VCFRecord record, FeatureVector features,
      AbstractVCFReader reader, int windowSize, int nBins, int snpCount[], int max)
  {
    int pos = record.getPos() + vcfView.getSequenceOffset(record.getChrom());
    if(!vcfView.showVariant(record, features, pos, reader, -1, vcfIndex))
      return max;
    
    int bin = (int)((pos-vcfView.getBaseAtStartOfView()) / windowSize);
    if(bin < 0 || bin > nBins-1)
      return max;

    if(type == 1)
    {
      int dp = 0;
      try
      {
        dp = Integer.parseInt(record.getInfoValue("DP"));
      }
      catch(Exception e){ e.printStackTrace(); }
      snpCount[bin]+=dp;
    }
    else
      snpCount[bin]+=1;
    
    if(snpCount[bin] > max)
      max = snpCount[bin];

    return max;
  }
  
  private void initPopupMenu(final JPanel graphPanel)
  {
    final JMenuItem setScale = new JMenuItem("Set the Window Size...");
    setScale.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent _)
      {
        GridBagLayout gridbag = new GridBagLayout();
        GridBagConstraints c = new GridBagConstraints();
        JPanel pane = new JPanel(gridbag);
        final JTextField newWinSize = new JTextField(Integer.toString(userWinSize), 10);
        final JLabel lab = new JLabel("Window size:");
        c.gridy = 0;
        pane.add(lab, c);
        pane.add(newWinSize, c);

        final JCheckBox autoSet = new JCheckBox("Automatically set window size", false);
        autoSet.addActionListener(new ActionListener()
        {
          public void actionPerformed(ActionEvent e)
          {
            lab.setEnabled(!autoSet.isSelected());
            newWinSize.setEnabled(!autoSet.isSelected());
          }
        });
        c.gridy = 1;
        c.gridwidth = GridBagConstraints.REMAINDER;
        pane.add(autoSet, c);

        String window_options[] = { "OK", "Cancel" };
        int select = JOptionPane.showOptionDialog(null, pane, "Window Size",
            JOptionPane.DEFAULT_OPTION, JOptionPane.QUESTION_MESSAGE, null,
            window_options, window_options[0]);

        if (select == 1)
          return;
        autoWinSize = autoSet.isSelected();
        try
        {
          userWinSize = Integer.parseInt(newWinSize.getText().trim());
        }
        catch (NumberFormatException nfe)
        {
          return;
        }
        graphPanel.repaint();
      }
    });
    popup.add(setScale);
    
    addMouseListener(new PopupListener());
  }
  
  /**
   * @param type the type to set
   */
  protected void setType(int type)
  {
    this.type = type;
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