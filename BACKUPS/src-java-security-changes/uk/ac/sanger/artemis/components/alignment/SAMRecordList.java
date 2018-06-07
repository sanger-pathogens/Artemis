/* SAMRecordList
 *
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
 */
package uk.ac.sanger.artemis.components.alignment;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Graphics;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.List;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollBar;

import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.components.DisplayAdjustmentEvent;
import uk.ac.sanger.artemis.components.DisplayAdjustmentListener;
import uk.ac.sanger.artemis.components.SwingWorker;
import uk.ac.sanger.artemis.components.alignment.BamViewRecord;


import htsjdk.samtools.SAMRecord;

public class SAMRecordList extends JPanel
                           implements DisplayAdjustmentListener
{
  private static final long serialVersionUID = 1L;
  private BamView bamView;
  private JScrollBar verticalScroll = new JScrollBar(JScrollBar.VERTICAL);
  private int font_size =12;
  private int highlight = -1;
  
  public SAMRecordList(final BamView bamView)
  {
    this.bamView = bamView;
    setBackground(Color.white);

    font_size = Options.getOptions().getFont().getSize();

    final JFrame f = new JFrame("Reads");
    JPanel mainPanel = (JPanel) f.getContentPane();
    JPanel header = new JPanel(new FlowLayout(FlowLayout.LEADING));
    header.add(new JLabel("NAME :: COORDINATES :: LENGTH :: QUALITY :: DIRECTION"));
    JButton close = new JButton("CLOSE");
    
    mainPanel.add(verticalScroll, BorderLayout.EAST);
    mainPanel.add(this, BorderLayout.CENTER);
    mainPanel.add(header, BorderLayout.NORTH);
    mainPanel.add(close, BorderLayout.SOUTH);
    
    setPreferredSize(new Dimension(450,500));
    f.pack();
    f.setVisible(true);
    f.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
    f.addWindowListener(new WindowAdapter()
    {
      public void windowClosed(WindowEvent e)
      {
        windowClose(f);
      }
    });
    
    verticalScroll.addAdjustmentListener(new AdjustmentListener(){
      public void adjustmentValueChanged(AdjustmentEvent arg0)
      {
         repaint();
      }
    });
    
    addMouseListener(new MouseAdapter()
    {
      public void mouseClicked(MouseEvent e)
      {
        highlight = e.getY()/font_size + verticalScroll.getValue();
        if(highlight >= bamView.getReadsInView().size())
          return;

        setHighlight();
        if(e.getClickCount() > 1)
        {
          setCursor(new Cursor(Cursor.WAIT_CURSOR));
          BamView.openFileViewer(bamView.getHighlightSAMRecord().sam,
            bamView.getMate(bamView.getHighlightSAMRecord()),
            bamView.bamList);
          setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
        }
      }
    });

    setFocusable(true);
    addKeyListener(new KeyAdapter()
    {      
      public void keyPressed(KeyEvent e)
      {
        int keyCode = e.getKeyCode();
        switch( keyCode ) { 
          case KeyEvent.VK_UP:
            if(highlight == 0)
              break;
            highlight--;
            setHighlight();
            
            if(highlight < verticalScroll.getValue())
              verticalScroll.setValue(highlight);
            break;
          case KeyEvent.VK_DOWN:
            if(highlight >= bamView.getReadsInView().size()-1)
              break;
            highlight++;
            setHighlight();
            
            if(highlight > verticalScroll.getValue()+verticalScroll.getVisibleAmount())
              verticalScroll.setValue(highlight);
            break;
        }

      }
    });

    
    close.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent arg0)
      {
        windowClose(f);
      }
    });
    
    if(bamView.getFeatureDisplay() != null)
      bamView.getFeatureDisplay().addDisplayAdjustmentListener(this);
  }
  
  
  private void setHighlight()
  {
    bamView.setHighlightSAMRecord( bamView.getReadsInView().get(highlight) );
    bamView.repaint();
    repaint();
  }

  protected void paintComponent(Graphics g)
  {
    super.paintComponent(g);

    int nrecordsShown = getHeight()/font_size;
    List<BamViewRecord> readsInView = bamView.getReadsInView();
    verticalScroll.setMaximum(readsInView.size());
    verticalScroll.setVisibleAmount(nrecordsShown);

    int fst = verticalScroll.getValue();
    int lst = nrecordsShown+fst+1;
    
    if(lst > readsInView.size())
      lst = readsInView.size();
    
    String highlightedSAMRecord = (bamView.getHighlightSAMRecord() == null ? 
        null : bamView.getHighlightSAMRecord().sam.getReadName());
    
    String fmt = getFormatString(fst, lst, readsInView);
    for(int i=fst; i<lst; i++)
    {
      SAMRecord thisRead = readsInView.get(i).sam;
      
      if(highlightedSAMRecord != null && highlightedSAMRecord.equals(thisRead.getReadName()))
      {
        g.setColor(Color.pink);
        g.fillRect(0, ((i-fst)*font_size)+1, getWidth(), font_size);
        g.setColor(Color.black);
      }
      
      String readStr = String.format(fmt,
                       thisRead.getReadName(),
                       thisRead.getAlignmentStart()+".."+thisRead.getAlignmentEnd(),
                       thisRead.getReadLength(),
                       thisRead.getMappingQuality(),
                       (thisRead.getReadNegativeStrandFlag() ? "-" : "+"));
      g.drawString(readStr, 5, (((i-fst)*font_size)+font_size));
    }
  }
  
  private String getFormatString(int fst, int lst, List<BamViewRecord> readsInView)
  {
    int nameWidth  = 10;
    int coordWidth = 10;
    int lenWidth   = 4;
    int qualWidth  = 4;
    
    for(int i=fst; i<lst; i++)
    {
      SAMRecord thisRead = readsInView.get(i).sam;
      int thisWidth = thisRead.getReadName().length();
      if(thisWidth > nameWidth)
        nameWidth = thisWidth;
      
      thisWidth = (thisRead.getAlignmentStart()+".."+thisRead.getAlignmentEnd()).length();
      if(thisWidth > coordWidth)
        coordWidth = thisWidth;
      
      thisWidth = (String.valueOf(thisRead.getReadLength()) ).length();
      if(thisWidth > lenWidth)
        lenWidth = thisWidth;
      
      thisWidth = (String.valueOf(thisRead.getMappingQuality())).length();
      if(thisWidth > qualWidth)
        qualWidth = thisWidth;
    }

    nameWidth++;
    coordWidth++;
    return "%1$-"+nameWidth+"s| %2$-"+coordWidth+"s| %3$-"+lenWidth+"s| %4$-"+qualWidth+"s| %5$-2s";
  }

  public void displayAdjustmentValueChanged(DisplayAdjustmentEvent event)
  {
    SwingWorker worker = new SwingWorker()
    {
      public Object construct()
      {
        try
        {
          Thread.sleep(1500);
        }
        catch (InterruptedException e){}
        
        repaint();
        return null;
      }
    };
    worker.start();
  }
  
  private void windowClose(JFrame f)
  {
    f.dispose();
    if(bamView.getFeatureDisplay() != null)
      bamView.getFeatureDisplay().removeDisplayAdjustmentListener(SAMRecordList.this);
  }
}
