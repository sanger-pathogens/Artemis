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
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Graphics;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.List;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollBar;


import net.sf.samtools.SAMRecord;

public class SAMRecordList extends JPanel
{
  private static final long serialVersionUID = 1L;
  private BamView bamView;
  private JScrollBar verticalScroll = new JScrollBar(JScrollBar.VERTICAL);

  public SAMRecordList(final BamView bamView)
  {
    this.bamView = bamView;
    setBackground(Color.white);

    final JFrame f = new JFrame("Reads");
    JPanel mainPanel = (JPanel) f.getContentPane();
    JPanel header = new JPanel(new FlowLayout(FlowLayout.LEADING));
    header.add(new JLabel("NAME :: COORDINATES :: LENGTH :: QUALITY"));
    JButton close = new JButton("CLOSE");
    
    mainPanel.add(verticalScroll, BorderLayout.EAST);
    mainPanel.add(this, BorderLayout.CENTER);
    mainPanel.add(header, BorderLayout.NORTH);
    mainPanel.add(close, BorderLayout.SOUTH);
    
    setPreferredSize(new Dimension(400,500));
    f.pack();
    f.setVisible(true);
    f.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
    
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
        int highlight = e.getY()/12 + verticalScroll.getValue();
        if(highlight > bamView.getReadsInView().size())
          return;

        bamView.setHighlightSAMRecord( bamView.getReadsInView().get(highlight) );
        bamView.repaint();
        repaint();
      }
    });
    
    close.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent arg0)
      {
        f.dispose();
      }
    });
  }

  protected void paintComponent(Graphics g)
  {
    super.paintComponent(g);

    int nrecordsShown = getHeight()/12;
    List<SAMRecord> readsInView = bamView.getReadsInView();
    verticalScroll.setMaximum(readsInView.size()-nrecordsShown);

    int fst = verticalScroll.getValue();
    int lst = nrecordsShown+fst+1;
    
    if(lst > readsInView.size())
      lst = readsInView.size();
    
    String highlightedSAMRecord = (bamView.getHighlightSAMRecord() == null ? 
        null : bamView.getHighlightSAMRecord().getReadName());
    
    String fmt = getFormatString(fst, lst, readsInView);
    for(int i=fst; i<lst; i++)
    {
      SAMRecord thisRead = readsInView.get(i);
      
      if(highlightedSAMRecord != null && highlightedSAMRecord.equals(thisRead.getReadName()))
      {
        g.setColor(Color.pink);
        g.fillRect(0, ((i-fst)*12)+1, getWidth(), 12);
        g.setColor(Color.black);
      }
      
      String readStr = String.format(fmt,
                       thisRead.getReadName(),
                       thisRead.getAlignmentStart()+".."+thisRead.getAlignmentEnd(),
                       thisRead.getReadLength(),
                       thisRead.getMappingQuality(),
                       (thisRead.getReadNegativeStrandFlag() ? "<--" : "-->"));
      g.drawString(readStr, 5, (((i-fst)*12)+12));
    }
  }
  
  private String getFormatString(int fst, int lst, List<SAMRecord> readsInView)
  {
    int nameWidth  = 10;
    int coordWidth = 10;
    int lenWidth   = 7;
    int qualWidth  = 7;
    
    for(int i=fst; i<lst; i++)
    {
      SAMRecord thisRead = readsInView.get(i);
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
    return "%1$-"+nameWidth+"s| %2$-"+coordWidth+"s| %3$-"+lenWidth+"s| %4$-"+qualWidth+"s| %5$-4s";
  }
}
