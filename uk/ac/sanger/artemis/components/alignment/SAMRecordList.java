package uk.ac.sanger.artemis.components.alignment;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollBar;


import net.sf.samtools.SAMRecord;

public class SAMRecordList extends JPanel
{
  private static final long serialVersionUID = 1L;
  private List<SAMRecord> readsInView;
  private JScrollBar verticalScroll = new JScrollBar(JScrollBar.VERTICAL);
  
  public SAMRecordList(List<SAMRecord> readsInView)
  {
    this.readsInView = readsInView;
    setBackground(Color.white);

    JFrame f = new JFrame("Reads");
    JPanel mainPanel = (JPanel) f.getContentPane();
    mainPanel.add(verticalScroll, BorderLayout.EAST);
    mainPanel.add(this, BorderLayout.CENTER);
    setPreferredSize(new Dimension(400,500));
    f.pack();
    f.setVisible(true);
  }

  protected void paintComponent(Graphics g)
  {
    super.paintComponent(g);
    for(int i=0; i<readsInView.size(); i++)
    {
      SAMRecord thisRead = readsInView.get(i);
      String readStr = String.format("%1$-60s|%2$-15s|%3$-15s|%4$-15s",
                       thisRead.getReadName()+"\t",
                       thisRead.getAlignmentStart()+".."+thisRead.getAlignmentEnd()+"\t",
                       thisRead.getReadLength()+"\t",
                       thisRead.getMappingQuality());
      g.drawString(readStr, 5, ((i*12)+12));
    }
  }
}
