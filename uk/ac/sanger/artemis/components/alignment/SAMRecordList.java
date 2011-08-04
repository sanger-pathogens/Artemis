package uk.ac.sanger.artemis.components.alignment;

import java.awt.BorderLayout;
import java.awt.Graphics;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollBar;

import net.sf.samtools.SAMRecord;

public class SAMRecordList extends JPanel
{
  private List<SAMRecord> readsInView;
  private JScrollBar verticalScroll = new JScrollBar(JScrollBar.VERTICAL);
  
  public SAMRecordList(List<SAMRecord> readsInView)
  {
    this.readsInView = readsInView;
    JFrame f = new JFrame("Reads");
    JPanel mainPanel = (JPanel) f.getContentPane();
    mainPanel.add(verticalScroll, BorderLayout.EAST);
    mainPanel.add(this, BorderLayout.CENTER);
    f.pack();
    f.setVisible(true);
  }

  protected void paintComponent(Graphics g)
  {
    super.paintComponent(g);
    for(int i=0; i<readsInView.size(); i++)
    {
      SAMRecord thisRead = readsInView.get(i);
      String readStr = thisRead.getReadName()+"\t"+
                       thisRead.getAlignmentStart()+".."+thisRead.getAlignmentEnd()+"\t"+
                       thisRead.getReadLength()+"\t"+
                       thisRead.getMappingQuality();
      g.drawString(readStr, 10, ((i+12)*12));
    }
  }
}
