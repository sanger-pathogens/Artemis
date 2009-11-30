package uk.ac.sanger.artemis.components.alignment;

import java.awt.GridLayout;

import javax.swing.JCheckBox;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;


class SAMRecordFilter extends JPanel
{
  private static final long serialVersionUID = 1L;

  public SAMRecordFilter(final BamView bamView)
  {
    super();
    int nflags = SAMRecordFlagPredicate.FLAGS.length;
    setLayout(new GridLayout(nflags, 2));

    final JCheckBox flagCheck[] = new JCheckBox[nflags];
    for(int j=0; j<nflags; j++)
    {
      flagCheck[j] = new JCheckBox(
          SAMRecordFlagPredicate.FLAGS_DESCRUIPTION[j], false);
      flagCheck[j].addChangeListener(new ChangeListener()
      {
        public void stateChanged(ChangeEvent e)
        {
          filterChange(bamView, flagCheck);
        }            
      });
      add(flagCheck[j]);
    }
    
    int status = JOptionPane.showConfirmDialog(bamView, 
        this, "Filter Out Reads Based on Flag", 
        JOptionPane.OK_CANCEL_OPTION);

    if(status != JOptionPane.OK_OPTION)
      return;
  } 
  
  private void filterChange(final BamView bamView,
                            final JCheckBox flagCheck[])
  {
    int nflags = SAMRecordFlagPredicate.FLAGS.length;
    int flagsChecked = 0;
    for(int j=0; j<nflags; j++)
    {
      if(flagCheck[j].isSelected())
        flagsChecked++;
    }
    
    bamView.setSamRecordFlagPredicate(null);
    
    if(flagsChecked == 0)
    {
      bamView.repaint();
      return;
    }
    
    int flagsOn[] = new int[flagsChecked];
    int num = 0;
    for(int j=0; j<nflags; j++)
    {
      if(flagCheck[j].isSelected())
        flagsOn[num++] = SAMRecordFlagPredicate.FLAGS[j];
    }
    
    bamView.setSamRecordFlagPredicate(new SAMRecordFlagPredicate(flagsOn));
    bamView.repaint();
  }
  
}