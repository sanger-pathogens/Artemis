package uk.ac.sanger.artemis.components.alignment;

import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import uk.ac.sanger.artemis.components.Utilities;


class SAMRecordFilter extends JPanel
{
  private static final long serialVersionUID = 1L;

  public SAMRecordFilter(final BamView bamView)
  {
    super();
    int nflags = SAMRecordFlagPredicate.FLAGS.length;
    setLayout(new GridLayout(nflags+1, 2));

    final JCheckBox flagCheck[] = new JCheckBox[nflags];
    final SAMRecordFlagPredicate predicate = bamView.getSamRecordFlagPredicate();
    
    for(int j=0; j<nflags; j++)
    {
      flagCheck[j] = new JCheckBox(
          SAMRecordFlagPredicate.FLAGS_DESCRUIPTION[j], false);
      
      if(predicate != null &&
         predicate.isFlagSet(SAMRecordFlagPredicate.FLAGS[j]))
        flagCheck[j].setSelected(true);
      
      flagCheck[j].addChangeListener(new ChangeListener()
      {
        public void stateChanged(ChangeEvent e)
        {
          filterChange(bamView, flagCheck);
        }            
      });
      add(flagCheck[j]);
    }

    final JFrame f = new JFrame("Filter Out Reads Based on Flag");
    JButton closeFrame = new JButton("Close");
    closeFrame.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        f.dispose();
      }
    });
    add(closeFrame);

    f.getContentPane().add(this);
    f.pack();
    Utilities.centreFrame(f);
    f.setVisible(true);
  } 
  
  private void filterChange(final BamView bamView,
                            final JCheckBox flagCheck[])
  {
    int nflags = SAMRecordFlagPredicate.FLAGS.length;
    int flagCombined = 0;
    for(int j=0; j<nflags; j++)
    {
      if(flagCheck[j].isSelected())
        flagCombined = flagCombined | SAMRecordFlagPredicate.FLAGS[j];
    }

    if(flagCombined == 0)
      bamView.setSamRecordFlagPredicate(null);
    else
      bamView.setSamRecordFlagPredicate(new SAMRecordFlagPredicate(flagCombined));
    bamView.repaint();
  }
  
}