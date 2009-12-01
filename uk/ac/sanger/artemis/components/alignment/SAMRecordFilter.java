package uk.ac.sanger.artemis.components.alignment;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Point;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSeparator;
import javax.swing.JTextField;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

/**
 * Filter reads bases on their mapping quality (mapq) or the
 * flags that are set in the BAM file.
 */
class SAMRecordFilter extends JPanel
{
  private static final long serialVersionUID = 1L;

  public SAMRecordFilter(final BamView bamView)
  {
    super(new GridBagLayout());
    int nflags = SAMRecordFlagPredicate.FLAGS.length;
    GridBagConstraints c = new GridBagConstraints();
    
    int nrows = 0;
    c.ipadx = 5;
    c.ipady = 2;
    c.gridx = 0;
    c.gridwidth = 2;
    c.anchor = GridBagConstraints.WEST;
    c.gridy = nrows++;
    
    // MAPQ
    add(new JLabel("Mappying quality (mapq) cut-off:"),c);
    
    c.gridy = nrows++;
    c.gridwidth = 1;
    final JTextField cutOff = new JTextField(12);
    if(bamView.getSamRecordMapQPredicate() != null)
      cutOff.setText(Integer.toString(bamView.getSamRecordMapQPredicate().cutOff));
    
    cutOff.addKeyListener(new KeyAdapter()
    {
      public void keyPressed(KeyEvent e)
      {
        if (e.getKeyCode() == KeyEvent.VK_ENTER)
          setQualityCutOff(cutOff, bamView);
      }
    });
    add(cutOff,c);
    
    c.gridx = 1;
    final JButton setCutOff = new JButton("SET");
    setCutOff.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        setQualityCutOff(cutOff, bamView);
      }
    });
    add(setCutOff,c);
    
    // FLAGS
    c.gridwidth = 2;
    c.gridx = 0;
    c.gridy = nrows++;
    add(new JSeparator(),c);
    
    c.gridy = nrows++;
    add(new JLabel("Reads can be filtered using the FLAG column in the BAM file."),c);
    c.gridy = nrows++;
    add(new JLabel("Selecting any of the options below hides the reads with that"),c);
    c.gridy = nrows++;
    add(new JLabel("flag set."),c);
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
      
      c.gridy = nrows++;
      add(flagCheck[j], c);
    }

    final JFrame f = new JFrame("Filter Settings");
    JButton closeFrame = new JButton("Close");
    closeFrame.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        f.dispose();
      }
    });
    c.gridx = 1;
    c.gridy = nrows++;
    c.fill = GridBagConstraints.NONE;
    add(closeFrame, c);

    f.getContentPane().add(this);
    f.pack();
    
    rightJustifyFrame(f);
    f.setVisible(true);
  } 
  
  private void rightJustifyFrame(JFrame frame)
  {
    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
    int yPos = (screen.height - frame.getSize().height) / 2;
    int xPos = screen.width - frame.getSize().width;
    if(yPos < 10)
      yPos = 10;

    frame.setLocation(new Point(xPos, yPos));
  }
  
  private void setQualityCutOff(final JTextField cutOff,
                                final BamView bamView)
  {
    String cutOffStr = cutOff.getText().trim();
    if(cutOffStr.equals(""))
      bamView.setSamRecordMapQPredicate(null);
    else
    {
      try
      {
        SAMRecordMapQPredicate predicate = 
          new SAMRecordMapQPredicate(Integer.parseInt(cutOffStr));
        bamView.setSamRecordMapQPredicate(predicate);
      }
      catch(NumberFormatException nfe)
      {
        bamView.setSamRecordMapQPredicate(null);
      }
    }
    bamView.repaint();
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