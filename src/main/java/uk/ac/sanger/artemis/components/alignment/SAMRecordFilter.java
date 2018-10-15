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
import java.util.Vector;

import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSeparator;
import javax.swing.JTextField;

/**
 * Filter reads bases on their mapping quality (mapq) or the
 * flags that are set in the BAM file.
 */
class SAMRecordFilter extends JFrame
{
  private static final long serialVersionUID = 1L;

  public SAMRecordFilter(final BamView bamView)
  {
    super("Filter Reads");
    
    JPanel pane = (JPanel) getContentPane();
    pane.setLayout(new GridBagLayout());
    
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
    pane.add(new JLabel(" By Mappying Quality (mapq) cut-off:"),c);
    
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
    pane.add(cutOff,c);
    
    c.gridx = 1;
    final JButton setCutOff = new JButton("SET");
    setCutOff.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        setQualityCutOff(cutOff, bamView);
      }
    });
    pane.add(setCutOff,c);
    
    // FLAGS
    c.gridwidth = 2;
    c.gridx = 0;
    c.gridy = nrows++;
    pane.add(new JSeparator(),c);
    
    c.gridy = nrows++;
    pane.add(new JLabel(" By SAM FLAG column:"),c);
    c.gridy = nrows++;
    pane.add(new JLabel(" Select below to show or hide only the reads with "),c);
    c.gridy = nrows++;
    pane.add(new JLabel(" the flag set."),c);
    
    final JComboBox flagCombo[] = new JComboBox[nflags];
    
    final String[] items = {"", "SHOW", "HIDE"};
    c.gridwidth = 1;
    for(int j=0; j<nflags; j++)
    {
      flagCombo[j] = new JComboBox(items);
      if(SAMRecordFlagPredicate.FLAGS_DESCRIPTION[j].equalsIgnoreCase("Read Unmapped"))
        flagCombo[j].setSelectedItem("HIDE");
      
      flagCombo[j].addActionListener(new ActionListener() 
      {
        public void actionPerformed(ActionEvent e)
        {
          filterChange(bamView, flagCombo);
        } 
      });
      
      c.gridy = nrows++;
      c.gridx = 0;
      pane.add(flagCombo[j], c);
      c.gridx = 1;
      pane.add(new JLabel(SAMRecordFlagPredicate.FLAGS_DESCRIPTION[j]), c);
    }

    JButton closeFrame = new JButton("Close");
    closeFrame.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        setVisible(false);
      }
    });
    c.gridx = 1;
    c.gridy = nrows++;
    c.fill = GridBagConstraints.NONE;
    pane.add(closeFrame, c);

    pack();
    
    rightJustifyFrame(this);
    setVisible(true);
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

  private void filterChange(final BamView bamView, final JComboBox flagCheck[])
  {
    int nflags = SAMRecordFlagPredicate.FLAGS.length;
    int flagCombined = 0;
    Vector <SAMRecordPredicate> predicates = new Vector<SAMRecordPredicate>();
    
    for (int j = 0; j < nflags; j++)
    {
      String opt = (String) flagCheck[j].getSelectedItem();
      if (opt.equals("HIDE"))
        flagCombined = flagCombined | SAMRecordFlagPredicate.FLAGS[j];
      else if(opt.equals("SHOW"))
        predicates.add(new SAMRecordFlagPredicate( SAMRecordFlagPredicate.FLAGS[j], false ));
    }

    if (flagCombined == 0 && predicates.size() == 0)
      bamView.setSamRecordFlagPredicate(null);
    else
    {
      final SAMRecordPredicate predicate;
      if(predicates.size()  == 0)
        predicate = new SAMRecordFlagPredicate(flagCombined);
      else if(flagCombined == 0)
      {
        predicate = 
          new SAMRecordFlagConjunctionPredicate(predicates, SAMRecordFlagConjunctionPredicate.OR);
      }
      else
      {
        SAMRecordFlagPredicate p1 = new SAMRecordFlagPredicate(flagCombined);
        SAMRecordFlagConjunctionPredicate p2 = 
          new SAMRecordFlagConjunctionPredicate(predicates, SAMRecordFlagConjunctionPredicate.OR);
        
        predicate = new SAMRecordFlagConjunctionPredicate(p1, p2, 
            SAMRecordFlagConjunctionPredicate.OR);
      }
      bamView.setSamRecordFlagPredicate(predicate);
    }
    bamView.repaint();
  }
}