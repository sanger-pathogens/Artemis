package uk.ac.sanger.artemis.components.alignment;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Point;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.List;
import java.util.Vector;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JColorChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.SwingUtilities;

import htsjdk.samtools.SAMReadGroupRecord;

class ReadGroupsFrame extends JFrame
{
  private static final long serialVersionUID = 1L;
  private LineAttributes[] colours;
  private List<SAMReadGroupRecord> hideReadGroups = new Vector<SAMReadGroupRecord>();
  
  ReadGroupsFrame(List<SAMReadGroupRecord> readGroups, BamView bamView)
  {
    super("Read Groups");
    
    colours = LineAttributes.init(readGroups.size()); 
    update(readGroups, bamView);
    
    JFrame f = (JFrame) SwingUtilities.getAncestorOfClass(JFrame.class, bamView);
    
    int hgt = 250;
    if(hgt < f.getHeight()/2)
      hgt = f.getHeight()/2;
    setPreferredSize(new Dimension(300, hgt));
    
    pack();
    
    final Point p = f.getLocationOnScreen();
    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
    int y = (screen.height-getSize().height)/2;
    if(y < 10) 
      y = 10;
    int x = p.x-getSize().width;
    if(x < 10)
      x = 10;
    setLocation(new Point(x, y));
    
    setDefaultCloseOperation(HIDE_ON_CLOSE);
  }
  
  private void update(final List<SAMReadGroupRecord> readGroups, 
                      final BamView bamView)
  {
    final JPanel mainPanel = new JPanel(new BorderLayout());
    final JPanel panel = new JPanel(new GridBagLayout());
    final JScrollPane jsp = new JScrollPane(panel);
    mainPanel.add(jsp, BorderLayout.CENTER);
    getContentPane().add(mainPanel);
    
    final JButton close = new JButton("CLOSE");
    close.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent arg0)
      {
        setVisible(false);
      }
    });
    mainPanel.add(close, BorderLayout.SOUTH);
    
    final GridBagConstraints c = new GridBagConstraints();
    c.gridy = 0;
    c.anchor = GridBagConstraints.WEST;

    if(readGroups.size() == 0)
    {
      panel.add(new JLabel("No read groups"), c);
      return;
    }
    
    c.gridx = 0;
    JLabel rgLab = new JLabel("Read Groups");
    rgLab.setFont(rgLab.getFont().deriveFont(Font.BOLD));
    panel.add(rgLab, c);
    
    c.gridx = 1;
    JLabel colLab = new JLabel("Colour");
    colLab.setFont(rgLab.getFont().deriveFont(Font.BOLD));
    panel.add(colLab, c);
    
    final JButton toggle = new JButton("Toggle");
    final List<JCheckBox> checkBoxes = new Vector<JCheckBox>();
    c.gridx = 2;
    panel.add(toggle, c);
    toggle.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent arg0)
      {
        for(JCheckBox cb: checkBoxes)
          cb.setSelected(!cb.isSelected());
      }      
    });
    
    c.gridy = 1;
    for(int i=0; i<readGroups.size(); i++)
    {
      final SAMReadGroupRecord grp = readGroups.get(i);

      c.gridx = 0;
      final JLabel lab = new JLabel(grp.getId());
      lab.setToolTipText(grp.getDescription());
      panel.add(lab, c);
      
      c.gridx = 1;
      final LineAttributes lnAttr = colours[i];
      final JButton butt = new JButton("Select", bamView.getImageIcon(lnAttr.getLineColour()));
      butt.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent actionEvent)
        {
          Color newColour = JColorChooser.showDialog(null, "Colour Chooser",
              lnAttr.getLineColour());
          
          lnAttr.setLineColour(newColour); 
          butt.setIcon(bamView.getImageIcon(lnAttr.getLineColour()));
          bamView.repaint();
        }
      });
      panel.add(butt, c);
      
      c.gridx = 2;
      final JCheckBox show = new JCheckBox("Show", true);
      checkBoxes.add(show);
      show.addItemListener(new ItemListener(){
        public void itemStateChanged(ItemEvent arg0)
        {
          if(show.isSelected())
            hideReadGroups.remove(grp);
          else
            hideReadGroups.add(grp);
          bamView.repaintBamView();
        }
      });
      panel.add(show, c);
      
      c.gridy+=1;
    }
  }
  
  protected Color getReadGroupColour(final List<SAMReadGroupRecord> readGroups, final SAMReadGroupRecord rg)
  {
    if(rg == null)
      return Color.BLACK;
    final String rgId = rg.getId();
    for(int i=0; i<readGroups.size(); i++)
    {
      final SAMReadGroupRecord grp = readGroups.get(i);
      if(grp.getId().equals(rgId))
        return colours[i].getLineColour();
    }
    return Color.BLACK;
  }
  
  protected boolean isReadGroupVisible(final SAMReadGroupRecord rg)
  {
    return !hideReadGroups.contains(rg);
  }
}