package uk.ac.sanger.artemis.components.variant;

import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;

import uk.ac.sanger.artemis.components.Utilities;

public class VCFFilter extends JFrame
{
  private static final long serialVersionUID = 1L;
  private static float MIN_QUALITY = 0;
  private static int MIN_DP = 0;
  private static float MIN_MQ = 0;
  private static float MIN_AF1 = 0;
  private static float MAX_CI95 = 10;
  
  /**
   * Filter VCF records by the variant type and/or by different values in 
   * the record, QUAL, DP, MQ and AF1.
   * @param vcfView
   */
  public VCFFilter(final VCFview vcfView)
  {
    super("Variant Filter");
    GridBagConstraints c = new GridBagConstraints();
    
    JPanel panel = (JPanel)getContentPane();
    panel.setLayout(new GridBagLayout());

    // Filter by type
    c.gridx = 0;
    c.gridy = 0;
    c.anchor = GridBagConstraints.WEST;
    JLabel typeLabel = new JLabel("TYPE:");
    typeLabel.setFont(typeLabel.getFont().deriveFont(Font.BOLD));
    panel.add(typeLabel, c);
    
    c.gridy = c.gridy + 1;
    final JCheckBox showSyn = new JCheckBox("Synonymous", vcfView.showSynonymous);
    panel.add(showSyn, c);
    showSyn.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        vcfView.showSynonymous = showSyn.isSelected();
        vcfView.repaint();
      }
    });

    c.gridy = c.gridy + 1;
    final JCheckBox showNonSyn = new JCheckBox("Non-synonymous", vcfView.showNonSynonymous);
    panel.add(showNonSyn, c);
    showNonSyn.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        vcfView.showNonSynonymous = showNonSyn.isSelected();
        vcfView.repaint();
      }
    });
    
    c.gridy = c.gridy + 1;
    final JCheckBox showDeletionsMenu = new JCheckBox("Deletions", vcfView.showDeletions);
    panel.add(showDeletionsMenu, c);
    showDeletionsMenu.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        vcfView.showDeletions = showDeletionsMenu.isSelected();
        vcfView.repaint();
      }
    });

    c.gridy = c.gridy + 1;
    final JCheckBox showInsertionsMenu = new JCheckBox("Insertions", vcfView.showInsertions);
    panel.add(showInsertionsMenu, c);
    showInsertionsMenu.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        vcfView.showInsertions = showInsertionsMenu.isSelected();
        vcfView.repaint();
      }
    });
    
    c.gridy = c.gridy + 1;
    final JCheckBox showMultiAllelesMenu = new JCheckBox("Multiple alleles", vcfView.showMultiAlleles);
    panel.add(showMultiAllelesMenu, c);
    showMultiAllelesMenu.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        vcfView.showMultiAlleles = showMultiAllelesMenu.isSelected();
        vcfView.repaint();
      }
    });
    
    c.gridy = c.gridy + 1;
    final JCheckBox showNonOverlappingsMenu = new JCheckBox("Varaints not overlapping CDS", vcfView.showNonOverlappings);
    panel.add(showNonOverlappingsMenu, c);
    showNonOverlappingsMenu.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        vcfView.showNonOverlappings = showNonOverlappingsMenu.isSelected();
        vcfView.repaint();
      }
    });
    
    if(vcfView.getEntryGroup() == null || vcfView.getEntryGroup().getAllFeaturesCount() == 0)
    {
      showSyn.setEnabled(false);
      showNonSyn.setEnabled(false);
      showNonOverlappingsMenu.setEnabled(false);
    }
    
    // Filter by property
    c.gridy = c.gridy+1;
    panel.add(new JLabel(" "), c);
    c.gridy = c.gridy+1;
    JLabel propLabel = new JLabel("PROPERTY:");
    propLabel.setFont(propLabel.getFont().deriveFont(Font.BOLD));
    panel.add(propLabel, c);
    
    // min quality
    c.gridy = c.gridy+1;
    c.gridx = 0;
    c.anchor = GridBagConstraints.WEST;
    panel.add(new JLabel("Minimum quality score (QUAL):"), c);
    final JTextField minQuality = new JTextField(Float.toString(MIN_QUALITY), 8);
    c.gridx = 1;
    panel.add(minQuality, c);
    
    // min DP
    c.gridy = c.gridy+1;
    c.gridx = 0;
    panel.add(new JLabel("Minimum combined depth across samples (DP):"), c);
    final JTextField minDP = new JTextField(Integer.toString(MIN_DP), 8);
    c.gridx = 1;
    panel.add(minDP, c);
    
    // min MQ
    c.gridy = c.gridy+1;
    c.gridx = 0;
    panel.add(new JLabel("Minimum RMS mapping quality (MQ):"), c);
    final JTextField minMQ = new JTextField(Float.toString(MIN_MQ),8);
    c.gridx = 1;
    panel.add(minMQ, c);
    
    // min AF1
    c.gridy = c.gridy+1;
    c.gridx = 0;
    panel.add(new JLabel("Minimum site frequency of strongest non-reference allele (AF1):"), c);
    final JTextField minAF1 = new JTextField(Float.toString(MIN_AF1),8);
    c.gridx = 1;
    panel.add(minAF1, c);
    
    // max CI95
    c.gridy = c.gridy+1;
    c.gridx = 0;
    panel.add(new JLabel("Maximum 95% confidence interval variation from AF (CI95):"), c);
    final JTextField maxCI95 = new JTextField(Float.toString(MAX_CI95),8);
    c.gridx = 1;
    panel.add(maxCI95, c);

    //
    c.gridy = c.gridy+1;
    c.gridx = 0;
    c.anchor = GridBagConstraints.EAST;
    JButton apply = new JButton("Apply");
    apply.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        try
        {
          MIN_QUALITY = Float.parseFloat(minQuality.getText());
          MIN_DP = Integer.parseInt(minDP.getText());
          MIN_MQ = Float.parseFloat(minMQ.getText());
          MIN_AF1 = Float.parseFloat(minAF1.getText());
          MAX_CI95 = Float.parseFloat(maxCI95.getText());
          vcfView.repaint();
        }
        catch(NumberFormatException ex)
        {
          JOptionPane.showMessageDialog(null, 
              ex.getMessage(), 
              "Format Error", JOptionPane.ERROR_MESSAGE);
        }
      }
    });
    panel.add(apply, c);

    c.gridx = 1;
    c.anchor = GridBagConstraints.WEST;
    JButton close = new JButton("Close");
    close.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        dispose();
      }
    });
    panel.add(close, c);
    
    pack();
    Utilities.centreFrame(this);
    setVisible(true);
  }
  
  /**
   * Test for a given VCF record to see if it passes the filters.
   * @param record
   * @return
   */
  protected static boolean passFilter(VCFRecord record)
  {
    try
    {
      if(record.getQuality() < VCFFilter.MIN_QUALITY)
        return false;
      
      try
      {
        if(VCFFilter.MIN_DP > 0 && Integer.parseInt(record.getInfoValue("DP")) < VCFFilter.MIN_DP)
          return false;
      }
      catch(NullPointerException npe){}
    
      try
      {
        if(VCFFilter.MIN_MQ > 0 && Float.parseFloat(record.getInfoValue("MQ")) < VCFFilter.MIN_MQ)
          return false;
      }
      catch(NullPointerException npe){}

      try
      {
        if(VCFFilter.MIN_AF1 > 0 && Float.parseFloat(record.getInfoValue("AF1")) < VCFFilter.MIN_AF1)
          return false;
      }
      catch(NullPointerException npe){}
      
      
      try
      {
        String vals[] = record.getInfoValue("CI95").split(",");
        for(int i=0; i<vals.length; i++)
        {
          if(VCFFilter.MAX_CI95 < 10 && Float.parseFloat(vals[i]) > VCFFilter.MAX_CI95)
            return false;
        }
      }
      catch(NullPointerException npe){}
    }
    catch(NumberFormatException e)
    {
      System.err.println(e.getMessage()); 
    }

    return true;
  }
}