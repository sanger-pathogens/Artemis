package uk.ac.sanger.artemis.components.alignment;


import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;

import uk.ac.sanger.artemis.components.StickyFileChooser;

/**
 * File selection panel to allow input of DNA sequences
 */
class FileSelectionPanel extends JPanel
{
  private static final long serialVersionUID = 1L;
  private GridBagConstraints c = new GridBagConstraints();

  private JTextField bamField = new JTextField(30);
  private JTextField referenceField = new JTextField(30);
  
  /**
   * Constructor to display any given input files and options provided and
   */
  public FileSelectionPanel()
  {
    super(new GridBagLayout());

    int row = 0;
    
    bamField.setPreferredSize(
        new Dimension(200,bamField.getPreferredSize().height));
    c.gridy = row;
    c.gridx = 0;
    c.anchor = GridBagConstraints.WEST;
    add(new JLabel(" BAM file: "), c);
    c.gridy = ++row;
    add(bamField, c);
    c.gridx = 1;
    JButton selectBam = new JButton("Select...");
    addActionListener(selectBam, bamField);
    add(selectBam, c);
    
    c.gridy = ++row;
    c.gridx = 0;
    add(new JLabel(" Reference sequence file (optional): "), c);
    c.gridy = ++row;
    add(referenceField, c);
    JButton selectReference = new JButton("Select...");
    addActionListener(selectReference, referenceField);
    c.gridx = 1;
    add(selectReference, c);
  }
  
  /**
   * Add action listener to a button.
   * @param fileSelectionButton
   * @param tfield
   */
  private void addActionListener(final JButton fileSelectionButton, 
                                 final JTextField tfield)
  {
    fileSelectionButton.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        StickyFileChooser fileChooser = new StickyFileChooser();
        int status = fileChooser.showOpenDialog(null);
        if(status == StickyFileChooser.CANCEL_OPTION)
          return;
        tfield.setText(fileChooser.getSelectedFile().getAbsolutePath());
      }
    });
  }
  
  /**
   * Open up the FileSelectionPanel in a JOptionPane
   * @param f
   * @param displayButtonListener
   * @return
   */
  protected void showPanel()
  {
    JOptionPane.showMessageDialog(null, this, 
        "BamView :: Select Files", JOptionPane.PLAIN_MESSAGE);
  }
  
  /**
   * Get the bam file
   * @return
   */
  protected String getBamFile()
  {
    return bamField.getText();
  }
  
  /**
   * Get the reference
   * @return
   */
  protected String getReferenceFile()
  {
    return referenceField.getText();
  }
}