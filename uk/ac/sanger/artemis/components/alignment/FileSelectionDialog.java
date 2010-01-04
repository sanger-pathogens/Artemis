package uk.ac.sanger.artemis.components.alignment;


import java.awt.Dimension;
import java.awt.Frame;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Point;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.List;
import java.util.Vector;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import uk.ac.sanger.artemis.components.StickyFileChooser;

/**
 * File selection panel to allow input of DNA sequences
 */
class FileSelectionDialog extends JDialog
{
  private static final long serialVersionUID = 1L;
  
  private JPanel dialog = new JPanel(new GridBagLayout());
  private GridBagConstraints c = new GridBagConstraints();
  private int row = 0;
  private List<JTextField> bamFields = new Vector<JTextField>(30);
  private JTextField referenceField = new JTextField(30);
  
  /**
   * Constructor to display any given input files and options provided and
   */
  public FileSelectionDialog(Frame f, boolean showReferenceOption)
  {
    super(f, "BamView :: Select Files", true);

    addBamField();
    
    JButton addMoreFiles = new JButton("Add More");
    addMoreFiles.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        addBamField();
      }
    });
    
    row += 100;
    c.gridx = 1;
    c.gridy = row;
    dialog.add(addMoreFiles, c);

    if(showReferenceOption)
    {
      c.gridy = ++row;
      c.gridx = 0;
      dialog.add(new JLabel(" Reference sequence file (optional): "), c);
      c.gridy = ++row;
      dialog.add(referenceField, c);
      JButton selectReference = new JButton("Select...");
      addActionListener(selectReference, referenceField);
      c.gridx = 1;
      dialog.add(selectReference, c);
    }
    
    JButton okButton = new JButton("OK");
    okButton.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        FileSelectionDialog.this.dispose();
      }
    });
    c.gridy = ++row;
    dialog.add(okButton, c);
    getContentPane ().add (dialog, "South");
    
    row = 1;
    
    pack();
    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
    setLocation (new Point ((screen.width - getSize ().width) / 2,
                            (screen.height - getSize ().height) / 2));
    setVisible (true);
  }
  
  /**
   * Add a text field to the dialog for adding in a path
   * to a BAM file.
   */
  private void addBamField()
  {
    JTextField bamField = new JTextField(30);
    bamFields.add(bamField);
    
    bamField.setPreferredSize(
        new Dimension(200,bamField.getPreferredSize().height));
    c.gridy = row;
    c.gridx = 0;
    c.anchor = GridBagConstraints.WEST;
    dialog.add(new JLabel(" BAM file: "), c);
    c.gridy = ++row;
    dialog.add(bamField, c);
    c.gridx = 1;
    JButton selectBam = new JButton("Select...");
    addActionListener(selectBam, bamField);
    dialog.add(selectBam, c);
    
    pack();
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
   * Get the BAM files as a <code>List</code> of <code>String</code>'s.
   * @return
   */
  protected List<String> getBamFiles()
  {
    List<String> bamFiles = new Vector<String>();
    for(int i=0; i<bamFields.size(); i++)
    {
      String file = bamFields.get(i).getText();
      if(file != null && !file.equals(""))
        bamFiles.add(file);
    }
    return bamFiles;
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