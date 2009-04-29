/*
 * Copyright (C) 2009  Genome Research Limited
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 */

package uk.ac.sanger.artemis.circular.digest;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.util.List;
import java.util.Vector;

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import uk.ac.sanger.artemis.components.StickyFileChooser;
import uk.ac.sanger.artemis.components.Utilities;

/**
 * File selection panel to allow input of DNA sequences
 */
class FileSelectionPanel extends JPanel
{
  private static final long serialVersionUID = 1L;
  private GridBagConstraints c = new GridBagConstraints();
  private int row = 0;
  private List<SelectionRow> selectionRows = new Vector<SelectionRow>();

  private List<File> restrictOutputs;
  private JTextField enzymeField;
  private JTextField embossRootField;
  private JCheckBox methylationCheckBox;
  private JFrame f;
  
  /**
   * Constructor to display any given input files and options provided and
   * provide a graphical input window for adding sequences.
   * that have been given on the command line.
   * @param enzymes  enzyme string, takes the value from the -enz parameter
   * @param sequenceFiles sequence file, takes the value(s) from the -seq parameter
   * @param restrictOutputs pre-computed EMBOSS restrict output, takes the value(s)
   *                        from the -restrict parameter
   * @param methylation if true takes into account methylation blocked sites, takes
   *                    the value of -methylation parameter
   */
  public FileSelectionPanel(
      String enzymes, 
      List<File> sequenceFiles,
      List<File> restrictOutputs,
      boolean methylation)
  {
    super(new GridBagLayout());

    this.restrictOutputs = restrictOutputs;
    if(enzymes == null || enzymes.equals(""))
      enzymes = "xbai";
    enzymeField = new JTextField(enzymes, 30);
    enzymeField.setPreferredSize(
        new Dimension(200,enzymeField.getPreferredSize().height));
    c.gridy = row;
    c.gridx = 0;
    c.gridwidth = 2;
    c.anchor = GridBagConstraints.WEST;
    add(new JLabel(" Comma separated list of digest enzymes: "), c);
    row++;
    c.gridx = 1;
    c.gridy = row;
    c.gridwidth = 1;
    add(enzymeField, c);
    row++;
    
    methylationCheckBox = new JCheckBox(
        "RE sites will not match methylated bases", methylation);
    c.gridy = row;
    add(methylationCheckBox, c);
    row++;
    
    c.gridy = row;
    add(Box.createVerticalStrut(10), c);
    row++;
    
    if(restrictOutputs == null || restrictOutputs.size() == 0)
    {
      embossRootField = new JTextField(System.getProperty("EMBOSS_ROOT"), 30);
      c.gridy = row;
      c.gridx = 0;
      add(new JLabel(" EMBOSS location: "), c);
      c.gridx = 1;
      add(embossRootField, c);
      row++;
      
      c.gridy = row;
      add(Box.createVerticalStrut(10), c);
      row++;
    }
    
    
    c.gridy = row;
    c.gridx = 0;
    c.anchor = GridBagConstraints.WEST;
    add(new JLabel(" Sequence file(s): "), c);
    row++;
    
    if(sequenceFiles == null || sequenceFiles.size() < 1)
      addSelectionRow(null);
    else
    {
      for(int i=0; i<sequenceFiles.size(); i++)
        addSelectionRow(sequenceFiles.get(i).getAbsolutePath());
    }
  
    JButton addMoreFiles = new JButton("Add More");
    addMoreFiles.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        addSelectionRow(null);
      }
    });
    c.gridx = 1;
    c.gridy = row+100;
    add(addMoreFiles, c);
  }

  /**
   * Add a file selection row to the input window.
   * @param pathToSequence
   */
  private void addSelectionRow(final String pathToSequence)
  {
    final SelectionRow sRow = new SelectionRow();
    selectionRows.add(sRow);
    
    JButton fileSelectionButton = new JButton("Select File...");
    fileSelectionButton.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        StickyFileChooser fileChooser = new StickyFileChooser();
        int status = fileChooser.showOpenDialog(null);
        if(status == StickyFileChooser.CANCEL_OPTION)
          return;
        sRow.fileNameField.setText(fileChooser.getSelectedFile().getAbsolutePath());
        FileSelectionPanel.this.repaint();
      }
    });
    
    if(pathToSequence != null)
      sRow.fileNameField.setText(pathToSequence);
    
    c.gridy = row;
    c.gridx = 0;
    add(fileSelectionButton, c);
    c.gridx = 1;
    add(sRow.fileNameField, c);
    row++;
    if(f == null)
      revalidate();
    else
      f.pack();
  }
  
  /**
   * Open up the FileSelectionPanel in a JFrame
   * @param f
   * @param displayButtonListener
   * @return
   */
  protected JButton showJFrame(final JFrame f,
                               ActionListener displayButtonListener)
  {
    this.f = f;
    f.getContentPane().add(this, BorderLayout.CENTER);
    JButton displayButton = new JButton("Display");
    displayButton.addActionListener(displayButtonListener);
    f.getContentPane().add(displayButton, BorderLayout.SOUTH);
    f.pack();
    Utilities.centreFrame(f);
    f.setVisible(true);
    
    return displayButton;
  }
  
  /**
   * Get the comma separated enzyme list
   * @return
   */
  protected String getEnzymes()
  {
    return enzymeField.getText();
  }

  protected List<File> getSequenceFiles()
  {
    List<File> sequenceFiles = new Vector<File>();
    for(int i=0; i<selectionRows.size(); i++)
    {
      SelectionRow r = selectionRows.get(i);
      if(!r.fileNameField.getText().equals(""))
        sequenceFiles.add(new File(r.fileNameField.getText()));
    }
    return sequenceFiles;
  }

  protected List<File> getRestrictOutputs()
  {
    return restrictOutputs;
  }
  
  protected JTextField getEmbossRootField()
  {
    return embossRootField;
  }

  protected boolean isMethylation()
  {
    return methylationCheckBox.isSelected();
  }
  
  class SelectionRow
  {
    final JTextField fileNameField = new JTextField(30);
    
    SelectionRow()
    {
      fileNameField.setPreferredSize(
          new Dimension(200,fileNameField.getPreferredSize().height));
    }
  }
}