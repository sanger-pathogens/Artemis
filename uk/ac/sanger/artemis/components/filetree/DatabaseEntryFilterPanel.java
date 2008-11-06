/********************************************************************
*
*  This library is free software; you can redistribute it and/or
*  modify it under the terms of the GNU Library General Public
*  License as published by the Free Software Foundation; either
*  version 2 of the License, or (at your option) any later version.
*
*  This library is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
*  Library General Public License for more details.
*
*  You should have received a copy of the GNU Library General Public
*  License along with this library; if not, write to the
*  Free Software Foundation, Inc., 59 Temple Place - Suite 330,
*  Boston, MA  02111-1307, USA.
*
*  Copyright (C) Genome Research Limited
*
********************************************************************/
package uk.ac.sanger.artemis.components.filetree;

import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import uk.ac.sanger.artemis.components.KeyChoice;
import uk.ac.sanger.artemis.io.GFFEntryInformation;
import uk.ac.sanger.artemis.util.DatabaseDocument;

/**
 * This displays the feature key information used to split the data read 
 * from the database into separate Artemis entries. The name of each entry
 * is displayed along with an editable list of feature keys that belong
 * to that entry.
 */
public class DatabaseEntryFilterPanel extends JPanel
{
  private static final long serialVersionUID = 1L;
  private JTextField[] nameField;
  private JComboBox[] keyList;
  
  public DatabaseEntryFilterPanel()
  {
    super(new GridBagLayout());

    final GridBagConstraints c = new GridBagConstraints();

    // current feature types in each entry
    String types[][][] = DatabaseDocument.getTYPES();

    c.gridx = 0;
    c.gridy = 0;
    add(new JLabel("Name"), c);
    c.gridx = 1;
    add(new JLabel("Feature Key List"), c);
    nameField = new JTextField[types.length];
    keyList = new JComboBox[types.length];
    for(int i = 0; i < types.length; i++)
    {
      final KeyChoice keyChoice = new KeyChoice(new GFFEntryInformation());
      c.gridx = 0;
      c.gridy = i + 1;
      
      nameField[i] = new JTextField(types[i][0][0], 20);
      add(nameField[i], c);
      
      keyList[i] = new JComboBox(types[i][1]);
      keyList[i].setPreferredSize(new Dimension(
          keyChoice.getPreferredSize().width, keyList[i].getPreferredSize().height));
      final int MAX_VISIBLE_ROWS = 30;
      keyList[i].setMaximumRowCount (MAX_VISIBLE_ROWS);
      keyList[i].setEditable(false);
      final JComboBox thiskeyList = keyList[i];
      
      c.gridx = 1;
      add(keyList[i], c);
      c.gridx = 2;
      JButton deleteKey = new JButton("Delete");
      add(deleteKey, c);
      deleteKey.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent arg0)
        {
          thiskeyList.removeItem(thiskeyList.getSelectedItem());
          thiskeyList.revalidate();
        }
      });
      
      c.gridx = 3;
      JButton addKey = new JButton("Add :");
      add(addKey, c);
      c.gridx = 4;
      
      add(keyChoice, c);
      
     
      addKey.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent arg0)
        {
          for(int i=0; i<thiskeyList.getItemCount(); i++)
            if(thiskeyList.getItemAt(i).equals(keyChoice.getSelectedItem()))
              return;
          thiskeyList.addItem(keyChoice.getSelectedItem());
          thiskeyList.setSelectedItem(keyChoice.getSelectedItem());
        }
      });
    }
  }

  /**
   * Method to set the Entry names (from the nameField text fields) and
   * set the feature types in each entry (based on keyList values).
   */
  protected void setTypesForEntries()
  {
    String newTypes[][][] = new String[nameField.length][2][];
    for(int i=0; i<nameField.length; i++)
    {
      newTypes[i][0] = new String[1];
      newTypes[i][0][0] = nameField[i].getText().trim() ;
      
      String[] keys = new String[keyList[i].getItemCount()];
      //System.out.print(newTypes[i][0][0]);
      
      for(int j=0; j<keyList[i].getItemCount(); j++)
      {
        keys[j] = (String) keyList[i].getItemAt(j);
        //System.out.print(" "+keys[j]);
      }
      //System.out.println();
      newTypes[i][1] = keys;
    }
    DatabaseDocument.setTYPES(newTypes);
  }

}