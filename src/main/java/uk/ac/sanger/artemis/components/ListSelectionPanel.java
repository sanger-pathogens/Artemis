/* ListSelectionPanel.java
 *
 * created: Fri Oct  9 1998
 *
 * This file is part of Artemis
 *
 * Copyright(C) 1998,1999,2000,2001  Genome Research Limited
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or(at your option) any later version.
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/ListSelectionPanel.java,v 1.6 2007-09-12 10:14:01 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.Box;
import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.io.GFFDocumentEntry;
import uk.ac.sanger.artemis.io.Key;

public class ListSelectionPanel extends JPanel
{

  /** */
  private static final long serialVersionUID = 1L;
  private DefaultListModel listModel;
  private JCheckBox save;
  
  public ListSelectionPanel(final EntryGroup entry_group,
                            final Object names[],
                            final String[] description,
                            final boolean saveOption)
  {
    this(entry_group, names, description, saveOption, null, null);
  }
  
  /**
   * Panel used to select from and order a list. 
   * @param entry_group
   * @param names
   * @param description
   */
  public ListSelectionPanel(final EntryGroup entry_group,
                            final Object names[],
                            final String[] description,
                            final boolean saveOption,
                            final Key key,
                            final String default_qualifier)
  {
    listModel = new DefaultListModel();

    for(int i = 0; i < names.length; i++)
      listModel.addElement(names[i]);

    final JLabel label = new JLabel("Qualifier order :");
    final JList name_list = new JList(listModel);
    final JScrollPane jsp = new JScrollPane(name_list);

    JButton remove_butt = new JButton("REMOVE");
    remove_butt.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        while(!name_list.isSelectionEmpty())
          listModel.remove(name_list.getSelectedIndex());
      }
    });

    Box bdown = Box.createVerticalBox();
    bdown = Box.createVerticalBox();
    JButton upButt = new JButton("UP");
    upButt.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        Object obj = name_list.getSelectedValue();
        int index = name_list.getSelectedIndex();

        if(index <= 0)
          return;
        listModel.removeElementAt(index);
        listModel.insertElementAt(obj, index - 1);
        name_list.setSelectedIndex(index - 1);
      }
    });
    bdown.add(upButt);

    JButton downButt = new JButton("DOWN");
    downButt.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        Object obj = name_list.getSelectedValue();
        int index = name_list.getSelectedIndex();

        if(index >= listModel.getSize() - 1)
          return;
        listModel.removeElementAt(index);
        listModel.insertElementAt(obj, index + 1);
        name_list.setSelectedIndex(index + 1);
      }
    });
    bdown.add(downButt);
    add(bdown, BorderLayout.WEST);

    bdown = Box.createVerticalBox();
  
    bdown.add(Box.createVerticalStrut(15));
    for(int i=0; i<description.length; i++)
    	bdown.add(new JLabel(description[i]));
    bdown.add(label);
    bdown.add(jsp);
    
    add(bdown, BorderLayout.CENTER);

    boolean isGFF = false;
    if(entry_group.getDefaultEntry().getEMBLEntry() instanceof GFFDocumentEntry)
      isGFF = true;
    
    final QualifierChoice qualifier_choice = new QualifierChoice(
        entry_group.elementAt(0).getEntryInformation(), key, default_qualifier, isGFF);
    
    JButton add_butt = new JButton("ADD");
    add_butt.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        listModel.addElement(qualifier_choice.getSelectedItem().toString());
      }
    });
    
    if(saveOption)
      save = new JCheckBox("Save between sessions", false);

    bdown = Box.createVerticalBox();
    bdown.add(Box.createVerticalGlue());
    bdown.add(qualifier_choice);
    bdown.add(add_butt);
    bdown.add(remove_butt);
    if(saveOption)
      bdown.add(save);
    
    add(bdown, BorderLayout.EAST);
  }
  
  public Object[] getResultArray()
  {
    return listModel.toArray();
  }
  
  public String getResultString()
  {
    Object listNames[] = listModel.toArray();
    String listNamesString = "";
    for(int i=0; i<listNames.length; i++)
    	listNamesString = listNamesString + listNames[i] + " ";
    return listNamesString;
  }
    
  public boolean isSaveOption()
  {
    return save.isSelected();
  }
}