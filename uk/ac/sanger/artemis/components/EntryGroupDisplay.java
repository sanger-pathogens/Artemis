/* EntryGroupDisplay.java
 *
 * created: Mon Dec  7 1998
 *
 * This file is part of Artemis
 * 
 * Copyright(C) 1998,1999,2000,2002  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/EntryGroupDisplay.java,v 1.5 2009-05-29 10:16:15 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.io.IndexFastaStream;

import java.awt.*;
import java.awt.event.*;
import java.util.Iterator;
import java.util.Vector;

import javax.swing.*;

import net.sf.picard.reference.FastaSequenceIndex;

/**
 *  This component allows the user to change the "active" setting of the
 *  objects in an EntryGroup.
 *
 *  @author Kim Rutherford
 *  @version $Id: EntryGroupDisplay.java,v 1.5 2009-05-29 10:16:15 tjc Exp $
 **/

public class EntryGroupDisplay extends JPanel
    implements EntryGroupChangeListener, EntryChangeListener 
{
  final protected static Color background_colour = new Color(200, 200, 200);

  /**
   *  This is a reference to the EntryEdit component that created this
   *  EntryGroupDisplay.  We need this reference so that we can watch the
   *  entry group.
   **/
  private EntryEdit owning_component;

  /**
   *  A vector containing the Entry objects that this EntryEdit object knows
   *  about.  This reference is obtained from owning_component.
   **/
  private EntryGroup entry_group;

  /**
   *  A vector containing one JCheckBox or Label for each Entry in the
   *  EntryGroup object.
   **/
  private Vector<JCheckBox> entry_components = new Vector<JCheckBox>();

  /**
   *  A label containing the message "Entry:".
   **/
  private JLabel label;
  
  private JComboBox indexFastaCombo;

  /**
   *  Create a new EntryGroupDisplay object.
   *  @param owning_component The EntryEdit object that this EntryGroupDisplay
   *    component is in.
   **/
  public EntryGroupDisplay(final EntryEdit owning_component) 
  {
    this.owning_component = owning_component;
    this.entry_group = owning_component.getEntryGroup();

    entry_group.addEntryGroupChangeListener(this);
    entry_group.addEntryChangeListener(this);

    final FlowLayout flow_layout = new FlowLayout(FlowLayout.LEFT,2,1); 

    label = new JLabel("Entry: ");

    setLayout(flow_layout);
    refreshButtons();
    
    setBackground(background_colour);
  }

 
  protected void printComponent(Graphics g)
  {
    super.paintComponent(g);
    super.printChildren(g);
  }


  /**
   *  Implementation of the EntryGroupChangeListener interface.  We listen to
   *  EntryGroupChange events so that we can update the display if entries
   *  are added or deleted.
   **/
  public void entryGroupChanged(final EntryGroupChangeEvent event) 
  {
    switch(event.getType()) 
    {
      case EntryGroupChangeEvent.ENTRY_ADDED:
      case EntryGroupChangeEvent.ENTRY_DELETED:
        refreshButtons();
        break;
      case EntryGroupChangeEvent.ENTRY_INACTIVE:
      case EntryGroupChangeEvent.ENTRY_ACTIVE:
        updateActive();
        break;
      case EntryGroupChangeEvent.NEW_DEFAULT_ENTRY:
        highlightDefaultEntry(event);
        break;
    }
  }

  /**
   *  Implementation of the EntryChangeListener interface.
   **/
  public void entryChanged(final EntryChangeEvent event) 
  {
    if(event.getType() == EntryChangeEvent.NAME_CHANGED) 
      refreshButtons();
  }

  /**
   *  Remove and then recreate the Buttons to the reflect the current contents
   *  of the EntryGroup.
   **/
  private void refreshButtons() 
  {
    removeAll();
    add(label);

    entry_components = new Vector<JCheckBox>();

    if(entry_group == null) 
      return;
    else 
    {
      for(int i = 0; i < entry_group.size(); ++i) 
        add(entry_group.elementAt(i));
    }

    validate();
  }

  /**
   *  Update the buttons to reflect the current state of the EntryGroup.
   **/
  private void updateActive()
  {
    for(int i = 0 ; i < entry_group.size() ; ++i) 
    {
      final JCheckBox menu_item = entry_components.elementAt(i);
      menu_item.setSelected(entry_group.isActive(entry_group.elementAt(i)));
    }
  }

  /**
   *  Add a JLabel or JCheckBox for the given Entry to this component.
   **/
  private void add(final Entry entry) 
  {
    final JCheckBox new_component;
    String entry_name = entry.getName();

    if(entry_name == null) 
      entry_name = "no name";
    
    new_component = new JCheckBox(entry_name, entry_group.isActive(entry));
    new_component.setOpaque(true);
    setEntryHighlight(entry, new_component);

    new_component.addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent event) 
      {
        final int button_index =
          entry_components.indexOf(event.getSource());

        if(event.getStateChange() == ItemEvent.SELECTED) 
          entry_group.setIsActive(button_index, true);
        else 
          entry_group.setIsActive(button_index, false);
      }
    });

    new_component.addMouseListener(new MouseAdapter() 
    {
      /**
       *  Listen for mouse press events so that we can change the default
       *  Entry when the popup trigger is pressed.
       **/
      public void mousePressed(MouseEvent event) 
      {
        if(event.isPopupTrigger()) 
          entry_group.setDefaultEntry(entry);
      }
    });

    entry_components.addElement(new_component);
    add(new_component);
    
    if(entry.getEMBLEntry().getSequence() instanceof IndexFastaStream)
    {
      if(indexFastaCombo == null)
      {
        FastaSequenceIndex indexFasta = 
          ((IndexFastaStream)entry.getEMBLEntry().getSequence()).getFastaIndex();
        Iterator it = indexFasta.iterator();
        Vector<String> contigs = new Vector<String>();
        while(it.hasNext())
        {
          String contig = it.next().toString().split(";")[0];
          if(contig.startsWith("contig "))
            contig = contig.substring(6);
          contigs.add(  contig );
        }

        indexFastaCombo = new JComboBox(contigs);
        indexFastaCombo.addActionListener(new ActionListener()
        {
          public void actionPerformed(ActionEvent e)
          {
            IndexFastaStream is = (IndexFastaStream)entry.getEMBLEntry().getSequence();
            is.setContigByIndex(indexFastaCombo.getSelectedIndex());
          
            owning_component.resetScrolls();
            owning_component.getFeatureDisplay().getBases().clearCodonCache();
            owning_component.repaint();
          }
        
        });
      }
      add(indexFastaCombo);
    }
  }

  /**
   *  Given a EntryGroupChangeEvent, this method will highlight the new
   *  default entry
   **/
  private void highlightDefaultEntry(final EntryGroupChangeEvent event) 
  {
    final EntryGroup entry_group = owning_component.getEntryGroup();

    for(int i = 0 ; i < entry_group.size() ; ++i) 
    {
      final JCheckBox check_box = entry_components.elementAt(i);
      setEntryHighlight(entry_group.elementAt(i), check_box);
    }
  }

  /**
   *  Highlight the given JCheckBox/Entry appropriately.  The default Entry
   *  will look different to the others
   **/
  private void setEntryHighlight(final Entry entry,
                                  final JCheckBox component)
  {
    //final String label = component.getText();

    if(entry_group.getDefaultEntry() == entry) 
      component.setBackground(Color.yellow);
    else 
      component.setBackground(background_colour);
  }
  
}
