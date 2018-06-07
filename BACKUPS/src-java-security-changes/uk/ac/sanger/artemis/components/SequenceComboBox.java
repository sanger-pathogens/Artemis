/* SequenceComboBox
 *
 * created: 2012
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2012  Genome Research Limited
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
 */
package uk.ac.sanger.artemis.components;

import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Vector;

import javax.swing.DefaultComboBoxModel;
import javax.swing.JComboBox;
import javax.swing.text.JTextComponent;

import uk.ac.sanger.artemis.components.genebuilder.AutoCompleteComboDocument;

public abstract class SequenceComboBox extends JComboBox implements IndexReferenceListener
{
  private static final long serialVersionUID = 1L;
  private ArrayList<IndexReferenceListener> indexReferenceListeners = new ArrayList<IndexReferenceListener>();

  public SequenceComboBox(Vector<String> sequenceNames)
  {
    super(sequenceNames);
    JTextComponent editor = (JTextComponent) getEditor().getEditorComponent();
    editor.setDocument(new AutoCompleteComboDocument(this));
    setEditable(true);
    setMaximumRowCount(20);
    
    addItemListener(new ItemListener()
    {
      public void itemStateChanged(ItemEvent e)
      {
        fireAction(new IndexReferenceEvent(SequenceComboBox.this));
      }
    });
  }
  
  public SequenceComboBox(String sequenceNames[])
  {
    this(new Vector<String>(Arrays.asList(sequenceNames)));
  }
  
  public void indexReferenceChanged(IndexReferenceEvent event)
  {
    final String contig = 
        ((String) ((SequenceComboBox)event.getSource()).getSelectedItem()).trim();
    if(!contains(contig))
      return;
    if(!getSelectedItem().equals(contig))
      setSelectedItem(contig);

    update(event);
  }
  
  private void fireAction(IndexReferenceEvent event)
  {
    update(event);
    for(IndexReferenceListener l: indexReferenceListeners)
      l.indexReferenceChanged(event);
  }
  
  private boolean contains(Object o) 
  {
    final DefaultComboBoxModel model = (DefaultComboBoxModel) getModel();
    for(int i = 0; i < model.getSize(); i++) 
    {
      Object obj = model.getElementAt(i);
      if(obj.equals(o))
        return true;
    }
    return false;
  }
  
  abstract public void update(IndexReferenceEvent event);
  
  /**
   *  Adds the specified index reference listener to receive
   *  change events from this object.
   *  @param l the event change listener.
   **/
  public void addIndexReferenceListener(IndexReferenceListener l)
  {
    indexReferenceListeners.add(l);
  }
}