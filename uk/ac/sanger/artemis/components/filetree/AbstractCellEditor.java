/*
 *
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2006  Genome Research Limited
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

package uk.ac.sanger.artemis.components.filetree;

import java.awt.Component;
import java.awt.event.*;
import java.awt.AWTEvent;
import javax.swing.*;
import javax.swing.event.*;
import java.util.EventObject;
import java.io.Serializable;

/**
 * 
 * A base class for CellEditors, providing default implementations for all 
 * methods in the CellEditor interface and support for managing a series 
 * of listeners. 
 *
 */
public class AbstractCellEditor implements CellEditor 
{

  protected EventListenerList listenerList = new EventListenerList();

  public Object getCellEditorValue() { return null; }
  public boolean isCellEditable(EventObject e) { return true; }
  public boolean shouldSelectCell(EventObject anEvent) { return false; }
  public boolean stopCellEditing() { return true; }
  public void cancelCellEditing() {}

  public void addCellEditorListener(CellEditorListener l) 
  {
    listenerList.add(CellEditorListener.class, l);
  }

  public void removeCellEditorListener(CellEditorListener l) 
  {
    listenerList.remove(CellEditorListener.class, l);
  }

  /**
   * Notify all listeners that have registered interest for
   * notification on this event type.  
   * @see EventListenerList
   */
  protected void fireEditingStopped() 
  {
    // Guaranteed to return a non-null array
    Object[] listeners = listenerList.getListenerList();
    // Process the listeners last to first, notifying
    // those that are interested in this event
    for(int i = listeners.length-2; i>=0; i-=2) 
    {
      if(listeners[i]==CellEditorListener.class) {
	((CellEditorListener)listeners[i+1]).editingStopped(new ChangeEvent(this));
      }	       
    }
  }

  /**
   * Notify all listeners that have registered interest for
   * notification on this event type.  
   * @see EventListenerList
   */
  protected void fireEditingCanceled() 
  {
    // Guaranteed to return a non-null array
    Object[] listeners = listenerList.getListenerList();
    // Process the listeners last to first, notifying
    // those that are interested in this event
    for (int i = listeners.length-2; i>=0; i-=2) 
    {
      if(listeners[i]==CellEditorListener.class) {
	((CellEditorListener)listeners[i+1]).editingCanceled(new ChangeEvent(this));
      }	       
    }
  }
}

