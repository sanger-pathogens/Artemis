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

import javax.swing.table.AbstractTableModel;
import javax.swing.JTree;
import javax.swing.tree.*;
import javax.swing.event.TreeExpansionEvent;
import javax.swing.event.TreeExpansionListener;

/**
 * This is a wrapper class takes a TreeTableModel and implements 
 * the table model interface. The implementation is trivial, with 
 * all of the event dispatching support provided by the superclass: 
 * the AbstractTableModel. 
 *
 * @version %I% %G%
 *
 * @author Philip Milne
 * @author Scott Violet
 */


public class TreeTableModelAdapter extends AbstractTableModel
{
    JTree tree;
    TreeModel treeTableModel;

    public TreeTableModelAdapter(TreeModel treeTableModel, JTree tree) {
        this.tree = tree;
        this.treeTableModel = treeTableModel;

	tree.addTreeExpansionListener(new TreeExpansionListener() {
	    // Don't use fireTableRowsInserted() here; 
	    // the selection model would get  updated twice. 
	    public void treeExpanded(TreeExpansionEvent event) {  
	      fireTableDataChanged(); 
	    }
            public void treeCollapsed(TreeExpansionEvent event) {  
	      fireTableDataChanged(); 
	    }
	});
    }

  // Wrappers, implementing TableModel interface. 

    public int getColumnCount() {
	return ((FileSystemModel)treeTableModel).getColumnCount();
    }

    public String getColumnName(int column) {
	return ((FileSystemModel)treeTableModel).getColumnName(column);
    }

    public Class getColumnClass(int column) {
	return ((FileSystemModel)treeTableModel).getColumnClass(column);
    }

    public int getRowCount() {
	return tree.getRowCount();
    }

    protected Object nodeForRow(int row) {
	TreePath treePath = tree.getPathForRow(row);
	return treePath.getLastPathComponent();         
    }

    public Object getValueAt(int row, int column) {
	return ((FileSystemModel)treeTableModel).getValueAt(nodeForRow(row), column);
    }

    public boolean isCellEditable(int row, int column) {
        return true;
//       return ((FileSystemModel)treeTableModel).isCellEditable(nodeForRow(row), column); 
    }

    public void setValueAt(Object value, int row, int column) {
//	((FileSystemModel)treeTableModel).setValueAt(value, nodeForRow(row), column);
    }
}


