/* DatabaseJTree.java
 *
 * created: October 2006
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

package uk.ac.sanger.artemis.components.database;


import java.awt.Cursor;
import java.awt.Point;
import java.awt.datatransfer.Transferable;
import java.awt.dnd.DnDConstants;
import java.awt.dnd.DragGestureEvent;
import java.awt.dnd.DragGestureListener;
import java.awt.dnd.DragSource;
import java.awt.dnd.DragSourceDragEvent;
import java.awt.dnd.DragSourceDropEvent;
import java.awt.dnd.DragSourceEvent;
import java.awt.dnd.DragSourceListener;
import java.awt.event.InputEvent;
import java.awt.event.MouseEvent;

import javax.swing.JTree;
import javax.swing.event.TreeExpansionEvent;
import javax.swing.event.TreeExpansionListener;
import javax.swing.tree.DefaultTreeModel;
import javax.swing.tree.TreePath;
import javax.swing.tree.TreeSelectionModel;


/**
*
* File node for local file tree manager
*
*/
public class DatabaseJTree extends JTree
               implements DragGestureListener,
                          DragSourceListener
{
   /** */
  private static final long serialVersionUID = 1L;

    /**
    *
    * @param file file node file
    *
    */
    public DatabaseJTree(DatabaseTreeNode rootNode)
    { 
      super(rootNode);
      getSelectionModel().setSelectionMode(
          TreeSelectionModel.SINGLE_TREE_SELECTION);
      
      
      DragSource dragSource = DragSource.getDefaultDragSource();

      dragSource.createDefaultDragGestureRecognizer(
         this,                             // component where drag originates
         DnDConstants.ACTION_COPY_OR_MOVE, // actions
         this);                            // drag gesture recognizer
      
      addTreeExpansionListener(new TreeExpansionListener() 
      {
        public void treeExpanded(TreeExpansionEvent e)
        {
          TreePath path = e.getPath();
          if(path != null) 
          {
            setCursor(new Cursor(Cursor.WAIT_CURSOR));
            DatabaseTreeNode node = (DatabaseTreeNode)path.getLastPathComponent();

            if(!node.isExplored()) 
            {
              node.explore();
              DefaultTreeModel model = (DefaultTreeModel)getModel();
              model.nodeStructureChanged(node);
            }
            setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
          }
        }
        public void treeCollapsed(TreeExpansionEvent e){}
      });


    }

    private boolean isLeafSelection()
    {
      TreePath path = getLeadSelectionPath();
      if(path == null)
        return false;

      DatabaseTreeNode node = (DatabaseTreeNode)path.getLastPathComponent();
      return node.isLeaf();
    }
    
    private DatabaseTreeNode getSelectedNode()
    {
      TreePath path = getLeadSelectionPath();
      if(path == null)
        return null;

      return (DatabaseTreeNode)path.getLastPathComponent();
    }
    
    public void dragGestureRecognized(DragGestureEvent dge)
    {
      // ignore if mouse popup trigger
      InputEvent ie = dge.getTriggerEvent();
      if(ie instanceof MouseEvent)
        if(((MouseEvent)ie).isPopupTrigger())
          return;
     
      // drag only files 
      if(isLeafSelection())
      {
        final int nlist = getSelectionCount();
        if(nlist > 1)
        {

        }
        else
        {
          dge.startDrag(DragSource.DefaultCopyDrop,     // cursor
                      null,
                      new Point(0,0),
                      (Transferable)getSelectedNode(), // transferable data
                                            this);     // drag source listener
        }
      }
      
    }

    public void dragDropEnd(DragSourceDropEvent dsde)
    {
    }

    public void dragEnter(DragSourceDragEvent dsde)
    {
    }

    public void dragExit(DragSourceEvent dse)
    {
    }

    public void dragOver(DragSourceDragEvent dsde)
    {  
    }

    public void dropActionChanged(DragSourceDragEvent dsde)
    { 
    }
}
