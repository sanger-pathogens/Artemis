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

import javax.swing.*;
import javax.swing.tree.*;
import javax.swing.table.*;
import javax.swing.filechooser.FileSystemView;

import java.io.File;
import java.io.IOException;
import java.io.FileOutputStream;
import java.io.FileNotFoundException;

import java.awt.*;
import java.awt.event.*;
import java.awt.datatransfer.*;
import java.awt.dnd.*;
import java.awt.image.BufferedImage;
import java.util.Vector;

import uk.ac.sanger.artemis.j2ssh.FileTransferProgressMonitor;
import uk.ac.sanger.artemis.j2ssh.FTProgress;

import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.SimpleEntryGroup;
import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.util.FileDocument;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.SimpleEntryInformation;
import uk.ac.sanger.artemis.components.EntryEdit;
import uk.ac.sanger.artemis.components.EntryFileDialog;
import uk.ac.sanger.artemis.components.SwingWorker;
import uk.ac.sanger.artemis.components.MessageDialog;
import uk.ac.sanger.artemis.sequence.NoSequenceException;

/**
*
* This example shows how to create a simple JTreeTable component, 
* by using a JTree as a renderer (and editor) for the cells in a 
* particular column in the JTable.  
*
* modified from the example at:
* http://java.sun.com/products/jfc/tsc/articles/treetable1/
* 
*/
public class JTreeTable extends JTable 
                 implements DragGestureListener,
                 DragSourceListener, DropTargetListener, ActionListener,
                 Autoscroll
{
  /** */
  private static final long serialVersionUID = 1L;
  /** popup menu */
  private JPopupMenu popup;
  /** busy cursor */
  private Cursor cbusy = new Cursor(Cursor.WAIT_CURSOR);
  /** done cursor */
  private Cursor cdone = new Cursor(Cursor.DEFAULT_CURSOR);
  /** file separator */
  private String fs = new String(System.getProperty("file.separator"));
  /** AutoScroll margin */
  private static final int AUTOSCROLL_MARGIN = 45;
  /** used by AutoScroll method */
  private Insets autoscrollInsets = new Insets( 0, 0, 0, 0 );

  protected TreeTableCellRenderer tree;

  public JTreeTable(TreeModel treeTableModel) 
  {
    super();
    setName("File Tree");

    DragSource dragSource = DragSource.getDefaultDragSource();

    dragSource.createDefaultDragGestureRecognizer(
       this,                             // component where drag originates
       DnDConstants.ACTION_COPY_OR_MOVE, // actions
       this);                            // drag gesture recognizer

    setDropTarget(new DropTarget(this,this));

    // Create the tree. It will be used as a renderer and editor. 
    tree = new TreeTableCellRenderer(treeTableModel); 

    // Install a tableModel representing the visible rows in the tree. 
    super.setModel(new TreeTableModelAdapter(treeTableModel, tree));

    // Force the JTable and JTree to share their row selection models. 
    tree.setSelectionModel(new DefaultTreeSelectionModel() 
    { 
      /***/
      private static final long serialVersionUID = 1L;

      // Extend the implementation of the constructor, as if: 
      /* public this() */
      {
	setSelectionModel(listSelectionModel); 
      } 
    }); 
    // Make the tree and table row heights the same. 
    tree.setRowHeight(getRowHeight());

    // Install the tree editor renderer and editor. 
    setDefaultRenderer(TreeTableModel.class, tree); 
    setDefaultEditor(TreeTableModel.class, new TreeTableCellEditor());  

    setShowGrid(false);
    setIntercellSpacing(new Dimension(3, 0)); 	        

    //Listen for when a file is selected
    MouseListener mouseListener = new MouseAdapter()
    {
      public void mouseClicked(MouseEvent me)
      {
        if(me.getClickCount() == 2 && isFileSelection() &&
           !me.isPopupTrigger())
        {
          setCursor(cbusy);
          FileNode node = getSelectedNode();
          String selected = node.getFile().getAbsolutePath();
          showFilePane(selected);
          setCursor(cdone);
        }
      }
    };
    this.addMouseListener(mouseListener);

    // Popup menu
    addMouseListener(new PopupListener());
    popup = new JPopupMenu();

    JMenuItem menuItem = new JMenuItem("Refresh");
    menuItem.addActionListener(this);
    popup.add(menuItem);
    
    JMenuItem gotoItem = new JMenuItem("GoTo Directory...");
    gotoItem.addActionListener(this);
    popup.add(gotoItem);
    popup.add(new JSeparator());
    //open menu
    JMenu openMenu = new JMenu("Open With");
    popup.add(openMenu);
    menuItem = new JMenuItem("Jemboss Alignment Editor");
    menuItem.addActionListener(this);
    openMenu.add(menuItem);
    menuItem = new JMenuItem("Artemis");
    menuItem.addActionListener(this);
    openMenu.add(menuItem);

    menuItem = new JMenuItem("Rename...");
    menuItem.addActionListener(this);
    popup.add(menuItem);
    menuItem = new JMenuItem("New Folder...");
    menuItem.addActionListener(this);
    popup.add(menuItem);
    menuItem = new JMenuItem("Delete...");
    menuItem.addActionListener(this);
    popup.add(menuItem);
    popup.add(new JSeparator());
    menuItem = new JMenuItem("De-select All");
    menuItem.addActionListener(this);
    popup.add(menuItem);

// Set the first visible column to 10 pixels wide
    TableColumn col = getColumnModel().getColumn(1);
    col.setPreferredWidth(10);
  }

  public TreeTableCellRenderer getTree()
  {
    return tree;
  }

  /**
  *
  * Get FileNode of selected node
  * @return     node that is currently selected
  *
  */
  public FileNode getSelectedNode()
  {
    TreePath path = tree.getLeadSelectionPath();
    if(path == null)
      return null;
    FileNode node = (FileNode)path.getLastPathComponent();
    return node;
  }

  /* Workaround for BasicTableUI anomaly. Make sure the UI never tries to 
   * paint the editor. The UI currently uses different techniques to 
   * paint the renderers and editors and overriding setBounds() below 
   * is not the right thing to do for an editor. Returning -1 for the 
   * editing row in this case, ensures the editor is never painted. 
   */
  public int getEditingRow()
  {
    return (getColumnClass(editingColumn) == TreeTableModel.class) ? -1 : editingRow;  
  }


  private void refresh(FileNode node)
  {
    node.removeAllChildren();
    node.reset();

    node.getChildren( ((FileSystemModel)tree.getModel()).getFilter() );
    ((DefaultTreeModel)tree.getModel()).nodeStructureChanged(node);
    tree.revalidate();
    tree.repaint(); 
    revalidate();
  }

  protected void refreshAll()
  {
    FileSystemModel model = (FileSystemModel)tree.getModel();
    Object root = model.getRoot();

    Vector vnodes = new Vector();
    addChildren(vnodes, (FileNode)root, model);

    for(int i = 0; i<vnodes.size(); i++)
      refresh((FileNode)vnodes.get(i));
    
    tree.revalidate();  
    repaint(); 
    revalidate();
  }

  private void addChildren(Vector v, FileNode node, FileSystemModel model)
  {
    int nchild = model.getChildCount(node);

    for(int i = 0; i<nchild; i++)
    {
      FileNode fn = (FileNode)model.getChild(node, i);

      if(fn.isDirectory())
        v.add(fn);
    }
  }

  /**
  *
  * Get FileNodes of selected nodes
  * @return     node that is currently selected
  *
  */
  private FileNode[] getSelectedNodes()
  {
    TreePath path[] = tree.getSelectionPaths();
    if(path == null)
      return null;

    int numberSelected = path.length;
    FileNode nodes[] = new FileNode[numberSelected];
    for(int i=0;i<numberSelected;i++)
       nodes[i] = (FileNode)path[i].getLastPathComponent();

    return nodes;
  }


  /**
  *
  * Get selected files
  * @return     node that is currently selected
  *
  */
  private File[] getSelectedFiles()
  {
    FileNode[] fn = getSelectedNodes();
    int numberSelected = fn.length;
    File files[] = new File[numberSelected];
    for(int i=0;i<numberSelected;i++)
       files[i] = fn[i].getFile();

    return files;
  }

  /**
  *
  * Popup menu actions
  * @param e    action event
  *
  */
  public void actionPerformed(ActionEvent e)
  {
    JMenuItem source = (JMenuItem)(e.getSource());
    FileNode node = getSelectedNode();

    if(source.getText().equals("Refresh"))
    {
      if(node == null)
        return;
      else if(node.isLeaf())
        node = (FileNode)node.getParent();

      node.removeAllChildren();
      node.reset();
      node.getChildren( ((FileSystemModel)tree.getModel()).getFilter() );
      ((DefaultTreeModel)tree.getModel()).nodeStructureChanged(node);
      tree.repaint();
      revalidate();
      repaint();
      return;
    }
   
    if(source.getText().startsWith("GoTo Directory"))
    {
      String dir = JOptionPane.showInputDialog(this, "GoTo Directory", 
          System.getProperty("home.dir"));
      
      File fileDir = new File(dir);
      if(!fileDir.exists())
      {
        JOptionPane.showMessageDialog(this, dir + " not found.");
        return;
      }
      else if(fileDir.isFile())
      {
        JOptionPane.showMessageDialog(this, dir + " is a file.");
        return;
      }
      FileNode rootNode = (FileNode) ((DefaultTreeModel)tree.getModel()).getRoot();
      FileNode newNode = new FileNode(fileDir);
      rootNode.add(newNode);
      
      ((DefaultTreeModel)tree.getModel()).nodeStructureChanged(rootNode);
      revalidate();
      return;
    }

    if(node == null)
    {
      JOptionPane.showMessageDialog(null,"No file selected.",
                        "Warning",
                        JOptionPane.WARNING_MESSAGE);
      return;
    }
   
    final File f = node.getFile();

    if(source.getText().equals("Jemboss Alignment Editor"))
    {
      org.emboss.jemboss.editor.AlignJFrame ajFrame =
               new org.emboss.jemboss.editor.AlignJFrame(f);
      ajFrame.setVisible(true);
    }
    else if(source.getText().equals("Artemis"))
    {
      setCursor(cbusy);
      String selected = node.getFile().getAbsolutePath();
      showFilePane(selected);
      setCursor(cdone);
    }
    else if(source.getText().equals("New Folder..."))
    {
      if(node.isLeaf())
        node = (FileNode)node.getParent();

      String path = node.getFile().getAbsolutePath();

      String inputValue = JOptionPane.showInputDialog(null,
                    "Folder Name","Create New Folder in "+path,
                    JOptionPane.QUESTION_MESSAGE);

      if(inputValue != null && !inputValue.equals("") )
      {
        String fullname = path+fs+inputValue;
        File dir = new File(fullname);

        if(dir.exists())
          JOptionPane.showMessageDialog(null, fullname+" alread exists!",
                                   "Error", JOptionPane.ERROR_MESSAGE);
        else
        {
          if(dir.mkdir())
            refresh(node);
          else
            JOptionPane.showMessageDialog(null,
                       "Cannot make the folder\n"+fullname,
                       "Error", JOptionPane.ERROR_MESSAGE);
        }
      }
    }
    else if(source.getText().equals("Delete..."))
    {
      File fn[] = getSelectedFiles();
      String[] names = new String[fn.length];
      for(int i=0; i<fn.length;i++)
        names[i] = fn[i].getAbsolutePath();  

      JList list = new JList(names);
      JScrollPane jsp = new JScrollPane(list);
      int n = JOptionPane.showConfirmDialog(null,
                                 jsp,
                                 "Delete "+fn.length+" Files",
                                 JOptionPane.YES_NO_OPTION);

      FileNode nodes[] = getSelectedNodes();
      if(n == JOptionPane.YES_OPTION)
        for(int i=0; i<nodes.length;i++)
          deleteFile(nodes[i]);
    }
    else if(source.getText().equals("Rename..."))
    {
      String inputValue = (String)JOptionPane.showInputDialog(null,
                              "New File Name","Rename "+f.getName(),
                              JOptionPane.QUESTION_MESSAGE,null,null,f.getName());

      if(inputValue != null && !inputValue.equals("") )
      {
        String path = f.getParent();
        String fullname   = path+fs+inputValue;
        File newFile = new File(fullname);

        try
        {
          renameFile(f,node,newFile.getCanonicalPath());
        }
        catch(IOException ioe){}
      }
    }
    else if(source.getText().equals("De-select All"))
      clearSelection(); 
  }

  /**
  *
  * Method to rename a file and update the filenode's.
  * @param oldFile      file to rename
  * @param oldNode      filenode to be removed
  * @param newFullName  name of the new file
  *
  */
  private void renameFile(final File oldFile, final FileNode oldNode,
                          String newFullName)
  {
    final File fnew = new File(newFullName);
    if(fnew.exists())
      JOptionPane.showMessageDialog(null, newFullName+" alread exists!",
                               "Warning", JOptionPane.ERROR_MESSAGE);
    else
    {
      if(oldFile.renameTo(fnew))
      {
        Runnable renameFileInTree = new Runnable()
        {
          public void run ()
          {
            refresh((FileNode)oldNode.getParent());
          };
        };
        SwingUtilities.invokeLater(renameFileInTree);
      }
      else
        JOptionPane.showMessageDialog(null,
                   "Cannot rename \n"+oldFile.getAbsolutePath()+
                   "\nto\n"+fnew.getAbsolutePath(), "Rename Error",
                   JOptionPane.ERROR_MESSAGE);
    }
    return;
  }


  /**
  *
  * Delete a file from the tree
  * @param node         node to delete
  *
  */
  public void deleteFile(final FileNode node)
  {
    File f = node.getFile();
    if(f.delete())
    {
      Runnable deleteFileFromTree = new Runnable()
      {
        public void run () { refresh((FileNode)node.getParent()); };
      };
      SwingUtilities.invokeLater(deleteFileFromTree);
    }
    else
      JOptionPane.showMessageDialog(null,"Cannot delete\n"+
                         f.getAbsolutePath(),"Warning",
                         JOptionPane.ERROR_MESSAGE);
  }

  /**
  *
  * Opens a JFrame with the file contents displayed.
  * @param filename     file name to display
  *
  */
  public void showFilePane(final String filename)
  {
    SwingWorker entryWorker = new SwingWorker()
    {
      EntryEdit entry_edit;
      public Object construct()
      {
        try
        {
          EntryInformation new_entry_information =
             new SimpleEntryInformation(Options.getArtemisEntryInformation());

          final Entry entry =  new Entry(EntryFileDialog.getEntryFromFile(
                         null, new FileDocument(new File(filename)),
                         new_entry_information, true));
          if(entry == null)
            return null;

          final EntryGroup entry_group =
              new SimpleEntryGroup(entry.getBases());

          entry_group.add(entry);
          entry_edit = new EntryEdit(entry_group);
          return null;
        }
        catch(NoSequenceException e)
        {
          new MessageDialog(null, "read failed: entry contains no sequence");
        }
        catch(OutOfRangeException e)
        {
          new MessageDialog(null, "read failed: one of the features in " +
                     " the entry has an out of range " +
                     "location: " + e.getMessage());

        }
        catch(NullPointerException npe)
        {
          npe.printStackTrace();
        }

        return null;
      }

      public void finished()
      {
        if(entry_edit != null)
          entry_edit.setVisible(true);
      }
    };
    entryWorker.start();

  }

  /**
  *
  * Return true if selected node is a file
  * @return true is a file is selected, false if
  *         a directory is selected
  *
  */
  public boolean isFileSelection()
  {
    TreePath path = tree.getLeadSelectionPath();
    if(path == null)
      return false;

    FileNode node = (FileNode)path.getLastPathComponent();
    return node.isLeaf();
  }


  // 
  // The renderer used to display the tree nodes, a JTree.  
  //
  public class TreeTableCellRenderer extends JTree implements TableCellRenderer 
  {
    /** */
    private static final long serialVersionUID = 1L;
    protected int visibleRow;
   
    public TreeTableCellRenderer(TreeModel model) 
    { 
      super(model); 
    }

    public void setBounds(int x, int y, int w, int h)
    {
      super.setBounds(x, 0, w, JTreeTable.this.getHeight());
    }

    public void paint(Graphics g)
    {
      g.translate(0, -visibleRow * getRowHeight());
      super.paint(g);
    }

    public Component getTableCellRendererComponent(JTable table,
					       Object value,
					       boolean isSelected,
					       boolean hasFocus,
					       int row, int column) 
    {
      if(isSelected)
        setBackground(table.getSelectionBackground());
      else
	setBackground(table.getBackground());
       
      visibleRow = row;
      return this;
    }
  }

  // 
  // The editor used to interact with tree nodes, a JTree.  
  //

  public class TreeTableCellEditor extends AbstractCellEditor implements TableCellEditor 
  {
    public Component getTableCellEditorComponent(JTable table, Object value,
					     boolean isSelected, int r, int c) 
    {
      return tree;
    }
  }

  /**
  *
  * Popup menu listener
  *
  */
  class PopupListener extends MouseAdapter
  {
    public void mousePressed(MouseEvent e)
    {
      maybeShowPopup(e);
    }

    public void mouseReleased(MouseEvent e)
    {
      maybeShowPopup(e);
    }

    private void maybeShowPopup(MouseEvent e)
    {
      if(e.isPopupTrigger())
        popup.show(e.getComponent(),
                e.getX(), e.getY());
    }
  }

  private boolean writeByteFile(byte[] contents, File fn)
  {
    if(fn.exists())
    {
      int n = JOptionPane.showConfirmDialog(null,
                                 "Overwrite \n"+fn.getName()+"?",
                                 "Overwrite File",
                                 JOptionPane.YES_NO_OPTION);
      if(n == JOptionPane.NO_OPTION)
        return false;
    }
    else if(!fn.getParentFile().canWrite())
      JOptionPane.showMessageDialog(null,"Cannot write to "+fn.getName(),
                        "Write Permission Denied",
                        JOptionPane.WARNING_MESSAGE);

    try
    {
      FileOutputStream out = new FileOutputStream(fn);
      out.write(contents);
      out.close(); 
    }
    catch(FileNotFoundException fnfe) {return false;}
    catch(IOException ioe) {return false;}
  
    return true;
  }


  private void localDrop(DropTargetDropEvent e, Vector vnode, FileNode dropNode)
  {
    try
    {
      for(int i=0; i<vnode.size(); i++)
      {
        FileNode fn = (FileNode)vnode.get(i);
//      fn = getNode(fn.getFile().getAbsolutePath());

        if (dropNode.isLeaf())
        {
          e.rejectDrop();
          return;
        }

        String dropDir = dropNode.getFile().getAbsolutePath();
        String newFullName = dropDir+fs+fn.toString();
        renameFile(fn.getFile(),fn,newFullName);
      }
    }
    catch(Exception ufe){}
  }

  private void remoteDrop(final DropTargetDropEvent e,
                          final Vector vnode, final FileNode dropNode)
  {
    SwingWorker getFileWorker = new SwingWorker()
    {
      FileTransferProgressMonitor monitor;

      public Object construct()
      {
        try
        {
          monitor = new FileTransferProgressMonitor(JTreeTable.this);
          for(int i=0; i<vnode.size(); i++)
          {
            final RemoteFileNode fn = (RemoteFileNode)vnode.get(i);
            final File dropDest;
            String dropDir = null;
            if (dropNode.isLeaf())
            {
              FileNode pn = (FileNode)dropNode.getParent();
              dropDir = pn.getFile().getAbsolutePath();
              dropDest = new File(dropDir,fn.getFile());
            }
            else
            {
              dropDir = dropNode.getFile().getAbsolutePath();
              dropDest = new File(dropDir,fn.getFile());
            }

            try
            {
              FTProgress progress = monitor.add(fn.getFile());

              final byte[] contents = fn.getFileContents(progress);
              //final String ndropDir = dropDir;

              Runnable updateTheTree = new Runnable()
              {
                public void run ()
                {
                  if(writeByteFile(contents, dropDest))
                  {
                    if(dropNode.isLeaf())
                      refresh((FileNode)dropNode.getParent());
                    else
                      refresh(dropNode);
                  }
                };
              };
              SwingUtilities.invokeLater(updateTheTree);
            }
            catch (Exception exp)
            {
              System.out.println("FileTree: caught exception");
              exp.printStackTrace();
            }
          }
          e.getDropTargetContext().dropComplete(true);
        }
        catch (Exception exp)
        {
          e.rejectDrop();
        }
        return null;
      }

      public void finished()
      {
        if(monitor != null)
          monitor.close();
      }
    };
    getFileWorker.start();
  }


////////////////////
// DRAG AND DROP
////////////////////
// drag source
  public void dragGestureRecognized(DragGestureEvent e)
  {
    // ignore if mouse popup trigger
    InputEvent ie = e.getTriggerEvent();
    if(ie instanceof MouseEvent)
      if(((MouseEvent)ie).isPopupTrigger())
        return;
   
    // drag only files 
    if(isFileSelection())
    {
      final int nlist = tree.getSelectionCount();
      if(nlist > 1)
      {
        TransferableFileNodeList list = new TransferableFileNodeList(nlist);
        FileNode nodes[] = getSelectedNodes();
        for(int i=0; i<nodes.length; i++)
          list.add(nodes[i]);

        BufferedImage buff = getImage(nodes[0].getFile());
        e.startDrag(DragSource.DefaultCopyDrop,     // cursor
                    buff,
                    new Point(0,0),
                   (Transferable)list,              // transferable data
                                         this);     // drag source listener
      }
      else
      {
        BufferedImage buff = getImage(getSelectedNode().getFile());
        e.startDrag(DragSource.DefaultCopyDrop,     // cursor
                    buff,
                    new Point(0,0),
                    (Transferable)getSelectedNode(), // transferable data
                                          this);     // drag source listener
      }
    }
  }

  /**
  *
  * FileSystemView provides a platform-independent way to get the 
  * appropriate icon 
  *
  */
  private BufferedImage getImage(File temp)
  {
    // get the right icon
    FileSystemView fsv = FileSystemView.getFileSystemView( );
    Icon icn = fsv.getSystemIcon(temp);

    Toolkit tk = Toolkit.getDefaultToolkit( );
    Dimension dim = tk.getBestCursorSize(
	icn.getIconWidth( ),icn.getIconHeight( ));
    BufferedImage buff = new BufferedImage(dim.width,dim.height,
			 BufferedImage.TYPE_INT_ARGB);
    icn.paintIcon(this,buff.getGraphics( ),0,0);
    return buff;
  }

  public void dragDropEnd(DragSourceDropEvent e) {}
  public void dragEnter(DragSourceDragEvent e) {}
  public void dragExit(DragSourceEvent e) {}
  public void dragOver(DragSourceDragEvent e) {}
  public void dropActionChanged(DragSourceDragEvent e) {}

// drop sink
  public void dragEnter(DropTargetDragEvent e)
  {
    if(e.isDataFlavorSupported(FileNode.FILENODE) ||
       e.isDataFlavorSupported(RemoteFileNode.REMOTEFILENODE) ||
       e.isDataFlavorSupported(TransferableFileNodeList.TRANSFERABLEFILENODELIST))
      e.acceptDrag(DnDConstants.ACTION_COPY_OR_MOVE);
  }

  public void drop(final DropTargetDropEvent e)
  {
    final Transferable t = e.getTransferable();
    final FileNode dropNode = getSelectedNode();
    if(dropNode == null)
    {
      e.rejectDrop();
      return;
    }

    if(t.isDataFlavorSupported(TransferableFileNodeList.TRANSFERABLEFILENODELIST))
    {
      try
      {
        TransferableFileNodeList filelist = (TransferableFileNodeList)
                  t.getTransferData(TransferableFileNodeList.TRANSFERABLEFILENODELIST);

        if(filelist.get(0) instanceof RemoteFileNode)
          remoteDrop(e, filelist, dropNode);
        else
          localDrop(e, filelist, dropNode);
      }
      catch(UnsupportedFlavorException exp){}
      catch(IOException ioe){}
    }
    else if(t.isDataFlavorSupported(FileNode.FILENODE))
    {
      try
      {
        Vector v = new Vector();
        FileNode fn = (FileNode)t.getTransferData(FileNode.FILENODE);
        v.add(fn);
        localDrop(e, v, dropNode);
      }
      catch(Exception ufe){}
    }
    else if(t.isDataFlavorSupported(RemoteFileNode.REMOTEFILENODE))
    {
      try
      {
        Vector v = new Vector();
        final RemoteFileNode fn =
             (RemoteFileNode)t.getTransferData(RemoteFileNode.REMOTEFILENODE);
        v.add(fn);
        remoteDrop(e, v, dropNode);
      }
      catch(Exception ufe){}
    }
    else
    {
      e.rejectDrop();
      return;
    }
  }

  /**
  *
  * When a suitable DataFlavor is offered over a remote file
  * node the node is highlighted/selected and the drag
  * accepted. Otherwise the drag is rejected.
  *
  */
  public void dragOver(DropTargetDragEvent e)
  {
    if(e.isDataFlavorSupported(FileNode.FILENODE))
    {
      Point ploc = e.getLocation();
      TreePath ePath = tree.getPathForLocation(ploc.x,ploc.y);
      if (ePath == null)
      {
        e.rejectDrag();
        return;
      }
      FileNode node = (FileNode)ePath.getLastPathComponent();
      if(!node.isDirectory())
        e.rejectDrag();
      else
      {
        tree.setSelectionPath(ePath);
        e.acceptDrag(DnDConstants.ACTION_COPY_OR_MOVE);
      }
    }
    else if(e.isDataFlavorSupported(RemoteFileNode.REMOTEFILENODE) ||
            e.isDataFlavorSupported(TransferableFileNodeList.TRANSFERABLEFILENODELIST))
    {
      Point ploc = e.getLocation();
      TreePath ePath = tree.getPathForLocation(ploc.x,ploc.y);
      if (ePath == null)
        e.rejectDrag();
      else
      {
        tree.setSelectionPath(ePath);
        e.acceptDrag(DnDConstants.ACTION_COPY_OR_MOVE);
      }
    }
    else
      e.rejectDrag();

    return;
  }

  public void dropActionChanged(DropTargetDragEvent e) {}
  public void dragExit(DropTargetEvent e){}

////////////////////
// AUTO SCROLLING //
////////////////////
  /**
  *
  * Handles the auto scrolling of the JTree.
  * @param location The location of the mouse.
  *
  */
  public void autoscroll( Point location )
  {
    int top = 0, left = 0, bottom = 0, right = 0;
    Dimension size = getSize();
    Rectangle rect = getVisibleRect();
    int bottomEdge = rect.y + rect.height;
    int rightEdge = rect.x + rect.width;
    if( location.y - rect.y < AUTOSCROLL_MARGIN && rect.y > 0 )
      top = AUTOSCROLL_MARGIN;
    if( location.x - rect.x < AUTOSCROLL_MARGIN && rect.x > 0 )
      left = AUTOSCROLL_MARGIN;
    if( bottomEdge - location.y < AUTOSCROLL_MARGIN && bottomEdge < size.height )
      bottom = AUTOSCROLL_MARGIN;
    if( rightEdge - location.x < AUTOSCROLL_MARGIN && rightEdge < size.width )
      right = AUTOSCROLL_MARGIN;
    rect.x += right - left;
    rect.y += bottom - top;
    scrollRectToVisible( rect );
  }


  /**
  *
  * Gets the insets used for the autoscroll.
  * @return The insets.
  *
  */
  public Insets getAutoscrollInsets()
  {
    Dimension size = getSize();
    Rectangle rect = getVisibleRect();
    autoscrollInsets.top = rect.y + AUTOSCROLL_MARGIN;
    autoscrollInsets.left = rect.x + AUTOSCROLL_MARGIN;
    autoscrollInsets.bottom = size.height - (rect.y+rect.height) + AUTOSCROLL_MARGIN;
    autoscrollInsets.right  = size.width - (rect.x+rect.width) + AUTOSCROLL_MARGIN;
    return autoscrollInsets;
  }

}

