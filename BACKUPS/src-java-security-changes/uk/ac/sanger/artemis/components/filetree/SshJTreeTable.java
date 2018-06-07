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

import java.awt.*;
import java.awt.event.*;
import java.awt.datatransfer.*;
import java.awt.dnd.*;
import java.awt.image.BufferedImage;
import java.util.Vector;
import java.util.Enumeration;

import uk.ac.sanger.artemis.j2ssh.FileTransferProgressMonitor;
import uk.ac.sanger.artemis.j2ssh.FTProgress;

import uk.ac.sanger.artemis.util.RemoteFileDocument;
import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.SimpleEntryGroup;
import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.SimpleEntryInformation;
import uk.ac.sanger.artemis.components.EntryEdit;
import uk.ac.sanger.artemis.components.EntryFileDialog;
import uk.ac.sanger.artemis.components.SwingWorker;
import uk.ac.sanger.artemis.components.MessageDialog;
import uk.ac.sanger.artemis.sequence.NoSequenceException;

/**
*
* This example shows how to create a simple SshJTreeTable component, 
* by using a JTree as a renderer (and editor) for the cells in a 
* particular column in the JTable.  
*
* modified from the example at:
* http://java.sun.com/products/jfc/tsc/articles/treetable1/
* 
*/
public class SshJTreeTable extends JTable 
                 implements DragGestureListener,
                 DragSourceListener, DropTargetListener, ActionListener,
                 Autoscroll
{
  /**  */
  private static final long serialVersionUID = 1L;
  /** popup menu */
  private JPopupMenu popup;
  /** busy cursor */
  private Cursor cbusy = new Cursor(Cursor.WAIT_CURSOR);
  /** done cursor */
  private Cursor cdone = new Cursor(Cursor.DEFAULT_CURSOR);
  /** line separator */
  private String ls = new String(System.getProperty("line.separator"));
  /** AutoScroll margin */
  private static final int AUTOSCROLL_MARGIN = 45;
  /** used by AutoScroll method */
  private Insets autoscrollInsets = new Insets( 0, 0, 0, 0 );
  protected TreeTableCellRenderer tree;
  private JFrame frame;

  public SshJTreeTable(TreeModel treeTableModel, final JFrame frame) 
  {
    super();

    this.frame = frame;
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
      /** */
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

// popup menu
    addMouseListener(new PopupListener());
    popup = new JPopupMenu();
    JMenuItem menuItem = new JMenuItem("Refresh");
    menuItem.addActionListener(this);
    popup.add(menuItem);
    popup.add(new JSeparator());
//open menu
    JMenu openMenu = new JMenu("Open With");
    popup.add(openMenu);
    menuItem = new JMenuItem("Jemboss Aligmnment Editor");
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


    //Listen for when a file is selected
    addMouseListener(new MouseListener()
    {
      public void mouseClicked(MouseEvent me)
      {
        if(me.getClickCount() == 2 && isFileSelection() &&
           !me.isPopupTrigger())
        {
          RemoteFileNode node = (RemoteFileNode)tree.getLastSelectedPathComponent();
          if(node==null)
            return;
          frame.setCursor(cbusy);
          if(node.isLeaf())
            showFilePane(node);
          frame.setCursor(cdone);
        }
      }
      public void mousePressed(MouseEvent me){}
      public void mouseEntered(MouseEvent me){}
      public void mouseExited(MouseEvent me){}
      public void mouseReleased(MouseEvent me){}
    });

  }

  /**
  *
  * Determine if selected node is a file
  * @return     true if the selected node is a file
  *
  */
  public boolean isFileSelection()
  {
    TreePath path = tree.getLeadSelectionPath();
    if(path == null)
      return false;
    RemoteFileNode node = (RemoteFileNode)path.getLastPathComponent();
    return !node.isDirectory();
  }

  /**
  *
  * Get FileNode of selected node
  * @return     node that is currently selected
  *
  */
  public RemoteFileNode getSelectedNode()
  {
    TreePath path = tree.getLeadSelectionPath();
    if(path == null)
      return null;
    RemoteFileNode node = (RemoteFileNode)path.getLastPathComponent();
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


  protected void refresh(RemoteFileNode node)
  {
    node.reset();
    node.removeAllChildren();
    node.getChildren();
    ((DefaultTreeModel)tree.getModel()).nodeStructureChanged(node);
  }

  /**
  *
  * Get FileNodes of selected nodes
  * @return     node that is currently selected
  *
  */
  private RemoteFileNode[] getSelectedNodes()
  {
    TreePath path[] = tree.getSelectionPaths();
    if(path == null)
      return null;

    int numberSelected = path.length;
    RemoteFileNode nodes[] = new RemoteFileNode[numberSelected];
    for(int i=0;i<numberSelected;i++)
       nodes[i] = (RemoteFileNode)path[i].getLastPathComponent();

    return nodes;
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
    final RemoteFileNode node = getSelectedNode();
    if(node == null)
    {
      JOptionPane.showMessageDialog(null,"No file selected.",
                                    "Warning",
                                    JOptionPane.WARNING_MESSAGE);
      return;
    }

    final String fn = node.getFullName();
    final String parent = node.getPathName();
    String rootPath = node.getRootDir();
    RemoteFileNode pn = node;

    if(source.getText().equals("Refresh"))
    {
      if(node.isLeaf() && node.getParentNode() != null)
        refresh(node.getParentNode());
      else
        refresh(node);
    }
    else if(source.getText().equals("Jemboss Aligmnment Editor"))
    {
      FileTransferProgressMonitor monitor =
                     new FileTransferProgressMonitor(SshJTreeTable.this);
      FTProgress progress = monitor.add(node.getFile());
  
      final byte[] contents = node.getFileContents(progress);
      monitor.close();

      org.emboss.jemboss.editor.AlignJFrame ajFrame =
                    new org.emboss.jemboss.editor.AlignJFrame(new String(contents), fn);
      ajFrame.setVisible(true);
    }
    else if(source.getText().equals("Artemis"))
      showFilePane(node);
    else if(source.getText().equals("New Folder..."))
    {
      final String inputValue = JOptionPane.showInputDialog(null,
                          "Folder Name","Create New Folder in",
                          JOptionPane.QUESTION_MESSAGE);

      if(node.isLeaf())
        pn = (RemoteFileNode)node.getParent();


      String newNode = pn.getServerName();
      if(!newNode.endsWith("/"))
        newNode = newNode.concat("/");
      newNode = newNode.concat(inputValue);

      if(nodeExists(pn,newNode))
        return;

      if(inputValue != null && !inputValue.equals("") )
      {
        node.mkdir(newNode);

        Runnable addDirToTree = new Runnable()
        {
          public void run () { refresh(node); };
        };
        SwingUtilities.invokeLater(addDirToTree);
      }
    }
    else if(source.getText().equals("Delete..."))
    {
      RemoteFileNode nodes[] = getSelectedNodes();
      String sname[] = new String[nodes.length];
      for(int i=0;i<nodes.length;i++)
        sname[i] = nodes[i].getServerName();

      JList list = new JList(sname);
      JScrollPane jsp = new JScrollPane(list);
      int n = JOptionPane.showConfirmDialog(null,
             jsp, "Delete "+nodes.length+" File(s)",
             JOptionPane.YES_NO_OPTION);

      if(n == JOptionPane.YES_OPTION)
      {
        frame.setCursor(cbusy);
        for(int i=0;i<nodes.length;i++)
          deleteNode(nodes[i]);
        frame.setCursor(cdone);
      }
    }
    else if(source.getText().equals("De-select All"))
      clearSelection();
    else if(source.getText().equals("Rename..."))
    {
      if(node.isLeaf())
      {
        String inputValue = (String)JOptionPane.showInputDialog(null,
                              "New File Name","Rename "+fn,
                              JOptionPane.QUESTION_MESSAGE,null,null,fn);

        pn = (RemoteFileNode)node.getParent();

        if(inputValue != null && !inputValue.equals("") )
        {
          String newfile = null;
          if(parent.endsWith("/"))
            newfile = parent+inputValue;
          else
            newfile = parent+"/"+inputValue;

          if(!nodeExists(pn,newfile))
          {
            rename(node,rootPath+"/"+inputValue);

            Runnable addDirToTree = new Runnable()
            {
              public void run () 
              {
                refresh(node.getParentNode()); 
              };
            };
            SwingUtilities.invokeLater(addDirToTree);
          }
        }
      }
    }

  }


  /**
  *
  * Rename a node from the tree
  * @param node                 node to rename
  * @param newfile              new file name
  *
  */
  private void rename(final RemoteFileNode node,
                      final String newfile)
  {
    frame.setCursor(cbusy);

    boolean lrename = node.rename(newfile);
    frame.setCursor(cdone);

    if(!lrename)
      return;
  }


  /**
  *
  * Delete a node (file or directory) from the tree
  * and from the server
  * @param node node to remove
  *
  */
  private void deleteNode(final RemoteFileNode node)
  {
    boolean deleted = false;
    deleted = node.delete();

    if(!deleted && !node.isLeaf())
      JOptionPane.showMessageDialog(null,"Cannot delete"+ls+
                 node.getServerName()+
                 ls+"this directory is not empty","Warning",
                 JOptionPane.ERROR_MESSAGE);
    else if(deleted)
    {
      Runnable deleteFileFromTree = new Runnable()
      {
        public void run () 
        {
          if(node.isLeaf() && node.getParentNode() != null)
            refresh(node.getParentNode());
          else
            refresh(node); 
        };
      };
      SwingUtilities.invokeLater(deleteFileFromTree);
    }

  }


  /**
  *
  * Test if a child node exists
  * @param parentNode   parent node
  * @param child        child to test for
  *
  */
  public boolean nodeExists(RemoteFileNode parentNode, final String child)
  {
    if(!parentNode.isDirectory())
      parentNode = (RemoteFileNode)parentNode.getParent();

    RemoteFileNode childNode = getChildNode(parentNode,child);

    if(childNode != null)
    {
      Runnable warnMsg = new Runnable()
      {
        public void run ()
        {
          String ls = System.getProperty("line.separator");
          JOptionPane.showMessageDialog(null, child+ls+" already exists!",
                                    "File Exists",
                                    JOptionPane.ERROR_MESSAGE);
        };
      };
      SwingUtilities.invokeLater(warnMsg);
      return true;
    }

    return false;
  }


  /**
  *
  * Gets the child node of a parent node
  * @param parent       parent node
  * @param childName    name of child
  * @return the child node
  *
  */
  private RemoteFileNode getChildNode(RemoteFileNode parent, String childName)
  {
    for(Enumeration children = parent.children();
                    children.hasMoreElements() ;)
    {
      RemoteFileNode childNode = (RemoteFileNode)children.nextElement();
      String nodeName = childNode.getServerName();
//    System.out.println(nodeName+" childName= "+childName);
      if(childName.equals(nodeName))
        return childNode;
    }

    return null;
  }


  /**
  *
  * Opens a JFrame with the file contents displayed.
  * @param filename     file name
  *
  */
  public static void showFilePane(final RemoteFileNode node)
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
                         null, new RemoteFileDocument(node),
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
        catch(NullPointerException npe){}

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
      super.setBounds(x, 0, w, SshJTreeTable.this.getHeight());
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

  private void remoteDrop(DropTargetDropEvent e, Vector vnode)
  {
    try
    {
      Point ploc = e.getLocation();
      TreePath dropPath = tree.getPathForLocation(ploc.x,ploc.y);
      if(dropPath != null)
      {
        RemoteFileNode drop = (RemoteFileNode)dropPath.getLastPathComponent();
        final RemoteFileNode fdropPath;
        if(!drop.isDirectory())
          fdropPath = (RemoteFileNode)drop.getParent();
        else
          fdropPath = drop;

        for(int i=0; i<vnode.size();i++)
        {
          final RemoteFileNode fn = (RemoteFileNode)vnode.get(i);
          
          String serverName = fdropPath.getServerName()+"/"+fn.getFile();
          if(!nodeExists(fdropPath,serverName))
          {
            rename(fn, serverName);

            Runnable addDirToTree = new Runnable()
            {
              public void run ()
              {
                refresh(fdropPath);
                repaint();
              };
            };
            SwingUtilities.invokeLater(addDirToTree);
          }
        }
      }
    }
    catch(Exception ex){ ex.printStackTrace(); }
  }

  private void localDrop(final DropTargetDropEvent e, final Vector fn)
  {
    //put this in separate thread for progress bar
    SwingWorker putWorker = new SwingWorker()
    {
      FileTransferProgressMonitor monitor;
      public Object construct()
      {
        try
        {
          Point ploc = e.getLocation();
          TreePath dropPath = tree.getPathForLocation(ploc.x,ploc.y);

          if(dropPath != null)
          {
            monitor = new FileTransferProgressMonitor(SshJTreeTable.this);
            for(int i=0; i<fn.size(); i++)
            {
              final File lfn = ((FileNode)fn.get(i)).getFile();

              RemoteFileNode pn = (RemoteFileNode)dropPath.getLastPathComponent();
              if(!pn.isDirectory())
                pn = (RemoteFileNode)pn.getParent();

              String serverName;
              if(pn.getServerName().endsWith("/"))
                serverName = pn.getServerName()+lfn.getName();
              else
                serverName = pn.getServerName()+"/"+lfn.getName();

              if(!nodeExists(pn,serverName))
              {
                FTProgress progress = monitor.add(lfn.getName());
                pn.put(lfn, progress);

                try
                {
                  //add file to remote file tree
                  if(pn.isLeaf())
                    pn = (RemoteFileNode)pn.getParent();

                  refresh(pn);
                  Object[] children = pn.getChildren();
                  for(int j = 0; j<children.length; j++)
                  {
                    String childNode = ((RemoteFileNode)children[j]).getFile();
                    if(childNode.equals(lfn.getName()))
                      tree.scrollPathToVisible(new TreePath( ((RemoteFileNode)children[j]).getPath()));
                  }
                }
                catch (Exception exp) {}
              }
              else
                e.rejectDrop();
            }
          }
        }
        catch (Exception ex) {}

        return null;
      }

      public void finished()
      {
        if(monitor != null)
          monitor.close();
      }
    };
    putWorker.start();
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
        RemoteFileNode nodes[] = getSelectedNodes();
        for(int i=0; i<nodes.length; i++)
          list.add(nodes[i]);

        BufferedImage buff = getImage();
        e.startDrag(DragSource.DefaultCopyDrop,     // cursor
                    buff,
                    new Point(0,0),
                   (Transferable)list,              // transferable data
                                         this);     // drag source listener
      }
      else
      {
        BufferedImage buff = getImage();
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
  private BufferedImage getImage()
  {
    File temp = new File(System.getProperty("user.dir"));
    File[] files = temp.listFiles();
    for(int i=0; i<files.length; i++)
    {
      if(files[i].isFile())
      {
        temp = files[i];
        break;
      }
    }

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

  /**
  *
  * FileSystemView provides a platform-independent way to get the 
  * appropriate icon 
  *
  */
  /*private BufferedImage getImage(File temp)
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
  }*/

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

    if(t.isDataFlavorSupported(TransferableFileNodeList.TRANSFERABLEFILENODELIST))
    {
      try
      {
        TransferableFileNodeList filelist = (TransferableFileNodeList)
                  t.getTransferData(TransferableFileNodeList.TRANSFERABLEFILENODELIST);

        if(filelist.get(0) instanceof RemoteFileNode)
          remoteDrop(e, filelist);
        else
          localDrop(e, filelist);
      }
      catch(UnsupportedFlavorException exp){}
      catch(IOException ioe){}
    }
    else if(t.isDataFlavorSupported(RemoteFileNode.REMOTEFILENODE))
    {
      try
      {
        Vector v = new Vector();
        RemoteFileNode node =
            (RemoteFileNode)t.getTransferData(RemoteFileNode.REMOTEFILENODE);
        v.add(node);
        remoteDrop(e, v);
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
        localDrop(e, v);
      }
      catch(UnsupportedFlavorException exp){}
      catch(IOException ioe){}
    }
    else
      e.rejectDrop();
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
    if (e.isDataFlavorSupported(FileNode.FILENODE))
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
    else if(e.isDataFlavorSupported(RemoteFileNode.REMOTEFILENODE) ||
            e.isDataFlavorSupported(TransferableFileNodeList.TRANSFERABLEFILENODELIST))
    {
      Point ploc = e.getLocation();
      TreePath ePath = tree.getPathForLocation(ploc.x,ploc.y);
      if (ePath == null)
      {
        e.rejectDrag();
        return;
      }

      RemoteFileNode node = (RemoteFileNode)ePath.getLastPathComponent();
      if(!node.isDirectory())
      {
        e.rejectDrag();
        return;
      }
      tree.setSelectionPath(ePath);
      e.acceptDrag(DnDConstants.ACTION_COPY_OR_MOVE);
    }
    else
      e.rejectDrag();
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

