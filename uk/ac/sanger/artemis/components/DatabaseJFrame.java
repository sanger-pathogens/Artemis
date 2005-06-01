/* DatabaseJFrame.java
 *
 * created: June 2005
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2000  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/DatabaseJFrame.java,v 1.2 2005-06-01 11:37:04 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.sequence.*;
import uk.ac.sanger.artemis.util.InputStreamProgressListener;
import uk.ac.sanger.artemis.util.OutOfRangeException;

import javax.swing.JFrame;
import javax.swing.JTree;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JMenuBar;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.Dimension;
import java.awt.Toolkit;
import java.awt.Cursor;
import java.io.*;
import javax.swing.tree.DefaultMutableTreeNode;

public class DatabaseJFrame extends JFrame
{

  public DatabaseJFrame(final DatabaseEntrySource entry_source,
                        final ArtemisMain art_main)
  {
    super("PSU Organism");
    final JTree tree = entry_source.getDatabaseTree();

    //Listen for when the selection changes.
    MouseListener mouseListener = new MouseAdapter()
    {
      public void mouseClicked(MouseEvent e)
      {
        if(e.getClickCount() == 2 && !e.isPopupTrigger())
          showSelected(entry_source,tree,art_main);
      }
    };
    tree.addMouseListener(mouseListener);

    JScrollPane scroll = new JScrollPane(tree);

    Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
    Dimension dim_frame = new Dimension(screen.width*2 / 10, screen.height*6 / 10);
    scroll.setPreferredSize(dim_frame);

    setJMenuBar(makeMenuBar(entry_source,tree,art_main));

    getContentPane().add(scroll);
    pack();
    Utilities.rightJustifyFrame(this);
  }


  /**
  *
  * Show the selected sequence in the tree
  *
  */
  private void showSelected(final DatabaseEntrySource entry_source,
                            final JTree tree, final ArtemisMain art_main)
  {
    Cursor cbusy = new Cursor(Cursor.WAIT_CURSOR);
    Cursor cdone = new Cursor(Cursor.DEFAULT_CURSOR);

    tree.setCursor(cbusy);
    DefaultMutableTreeNode node = entry_source.getSelectedNode(tree);
    if(!node.isLeaf())
      return;
    String id =  entry_source.getEntryID((String)node.getUserObject());
    if(id != null)
      getEntryEditFromDatabase(id, entry_source, art_main);

    tree.setCursor(cdone);
  }


  /**
  *
  * Retrieve a database entry.
  *
  */
  private void getEntryEditFromDatabase(final String id,
                                        final DatabaseEntrySource entry_source,
                                        final ArtemisMain art_main)
  {
    SwingWorker entryWorker = new SwingWorker()
    {
      public Object construct()
      {
        try
        {
          final InputStreamProgressListener progress_listener =
                                     art_main.getInputStreamProgressListener();

          final Entry entry = entry_source.getEntry(id, progress_listener);
          if(entry == null)
            return null;

          final EntryEdit new_entry_edit = art_main.makeEntryEdit(entry);
          new_entry_edit.setVisible(true);
        }
        catch(OutOfRangeException e)
        {
          new MessageDialog(art_main, "read failed: one of the features in " +
                 " the entry has an out of range " +
                 "location: " + e.getMessage());
        }
        catch(NoSequenceException e)
        {
          new MessageDialog(art_main, "read failed: entry contains no sequence");
        }
        catch(IOException e)
        {
          new MessageDialog(art_main, "read failed due to IO error: " + e);
        }
        return null;
      }

    };
    entryWorker.start();

  }

  private JMenuBar makeMenuBar(final DatabaseEntrySource entry_source,
                         final JTree tree, final ArtemisMain art_main)
  {
    JMenuBar mBar = new JMenuBar();
    JMenu fileMenu = new JMenu("File");
    mBar.add(fileMenu);

    JMenuItem fileShow = new JMenuItem("Open Selected Sequence ...");
    fileShow.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        showSelected(entry_source,tree,art_main);
      }
    });
    fileMenu.add(fileShow);
    fileMenu.add(new JSeparator());

    JMenuItem fileMenuClose = new JMenuItem("Close");
    fileMenuClose.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        setVisible(false);
      }
    });
    fileMenu.add(fileMenuClose);

    return mBar;
  }
}
