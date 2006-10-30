/* DatabaseJPanel.java
 *
 * created: June 2005
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2005,2006  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/database/DatabaseJPanel.java,v 1.2 2006-10-30 10:45:58 tjc Exp $
 */

package uk.ac.sanger.artemis.components.database;

import uk.ac.sanger.artemis.components.*;

import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.sequence.*;
import uk.ac.sanger.artemis.util.InputStreamProgressEvent;
import uk.ac.sanger.artemis.util.InputStreamProgressListener;
import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.io.DatabaseDocumentEntry;

import javax.swing.JTree;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JMenuBar;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JPanel;
import javax.swing.JLabel;
import javax.swing.BorderFactory;
import javax.swing.border.Border;
import javax.swing.tree.TreePath;

import org.gmod.schema.organism.Organism;

import java.awt.BorderLayout;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.Dimension;
import java.awt.Toolkit;
import java.awt.Cursor;
import java.awt.FontMetrics;
import java.io.*;
import java.net.ConnectException;
import java.sql.SQLException;
import java.util.Arrays;
import java.util.Hashtable;
import java.util.List;

public class DatabaseJPanel extends JPanel
{
  private JLabel status_line = new JLabel("");
  private boolean splitGFFEntry = false;
  private JTree tree;
  
  public DatabaseJPanel(final DatabaseEntrySource entry_source,
                        final Splash splash_main)
  {
    setLayout(new BorderLayout());
    tree = getDatabaseTree(entry_source);

    // Listen for when the selection changes.
    MouseListener mouseListener = new MouseAdapter()
    {
      public void mouseClicked(MouseEvent e)
      {
        if(e.getClickCount() == 2 && !e.isPopupTrigger())
          showSelected(entry_source, tree, splash_main);
      }
    };
    tree.addMouseListener(mouseListener);

    JScrollPane scroll = new JScrollPane(tree);

    Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
    Dimension dim_frame = new Dimension(screen.width * 2 / 10,
                                        screen.height * 6 / 10);
    scroll.setPreferredSize(dim_frame);

    add(scroll, BorderLayout.CENTER);

    final FontMetrics fm = this.getFontMetrics(status_line.getFont());

    final int font_height = fm.getHeight() + 10;
    status_line.setMinimumSize(new Dimension(100, font_height));
    status_line.setPreferredSize(new Dimension(100, font_height));

    Border loweredbevel = BorderFactory.createLoweredBevelBorder();
    Border raisedbevel = BorderFactory.createRaisedBevelBorder();
    Border compound = BorderFactory.createCompoundBorder(raisedbevel,
                                                         loweredbevel);
    status_line.setBorder(compound);
    add(status_line, BorderLayout.SOUTH);
    
    JLabel title_line = new JLabel("Database List");
    title_line.setMinimumSize(new Dimension(100, font_height));
    title_line.setPreferredSize(new Dimension(100, font_height));
    title_line.setBorder(compound);
    add(title_line, BorderLayout.NORTH);
  }

  /**
   * Show the selected sequence in the tree
   * @param entry_source
   * @param tree
   * @param splash_main
   */
  private void showSelected(final DatabaseEntrySource entry_source,
      final JTree tree, final Splash splash_main)
  {
    try
    {
      TreePath path = tree.getLeadSelectionPath();
      if(path == null)
        return;
      DatabaseTreeNode seq_node = 
        (DatabaseTreeNode)path.getLastPathComponent();
      String node_name = (String)seq_node.getUserObject();
      String schema = seq_node.getSchema();
      
      String id = seq_node.getFeatureId(); 
      if(id != null)
        getEntryEditFromDatabase(id, entry_source, tree,
                                 splash_main, node_name,
                                 schema);
    }
    catch(NullPointerException npe)
    {
    }
  }

  /**
   * Retrieve a database entry.
   * @param id
   * @param entry_source
   * @param tree
   * @param splash_main
   * @param node_name
   * @param schema
   */
  private void getEntryEditFromDatabase(final String id,
      final DatabaseEntrySource entry_source, final JTree tree,
      final Splash splash_main, final String node_name,
      final String schema)
  {
    SwingWorker entryWorker = new SwingWorker()
    {
      public Object construct()
      {
        Cursor cbusy = new Cursor(Cursor.WAIT_CURSOR);
        Cursor cdone = new Cursor(Cursor.DEFAULT_CURSOR);

        status_line.setText("Retrieving sequence....");
        tree.setCursor(cbusy);
        try
        {
          entry_source.setSplitGFF(splitGFFEntry);

          final Entry entry = entry_source.getEntry(id, schema,
              stream_progress_listener);


          DatabaseDocumentEntry db_entry =
                         (DatabaseDocumentEntry)entry.getEMBLEntry();
          ((DatabaseDocument)db_entry.getDocument()).setName(node_name);

          if(entry == null)
          {
            tree.setCursor(cdone);
            status_line.setText("No entry.");
            return null;
          }

          final EntryEdit new_entry_edit = ArtemisMain.makeEntryEdit(entry);

          // add gff entries
          if(splitGFFEntry)
          {
//          DatabaseDocumentEntry db_entry = 
//                       (DatabaseDocumentEntry)entry.getEMBLEntry();

            final DatabaseDocumentEntry[] entries = entry_source.makeFromGff(
                        (DatabaseDocument)db_entry.getDocument(), id, schema);
            for(int i = 0; i < entries.length; i++)
            {
              if(entries[i] == null)
                continue;

              final Entry new_entry = new Entry(new_entry_edit.getEntryGroup().getBases(),
                                                entries[i]);
              new_entry_edit.getEntryGroup().add(new_entry);
            }
          }

          new_entry_edit.setVisible(true);
          status_line.setText("Sequence loaded.");
        }
        catch(OutOfRangeException e)
        {
          new MessageDialog(splash_main, "read failed: one of the features in "
              + " the entry has an out of range " + "location: "
              + e.getMessage());
        }
        catch(NoSequenceException e)
        {
          new MessageDialog(splash_main, "read failed: entry contains no sequence");
        }
        catch(IOException e)
        {
          new MessageDialog(splash_main, "read failed due to IO error: " + e);
        }
        tree.setCursor(cdone);
        return null;
      }

    };
    entryWorker.start();

  }

  /**
   * Create a menu bar 
   * @param entry_source
   * @param splash_main
   * @return
   */
  public JMenuBar makeMenuBar(final DatabaseEntrySource entry_source,
                              final Splash splash_main)
  {
    JMenuBar mBar = new JMenuBar();
    JMenu fileMenu = new JMenu("File");
    mBar.add(fileMenu);

    JMenuItem fileShow = new JMenuItem("Open Selected Sequence ...");
    fileShow.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        showSelected(entry_source, tree, splash_main);
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

    JMenu optionMenu = new JMenu("Options");
    mBar.add(optionMenu);

    final JCheckBoxMenuItem splitGFF = new JCheckBoxMenuItem(
                    "Split GFF into entries", splitGFFEntry);
    splitGFF.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        if(splitGFF.isSelected())
          splitGFFEntry = true;
        else
          splitGFFEntry = false;
      }
    });
    optionMenu.add(splitGFF);

    return mBar;
  }
  
  /**
   * Create database organism JTree.
   */
  private JTree getDatabaseTree(final DatabaseEntrySource entry_source)
  {
    Hashtable entries = null;
    
    while(entries == null)
    {
      try
      {
        DatabaseDocument doc = entry_source.getDatabaseDocument();
        entries = doc.getDatabaseEntries();
        final DatabaseTreeNode top = new DatabaseTreeNode("");
        createNodes(top, doc.getSchema(), entries);
        return new DatabaseJTree(top);
      }
      catch(ConnectException e)
      {
        entry_source.setLocation(true);
      }
      catch(SQLException e)
      {
        entry_source.setLocation(true);
      }
      catch(RuntimeException p)
      {
        entry_source.setLocation(true);
      }
    }
    
    return null;
  }


  /**
   * Create the nodes of the organism JTree
   * @param top       root node
   * @param schema    <code>List</code>
   * @param organism  sequences collection
   */
  private void createNodes(DatabaseTreeNode top, List schema,
                           Hashtable entries)
  {
    DatabaseTreeNode schema_node;
    DatabaseTreeNode seq_node;
    DatabaseTreeNode typ_node;

    final Object v_organism[] = entries.keySet().toArray();
    final int v_organism_size = v_organism.length;
    Arrays.sort(v_organism);  
    
    for(int i=0; i<schema.size(); i++)
    {
      int nchild = 0;
      String name;
      
      if(schema.get(i) instanceof String)
        name = (String)schema.get(i);
      else
        name = ((Organism)schema.get(i)).getCommonName();
      
      schema_node = new DatabaseTreeNode(name);
      final Hashtable seq_type_node = new Hashtable();

      for(int j = 0; j < v_organism_size; j++)
      {
        String seq_name  = (String)v_organism[j];
        String featureId = (String)entries.get(seq_name);
        
        if(seq_name.startsWith(name))
        {
          int ind1 = seq_name.indexOf("- ");
          int ind2 = seq_name.lastIndexOf("- ");

          String type = seq_name.substring(ind1 + 2, ind2 - 1);
          seq_name = seq_name.substring(ind2 + 2);

          if(!seq_type_node.containsKey(type))
          {
            typ_node = new DatabaseTreeNode(type);
            seq_type_node.put(type, typ_node);
            schema_node.add(typ_node);
          }
          else
            typ_node = (DatabaseTreeNode) seq_type_node.get(type);

          seq_node = new DatabaseTreeNode(seq_name,
                                          featureId, name);
          typ_node.add(seq_node);
          nchild++;
        }
      }
      if(nchild > 0)
        top.add(schema_node);
    }
  }
  
  /**
   *  An InputStreamProgressListener used to update the error label with the
   *  current number of chars read.
   **/
  private final InputStreamProgressListener stream_progress_listener =
    new InputStreamProgressListener() 
  {
    public void progressMade(final InputStreamProgressEvent event) 
    {
      final int char_count = event.getCharCount();
      if(char_count == -1) 
        status_line.setText("");
      else 
        status_line.setText("chars read so far: " + char_count);
    }

    public void progressMade(String progress)
    {
      status_line.setText(progress);
    }
  };
}
