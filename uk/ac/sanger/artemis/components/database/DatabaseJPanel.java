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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/database/DatabaseJPanel.java,v 1.21 2008-12-10 16:43:38 tjc Exp $
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
import uk.ac.sanger.artemis.io.Range;

import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JTree;
import javax.swing.JScrollPane;
import javax.swing.JPanel;
import javax.swing.JLabel;
import javax.swing.BorderFactory;
import javax.swing.border.Border;
import javax.swing.tree.TreePath;

import org.gmod.schema.organism.Organism;
import org.gmod.schema.organism.OrganismProp;
import org.gmod.schema.sequence.Feature;

import java.awt.BorderLayout;
import java.awt.event.MouseListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.Dimension;
import java.awt.Frame;
import java.awt.Toolkit;
import java.awt.Cursor;
import java.awt.FontMetrics;
import java.io.*;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

public class DatabaseJPanel extends JPanel
{
  /** */
  private static final long serialVersionUID = 1L;
  private JLabel status_line = new JLabel("");
  private boolean splitGFFEntry = false;
  private JTree tree;
  private DatabaseDocument doc;
  private static org.apache.log4j.Logger logger4j = 
    org.apache.log4j.Logger.getLogger(DatabaseJPanel.class);
  
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
          showSelected(entry_source, splash_main);
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
  public void showSelected(final DatabaseEntrySource entry_source,
                           final Splash splash_main)
  {
    try
    {
      TreePath path = tree.getLeadSelectionPath();
      if(path == null)
        return;
      DatabaseTreeNode seq_node = 
        (DatabaseTreeNode)path.getLastPathComponent();
      String node_name = (String)seq_node.getUserObject();
      String userName = doc.getUserName();

      final String id = seq_node.getFeatureId();
      if(id != null)
      {
    	getEntryEditFromDatabase(id, entry_source, tree, 
            status_line, stream_progress_listener, 
            splitGFFEntry, splash_main, 
            node_name, userName);
    	
    	Frame[] frames = JFrame.getFrames();
    	Splash splash = null;
        for(int i=0;i<frames.length;i++)
        {
          if(frames[i] instanceof Splash)
          {
        	splash = (Splash)frames[i];
        	break;
          }
        }
        
        if(splash != null)
        {
          final Iterator it = seq_node.getOrganism().getOrganismProps().iterator();
          while (it.hasNext()) 
          {
	        OrganismProp organismProp = (OrganismProp) it.next();
     	    if (organismProp.getCvTerm().getName().equals("translationTable"))
		   	  splash.setTranslationTable(organismProp.getValue());
		  }
        }
      }
    }
    catch(NullPointerException npe)
    {
    	npe.printStackTrace();
    }
  }
  
  /**
   * Open an Artemis EntryEdit window
   * @param entry_source
   * @param splash_main
   * @param stream_progress_listener
   * @param entryName                 e.g. Pfalciparum:Pf3D7_09 or 
   *                                  Pfalciparum:Pf3D7_09:20050..40000
   * @return
   */
  public static EntryEdit show(final DatabaseEntrySource entry_source,
                          final Splash splash_main,
                          final InputStreamProgressListener stream_progress_listener,
                          final String entryName)
  {
    final String entry[] = entryName.split(":");
    String url = (String)entry_source.getLocation();
    int index  = url.indexOf("?");
    
    String userName = url.substring(index+1).trim();
    if(userName.startsWith("user="))
      userName = userName.substring(5);
    
    DatabaseDocument doc = entry_source.getDatabaseDocument();
    
    Range range = null;
    if(entry.length>2)
    {
      if(entry[2].indexOf("..") > -1)
      {
        String ranges[] = entry[2].split("\\.\\.");
        if(ranges.length == 2)
        {
          try
          {
            range = new Range(Integer.parseInt(ranges[0]), 
                              Integer.parseInt(ranges[1]));
          }
          catch(Exception e){ e.printStackTrace(); }
        }
      }
    }
    
    Feature f = doc.getFeatureByUniquename(entry[1]);
    
    return openEntry(Integer.toString(f.getFeatureId()), entry_source, 
        splash_main.getCanvas(), splash_main.getStatusLabel(), stream_progress_listener,
        false, splash_main,  entry[1], userName, range);
  }

  /**
   * Retrieve a database entry.
   * @param srcfeatureId
   * @param entry_source
   * @param srcComponent
   * @param status_line
   * @param stream_progress_listener
   * @param splitGFFEntry
   * @param splash_main
   * @param dbDocumentName
   * @param userName
   */
  private static void getEntryEditFromDatabase(
      final String srcfeatureId,
      final DatabaseEntrySource entry_source, 
      final JComponent srcComponent,
      final JLabel status_line,
      final InputStreamProgressListener stream_progress_listener,
      final boolean splitGFFEntry,
      final Splash splash_main, 
      final String dbDocumentName,
      final String userName)
  {
    SwingWorker entryWorker = new SwingWorker()
    {
      public Object construct()
      {
        Cursor cbusy = new Cursor(Cursor.WAIT_CURSOR);
        Cursor cdone = new Cursor(Cursor.DEFAULT_CURSOR);

        status_line.setText("Retrieving sequence....");
        srcComponent.setCursor(cbusy);
        try
        {
          while(DatabaseDocument.isCvThreadAlive())
            Thread.sleep(5);
          openEntry(srcfeatureId, entry_source, srcComponent, status_line, stream_progress_listener,
              splitGFFEntry, splash_main, dbDocumentName, userName, null);
        }
        catch(RuntimeException re)
        {
          logger4j.warn(re.getMessage());

          re.printStackTrace();
          final DatabaseEntrySource entry_source = new DatabaseEntrySource();
          entry_source.setLocation(true);

          String url = entry_source.getLocation();
          int index  = url.indexOf("?");
            
          String userName = url.substring(index+1).trim();
          if(userName.startsWith("user="))
            userName = userName.substring(5);
            
          openEntry(srcfeatureId, entry_source, srcComponent, status_line, stream_progress_listener,
              splitGFFEntry, splash_main, dbDocumentName, userName, null);
        }
        catch(InterruptedException e)
        {
          e.printStackTrace();
        }
        srcComponent.setCursor(cdone);
        return null;
      }

    };
    entryWorker.start();

  }

  /**
   * 
   * @param srcfeatureId
   * @param entry_source
   * @param srcComponent
   * @param status_line
   * @param stream_progress_listener
   * @param splitGFFEntry
   * @param splash_main
   * @param dbDocumentName
   * @param userName
   * @param range           range for to retrieve features in
   * @return
   */
  private static EntryEdit openEntry(
      final String srcfeatureId,
      final DatabaseEntrySource entry_source, 
      final JComponent srcComponent,
      final JLabel status_line,
      final InputStreamProgressListener stream_progress_listener,
      final boolean splitGFFEntry,
      final Splash splash_main, 
      final String dbDocumentName,
      final String userName,
      final Range range) 
  {
    Cursor cdone = new Cursor(Cursor.DEFAULT_CURSOR);
    try
    {
      if(range != null)
        logger4j.info("LOAD FEATURES IN THE RANGE : "+range.toString());
      entry_source.setSplitGFF(splitGFFEntry);

      final Entry entry = entry_source.getEntry(srcfeatureId, userName,
          stream_progress_listener, range);

      DatabaseDocumentEntry db_entry = (DatabaseDocumentEntry) entry
          .getEMBLEntry();
      DatabaseDocument doc = (DatabaseDocument) db_entry.getDocument();
      doc.setName(dbDocumentName);

      if(entry == null)
      {
        srcComponent.setCursor(cdone);
        status_line.setText("No entry.");
        return null;
      }

      final EntryEdit new_entry_edit = ArtemisMain.makeEntryEdit(entry);

      // add gff entries
      if(splitGFFEntry)
      {
        final DatabaseDocumentEntry[] entries = entry_source.makeFromGff(
            (DatabaseDocument) db_entry.getDocument(), srcfeatureId, userName);

        for(int i = 0; i < entries.length; i++)
        {
          if(entries[i] == null)
            continue;

          final Entry new_entry = new Entry(new_entry_edit.getEntryGroup()
              .getBases(), entries[i]);
          new_entry_edit.getEntryGroup().add(new_entry);
        }
      }

      new_entry_edit.setVisible(true);
      status_line.setText("Sequence loaded.");
      return new_entry_edit;
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
    return null;
  }
  
  /**
   * Create database organism JTree.
   */
  private JTree getDatabaseTree(final DatabaseEntrySource entry_source)
  {
    HashMap entries = null;
    
    while(entries == null)
    {
      try
      {
        doc = entry_source.getDatabaseDocument();
        final DatabaseTreeNode top = new DatabaseTreeNode("");
        
        final List organisms = doc.getOrganismsContainingSrcFeatures();
        for(int i=0; i<organisms.size(); i++)
        {
          Organism org = (Organism)organisms.get(i);
          
          String name = org.getCommonName();
          if(name == null || name.equals(""))
            name = org.getGenus() + "." + org.getSpecies();
          DatabaseTreeNode orgNode = 
            new DatabaseTreeNode(name, false, org, doc.getUserName(), doc);
          top.add(orgNode); 
        }
        return new DatabaseJTree(top);
      }
      catch(Exception e)
      {
        if(!entry_source.setLocation(true))
          return null;
      }
    }

    return null;
  }
  
  /**
   * An InputStreamProgressListener used to update the error label with the
   * current number of chars read.
   */
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

  public void setSplitGFFEntry(boolean splitGFFEntry)
  {
    this.splitGFFEntry = splitGFFEntry;
  }
}
