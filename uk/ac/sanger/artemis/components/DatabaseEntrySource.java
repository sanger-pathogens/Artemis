/* DatabaseEntrySource.java
 *
 * created: Mar 2005
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2001  Genome Research Limited
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
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

import javax.swing.*;
import javax.swing.tree.DefaultMutableTreeNode;
import javax.swing.tree.TreeSelectionModel;
import javax.swing.tree.TreePath;
import java.awt.Container;
import java.awt.GridLayout;
import java.net.*;
import java.io.*;
import java.util.*;

import uk.ac.sanger.artemis.util.*;
import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.sequence.*;
import uk.ac.sanger.artemis.io.DatabaseDocumentEntry;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.SimpleEntryInformation;
import uk.ac.sanger.artemis.io.InvalidKeyException;
import uk.ac.sanger.artemis.io.EntryInformationException;

/**
 * 
 * This is an EntrySource that reads Entry objects from a relational database.
 * 
 */

public class DatabaseEntrySource implements EntrySource
{
  private String location;

  private JPasswordField pfield;

  private Hashtable entries;

  private Hashtable schemas;

  private boolean splitGFFEntry;

  /**
   * Create a new DatabaseEntrySource.
   */
  public DatabaseEntrySource()
  {
  }

  /**
   * @param bases
   *          The Bases object to pass to the Entry constructor.
   * @param show_progress
   *          If true a InputStreamProgressDialog will be shown while loading.
   *          (Not implemented)
   * @exception OutOfRangeException
   *              Thrown if one of the features in the Entry is out of range of
   *              the Bases object.
   * @return null if and only if the read is cancelled by the user or if the
   *         read fails.
   */
  public Entry getEntry(final Bases bases, final boolean show_progress)
      throws OutOfRangeException, IOException
  {
    return null;
  }

  public Entry getEntry(final boolean show_progress)
      throws OutOfRangeException, IOException
  {
    return null;
  }

  /**
   * Get an Entry object from this source.
   * 
   * @param id
   *          Feature ID to read in.
   * @exception OutOfRangeException
   *              Thrown if one of the features in the Entry is out of range of
   *              the Bases object.
   * @exception NoSequenceException
   *              Thrown if the entry that we read has no sequence.
   * @return null if and only if the user cancels the read or if the read fails.
   */
  public Entry getEntry(String id, String schema,
      InputStreamProgressListener progress_listener)
      throws OutOfRangeException, NoSequenceException, IOException
  {
    return makeFromID(null, id, schema, false, progress_listener);
  }

  /**
   * Returns true if and only if this EntrySource always returns "full" entries.
   * ie. entries that contain features and sequence.
   */
  public boolean isFullEntrySource()
  {
    return true;
  }

  /**
   * Return the name of this source (for display to the user in menus).
   */
  public String getSourceName()
  {
    return "Database";
  }

  /**
   * 
   * Set the database location as:
   * jdbc:postgresql://localhost:13001/chadoCVS?user=es2
   * 
   */
  protected boolean setLocation(final boolean prompt_user)
  {
    Container bacross = new Container();
    bacross.setLayout(new GridLayout(6, 2, 5, 5));

    JLabel lServer = new JLabel("Server : ");
    bacross.add(lServer);
    JTextField inServer = new JTextField("localhost");
    bacross.add(inServer);

    JLabel lPort = new JLabel("Port : ");
    bacross.add(lPort);
    JTextField inPort = new JTextField("5432");
    bacross.add(inPort);

    JLabel lDB = new JLabel("Database : ");
    bacross.add(lDB);
    JTextField inDB = new JTextField("chado");
    bacross.add(inDB);

    JLabel lUser = new JLabel("User : ");
    bacross.add(lUser);
    JTextField inUser = new JTextField("afumigatus");
    bacross.add(inUser);

    JLabel lpasswd = new JLabel("Password : ");
    bacross.add(lpasswd);
    pfield = new JPasswordField(16);
    bacross.add(pfield);

    // given -Dchado=localhost:port/dbname?username
    if(System.getProperty("chado") != null)
    {
      String db_url = System.getProperty("chado").trim();
      int index;
      if((index = db_url.indexOf(":")) > -1)
      {
        inServer.setText(db_url.substring(0, index));
        int index2;
        if((index2 = db_url.indexOf("/")) > -1)
        {
          inPort.setText(db_url.substring(index + 1, index2));
          int index3;
          if((index3 = db_url.indexOf("?")) > -1)
          {
            inDB.setText(db_url.substring(index2 + 1, index3));
            inUser.setText(db_url.substring(index3 + 1));

            /*
             * if(!prompt_user) { location = "jdbc:postgresql://"
             * +inServer.getText().trim()+ ":" +inPort.getText().trim()+ "/"
             * +inDB.getText().trim()+ "?user=" +inUser.getText().trim(); return
             * true; }
             */
          }
        }
      }
    }

    int n = JOptionPane.showConfirmDialog(null, bacross,
                                          "Enter Database Address",
                                          JOptionPane.OK_CANCEL_OPTION,
                                          JOptionPane.QUESTION_MESSAGE);
    if(n == JOptionPane.CANCEL_OPTION)
      return false;

    location = "jdbc:postgresql://" + 
               inServer.getText().trim() + ":" +
               inPort.getText().trim() + "/" +
               inDB.getText().trim() + "?user=" +
               inUser.getText().trim();

    return true;
  }

  /**
   * 
   * Given the entry name get the seq feature_id
   * 
   */
  protected String getEntryID(String name)
  {
    return (String)entries.get(name);
  }

  /**
   * 
   * Create database organism JTree.
   * 
   */
  protected JTree getDatabaseTree()
  {
    DatabaseDocument doc = new DatabaseDocument(location, pfield);

    entries = doc.getDatabaseEntries();
    Vector organism = doc.getOrganism();
    schemas = doc.getSchemaEntries();

    DefaultMutableTreeNode top = new DefaultMutableTreeNode("PSU Organism List");
    createNodes(top, doc.getSchema(), entries);
    final JTree tree = new JTree(top);
    tree.getSelectionModel().setSelectionMode(
        TreeSelectionModel.SINGLE_TREE_SELECTION);

    return tree;
  }

  /**
   * 
   * Get DefaultMutableTreeNode of selected node
   * 
   * @return node that is currently selected
   * 
   */
  protected String getSelectedNode(JTree tree)
  {
    TreePath path = tree.getLeadSelectionPath();
    if (path == null)
      return null;

    DefaultMutableTreeNode seq_node = 
                            (DefaultMutableTreeNode)path.getLastPathComponent();
    DefaultMutableTreeNode type_node =
                            (DefaultMutableTreeNode)seq_node.getParent();
    DefaultMutableTreeNode org_node = 
                            (DefaultMutableTreeNode)type_node.getParent();

    return (String)org_node.getUserObject() + " - " +
           (String)type_node.getUserObject() + " - " +
           (String)seq_node.getUserObject();
  }

  /**
   * 
   * Get Organism of selected node
   * 
   * @return name of Organism that is top level of selected node
   * 
   */
  protected String getSchemaOfSelectedNode(JTree tree)
  {
    TreePath path = tree.getLeadSelectionPath();
    if(path == null)
      return null;

    DefaultMutableTreeNode seq_node = 
                            (DefaultMutableTreeNode)path.getLastPathComponent();
    DefaultMutableTreeNode type_node = 
                            (DefaultMutableTreeNode)seq_node.getParent();
    DefaultMutableTreeNode org_node = 
                            (DefaultMutableTreeNode)type_node.getParent();

    String org = (String)org_node.getUserObject();

    System.out.println(org);
    return org;
//  if(schemas.containsKey(org))
//    return (String) schemas.get(org);
//  else
//    return null;
  }

  /**
   * 
   * Create the nodes of the organism JTree
   * 
   * @param top
   *          root node
   * @param org
   *          organism collection
   * @param organism
   *          sequences collection
   * 
   */
  private void createNodes(DefaultMutableTreeNode top, List org,
                           Hashtable organism)
  {
    DefaultMutableTreeNode org_node;
    DefaultMutableTreeNode seq_node;
    DefaultMutableTreeNode typ_node;

    final Object v_organism[] = organism.keySet().toArray();
    final int v_organism_size = v_organism.length;
    Arrays.sort(v_organism);

    for(int i=0; i<org.size(); i++)
    {
      String name = (String)org.get(i);
      org_node = new DefaultMutableTreeNode(name);
      top.add(org_node);

      Hashtable seq_type_node = new Hashtable();

      for(int j = 0; j < v_organism_size; j++)
      {
        String seq_name = (String)v_organism[j];
        if(seq_name.startsWith(name))
        {
          int ind1 = seq_name.indexOf("- ");
          int ind2 = seq_name.lastIndexOf("- ");

          String type = seq_name.substring(ind1 + 2, ind2 - 1);
          seq_name = seq_name.substring(ind2 + 2);

          if(!seq_type_node.containsKey(type))
          {
            typ_node = new DefaultMutableTreeNode(type);
            seq_type_node.put(type, typ_node);
            org_node.add(typ_node);
          }
          else
            typ_node = (DefaultMutableTreeNode) seq_type_node.get(type);

          seq_node = new DefaultMutableTreeNode(seq_name);
          typ_node.add(seq_node);
        }
      }
    }
  }

  protected void setSplitGFF(final boolean splitGFFEntry)
  {
    this.splitGFFEntry = splitGFFEntry;
  }

  /**
   * 
   * Given an accession number and the handle of an EMBL corba server, this
   * method will ask the user (using a TextRequester) for the id of a entry in
   * the server and will then attempt to get it.
   * 
   * @param bases
   *          If this is null a new Bases object will be created for the Entry
   *          once it has been read from the server. If not null then it will be
   *          passed to the Entry constructor.
   * @param id
   *          The id of the entry in the database
   * @param schema
   *          The schema of the entry in the database
   * @param read_only
   *          true if and only if a read-only Entry should be created (some are
   *          always read only).
   * @exception OutOfRangeException
   *              Thrown if one of the features in the Entry is out of range of
   *              the Bases object.
   */
  protected Entry makeFromID(final Bases bases, final String id,
                             final String schema, final boolean read_only,
                             InputStreamProgressListener progress_listener)
             throws OutOfRangeException, IOException
  {
    try
    {
      DatabaseDocumentEntry db_entry = null;

      if(read_only)
      {
      }
      else
      {
        DatabaseDocument doc = new DatabaseDocument(location, pfield, id,
                                                    schema, splitGFFEntry, progress_listener);
        db_entry = new DatabaseDocumentEntry(doc);
      }

      final Bases real_bases;

      if(bases == null)
      {
        if (db_entry.getSequence() == null)
        {
          JOptionPane.showMessageDialog(null,
              "The selected entry contains no sequence: " + id, "No Sequence",
              JOptionPane.ERROR_MESSAGE);

          return null;
        }

        real_bases = new Bases(db_entry.getSequence());
      }
      else
        real_bases = bases;

      return new Entry(real_bases, db_entry);
    }
    catch(InvalidKeyException e)
    {
      JOptionPane.showMessageDialog(null, 
                                    "Unexpected error while accessing " +
                                    id + ": " + e,
                                    "Invalid Key", JOptionPane.ERROR_MESSAGE);
    }
    catch(EntryInformationException e)
    {
      JOptionPane.showMessageDialog(null, 
                                    "Failed to get entry: " + e,
                                    "Entry Information Exception",
                                    JOptionPane.ERROR_MESSAGE);
    }

    return null;
  }

  protected DatabaseDocumentEntry[] makeFromGff(final DatabaseDocument doc,
                                                String id, String schema)
       throws OutOfRangeException, IOException
  {
    DatabaseDocumentEntry[] db_entry = null;

    try
    {
      DatabaseDocument[] new_docs = doc.getGffDocuments(location, id, schema);
      db_entry = new DatabaseDocumentEntry[new_docs.length];

      for(int i = 0; i < new_docs.length; i++)
        db_entry[i] = new DatabaseDocumentEntry(new_docs[i]);
    }
    catch (EntryInformationException e)
    {
      JOptionPane.showMessageDialog(null, 
                                    "Failed to get entry: " + e,
                                    "Entry Information Exception",
                                    JOptionPane.ERROR_MESSAGE);
    }

    return db_entry;
  }

}
