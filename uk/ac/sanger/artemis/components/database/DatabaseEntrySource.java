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

package uk.ac.sanger.artemis.components.database;

import javax.swing.*;
import java.awt.Container;
import java.awt.GridLayout;
import java.io.*;
import java.util.*;
import java.net.ConnectException;

import org.gmod.schema.organism.Organism;

import uk.ac.sanger.artemis.util.*;
import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.sequence.*;
import uk.ac.sanger.artemis.io.DatabaseDocumentEntry;
import uk.ac.sanger.artemis.io.InvalidKeyException;
import uk.ac.sanger.artemis.io.EntryInformationException;

/**
 * 
 * This is an EntrySource that reads Entry objects from a relational database.
 * 
 */
public class DatabaseEntrySource implements EntrySource, Serializable
{
  private String location;

  private JPasswordField pfield;

  private Hashtable entries;

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
   * jdbc:postgresql://host:port/database?user=username
   * 
   */
  public boolean setLocation(final boolean prompt_user)
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
   * Create database organism JTree.
   */
  protected JTree getDatabaseTree()
  {
    DatabaseDocument doc = null;

    while(entries == null)
    {
      try
      {
        doc = new DatabaseDocument(location, pfield);
        entries = doc.getDatabaseEntries();
      }
      catch(java.sql.SQLException sqle)
      {
        setLocation(true);
      }
      catch(ConnectException ce)
      {
        setLocation(true);
      }
    }
    
    DatabaseTreeNode top = new DatabaseTreeNode("", DatabaseEntrySource.this);
    createNodes(top, doc.getSchema(), entries);
    final DatabaseJTree tree = new DatabaseJTree(top);
 
    return tree;
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
      
      schema_node = new DatabaseTreeNode(name, DatabaseEntrySource.this);
      Hashtable seq_type_node = new Hashtable();

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
            typ_node = new DatabaseTreeNode(type, DatabaseEntrySource.this);
            seq_type_node.put(type, typ_node);
            schema_node.add(typ_node);
          }
          else
            typ_node = (DatabaseTreeNode) seq_type_node.get(type);

          seq_node = new DatabaseTreeNode(seq_name, 
                                          DatabaseEntrySource.this, 
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
   * Set whether to split into separate entries.
   * @param splitGFFEntry
   */
  protected void setSplitGFF(final boolean splitGFFEntry)
  {
    this.splitGFFEntry = splitGFFEntry;
  }

  /**
   * Given an database feature identifier, this makes an <code>Entry</code>.
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
        db_entry = new DatabaseDocumentEntry(doc, null);
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
        db_entry[i] = new DatabaseDocumentEntry(new_docs[i], null);
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

  public String getLocation()
  {
    return location;
  }

  public JPasswordField getPfield()
  {
    return pfield;
  }
}
