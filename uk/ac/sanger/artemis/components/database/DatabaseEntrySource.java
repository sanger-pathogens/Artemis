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

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.IOException;
import java.io.Serializable;
import javax.swing.JCheckBox;

import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPasswordField;
import javax.swing.JTextField;

import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.EntrySource;
import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.sequence.Bases;
import uk.ac.sanger.artemis.sequence.NoSequenceException;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.DatabaseLocationParser;
import uk.ac.sanger.artemis.util.InputStreamProgressListener;
import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.util.StringVector;
import uk.ac.sanger.artemis.components.filetree.LocalAndRemoteFileManager;
import uk.ac.sanger.artemis.components.genebuilder.JExtendedComboBox;
import uk.ac.sanger.artemis.io.UI;
import uk.ac.sanger.artemis.io.DatabaseDocumentEntry;
import uk.ac.sanger.artemis.io.InvalidKeyException;
import uk.ac.sanger.artemis.io.EntryInformationException;
import uk.ac.sanger.artemis.io.Range;

/**
 * 
 * This is an EntrySource that reads Entry objects from a relational database.
 * 
 */
public class DatabaseEntrySource implements EntrySource, Serializable
{
  /** */
  private static final long serialVersionUID = 1L;

  private String location;

  private JPasswordField pfield;

  private boolean splitGFFEntry;
  
  private boolean readOnly = false;

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
    return makeFromID(null, id, schema, progress_listener, null);
  }
  
  public Entry getEntry(String id, String schema,
      InputStreamProgressListener progress_listener, Range range)
      throws OutOfRangeException, NoSequenceException, IOException
  {
    return makeFromID(null, id, schema, progress_listener, range);
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
  public boolean setLocation(final boolean promptUser)
  {
    ILoginPrompt promptPanel; // = new DatabaseLoginPrompt();
    
    if (UI.mode == UI.UIMode.SWING)
	{
		promptPanel = new DatabaseLoginPrompt();
	}
    else
	{
    	promptPanel = new DatabaseLoginPromptConsole();
	}
	
    
    if(promptUser)
    {
    	
    	boolean userEntered = promptPanel.prompt();
    	if (userEntered == false)
    	{
    		return false;
    	}
    	
//    	
//      int n = JOptionPane.showConfirmDialog(null, promptPanel,
//                                          "Enter Database Address",
//                                          JOptionPane.OK_CANCEL_OPTION,
//                                          JOptionPane.QUESTION_MESSAGE);
//      if(n == JOptionPane.CANCEL_OPTION)
//        return false;
    }
    
    pfield = promptPanel.getPfield();
    DatabaseLocationParser dlp = new DatabaseLocationParser();
    dlp.setHost(promptPanel.getServer());
    dlp.setPort(promptPanel.getPort());
    dlp.setDatabase(promptPanel.getDb());
    dlp.setUsername(promptPanel.getUser());
    dlp.setSSL(promptPanel.getSSL());
    
    location = dlp.getCompleteURL();
    return true;
  }
  
  

  public DatabaseDocument getDatabaseDocument()
  {
//	ArtemisConsole.info(location, "location");
//	ArtemisConsole.info(pfield.getText(), "password");
    return new DatabaseDocument(location, pfield);
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
   * @exception OutOfRangeException
   *              Thrown if one of the features in the Entry is out of range of
   *              the Bases object.
   */
  private Entry makeFromID(final Bases bases, final String id,
                             final String schema, 
                             InputStreamProgressListener progress_listener, 
                             final Range range)
             throws OutOfRangeException, IOException
  {
    try
    {
      DatabaseDocumentEntry db_entry = null;

      DatabaseDocument doc = new DatabaseDocument(location, pfield, id,
                                                  schema, splitGFFEntry, 
                                                  progress_listener);
      doc.setLazyFeatureLoad(LocalAndRemoteFileManager.lazyLoad.isSelected());
      if(range != null)
        doc.setRange(range);
      db_entry = new DatabaseDocumentEntry(doc, null);
        
      db_entry.setReadOnly(isReadOnly());
      
      final Bases real_bases;

      if(bases == null)
      {
        if (db_entry.getSequence() == null)
        {
          UI.error("The selected entry contains no sequence: " + id, "No Sequence");
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
      UI.error("Unexpected error while accessing " + id + ": " + e, "Invalid Key");
    }
    catch(EntryInformationException e)
    {
      UI.error("Failed to get entry: " + e,"Entry Information Exception");
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
      {
        new_docs[i].setLazyFeatureLoad(doc.isLazyFeatureLoad());
        db_entry[i] = new DatabaseDocumentEntry(new_docs[i], null);
      }
    }
    catch (EntryInformationException e)
    {
//      JOptionPane.showMessageDialog(null, 
//                                    "Failed to get entry: " + e,
//                                    "Entry Information Exception",
//                                    JOptionPane.ERROR_MESSAGE);
      UI.error("Failed to get entry: " + e,"Entry Information Exception");
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

  public boolean isReadOnly()
  {
    return readOnly;
  }

  public void setReadOnly(boolean readOnly)
  {
    this.readOnly = readOnly;
  }
 
}

/*
 * A common interface for login prompting, irrespective of the implementation. 
 */
interface ILoginPrompt
{
	boolean prompt();
	JPasswordField getPfield();
	String getServer();
    String getPort();
    String getDb();
    String getUser();
    boolean getSSL();
}

/**
 * Database login panel
 */
class DatabaseLoginPrompt extends JPanel implements ILoginPrompt
{
  private static final long serialVersionUID = 1L;
  private JPasswordField pfield;
  private JTextField server;
  private JTextField port;
  private JTextField db;
  private JTextField user;
  private JCheckBox ssl;
  public DatabaseLoginPrompt()
  {
    super(new GridBagLayout());

    final GridBagConstraints c = new GridBagConstraints();
    
    // column 1
    int nrow = 1;
    c.anchor = GridBagConstraints.EAST;
    c.gridx = 0;
    c.gridy = nrow;
    add(new JLabel("Server : "), c);

    c.gridy = ++nrow;
    add(new JLabel("Port : "), c);

    c.gridy = ++nrow;
    add(new JLabel("Database : "), c);

    c.gridy = ++nrow;
    add(new JLabel("User : "), c);

    c.gridy = ++nrow;
    add(new JLabel("Password : "), c);

    c.gridy = ++nrow;
    add(new JLabel("SSL : "), c);
    
    // column 2
    nrow = 1;
    c.anchor = GridBagConstraints.WEST;
    c.fill = GridBagConstraints.HORIZONTAL;
    c.gridx = 1;
    c.gridy = nrow;
    server = new JTextField("localhost");
    add(server, c);
    
    c.gridy = ++nrow;
    port = new JTextField("5432");
    add(port, c);

    c.gridy = ++nrow;
    db = new JTextField("chado");
    add(db, c);

    c.gridy = ++nrow;
    user = new JTextField("");
    add(user, c);
    
    c.gridy = ++nrow;
    pfield = new JPasswordField(16);
    add(pfield, c);
    
    c.gridy = ++nrow;
    ssl = new JCheckBox();
    add(ssl, c);

    // given -Dchado=localhost:port/dbname?username
    if(System.getProperty("chado") != null)
      setFromURL(System.getProperty("chado").trim());
    
    if (System.getProperty("chadoPassword") != null)
    {
    	//System.out.println("detected chadoPassword: " + System.getProperty("chadoPassword"));
    	pfield.setText(System.getProperty("chadoPassword"));
    }
    
    c.gridx = 1;
    c.gridy = ++nrow;
    add(getServerComboBox(), c);
    
    c.gridx = 0;
    c.anchor = GridBagConstraints.EAST;
    add(new JLabel("Available databases : "), c);
  }
  
  public boolean prompt()
  {
	  int n = JOptionPane.showConfirmDialog(null, this,
              "Enter Database Address",
              JOptionPane.OK_CANCEL_OPTION,
              JOptionPane.QUESTION_MESSAGE);
	  
	  if(n == JOptionPane.CANCEL_OPTION)
		  return false;
	return true;
  }
  
  private void setFromURL(String db_url)
  {
      DatabaseLocationParser dlp = new DatabaseLocationParser(db_url);
      server.setText(dlp.getHost());
      port.setText("" + dlp.getPort());
      db.setText(dlp.getDatabase());
      user.setText(dlp.getUsername());
      ssl.setSelected(dlp.isSSLEnabled());
  }
  
  /**
   * Get a combo box of the available database servers.
   * @return
   */
  private JExtendedComboBox getServerComboBox()
  {
    // known servers
    StringVector serversVector = Options.getOptions().getOptionValues("chado_servers");
    if(serversVector == null)
    {
      serversVector = new StringVector();
      serversVector.add("GeneDB (read-only)");
      serversVector.add("db.genedb.org:5432/snapshot?genedb_ro");
    }
    
    
    // set the default to the chado property value
    if(System.getProperty("chado") != null &&
       !serversVector.contains(System.getProperty("chado").trim()))
    {
      serversVector.add(0, "Default");
      serversVector.add(1, System.getProperty("chado").trim());
    }
    
    final String servers[][] = new String[2][serversVector.size()/2];
    
    int nserver = 0;
    for(int i=0; i<serversVector.size(); i+=2)
    {
      servers[0][nserver] = (String) serversVector.get(i);
      servers[1][nserver] = (String) serversVector.get(i+1);
      nserver++;
    }
    
    final JExtendedComboBox serverSelection = new JExtendedComboBox(
        servers[0], false);
    
    serverSelection.addItemListener(new ItemListener()
    {
      public void itemStateChanged(ItemEvent arg0)
      {
        String selectedServerUrl = null;
        for(int i=0; i<servers[0].length; i++)
        {
          if(serverSelection.getSelectedItem().equals(servers[0][i]))
            selectedServerUrl = servers[1][i];
        }
        if(selectedServerUrl != null)
          setFromURL(selectedServerUrl);
      } 
    });
    return serverSelection;
  }
  
  public JPasswordField getPfield()
  {
    return pfield;
  }

  public String getDb()
  {
    return db.getText().trim();
  }

  public String getPort()
  {
    return port.getText().trim();
  }

  public String getServer()
  {
    return server.getText().trim();
  }

  public String getUser()
  {
    return user.getText().trim();
  }
  
  public boolean getSSL()
  {
    return ssl.isSelected();
  }
}

class DatabaseLoginPromptConsole implements ILoginPrompt
{
	private String password;
	private String server;
	private String port;
	private String db;
	private String user;
        private String ssl;

	public DatabaseLoginPromptConsole()
	{
		
		if(System.getProperty("chado") != null)
		{
			//System.out.println("detected chado : " + System.getProperty("chado"));
			setFromURL(System.getProperty("chado").trim());
		}
			
		
		if (System.getProperty("chadoPassword") != null)
	    {
	    	//System.out.println("detected chadoPassword: " + System.getProperty("chadoPassword"));
	    	password = System.getProperty("chadoPassword");
	    }
	}
	
	private void setFromURL(String db_url)
	  {
              DatabaseLocationParser dlp = new DatabaseLocationParser(db_url);
              server = dlp.getHost();
              port = "" + dlp.getPort();
              db = dlp.getDatabase();
              user = dlp.getUsername();
              ssl = dlp.isSSLEnabled() ? "y" : "";
	  }
	
	public JPasswordField getPfield()
	{
		JPasswordField pfield = new JPasswordField(16);
		pfield.setText(password);
		return pfield;
	}
	
	public String getDb()
	{
		return db;
	}
	
	public String getPort()
	{
		return port;
	}
	
	public String getServer()
	{
		return server;
	}
	
	public String getUser()
	{
		return user;
	}
        
        public boolean getSSL()
        {
                return ssl.toLowerCase().equals("y") || ssl.toLowerCase().equals("yes");
        }

	public boolean prompt() 
	{
		boolean userEntered = false;
		if (UI.mode == UI.UIMode.SCRIPT)
		{
			return false;
		}
        
		if (db == null)
		{
			db = UI.userInput("Enter Database", false);
			userEntered = true;
		}
		if (port == null)
		{
			port = UI.userInput("Enter Port", false);
			userEntered = true;
		}
		if (server == null)
		{
			server = UI.userInput("Enter Server",false);
			userEntered = true;
		}
		if (user == null)
		{
			user = UI.userInput("Enter User", false);
			userEntered = true;
		}
		if (password == null)
		{
			password = UI.userInput("Enter Password", true);
			userEntered = true;
		}
                if (ssl == null)
		{
			ssl = UI.userInput("Use SSL [y/n]", false);
			userEntered = true;
		}
		
		/**
		    GSV found that that this should always return true
		    or else the connection doesn't use some of the
		    parameters passed in on the command line and fails
		    @TODO must investigate why...
		*/
		return true;
	}
}


