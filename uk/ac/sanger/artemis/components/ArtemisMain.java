/* ArtemisMain.java
 *
 * created: Wed Feb 23 2000
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/ArtemisMain.java,v 1.33 2008-12-10 16:43:38 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.components.database.DatabaseEntrySource;
import uk.ac.sanger.artemis.components.database.DatabaseJPanel;
import uk.ac.sanger.artemis.components.filetree.FileManager;
import uk.ac.sanger.artemis.components.filetree.LocalAndRemoteFileManager;
import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.sequence.NoSequenceException;
import uk.ac.sanger.artemis.sequence.Bases;
import uk.ac.sanger.artemis.util.Document;
import uk.ac.sanger.artemis.util.TextDocument;
import uk.ac.sanger.artemis.util.DocumentFactory;
import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.util.InputStreamProgressListener;
import uk.ac.sanger.artemis.io.EntryInformation;

import org.biojava.bio.seq.io.SequenceFormat;

import java.awt.event.*;
import java.awt.Toolkit;
import java.io.*;
import java.util.Vector;
import java.awt.datatransfer.*;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.SwingUtilities;

/**
 *  The main window for the Artemis sequence editor.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: ArtemisMain.java,v 1.33 2008-12-10 16:43:38 tjc Exp $
 **/

public class ArtemisMain extends Splash 
{
  /** */
  private static final long serialVersionUID = 1L;

  /** A vector containing all EntryEdit object we have created. */
  private Vector<EntryEdit> entry_edit_objects = new Vector<EntryEdit>();

  protected static FileManager filemanager = null;
  
  private LocalAndRemoteFileManager fm;
  
  /**
   *  The constructor creates all the components for the main Artemis 
   *  window and sets up all the menu callbacks.
   **/
  public ArtemisMain(final String args[]) 
  {
    super("Artemis", "Artemis", "15");

    makeMenuItem(file_menu, "Open Project Manager ...", new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        new ProjectProperty(ArtemisMain.this);
      }
    });
    
    ActionListener menu_listener = new ActionListener()
    {
      public void actionPerformed(ActionEvent event)
      {
        if(filemanager == null)
          filemanager = new FileManager(ArtemisMain.this);
        else
          filemanager.setVisible(true);
      }
    };
    makeMenuItem(file_menu, "Open File Manager ...", menu_listener);

    ActionListener menu_listener_ssh = new ActionListener()
    {
      public void actionPerformed(ActionEvent event)
      {
        if(fm == null)
          fm = new LocalAndRemoteFileManager(ArtemisMain.this);
        else
          fm.setVisible(true);
      }
    };
    
    if(System.getProperty("chado") != null)
      makeMenuItem(file_menu, "Open Database and SSH File Manager ...", menu_listener_ssh);
    else 
      makeMenuItem(file_menu, "Open SSH File Manager ...", menu_listener_ssh);


    final EntrySourceVector entry_sources = getEntrySources(this);

    for(int source_index=0; source_index<entry_sources.size();
                            ++source_index) 
    {
      final EntrySource this_entry_source =
        entry_sources.elementAt(source_index);

      String entry_source_name = this_entry_source.getSourceName();
      String menu_name = null;

      if(entry_source_name.equals("Filesystem")) 
        menu_name = "Open ...";
      else 
        menu_name = "Open from " + entry_source_name + " ...";

      menu_listener = new ActionListener() 
      {
        public void actionPerformed(ActionEvent event) 
        {
          getEntryEditFromEntrySource(this_entry_source);
        }
      };
      makeMenuItem(file_menu, menu_name, menu_listener);
    }

    menu_listener = new ActionListener()
    {
      public void actionPerformed(ActionEvent event)
      {
        launchDatabaseJFrame();
      }
    };

    //final boolean sanger_options =
    //  Options.getOptions().getPropertyTruthValue("sanger_options");

    //makeMenuItem(file_menu, "Database Entry ...", menu_listener);

    menu_listener = new ActionListener() 
    {
      public void actionPerformed(ActionEvent event)
      {
        exit();
      }
    };
    makeMenuItem(file_menu, "Quit", menu_listener);

    getCanvas().addMouseListener(new MouseAdapter() 
    {
      /**
       *  Listen for mouse press events so that we can do popup menus and
       *  selection.
       **/
      public void mousePressed(MouseEvent event)  
      {
        handleCanvasMousePress(event);
      }
    });
  }


  /**
   *  Handle a mouse press event on the drawing canvas - select on click,
   *  select and broadcast it on double click.
   **/
  private void handleCanvasMousePress(MouseEvent event) 
  {
    if(event.getID() != MouseEvent.MOUSE_PRESSED)
      return;

    if((event.getModifiers() & InputEvent.BUTTON2_MASK) != 0)
    {
      openClipboardContents();
    }
  }

  /**
  * Get the String residing on the clipboard.
  *
  * @return any text found on the Clipboard; if none found, return an
  * empty String.
  */
  public void openClipboardContents() 
  {
    Clipboard clipboard = Toolkit.getDefaultToolkit().getSystemClipboard();
    //odd: the Object param of getContents is not currently used
    Transferable contents = clipboard.getContents(null);
    boolean hasTransferableText = (contents != null) &&
                                  contents.isDataFlavorSupported(DataFlavor.stringFlavor);
    if(hasTransferableText)
    {
      TextDocument entry_document = new TextDocument();
      final InputStreamProgressListener progress_listener =
                                     getInputStreamProgressListener();

      entry_document.addInputStreamProgressListener(progress_listener);

      final EntryInformation artemis_entry_information =
                          Options.getArtemisEntryInformation();

      final uk.ac.sanger.artemis.io.Entry new_embl_entry =
          EntryFileDialog.getEntryFromFile(this, entry_document,
                                           artemis_entry_information,
                                           false);

      if(new_embl_entry == null)  // the read failed
        return;

      try 
      {
        final Entry entry = new Entry(new_embl_entry);
        EntryEdit last_entry_edit = makeEntryEdit(entry);
        addEntryEdit(last_entry_edit);
        getStatusLabel().setText("");
        last_entry_edit.setVisible(true);
      }
      catch(OutOfRangeException e) 
      {
        new MessageDialog(this, "read failed: one of the features in " +
                           " cut and paste has an out of range " +
                           "location: " + e.getMessage());
      } 
      catch(NoSequenceException e) 
      {
        new MessageDialog(this, "read failed: " +
                           " cut and paste contains no sequence");
      }
    }
  }


  /**
  *
  * Launch database manager window
  *
  */
  private void launchDatabaseJFrame()
  {
    SwingWorker entryWorker = new SwingWorker()
    {
      public Object construct()
      {
        getStatusLabel().setText("Connecting ...");
        DatabaseEntrySource entry_source = new DatabaseEntrySource();
        boolean promptUser = true;
        if(System.getProperty("read_only") != null)
          promptUser = false;
        
        if(!entry_source.setLocation(promptUser))
          return null;

        JFrame frame = new JFrame("Organism List");
        final DatabaseJPanel pane = new DatabaseJPanel(entry_source,
                                               ArtemisMain.this);
        frame.getContentPane().add(pane);
        frame.pack();
        Utilities.rightJustifyFrame(frame);
        frame.setVisible(true);
        //frame.setJMenuBar(pane.makeMenuBar(entry_source, ArtemisMain.this));
        getStatusLabel().setText("");
        return null;
      }
    };
    entryWorker.start();
  }

  /**
   *  Read the entries named in args and in the diana.ini file.
   **/
  protected void readArgsAndOptions(final String [] args, final JFrame f)
  {
    // JNLP properties
    if(System.getProperty("jnlp.chado") != null)
      System.setProperty("chado", System.getProperty("jnlp.chado"));
    if(System.getProperty("jnlp.offset") != null)
      System.setProperty("offset", System.getProperty("jnlp.offset"));
    if(System.getProperty("jnlp.artemis.environment") != null)
      System.setProperty("artemis.environment", System.getProperty("jnlp.artemis.environment"));
    if(System.getProperty("jnlp.sanger_options") != null)
      System.setProperty("sanger_options", System.getProperty("jnlp.sanger_options"));
    
    if(args.length == 0) 
    {
      if(System.getProperty("chado") != null && 
         (args.length < 1 || args[0].indexOf(':') == -1))
        fm = new LocalAndRemoteFileManager(ArtemisMain.this);
        
      // open the entries given in the options file(diana.ini)
      readDefaultEntries();
      return;
    }

    if(args[0].equals("-biojava")) 
    {
      handleBioJava(args);
      return;
    }

    final EntryInformation artemis_entry_information =
                          Options.getArtemisEntryInformation();

    EntryEdit last_entry_edit = null;
    boolean seen_plus = false;

    for(int i = 0 ; i<args.length ; ++i) 
    {
      String new_entry_name = args[i];

      if(new_entry_name.length() == 0) 
        continue;

      if(new_entry_name.equals("+")) 
      {
        seen_plus = true;
        continue;
      }

      if(new_entry_name.startsWith("+") && last_entry_edit != null ||
         seen_plus) 
      {
        // new feature file

        final Document entry_document;

        if(seen_plus) 
          entry_document = DocumentFactory.makeDocument(new_entry_name);
        else 
          entry_document =
            DocumentFactory.makeDocument(new_entry_name.substring(1));

        final InputStreamProgressListener progress_listener =
                                     getInputStreamProgressListener();

        entry_document.addInputStreamProgressListener(progress_listener);

        final uk.ac.sanger.artemis.io.Entry new_embl_entry =
          EntryFileDialog.getEntryFromFile(f, entry_document,
                                           artemis_entry_information,
                                           false);

        if(new_embl_entry == null)  // the read failed
          break;

        try
        {
          final Entry new_entry =
            new Entry(last_entry_edit.getEntryGroup().getBases(),
                      new_embl_entry);

          last_entry_edit.getEntryGroup().add(new_entry);
        } 
        catch(OutOfRangeException e) 
        {
          new MessageDialog(this, "read failed: one of the features in " +
                             new_entry_name + " has an out of range " +
                             "location: " + e.getMessage());
        }
        seen_plus = false; // reset
      }
      else if(System.getProperty("chado") != null && new_entry_name.indexOf(':')>-1)
      {
        // open from database e.g. Pfalciparum:Pf3D7_09:95000..150000
        Splash.logger4j.info("OPEN ENTRY "+new_entry_name);
        getStatusLabel().setText("Connecting ...");
        DatabaseEntrySource entry_source = new DatabaseEntrySource();

        boolean promptUser = true;
        if(System.getProperty("read_only") != null)
        {
          promptUser = false;
          entry_source.setReadOnly(true);
        }
        
        last_entry_edit = dbLogin(entry_source, promptUser, new_entry_name);
        if(last_entry_edit == null)
          return;
      }
      else
      {
        // new sequence file

        if(last_entry_edit != null) 
        {
          last_entry_edit.setVisible(true);
          last_entry_edit = null;
        }

      
        if (new_entry_name.indexOf ("://") == -1) 
        {
          File file = new File (new_entry_name);
          if(!file.exists())
          {
            JOptionPane.showMessageDialog(null, "File "+ 
                new_entry_name +" not found.\n"+
                "Check the file name.", "File Not Found", 
                JOptionPane.WARNING_MESSAGE);
          }
        }  

        final Document entry_document =
          DocumentFactory.makeDocument(new_entry_name);

        entry_document.addInputStreamProgressListener(getInputStreamProgressListener());
 
        final uk.ac.sanger.artemis.io.Entry new_embl_entry =
          EntryFileDialog.getEntryFromFile(f, entry_document,
                                           artemis_entry_information,
                                           false);

        if(new_embl_entry == null)  // the read failed
          break;

        try 
        {
          final Entry entry = new Entry(new_embl_entry);
          last_entry_edit = makeEntryEdit(entry);
          addEntryEdit(last_entry_edit);
          getStatusLabel().setText("");
        }
        catch(OutOfRangeException e) 
        {
          new MessageDialog(this, "read failed: one of the features in " +
                             new_entry_name + " has an out of range " +
                             "location: " + e.getMessage());
          break;
        } 
        catch(NoSequenceException e) 
        {
          new MessageDialog(this, "read failed: " +
                             new_entry_name + " contains no sequence");
          break;
        }
      }
    }

    if(System.getProperty("offset") != null)
      last_entry_edit.getGotoEventSource().gotoBase(
          Integer.parseInt(System.getProperty("offset")));
    
    last_entry_edit.setVisible(true);
    for(int entry_index=0; entry_index<entry_edit_objects.size();
        ++entry_index)
    {
      if(System.getProperty("offset") != null)
        entry_edit_objects.elementAt(entry_index).getGotoEventSource().gotoBase(
            Integer.parseInt(System.getProperty("offset")));
    }
  }
  
  /**
   * Handle database connection and construction of EntryEdit
   * @param entry_source
   * @param promptUser
   * @param new_entry_name
   * @return
   */
  private EntryEdit dbLogin(final DatabaseEntrySource entry_source, 
                            final boolean promptUser, 
                            final String new_entry_name)
  {
    EntryEdit last_entry_edit = null;
    
    // allow 3 attempts to login
    for(int i = 0; i < 3; i++)
    {
      if (!entry_source.setLocation(promptUser))
        return null;

      try
      {
        last_entry_edit = DatabaseJPanel.show(entry_source, this,
            getInputStreamProgressListener(), new_entry_name);
        break;
      }
      catch (Exception e)
      {
        new MessageDialog(this, e.getMessage());
        entry_source.getDatabaseDocument().reset();
      }
    }
    
    return last_entry_edit;
  }
  
  /**
   *
   *  Handle the -biojava option
   * 
   *  Command line syntax:  
   *  art -biojava org.biojava.bio.seq.io.EmblLikeFormat foo.embl
   *
   *  BioJava formats:
   *  EmblLikeFormat, FastaFormat, GAMEFormat, GenbankFormat, PhredFormat
   *
   **/
  private void handleBioJava(final String [] args) 
  {
    if(args.length == 3) 
    {
      final String class_name = args[1];
      final String location   = args[2];

      final Document location_document =
        DocumentFactory.makeDocument(location);

      try 
      {
        final Object biojava_object =
          Class.forName(class_name).newInstance();

        final EntryInformation entry_information =
          Options.getArtemisEntryInformation();

        final uk.ac.sanger.artemis.io.BioJavaEntry emblEntry;

        if(biojava_object instanceof SequenceFormat)
        {
          final SequenceFormat sequence_format = (SequenceFormat)biojava_object;

          emblEntry =
            new uk.ac.sanger.artemis.io.BioJavaEntry(entry_information,
                                                     location_document,
                                                     sequence_format);

          final Entry new_entry = new Entry(emblEntry);
          final EntryEdit new_entry_edit = makeEntryEdit(new_entry);
          new_entry_edit.setVisible(true);
        } 
        else 
          new MessageDialog(this, "not a SequenceFormat: " + class_name);
      } 
      catch(IllegalAccessException e) 
      {
        new MessageDialog(this, "cannot create class: " + class_name +
                           " - IllegalAccessException");
      }
      catch(ClassNotFoundException e)
      {
        new MessageDialog(this, "cannot find class: " + class_name);
      } 
      catch(ClassCastException e)
      {
        new MessageDialog(this, class_name + " is not a sub-class of " +
                           "SequenceFormat");
      }
      catch(IOException e) 
      {
        new MessageDialog(this, "I/O error while reading from " +
                           location + ": " + e.getMessage());
      }
      catch(NoSequenceException e) 
      {
        new MessageDialog(this, location + " contained no sequence");
      } 
      catch(InstantiationException e) 
      {
        new MessageDialog(this, "cannot instantiate " + class_name);
      } 
      catch(OutOfRangeException e) 
      {
        new MessageDialog(this, "read failed: one of the features in " +
                           location +
                           " has an out of range location: " +
                           e.getMessage());
      }
    } 
    else 
      new MessageDialog(this, "the -biojava option needs two arguments");
  }

  /**
   *  Read the entries given in the uk.ac.sanger.artemis.ini file.
   **/
  private void readDefaultEntries() 
  {
    final EntryInformation artemis_entry_information =
                     Options.getArtemisEntryInformation();

    final String default_sequence_file_name =
                     Options.getOptions().getDefaultSequenceFileName();

    final String default_feature_file_name =
                     Options.getOptions().getDefaultFeatureFileName();

    if(default_sequence_file_name != null) 
    {
      final String default_sequence_file_name_embl =
                     default_sequence_file_name + "_embl";

      uk.ac.sanger.artemis.io.Entry new_embl_entry = null;

      // try opening the default sequence file with "_embl" added to the name
      // if that fails try the plain sequence file name

      final Document entry_document =
        DocumentFactory.makeDocument(default_sequence_file_name_embl);

      final InputStreamProgressListener progress_listener =
        getInputStreamProgressListener();

      entry_document.addInputStreamProgressListener(progress_listener);

      if(entry_document.readable()) 
      {
        new_embl_entry =
          EntryFileDialog.getEntryFromFile(this,
                                            entry_document,
                                            artemis_entry_information,
                                            false);
      }

      if(new_embl_entry == null || new_embl_entry.getSequence() == null ||
          new_embl_entry.getSequence().length() == 0) 
      {
        final File entry_file = new File(default_sequence_file_name);

        if(entry_file.exists())
        {
          new_embl_entry =
            EntryFileDialog.getEntryFromFile(this,
                                              entry_document,
                                              artemis_entry_information,
                                              false);
        }
        else
        {
          // read failed
          System.err.println("file does not exist: " +
                              default_sequence_file_name +
                              "(given in options files)");
          return;
        }
      }

      if(new_embl_entry == null || new_embl_entry.getSequence() == null ||
          new_embl_entry.getSequence().length() == 0) 
      {
        // read failed
        System.err.println("failed to read " + default_sequence_file_name +
                            "(given in options files)");
        return;
      }

      getStatusLabel().setText("");

      try 
      {
        final Entry entry = new Entry(new_embl_entry);

        final EntryEdit new_entry_edit = makeEntryEdit(entry);

        new_entry_edit.setVisible(true);

        if(default_feature_file_name != null) 
        {
          final Document feature_document =
            DocumentFactory.makeDocument(default_feature_file_name);

          final uk.ac.sanger.artemis.io.Entry new_embl_table_entry =
            EntryFileDialog.getEntryFromFile(this,
                                              feature_document,
                                              artemis_entry_information,
                                              false);

          if(new_embl_table_entry == null)  // open failed
            return;

          final EntryGroup entry_group = new_entry_edit.getEntryGroup();

          final Entry new_table_entry =
            new Entry(entry.getBases(), new_embl_table_entry);

          entry_group.add(new_table_entry);
        }
      } catch(OutOfRangeException e) {
        new MessageDialog(this, "read failed: one of the features in " +
                           default_feature_file_name +
                           " has an out of range location: " +
                           e.getMessage());
      } catch(NoSequenceException e) {
        new MessageDialog(this, "read failed: " +
                           new_embl_entry.getName() +
                           " contains no sequence");
      }
    }
  }

  /**
   *  Make an EntryEdit component from the given Entry.
   **/
  public static EntryEdit makeEntryEdit(final Entry entry) 
  {
    final Bases bases = entry.getBases();
    final EntryGroup entry_group = new SimpleEntryGroup(bases);
    entry_group.add(entry);
    final EntryEdit entry_edit = new EntryEdit(entry_group);

    return entry_edit;
  }

  /**
   *  Add an EntryEdit object to our list of objects.
   *  @param entry_edit The object to add.
   **/
  private synchronized void addEntryEdit(EntryEdit entry_edit) 
  {
    entry_edit_objects.addElement(entry_edit);
  }


  /**
   *  Read an Entry from the given EntrySource and make a new EntryEdit
   *  component for the Entry.
   **/
  private void getEntryEditFromEntrySource(final EntrySource entry_source) 
  {
    try
    {
      final Entry entry = entry_source.getEntry(true);
      if(entry == null)
        return ;

      final EntryGroup entry_group =
           new SimpleEntryGroup(entry.getBases());
      entry_group.add(entry);

      SwingUtilities.invokeLater(new Runnable() 
      {
        public void run() 
        {
          EntryEdit entry_edit = new EntryEdit(entry_group);
          entry_edit.setVisible(true);
        }
      });
    }
    catch(OutOfRangeException e)
    {
      new MessageDialog(ArtemisMain.this, "read failed: one of the features in " +
               " the entry has an out of range " +
               "location: " + e.getMessage());
    }
    catch(NoSequenceException e)
    {
      new MessageDialog(ArtemisMain.this, "read failed: entry contains no sequence");
    }
    catch(IOException e)
    {
      new MessageDialog(ArtemisMain.this, "read failed due to IO error: " + e);
    }
  }

  
  /**
   *  Close the main frame and all EntryEdit frames and then this frame,
   *  then exit.
   **/
  protected void exit() 
  {
    Splash.exitApp();
//  for(int i=0 ; i<entry_edit_objects.size() ;++i) 
//    entryEditFinished(entry_edit_objects.elementAt(i));
 
//  if(filemanager != null)
//    filemanager.setVisible(false);

//  setVisible(false);
//  dispose();
//  System.gc();
  }

  /**
   *  Main entry point for the stand-alone version of Artemis.
   **/
  public static void main(final String [] args) 
  { 
    SwingUtilities.invokeLater(new Runnable() 
    {
      public void run() 
      {
        final ArtemisMain main_window = new ArtemisMain(args);
        main_window.setVisible(true);
        main_window.readArgsAndOptions(args, main_window);
      }
    });
  }

}
