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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/ArtemisMain.java,v 1.12 2005-05-31 15:35:53 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.sequence.NoSequenceException;
import uk.ac.sanger.artemis.sequence.Bases;
import uk.ac.sanger.artemis.util.Document;
import uk.ac.sanger.artemis.util.DocumentFactory;
import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.util.InputStreamProgressListener;
import uk.ac.sanger.artemis.io.EntryInformation;

import org.biojava.bio.seq.io.SequenceFormat;

import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.io.*;
import java.util.Hashtable;
import java.util.Enumeration;

import javax.swing.ListSelectionModel;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JScrollPane;

/**
 *  The main window for the Artemis sequence editor.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: ArtemisMain.java,v 1.12 2005-05-31 15:35:53 tjc Exp $
 **/

public class ArtemisMain extends Splash 
{
  /** Version String use for banner messages and title bars. */
  public static final String version = "Release 8";

  /** A vector containing all EntryEdit object we have created. */
  private EntryEditVector entry_edit_objects = new EntryEditVector();

  protected static FileManager filemanager = null;
  /**
   *  The constructor creates all the components for the main Artemis 
   *  window and sets up all the menu callbacks.
   **/
  public ArtemisMain() 
  {
    super("Artemis", "Artemis", version);

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
        DatabaseEntrySource entry_source = new DatabaseEntrySource();
        if(!entry_source.setLocation())
          return;

        String id = displayDatabases(entry_source);

        if(id != null)
          getEntryEditFromDatabase(id, entry_source);
      }
    };

    final boolean sanger_options =
      Options.getOptions().getPropertyTruthValue("sanger_options");

    if(sanger_options)
      makeMenuItem(file_menu, "Database Entry ...", menu_listener);
   
    menu_listener = new ActionListener() 
    {
      public void actionPerformed(ActionEvent event)
      {
        exit();
      }
    };
    makeMenuItem(file_menu, "Quit", menu_listener);

//      getCanvas().addMouseListener(new MouseAdapter() {
//        /**
//         *  Listen for mouse press events so that we can do popup menus and
//         *  selection.
//         **/
//        public void mousePressed(MouseEvent event) {
//          handleCanvasMousePress(event);
//        }
//      });

//  java.util.Properties props = System.getProperties();
//  java.util.Enumeration en = props.propertyNames();
//  while(en.hasMoreElements())
//  {
//    String prop = (String)en.nextElement();
//    System.out.println(prop+":: "+props.getProperty(prop));
//  }

  }


// XXX add pasteClipboard() one day

//    /**
//     *  Handle a mouse press event on the drawing canvas - select on click,
//     *  select and broadcast it on double click.
//     **/
//    private void handleCanvasMousePress(MouseEvent event) {
//      if(event.getID() != MouseEvent.MOUSE_PRESSED) {
//        return;
//      }

//      if((event.getModifiers() & InputEvent.BUTTON2_MASK) != 0) {
//        pasteClipboard();
//      }
//    }

  /**
   *  Read the entries named in args and in the diana.ini file.
   **/
  public void readArgsAndOptions(final String [] args)
  {
    if(args.length == 0) 
    {
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
          EntryFileDialog.getEntryFromFile(this, entry_document,
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
      } 
      else
      {
        // new sequence file

        if(last_entry_edit != null) 
        {
          last_entry_edit.setVisible(true);
          last_entry_edit = null;
        }

        final Document entry_document =
          DocumentFactory.makeDocument(new_entry_name);

        entry_document.addInputStreamProgressListener(getInputStreamProgressListener());
 
        final uk.ac.sanger.artemis.io.Entry new_embl_entry =
          EntryFileDialog.getEntryFromFile(this, entry_document,
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

    for(int entry_index=0; entry_index<entry_edit_objects.size();
        ++entry_index) 
      entry_edit_objects.elementAt(entry_index).setVisible(true);
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
  private EntryEdit makeEntryEdit(final Entry entry) 
  {
    final Bases bases = entry.getBases();
    final EntryGroup entry_group = new SimpleEntryGroup(bases);
    entry_group.add(entry);
    final EntryEdit entry_edit = new EntryEdit(entry_group);

    return entry_edit;
  }

  /**
   *  This method gets rid of an EntryEdit object and it's frame.  Each object
   *  removed with this method must have been added previously with
   *  addEntryEdit().
   *  @param entry_edit The object to get rid of.
   **/
  public void entryEditFinished(EntryEdit entry_edit) 
  {
    if(null == entry_edit) 
      throw new Error("entryEditFinished() was passed a null object");

    if(!entry_edit_objects.removeElement(entry_edit)) 
      throw new Error("entryEditFinished() - could not remove a " +
                       "object from an empty vector");

    entry_edit.setVisible(false);
    entry_edit.dispose();
  }

  /**
   *  Add an EntryEdit object to our list of objects.
   *  @param entry_edit The object to add.
   **/
  public synchronized void addEntryEdit(EntryEdit entry_edit) 
  {
    entry_edit_objects.addElement(entry_edit);
  }


  /**
  *
  * Display a list of the available relational database entries.
  *
  */
  private String displayDatabases(DatabaseEntrySource entry_source)
  {
    Hashtable db = entry_source.getDatabaseEntries();

    String db_array[] = new String[db.size()];

    Enumeration enum_db = db.keys();
    for(int i=0; enum_db.hasMoreElements(); i++)
      db_array[i] = (String)enum_db.nextElement();

    java.util.Arrays.sort(db_array);
    JList list_db = new JList(db_array);
    list_db.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
    JScrollPane scroll_list = new JScrollPane(list_db);

    Object[] options = { "OPEN IN ARTEMIS >", "CANCEL"};
    boolean selecting = true;

    int select = JOptionPane.showOptionDialog(null, scroll_list,
                                "Database Selection",
                                 JOptionPane.YES_NO_CANCEL_OPTION,
                                 JOptionPane.QUESTION_MESSAGE,
                                 null,
                                 options,
                                 options[0]);
    if(select == 1 || list_db.getSelectedValue() == null)
      return null;

    return (String)db.get((String)list_db.getSelectedValue());
  }

  /**
  *
  * Retrieve a database entry.
  *
  */
  private void getEntryEditFromDatabase(final String id,
                                        final DatabaseEntrySource entry_source)
  {
    SwingWorker entryWorker = new SwingWorker()
    {
      public Object construct()
      {
        try
        {
          final InputStreamProgressListener progress_listener =
                                     getInputStreamProgressListener();

//        DatabaseEntrySource entry_source = new DatabaseEntrySource();

          final Entry entry = entry_source.getEntry(id, progress_listener);
          if(entry == null)
            return null;

          final EntryEdit new_entry_edit = makeEntryEdit(entry);
          new_entry_edit.setVisible(true);
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
        return null;
      }

    };
    entryWorker.start();

  }


  /**
   *  Read an Entry from the given EntrySource and make a new EntryEdit
   *  component for the Entry.
   **/
  private void getEntryEditFromEntrySource(final EntrySource entry_source) 
  {
    SwingWorker entryWorker = new SwingWorker()
    { 
      EntryEdit entry_edit;
      public Object construct()
      {
        try
        {
          final Entry entry = entry_source.getEntry(true);
          if(entry == null)
            return null;

          final EntryGroup entry_group =
              new SimpleEntryGroup(entry.getBases());

          entry_group.add(entry);
          entry_edit = new EntryEdit(entry_group);
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
   *  Force all the EntryEdit components to be redisplayed.
   **/
  private void redisplayAll() 
  {
    for(int i=0 ; i<entry_edit_objects.size() ; ++i) 
      entry_edit_objects.elementAt(i).redisplay();
  }

  /**
   *  Close the main frame and all EntryEdit frames and then this frame,
   *  then exit.
   **/
  protected void exit() 
  {
//  for(int i=0 ; i<entry_edit_objects.size() ;++i) 
//    entryEditFinished(entry_edit_objects.elementAt(i));
 
//  if(filemanager != null)
//    filemanager.setVisible(false);

//  setVisible(false);
//  dispose();
//  System.gc();
    System.exit(0);
  }

  /**
   *  Main entry point for the stand-alone version of Artemis.
   **/
  public static void main(final String [] args) 
  {
    final ArtemisMain main_window = new ArtemisMain();
    main_window.setVisible(true);

    SwingWorker entryWorker = new SwingWorker()
    {
      public Object construct()
      {
        // read the entries given on the command line and in the diana.ini file
        main_window.readArgsAndOptions(args);
        return null;
      }
    };
    entryWorker.start();
  }

}
