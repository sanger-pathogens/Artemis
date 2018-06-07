/* ActMain.java
 *
 * created: Wed May 10 2000
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/ActMain.java,v 1.17 2008-11-17 13:52:34 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.components.filetree.FileManager;
import uk.ac.sanger.artemis.components.filetree.LocalAndRemoteFileManager;
import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.sequence.Bases;
import uk.ac.sanger.artemis.sequence.NoSequenceException;
import uk.ac.sanger.artemis.components.database.DatabaseEntrySource;
import uk.ac.sanger.artemis.components.database.DatabaseTreeNode;

import uk.ac.sanger.artemis.util.*;
import uk.ac.sanger.artemis.io.DatabaseDocumentEntry;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.SimpleEntryInformation;

import java.awt.event.*;
import java.io.File;
import java.io.IOException;
import java.net.URL;

import javax.swing.JFrame;

/**
 *  The main window for the Artemis Comparison Tool.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: ActMain.java,v 1.17 2008-11-17 13:52:34 tjc Exp $
 **/

public class ActMain extends Splash 
{
  /**  */
  private static final long serialVersionUID = 1L;
  /** Version String use for banner messages and title bars. */
  public static final String version = "Release 6";
  /** File manager */
  protected static FileManager filemanager = null;
  private static DatabaseEntrySource dbEntrySource;

  /**
   *  The constructor creates all the components for the main ACT window
   *  and sets up all the menu callbacks.
   **/
  public ActMain() 
  {
    super("Artemis Comparison Tool", "ACT", version);

    ActionListener open_listener = new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        makeOpenDialog();
      }
    };

    makeMenuItem(file_menu, "Open ...", open_listener);

    ActionListener quit_listener = new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        exit();
      }
    };

    ActionListener menu_listener_ssh = new ActionListener()
    {
      private LocalAndRemoteFileManager fm;
      public void actionPerformed(ActionEvent event)
      {
        if(fm == null)
          fm = new LocalAndRemoteFileManager(ActMain.this);
        else
          fm.setVisible(true);
        dbEntrySource = fm.getDatabaseEntrySource();
        new ComparatorDialog(ActMain.this).setVisible(true);
      }
    }; 
    
    if(System.getProperty("chado") != null)
      makeMenuItem(file_menu, "Open Database and SSH File Manager ...", menu_listener_ssh);
    else 
      makeMenuItem(file_menu, "Open SSH File Manager ...", menu_listener_ssh);

/*    final boolean sanger_options =
      Options.getOptions().getPropertyTruthValue("sanger_options");

    if(sanger_options)
    {
      ActionListener menu_listener = new ActionListener()
      {
        public void actionPerformed(ActionEvent event)
        {
          launchDatabaseJFrame(true);
        }
      };
      
      makeMenuItem(file_menu, "Database Entry ...", menu_listener);
      if(System.getProperty("chado") != null)
        launchDatabaseJFrame(false);
    }*/
    
    makeMenuItem(file_menu, "Quit", quit_listener);
  }

  /**
   *  Make a new Comparator component from the given files.
   *  @param frame The JFrame used when making a new MessageDialog.
   *  @param progress_listener The object to which InputStreamProgressEvents
   *    will be send while reading.  Can be null.
   *  @param file_names Alternating sequence and comparison data file names.
   *    Must be >= 3.  I there are an even number of file names the first
   *    file/sequence object will be added to the send of the display and the
   *    last comparison file will be assumed to be a a comparison between the
   *    last and first sequence files.
   **/
  public static boolean makeMultiComparator(final JFrame frame,
                          final InputStreamProgressListener progress_listener,
                          final Object[] file_names) 
  {
    processJnlpAttributes();
    final ProgressThread progress_thread = new ProgressThread(null,
                                                "Loading Entry...");

    SwingWorker entryWorker = new SwingWorker()
    {
      public Object construct()
      {
        progress_thread.start();
        final EntryGroup[] entry_group_array =
          new EntryGroup[file_names.length / 2 + 1];

        final ComparisonData[] comparison_data_array =
          new ComparisonData[file_names.length / 2];

        for(int i = 0; i<file_names.length; i += 2) 
        {
          final EntryInformation entry_information =
            new SimpleEntryInformation(Options.getArtemisEntryInformation());

          final Object this_file_name = file_names[i];
          //File this_file = new File(this_file_name);
 
          try 
          {
            if(!openEntry(this_file_name, entry_group_array, 
                               entry_information, i))
              return null;
          } 
          catch(OutOfRangeException e) 
          {
            new MessageDialog(frame, "read failed: one of the features has an " +
                               "out of range location: " + e.getMessage());
            return null;
          }
        }

        // add the first entry at the end to make the MultiComparator
        // circular(-ish)
        if(file_names.length % 2 == 0) 
          entry_group_array[entry_group_array.length - 1] = entry_group_array[0];

        try
        {
          for(int i = 1 ; i < file_names.length ; i += 2) 
          {
            final String comparison_data_file_name = (String)file_names[i];
            final Document comparison_data_document =
              DocumentFactory.makeDocument(comparison_data_file_name);

            comparison_data_array[i / 2] =
              ComparisonDataFactory.readComparisonData(comparison_data_document);

            final Bases prev_bases = entry_group_array[i/2].getBases();
            final Bases next_bases = entry_group_array[i/2 + 1].getBases();

            final ComparisonData swapped_comparison_data =
              comparison_data_array[i / 2].flipMatchesIfNeeded(prev_bases,
                                                            next_bases);

            if(swapped_comparison_data != null) 
              comparison_data_array[i / 2] = swapped_comparison_data;

            if(swapped_comparison_data != null)
            {
              final MessageFrame message_frame =
                new MessageFrame("note: hits from " + comparison_data_file_name +
                                  " have been flipped to match the " +
                                  "sequences");

              message_frame.setVisible(true);
            }
          }
        } 
        catch(IOException e) 
        {
          new MessageDialog(frame, "error while reading: " + e.getMessage());
          return null;
        }
        catch(OutOfRangeException e) 
        {
          new MessageDialog(frame, "comparison file read failed:  " +
                             "out of range error: " + e.getMessage());
          return null;
        }
        
        final MultiComparator comparator =
          new MultiComparator(entry_group_array,
                              comparison_data_array,
                              progress_listener);

        comparator.setVisible(true);
        return null;
      }

      public void finished()
      {
        if(progress_thread !=null)
          progress_thread.finished();
      }

      private boolean openEntry(Object this_file_name, EntryGroup[] entry_group_array,
                                final EntryInformation entry_information, int i)
                      throws OutOfRangeException
      {
        uk.ac.sanger.artemis.io.Entry embl_entry = null;
        Entry entry = null;
        // test if this is a database entry rather than a file 
        if(this_file_name instanceof DatabaseTreeNode)
        {
          DatabaseTreeNode dbNode = (DatabaseTreeNode)this_file_name;
          try
          {
            entry = dbEntrySource.getEntry(dbNode.getFeatureId(), 
                dbNode.getUserName(), progress_listener);
            
            boolean isMitochondrial = false;
            if(dbNode.getFeatureType() != null &&
               dbNode.getFeatureType().startsWith("mitochondrial_"))
              isMitochondrial = true;
            boolean readOnly = DatabaseTreeNode.setOrganismProps(
                dbNode.getOrganism().getOrganismProps(), isMitochondrial);
            embl_entry = (DatabaseDocumentEntry)entry.getEMBLEntry();
            ((DatabaseDocumentEntry)embl_entry).setReadOnly(readOnly);
          }
          catch(NoSequenceException e)
          {
            e.printStackTrace();
          }
          catch(IOException e)
          {
            e.printStackTrace();
          }
          
          DatabaseDocument doc = 
            (DatabaseDocument)((DatabaseDocumentEntry)embl_entry).getDocument();
          doc.setName((String)dbNode.getUserObject());
        }
        else
        {
          final Document entry_document =
              DocumentFactory.makeDocument((String)this_file_name);

          if(progress_listener != null)
            entry_document.addInputStreamProgressListener(progress_listener);

           embl_entry =
              EntryFileDialog.getEntryFromFile(frame, entry_document,
                                         entry_information,
                                         false);
        }
        
        // getEntryFromFile() has alerted the user so we just need to quit
        if(embl_entry == null)
          return false;

        final uk.ac.sanger.artemis.io.Sequence sequence =
                                       embl_entry.getSequence();

        if(sequence == null)
        {
          new MessageDialog(frame, "This file contains no sequence: " +
                             (String)this_file_name);
          return false;
        }
        
        if(entry == null)
        {
          final Bases embl_bases = new Bases(sequence);
          entry = new Entry(embl_bases, embl_entry);
        }
        
        EntryGroup entry_group = new SimpleEntryGroup(entry.getBases());
        
        entry_group.add(entry);
        entry_group_array[i / 2] = entry_group;
        return true;
      }
    };
    entryWorker.start();
    return true;
  }


  /**
   *  Create a dialog that allow the user to the choose two files to compare
   *  and a file containing comparison data.
   **/
  private void makeOpenDialog() 
  {
    if(filemanager == null)
      filemanager = new FileManager(this,null);
    else
      filemanager.setVisible(true);
    new ComparatorDialog(this).setVisible(true);
  }

  /**
   *  Exit from ACT.
   **/
  protected void exit() 
  {
    System.exit(0);
  }
  
  /**
   * Validate the program input arguments and exit
   * if they are not valid, displaying an error message(s).
   * @param args String array
   */
  protected static void validateStartupArguments(final String[] args)
  {
	boolean valid = true;
	  
	if(args.length >= 3) 
    {
	  // Make sure the files provided are actually valid as far as possible...
      for (String file : args)
      {
    	    if (file.startsWith("ftp") || file.startsWith("http"))
    	    {
    	    	  // web resource
    	    	  try 
    	    	  {
    	    	    URL url = new URL(file);
    	    	  }
    	    	  catch (Exception e)
    	    	  {
    	    		valid = false;
    			System.err.println("\nError - " + file + " is not a valid URL.");
    	    	  }
    	    }
    	    	else
    	    	{
    	    	  // normal file
		  File argFile = new File(file);
		  if (!argFile.exists() || !argFile.isFile())
		  {
		    	valid = false;
		    	System.err.println("\nError - " + argFile + " is not a valid file.");
		  }
    	    }
      }
	 
    }
    else 
    {
      if(args.length != 0) 
      {
        System.err.println("\nError - this program needs either no" +
                           " arguments or an odd number\n" +
                           "(3 or more):");
        System.err.println("   act sequence_1 comparison_data sequence_2");
        System.err.println("or");
        System.err.println("   act seq_1 comparison_data_2_v_1 seq_2 comparison_data_3_v_2 seq_3");
        System.err.println("or");
        System.err.println("   act");
        valid = false;
      }
    }  
	
	if (!valid)
	{
	  System.exit(1);
	}
	
  }
  
  /**
   *  Main entry point for ACT
   **/
  public static void main(final String [] args) 
  {
	validateStartupArguments(args);
    
    final ActMain main_window = new ActMain();
    main_window.setVisible(true);

    final InputStreamProgressListener progress_listener =
      main_window.getInputStreamProgressListener();
    
    if(args.length >= 3) 
    {
    	  ActMain.makeMultiComparator(main_window, progress_listener,
                                 args);
    }
  }

}
