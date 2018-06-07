/* DocumentEntryAutosaveThread.java
 *
 * created: Wed Aug 16 2000
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2000  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/DocumentEntryAutosaveThread.java,v 1.8 2008-08-01 12:50:16 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.*;
import java.io.*;

/**
 *  This is a Thread that automatically saves a DocumentEntry to a backup file
 *  called #entry_name# every 120 seconds.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: DocumentEntryAutosaveThread.java,v 1.8 2008-08-01 12:50:16 tjc Exp $
 **/

public class DocumentEntryAutosaveThread extends Thread {
	
  private static boolean disabled = false;
  
  /**
   *  Create a new DocumentEntryAutosaveThread with MIN_PRIORITY.
   **/
  public DocumentEntryAutosaveThread (final DocumentEntry document_entry) {
    this.document_entry = document_entry;

    setPriority (Thread.MIN_PRIORITY);
  }

  /**
   *  The length of time (in seconds) to sleep for.
   **/
  private final static int SLEEP_TIME = 2 * 60 * 1000;

  /**
   *  The main code for DocumentEntryAutosaveThread.  Attempts to autosave the
   *  DocumentEntry every 2 minutes if it has unsaved changes.
   **/
  public void run () {
    java.util.Date last_save_time = null;
    
    if (disabled)
    	  return;

    // Set to true the first time we save.
    boolean have_saved = false;

    try {
      // sleep for 240 seconds before starting - there is no point in attempting
      // to save straight away
      Thread.sleep (240000);
    } catch (InterruptedException ie) {
    }

    while (true) {
      final String entry_name = document_entry.getName ();

      if (entry_name == null) {
        continue;
      }

      final File save_file;

      if(document_entry.getDocument() instanceof RemoteFileDocument ||
          isMac())
        save_file = new File(System.getProperty("user.dir")+
                             System.getProperty("file.separator")+
                             "#" + entry_name + "#");
      else
        save_file = new File ("#" + entry_name + "#");

      final java.util.Date last_change_time =
        document_entry.getLastChangeTime ();

      // only save if there are unsaved changes and we haven't autosaved
      // since the last change
      if (document_entry.hasUnsavedChanges () &&
          (last_save_time == null ||
           last_change_time != null &&
           last_change_time.after (last_save_time))) {
        final Document save_document = new FileDocument (save_file);
        try {
          final Writer out_file = save_document.getWriter ();
          document_entry.writeToStream (out_file);
          out_file.close ();
          have_saved = true;
        } catch (IOException e) {
          System.err.println ("warning: could not auto save to: " +
                              save_file.getAbsolutePath() +
                              " (will try again later)");
        } catch (java.util.ConcurrentModificationException e) {
          // this Exception means that the tree that stores the features has
          // changed since we started writing.  since this thread doesn't
          // change the tree, this exception is harmless, so we ignore it
          // and try again later
        }
        catch(NullPointerException npe)
        {
          break;
        }
      } else {
        if (have_saved) {
          // auto save file isn't needed now so turn it into a backup file
          final File new_name;

          if(document_entry.getDocument() instanceof RemoteFileDocument ||
             isMac())
            new_name = new File (System.getProperty("user.dir")+
                                 System.getProperty("file.separator")+
                                 entry_name + "~");
          else
            new_name = new File (entry_name + "~");

          save_file.renameTo (new_name);
        } else {
          // the file wasn't made by us so leave it
        }
      }

      last_save_time = last_change_time;

      try {
        Thread.sleep (SLEEP_TIME);
      } catch (InterruptedException ie) {
      }
    }
  }

  private boolean isMac() 
  {
    return System.getProperty("mrj.version") != null ||
           System.getProperty("os.name").toLowerCase().indexOf("mac") >= 0;
  }
  
  public static boolean isDisabled()
  {
	return disabled;
  }

  public static void setDisabled(boolean disabled)
  {
	DocumentEntryAutosaveThread.disabled = disabled;
  }

/**
   *  The DocumentEntry we will save.
   **/
  private DocumentEntry document_entry;
}
