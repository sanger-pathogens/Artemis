/* DbfetchEntrySource.java
 *
 * created: Fri Nov 28 2003
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2018  Genome Research Limited
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

import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.sequence.*;

import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.util.Document;
import uk.ac.sanger.artemis.util.DocumentFactory;
import uk.ac.sanger.artemis.io.DocumentEntryFactory;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.EntryInformationException;
import uk.ac.sanger.artemis.io.SimpleEntryInformation;

import java.io.*;
import java.util.regex.Pattern;

import javax.swing.*;

/**
 *  This is an EntrySource that reads Entry objects from the EMBL Dbfetch
 *  server.
 *
 *  @author Kim Rutherford
 **/

public class DbfetchEntrySource
    implements EntrySource 
{
  // The number of databases may be expanded in the future and if so these
  // regular expressions would need to be made more specific, or provide a 
  // means of selecting a database.
  private static final Pattern ENASEQ_PATTERN  = Pattern.compile("^[a-zA-Z]{1,4}\\d+([.]\\d+)?$");
  private static final Pattern REFSEQN_PATTERN = Pattern.compile("^[a-zA-Z]{2}_\\d+([.]\\d+)?$");
  
  static final String EBI_URL  			= "https://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=%s&id=%s&format=%s&style=raw";
  static final String EBI_ENA_SEQ_DB 	= "ena_sequence";
  static final String EBI_REF_SEQ_DB 	= "refseqn";
  static final String EBI_FORMAT_FASTA 	= "fasta";
  static final String EBI_FORMAT_EMBL 	= "embl";
  
  /**
   *  The JFrame that was passed to the constructor.
   **/
  private JFrame frame = null;
  
  
  /**
   *  Create a new DbfetchEntrySource.
   *  @param frame The component that created this EntrySource.  (Used for
   *    requesters.)
   **/
  public DbfetchEntrySource (final JFrame frame) 
  {
	  this.frame = frame;
  }

  /**
   *  Get an Entry object from the Ensembl Dbfetch server.
   *  @param bases The Bases object to pass to the Entry constructor.
   *  @param show_progress If true a InputStreamProgressDialog will be shown
   *    while loading.  (Not implemented)
   *  @exception OutOfRangeException Thrown if one of the features in
   *    the Entry is out of range of the Bases object.
   *  @return null if and only if the user cancels the read or if the read
   *    fails.
   **/
  public Entry getEntry (final Bases bases, final boolean show_progress)
      throws OutOfRangeException, IOException 
  {
	MessageDialog waitingDialog = null;
	  
    final TextDialog text_dialog =
      new TextDialog (getFrame (), "Enter an accession number:", 10, "");

    final String text = text_dialog.getText ();

    if (text == null) {
      // user cancel
      return null;
    }

    final String embl_id = text.trim ();

    if(embl_id.length () > 0) 
    {
    	  Document url_document = null;
    	
      final LogReadListener read_event_logger = new LogReadListener (embl_id);

      final EntryInformation entry_information =
        new SimpleEntryInformation(Options.getArtemisEntryInformation ());
      
      String url = constructEbiUrl(embl_id);
      if (url != null) 
      {
    	  	url_document = DocumentFactory.makeDocument(url);
      }
      else 
      {
    	  	// Accession number format is not recognised
    	  	new MessageDialog (getFrame(), "invalid accession number format");
        return null;
      }

      try
      {
    	    uk.ac.sanger.artemis.io.Entry new_embl_entry = null;
    	  
    	  	try 
    	  	{
    	      waitingDialog = new MessageDialog (getFrame(), "downloading", "fetching from EBI, please wait", false);

          new_embl_entry = DocumentEntryFactory
        		  					.makeDocumentEntry(entry_information,
                                                 url_document,
                                                 read_event_logger);
    	  	}
        finally 
    	  	{
  		  if (waitingDialog != null)
  		  {
  		    waitingDialog.dispose();
  		  }
    	  	}

        if (read_event_logger.seenMessage()) {
          final YesNoDialog yes_no_dialog =
            new YesNoDialog (frame,
                             "there were warnings while reading - view now?");

          if (yes_no_dialog.getResult ()) {
            Splash.showLog ();
          }
        }

        final Bases real_bases;

        if (bases == null) {
          if (new_embl_entry.getSequence () == null) {
            final String message =
              "the entry contains no sequence: " + embl_id;
            new MessageDialog (getFrame (), message);
            return null;
          }

          real_bases = new Bases (new_embl_entry.getSequence ());
        } else {
          real_bases = bases;
        }

        return new Entry (real_bases, new_embl_entry);
      }
      catch (EntryInformationException e) 
      {
        throw new Error ("internal error - unexpected exception: " + e);
      }
    }

    return null;
  }

  /**
   *  Get an Entry object from the Ensembl Dbfetch server.
   *  @param show_progress If true a InputStreamProgressDialog will be shown
   *    while loading.  (Not implemented)
   *  @exception OutOfRangeException Thrown if one of the features in
   *    the Entry is out of range of the Bases object.
   *  @exception NoSequenceException Thrown if the entry that we read has no
   *    sequence.
   *  @return null if and only if the user cancels the read or if the read
   *    fails.
   **/
  public Entry getEntry (final boolean show_progress)
      throws OutOfRangeException, NoSequenceException, IOException 
  {
    return getEntry (null, show_progress);
  }

  /**
   *  Return the name of this source (for display to the user in menus).
   **/
  public String getSourceName () 
  {
    return "EBI - Dbfetch";
  }

  /**
   *  Returns true if and only if this EntrySource always returns "full"
   *  entries.  ie. entries that contain features and sequence.  Entries that
   *  are read from EMBL always contain sequence so in this class this method
   *  returns false.
   **/
  public boolean isFullEntrySource () 
  {
    return true;
  }

  /**
   *  Return the JFrame that was passed to the constructor.
   **/
  public JFrame getFrame () 
  {
    return frame;
  }
  
  /**
   * Determine a URL to use for an EBI database query, based on the given
   * accession number. Returns null if format is not recognised.
   * @param accessionNumber String
   */
  String constructEbiUrl(String accessionNumber) 
  {
	String result = null;
	
	if( ENASEQ_PATTERN.matcher(accessionNumber).matches() )
	{
	  result = String.format(EBI_URL, EBI_ENA_SEQ_DB, accessionNumber, EBI_FORMAT_EMBL);
	}
    else if (REFSEQN_PATTERN.matcher(accessionNumber).matches() )
    {
    	  result = String.format(EBI_URL, EBI_REF_SEQ_DB, accessionNumber, EBI_FORMAT_FASTA);
    }
    else 
    {
      result = null;
    }
	
	return result;
  }

}
