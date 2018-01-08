/* EntryFileDialog.java
 *
 * created: Mon Dec  7 1998
 *
 * This file is part of Artemis
 *
 * Copyright(C) 1998-2003  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/EntryFileDialog.java,v 1.15 2009-04-01 12:23:20 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.util.*;
import uk.ac.sanger.artemis.io.DatabaseDocumentEntry;
import uk.ac.sanger.artemis.io.Entry;
import uk.ac.sanger.artemis.io.DocumentEntryFactory;
import uk.ac.sanger.artemis.io.GFFDocumentEntry;
import uk.ac.sanger.artemis.io.ReadFormatException;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.EntryInformationException;
import uk.ac.sanger.artemis.io.DocumentEntry;
import uk.ac.sanger.artemis.io.ReadAndWriteEntry;

import java.io.*;
import javax.swing.*;

/**
 *  This class is a JFileChooser that can read EMBL Entry objects.
 *
 *  @author Kim Rutherford
 *  @version $Id: EntryFileDialog.java,v 1.15 2009-04-01 12:23:20 tjc Exp $
 **/

public class EntryFileDialog extends StickyFileChooser 
{

  /** */
  private static final long serialVersionUID = 1L;
  /** JFrame reference that was passed to the constructor. */
  private JFrame owner = null;

  /**
   *  Create a new EntryFileDialog component.
   *  @param owner The component where this dialog was created.
   *  @param sequence_only If true default to showing only file that have
   *    suffixes that suggest that the files contain sequence(eg. .embl,
   *    .seq, .dna).  Of false show all files that can be read.
   **/
  public EntryFileDialog(final JFrame owner,
                         final boolean show_sequence_only) 
  {
    super();
    this.owner = owner;

    setFileSelectionMode(JFileChooser.FILES_ONLY);
    setMultiSelectionEnabled(false);

    final StringVector sequence_suffixes =
      Options.getOptions().getOptionValues("sequence_file_suffixes");

    final StringVector feature_suffixes =
      Options.getOptions().getOptionValues("feature_file_suffixes");

    final javax.swing.filechooser.FileFilter artemis_filter =
      new javax.swing.filechooser.FileFilter()
    {
        public boolean accept(final File file) 
        {
          if(file.isDirectory()) 
            return true;

          for(int i = 0; i<sequence_suffixes.size(); ++i) 
          {
            final String this_suffix = (String)sequence_suffixes.elementAt(i);

            if(file.getName().endsWith("." + this_suffix) ||
                file.getName().endsWith("." + this_suffix + ".gz")) 
              return true;
          }

          for(int i = 0; i<feature_suffixes.size(); ++i) 
          {
            final String this_suffix = (String)feature_suffixes.elementAt(i);

            if(file.getName().endsWith("." + this_suffix) ||
               file.getName().endsWith("." + this_suffix + ".gz")) 
              return true;
          }
          return false;
        }

        public String getDescription() 
        {
          return "Artemis files";
        }
      };

    final javax.swing.filechooser.FileFilter feature_filter =
      new javax.swing.filechooser.FileFilter() 
    {
        public boolean accept(final File file) 
        {
          if(file.isDirectory()) 
            return true;

          for(int i = 0 ; i < feature_suffixes.size() ; ++i) 
          {
            final String this_suffix = (String)feature_suffixes.elementAt(i);

            if(file.getName().endsWith("." + this_suffix) ||
               file.getName().endsWith("." + this_suffix + ".gz")) 
              return true;
          }

          return false;
        }

        public String getDescription() 
        {
          return "Feature files";
        }
      };

    final javax.swing.filechooser.FileFilter sequence_filter =
      new javax.swing.filechooser.FileFilter() 
      {
        public boolean accept(final File file) 
        {
          if(file.isDirectory()) 
            return true;

          for(int i = 0 ; i<sequence_suffixes.size() ; ++i) 
          {
            final String this_suffix = (String)sequence_suffixes.elementAt(i);

            if(file.getName().endsWith("." + this_suffix) ||
               file.getName().endsWith("." + this_suffix + ".gz")) 
              return true;
          }

          return false;
        }

        public String getDescription() 
        {
          return "Sequence files";
        }
      };

    addChoosableFileFilter(artemis_filter);
    addChoosableFileFilter(feature_filter);
    addChoosableFileFilter(sequence_filter); 

    if(show_sequence_only) 
      setFileFilter(sequence_filter);
    else 
      setFileFilter(artemis_filter);
  }

  /**
   *  Return an uk.ac.sanger.artemis.io.Entry object representing the file
   *  the user has selected with the dialog or null if the read failed for any
   *  reason.
   *  @param entry_information The EntryInformation to use when reading.  This
   *    supplies the list of valid keys and qualifiers.  If a key or qualifier
   *    is read that is incompatible with this EntryInformation object the
   *    EntryInformation will be changed to cope.
   *  @param listener The object to which InputStreamProgressEvents will be
   *    send while reading.
   *  @exception EntryInformationException Thrown if an Entry using the given
   *    EntryInformation object cannot contain the Key, Qualifier or
   *    Key/Qualifier combination of one of the features in the Document.
   *  @param show_progress If true a InputStreamProgressDialog will be shown
   *    while loading.
   **/
  public Entry getEntry(final EntryInformation entry_information,
                        final InputStreamProgressListener listener,
                        final boolean show_progress) 
  {
    setDialogTitle("Select a file ...");
    setDialogType(JFileChooser.OPEN_DIALOG);
    SecurityManager sm = System.getSecurityManager();
    System.setSecurityManager(null);

    final int status = showOpenDialog(owner);

    System.setSecurityManager(sm);

    if(status != JFileChooser.APPROVE_OPTION ||
       getSelectedFile() == null) 
      return null;

    final File file = new File(getCurrentDirectory(),
                               getSelectedFile().getName());

    return getEntryFromFile(owner, new FileDocument(file),
                            entry_information, show_progress); 
  }

  /**
   *  This exists only to get around a bug in the 1.1/1.2 code generator,
   *  which generates unverifiable code.
   *  @param frame Used when creating MessageDialog components.
   *  @param file_document The Document to read the Entry from(should be made
   *    from entry_file).
   *  @param entry_information The EntryInformation to use when reading.  This
   *    supplies the list of valid keys and qualifiers.  If a key or qualifier
   *    is read that is incompatible with this EntryInformation object the
   *    EntryInformation will be changed to cope.
   **/
  private static Entry getEntryFromFileHelper(final JFrame frame,
                            final Document file_document,
                            final EntryInformation entry_information)
      throws ReadFormatException, IOException 
  {

    final LogReadListener read_event_logger =
      new LogReadListener(file_document.getName());

    final Entry new_entry;

    try 
    {
      new_entry =
        DocumentEntryFactory.makeDocumentEntry(entry_information,
                                 file_document,read_event_logger);
    }
    catch(EntryInformationException e) 
    {
      throw new Error("internal error - unexpected exception: " + e);
    }

    if(read_event_logger.seenMessage() && !Options.isBlackBeltMode())
    {
      final YesNoDialog yes_no_dialog = new YesNoDialog(frame,
                         "there were warnings while reading - view now?");

      if(yes_no_dialog.getResult()) 
        Splash.showLog();
    }

    return new_entry;
  }

  /**
   *  Read and return an Entry from the given File.
   *  @param frame Used when creating MessageDialog components.
   *  @param entry_file The file to read the Entry from.
   *  @param entry_information The EntryInformation to use when reading.  This
   *    supplies the list of valid keys and qualifiers.  If a key or qualifier
   *    is read that is incompatible with this EntryInformation object the
   *    EntryInformation will be changed to cope.
   *  @param listener The object to which InputStreamProgressEvents will be
   *    send while reading.
   *  @param show_progress If true a InputStreamProgressDialog will be shown
   *    while loading.
   **/
  public static Entry getEntryFromFile(final JFrame frame,
                     final Document entry_document,
                     final EntryInformation entry_information,
                     final boolean show_progress) 
  {
    InputStreamProgressDialog progress_dialog = null;

    if(show_progress) 
    {

    	   // TODO This doesn't work because getEntryFromFile() is called from the Swing thread so the Dialog never gets updated

       progress_dialog =
         new InputStreamProgressDialog(frame, "Reading ...",
                                        "Reading from " +
                                        entry_document.getName(), false);
       final InputStreamProgressListener listener =
         progress_dialog.getInputStreamProgressListener();

       entry_document.addInputStreamProgressListener(listener);
    }

    try 
    {
      return getEntryFromFileHelper(frame, entry_document,
                                    entry_information);
    }
    catch(ReadFormatException e) 
    {
      final String message =
        "failed to open " + entry_document.getName() + ": " +
        e.getMessage() +(e.getLineNumber() > 1 ?
                           " at line: " + e.getLineNumber() :
                           "");
      System.out.println(message);
      new MessageDialog(frame, message);
    }
    catch(FileNotFoundException e) 
    {
      final String message =
        "failed to open " + entry_document.getName() + ": file not found";
      new MessageDialog(frame, message);
    }
    catch(IOException e) 
    {
      final String message =
        "failed to open " + entry_document.getName() + ": " + e.getMessage();
      new MessageDialog(frame, message);
    }
    finally
    {
      if(progress_dialog != null) 
        progress_dialog.dispose();
    }
    return null;
  }

  /**
   *  Save the given entry, prompting for a file name if necessary.
   *  @param include_diana_extensions If true then the any diana additions to
   *    the embl file format will be included in the output, otherwise they
   *    will be removed.  Also possible problems that would cause an entry to
   *    bounce from the EMBL submission process will be flagged if this is
   *    true.
   *  @param ask_for_name If true then always prompt for a new filename,
   *    otherwise prompt only when the entry name is not set.
   *  @param keep_new_name If ask_for_name is true a file will be written with
   *    the new name the user selects - if keep_new_name is true as well, then
   *    the entry will have it's name set to the new name, otherwise it will
   *    be used for this save and then discarded.
   *  @param destination_type Should be one of EMBL_FORMAT, GENBANK_FORMAT or
   *    ANY_FORMAT.  If ANY_FORMAT then the Entry will be saved in the
   *    same format it was created, otherwise it will be saved in the given
   *    format.
   **/
  void saveEntry(final uk.ac.sanger.artemis.Entry entry,
                  final boolean include_diana_extensions,
                  final boolean ask_for_name,
                  final boolean keep_new_name,
                  final int destination_type) 
  {

    JCheckBox remoteSave = new JCheckBox("Ssh/Remote Save",
                                         true);
    File file = null;

    try 
    {
      if(ask_for_name || entry.getName() == null) 
      {
        Box yBox = Box.createVerticalBox();
        boolean useAccessory = false;

        JCheckBox emblHeader = new JCheckBox("Add EMBL Header", false);
        
        JCheckBox removeProductForPseudo = new JCheckBox(
            "Remove products from pseudogenes", false);

        setDialogTitle("Save to ...");
        setDialogType(JFileChooser.SAVE_DIALOG);

        if( destination_type == DocumentEntryFactory.EMBL_FORMAT )
        {
          if((entry.getHeaderText() == null || !isHeaderEMBL(entry.getHeaderText())))
            yBox.add(emblHeader);
          
          if(!include_diana_extensions)
            yBox.add(removeProductForPseudo);
          useAccessory = true;
        }

        final JCheckBox flattenGeneModel = new JCheckBox("Flatten Gene Model",
                                                          false);
        final JCheckBox ignoreObsoleteFeatures = new JCheckBox(
                              "Ignore obsolete features", true);
        if(((DocumentEntry)entry.getEMBLEntry()).getDocument() 
                                       instanceof RemoteFileDocument)
        {
          yBox.add(remoteSave);
          useAccessory = true;
        }
        else if(entry.getEMBLEntry() instanceof DatabaseDocumentEntry ||
                entry.getEMBLEntry() instanceof GFFDocumentEntry)
        {
          yBox.add(flattenGeneModel);
          if(entry.getEMBLEntry() instanceof DatabaseDocumentEntry)
            yBox.add(ignoreObsoleteFeatures);
          useAccessory = true;
        }

        if(useAccessory)
          setAccessory(yBox);

        switch(destination_type) // provide a suffix
        {
          case DocumentEntryFactory.EMBL_FORMAT:
            super.setSelectedFile(new File(".embl"));
            flattenGeneModel.setSelected(true);
            break;
          case DocumentEntryFactory.GENBANK_FORMAT:
            super.setSelectedFile(new File(".gbk"));
            flattenGeneModel.setSelected(true);
            break;
          case DocumentEntryFactory.GFF_FORMAT:
            super.setSelectedFile(new File(".gff"));
            break;
          default:
            break;
        }
        int status = showSaveDialog(owner);

        if(status != JFileChooser.APPROVE_OPTION ||
           getSelectedFile() == null) 
          return;

        if(emblHeader.isSelected())
        {
          Box bdown = Box.createVerticalBox();
          JTextField idField = new JTextField("");
          bdown.add(idField);

          int n = JOptionPane.showConfirmDialog(null, bdown,
                            "Enter the entry ID",
                            JOptionPane.OK_CANCEL_OPTION,
                            JOptionPane.QUESTION_MESSAGE);

          if(n != JOptionPane.CANCEL_OPTION &&
             !idField.getText().trim().equals(""))
          {
            int length = entry.getBases().getLength();
            String header = "ID   "+idField.getText().trim()+"; SV ; ; ; ; ; "+length+" BP.";
            //String header = "ID   "+idField.getText().trim();
            if(entry.getFeatureCount() > 0)
              header = header.concat("\nFH   Key             "+
                                     "Location/Qualifiers\nFH\n");
            entry.setHeaderText(header);
          }
        }

        file = new File(getCurrentDirectory(),
                        getSelectedFile().getName());

        if(file.exists()) 
        {
          final YesNoDialog yes_no_dialog = new YesNoDialog(owner,
                           "this file exists: " + file.getName() +
                           " overwrite it?");

          if(!yes_no_dialog.getResult()) 
            return;
        }

        final MessageDialog message = new MessageDialog(owner,
                             "saving to " + file.getName() + " ...",
                             false);
        try 
        {
          DocumentEntryFactory.REMOVE_PRODUCT_FROM_PSEUDOGENE = removeProductForPseudo.isSelected();
          if(entry.getEMBLEntry() instanceof DatabaseDocumentEntry ||
             entry.getEMBLEntry() instanceof GFFDocumentEntry)
            ReadAndWriteEntry.writeDatabaseEntryToFile(entry, file, 
                flattenGeneModel.isSelected(), 
                ignoreObsoleteFeatures.isSelected(), false, 
                include_diana_extensions, destination_type, owner);
          else if(include_diana_extensions) 
            entry.save(file, destination_type, false);
          else 
            entry.saveStandardOnly(file, destination_type, true);
        }
        catch(EntryInformationException e) 
        {
          final YesNoDialog yes_no_dialog = new YesNoDialog(owner,
                             "destination format can't handle all " +
                             "keys/qualifiers - continue?");
          
          if(yes_no_dialog.getResult()) 
          {
            try 
            {
              if(entry.getEMBLEntry() instanceof DatabaseDocumentEntry)
                ReadAndWriteEntry.writeDatabaseEntryToFile(entry, file, 
                    flattenGeneModel.isSelected(), 
                    ignoreObsoleteFeatures.isSelected(), true, 
                    include_diana_extensions, destination_type, null);
              else if(include_diana_extensions) 
                entry.save(file, destination_type, true);
              else 
                entry.saveStandardOnly(file, destination_type, true);
            }
            catch(EntryInformationException e2) 
            {
              throw new Error("internal error - unexpected exception: "+ e);
            }
          } 
          else 
            return;
        }
        finally 
        {
          DocumentEntryFactory.REMOVE_PRODUCT_FROM_PSEUDOGENE = false;
          if(message != null) 
            message.dispose();
        }

        if(keep_new_name) 
          entry.setName(file.getName());
      }
      else 
      {
        final MessageDialog message = new MessageDialog(owner,
                             "saving to " + entry.getName() + " ...",
                             false);
        try 
        {
          if(include_diana_extensions) 
            entry.save(destination_type);
          else 
            entry.saveStandardOnly(destination_type);
        }
        finally 
        {
          message.dispose();
        }
      }
    } 
    catch(ReadOnlyException e) 
    {
      new MessageDialog(owner, "this entry is read only");
      return;
    }
    catch(IOException e) 
    {
      new MessageDialog(owner, "error while saving: " + e);
      return;
    } 
    catch(EntryInformationException e) 
    {
      new MessageDialog(owner, "error while saving: " + e);
      return;
    }

    // save it back to ssh server
    if(((DocumentEntry)entry.getEMBLEntry()).getDocument()
                                   instanceof RemoteFileDocument &&
        remoteSave.isSelected())
    {
       RemoteFileDocument node =
           (RemoteFileDocument)(((DocumentEntry)entry.getEMBLEntry()).getDocument());

       if(file == null)
         file = new File( ((DocumentEntry)entry.getEMBLEntry()).getDocument().getLocation().toString() );
       node.saveEntry(file);
    }

  }

  /**
  *
  * Test for that the header of an EMBL entry begins
  * with an ID line.
  * @param header embl header
  *
  */
  private boolean isHeaderEMBL(String header)
  {
    StringReader reader = new StringReader(header);
    BufferedReader buff_reader = new BufferedReader(reader);

    try
    {  
      if(!buff_reader.readLine().startsWith("ID"))
        return false;
    }
    catch(IOException ioe){}
    return true;
  }


}

