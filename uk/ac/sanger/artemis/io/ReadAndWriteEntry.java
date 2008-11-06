/* ReadAndWriteEntry.java

 * This file is part of Artemis
 *
 * Copyright(C) 2008  Genome Research Limited
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
 */

package uk.ac.sanger.artemis.io;

import java.awt.BorderLayout;
import java.io.File;
import java.io.IOException;

import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;

import uk.ac.sanger.artemis.components.database.DatabaseEntrySource;
import uk.ac.sanger.artemis.components.genebuilder.GeneUtils;
import uk.ac.sanger.artemis.sequence.NoSequenceException;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.InputStreamProgressEvent;
import uk.ac.sanger.artemis.util.InputStreamProgressListener;
import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.Options;

public class ReadAndWriteEntry
{

  private static org.apache.log4j.Logger logger4j = 
    org.apache.log4j.Logger.getLogger(ReadAndWriteEntry.class);

  private static DatabaseEntrySource ENTRY_SOURCE;
  
  /**
   * Read from the database, given a srcFeature uniquename
   * @param uniqueName
   * @return
   * @throws OutOfRangeException
   * @throws NoSequenceException
   * @throws IOException
   */
  public static Entry readEntryFromDatabase(final String uniqueName,
                                            DatabaseEntrySource entry_source) 
         throws OutOfRangeException, NoSequenceException, IOException
  {
    if(entry_source == null)
    {
      ReadAndWriteEntry.ENTRY_SOURCE = new DatabaseEntrySource();
      entry_source = ENTRY_SOURCE;
      boolean promptUser = true;
      if(System.getProperty("read_only") != null)
        promptUser = false;

      if(!entry_source.setLocation(promptUser))
        return null;
    }
    
    String url = (String)entry_source.getLocation();
    int index  = url.indexOf("?");
    
    String userName = url.substring(index+1).trim();
    if(userName.startsWith("user="))
      userName = userName.substring(5);
    
    final String srcFeatureId = getFeatureId(entry_source, uniqueName);
    
    final InputStreamProgressListener stream_progress_listener =
        new InputStreamProgressListener() 
    {
      public void progressMade(final InputStreamProgressEvent event) 
      {
        final int char_count = event.getCharCount();
        if(char_count != -1) 
          logger4j.debug("chars read so far: " + char_count);
      }
      public void progressMade(String progress)
      {
        logger4j.debug(progress);
      }
    };
    return entry_source.getEntry(srcFeatureId, userName, 
                                 stream_progress_listener);
  }
  
  /**
   * Read from the database, given a srcFeature uniquename
   * @param uniqueName
   * @return
   * @throws OutOfRangeException
   * @throws NoSequenceException
   * @throws IOException
   */
  public static Entry readEntryFromDatabase(final String uniqueName) 
         throws OutOfRangeException, NoSequenceException, IOException
  {
    return readEntryFromDatabase(uniqueName, null);
  }
  
  /**
   * Write entry to a file
   * @param entry
   * @param file
   * @param flatten Flatten the gene model and combine the qualifiers if true.
   *    If false it will write all features and qualifiers out.
   * @param ignore obsolete features if true
   * @param force invalid qualifiers and any features with invalid keys will 
   *    be quietly thrown away when saving.
   * @param destination_type Should be one of EMBL_FORMAT, GENBANK_FORMAT,
   *    GFF_FORMAT or ANY_FORMAT.  If ANY_FORMAT then the Entry will
   *    be saved in the same format it was created, otherwise it will be saved
   *    in the given format.
   * @throws IOException
   * @throws EntryInformationException
   */
  public static void writeDatabaseEntryToFile(final Entry entry, final File file,
                                              final boolean flatten,
                                              final boolean ignoreObsolete,
                                              final boolean force,
                                              final int destination_type,
                                              final JFrame parent) 
         throws IOException, EntryInformationException
  {
    GeneUtils.lazyLoadAll(entry, parent);

    EntryInformation artemis_entry_information = Options.getArtemisEntryInformation();
    if(!flatten)
    {
      final FeatureVector features = entry.getAllFeatures();
      for(int i=0; i<features.size(); i++)
        addAllKeysQualifiers(artemis_entry_information, features.elementAt(i).getEmblFeature());
    }
    PublicDBDocumentEntry.IGNORE_OBSOLETE_FEATURES = ignoreObsolete;
    entry.save(file, destination_type, force, artemis_entry_information);
  }
  
  /**
   * Add all keys and qualifiers for a given feature to the EntryInformation 
   * @param entry_information
   * @param feature
   */
  private static void addAllKeysQualifiers(final EntryInformation entry_information,
                                           final Feature feature)
  {
    Key new_key = feature.getKey();
    
    new_key = PublicDBDocumentEntry.mapKeys(new_key);
    
    boolean keyAdded = false;
    if(!entry_information.isValidKey(new_key))
    {
      entry_information.addKey(new_key);
      keyAdded = true;
    }
    
    final QualifierVector feature_qualifiers = feature.getQualifiers();
    
    // check the qualifiers
    for(int i = 0 ; i < feature_qualifiers.size() ; ++i)
    {
      final Qualifier this_qualifier = (Qualifier)feature_qualifiers.elementAt(i);
      final String this_qualifier_name = this_qualifier.getName();

      if(!entry_information.isValidQualifier(this_qualifier_name) ||
         !entry_information.isValidQualifier(new_key, this_qualifier_name) ||
         keyAdded) 
      {
        QualifierInfo qualifierInfo = entry_information.getQualifierInfo(this_qualifier_name);
        
        if(qualifierInfo == null)
        {
          KeyVector keys = new KeyVector();
          qualifierInfo = new QualifierInfo(this_qualifier_name, QualifierInfo.QUOTED_TEXT,
                                            keys, null, false);
          try
          {
            entry_information.addQualifierInfo(qualifierInfo);
          }
          catch(QualifierInfoException e)
          {
            e.printStackTrace();
          }
        }
        
        if(qualifierInfo.getValidKeys() != null)
          qualifierInfo.getValidKeys().add(new_key);
      }
    }

  }
  
  /**
   * Get feature id
   * @param entry_source
   * @param srcUniqueName
   * @return
   */
  public static String getFeatureId(
      final DatabaseEntrySource entry_source, final String srcUniqueName)
  {
    final DatabaseDocument doc = entry_source.getDatabaseDocument();
    org.gmod.schema.sequence.Feature feature = doc.getFeatureByUniquename(srcUniqueName);
    return Integer.toString(feature.getFeatureId());
  }
  
  public static void main(final String args[])
  {
    try
    {
      String names[];
      boolean flatten = true;
      boolean ignoreObsolete = true;
      
      if( (args != null && args.length == 1 && args[0].startsWith("-h")) ||
          (args == null || args.length < 1))
      {
        System.out.println("-h\tshow help");
        System.out.println("-f\t[y|n] flatten the gene model, default is y");
        System.out.println("-i\t[y|n] ignore obsolete features, default is y");
        System.out.println("-s\tspace separated list of sequences to read and write out");
        System.out.println("-o\t[EMBL|GFF] output format, default is EMBL");
        System.exit(0);
      }
      

      names = args;
      int format = DocumentEntryFactory.EMBL_FORMAT;
      String suffix = ".embl";
      
      for(int i = 0; i < args.length; i++)
      {
        if(args[i].toLowerCase().equals("-f"))
        {
          if(i + 1 < args.length && args[i + 1].toLowerCase().equals("n"))
            flatten = false;
        }

        if(args[i].toLowerCase().equals("-i"))
        {
          if(i + 1 < args.length && args[i + 1].toLowerCase().equals("n"))
            ignoreObsolete = false;
        }
        
        if(args[i].toLowerCase().equals("-o"))
        {
          if(i + 1 < args.length && args[i + 1].toLowerCase().equals("gff"))
          {
            format = DocumentEntryFactory.GFF_FORMAT;
            suffix = ".gff";
          }
        }
      }
      
      java.util.Vector files = null;
      for(int i = 0; i < args.length; i++)
      {
        if(args[i].toLowerCase().equals("-s"))
        {
          if(files == null)
            files = new java.util.Vector();
          for(int j = i + 1; j < args.length; j++)
          {
            if(args[j].startsWith("-"))
              break;
            files.add(args[j]);
          }
        }
        else if(args[i].startsWith("-"))
        {
          i++;
        }
        else
        {
          if(files == null)
            files = new java.util.Vector();
          if(!files.contains(args[i]))
            files.add(args[i]);
        }
      }
      if(files != null && files.size() > 0)
      {
        names = new String[files.size()];
        files.toArray(names);
      }
      
      for(int i=0;i < names.length; i++)
      {
        System.out.println("read :: "+names[i]+" write :: "+names[i]+suffix);
        Entry entry = ReadAndWriteEntry.readEntryFromDatabase(names[i], ENTRY_SOURCE);
        
        try
        {
          ReadAndWriteEntry.writeDatabaseEntryToFile(
            entry, new File(names[i]+suffix), flatten, ignoreObsolete, false, 
            format, null);
        }
        catch(EntryInformationException eie)
        {
          final JPanel msgPanel = new JPanel(new BorderLayout());
          msgPanel.add(new JLabel("Destination format can't handle all " +
              "keys/qualifiers - continue?"), BorderLayout.NORTH);
          JTextArea msgError = new JTextArea(eie.getMessage());
          msgError.setLineWrap(true);
          msgError.setEditable(false);
          JScrollPane scollMsg = new JScrollPane(msgError);
          msgPanel.add(scollMsg, BorderLayout.CENTER);
          
          int val = JOptionPane.showConfirmDialog(null,
              msgPanel, 
              "Keys/Qualifier", JOptionPane.OK_CANCEL_OPTION, 
              JOptionPane.QUESTION_MESSAGE);

          if(val == JOptionPane.OK_OPTION) 
          {
            ReadAndWriteEntry.writeDatabaseEntryToFile(
              entry, new File(names[i]+suffix), flatten, ignoreObsolete, true, 
              format, null);
          }
        }
      }
    }
    catch(Exception e)
    {
      e.printStackTrace();
    }
    System.exit(0);
  }
  
}