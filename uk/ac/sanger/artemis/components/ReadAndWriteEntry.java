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

package uk.ac.sanger.artemis.components;

import java.io.File;
import java.io.IOException;

import uk.ac.sanger.artemis.components.database.DatabaseEntrySource;
import uk.ac.sanger.artemis.components.genebuilder.GeneUtils;
import uk.ac.sanger.artemis.io.DocumentEntryFactory;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.EntryInformationException;
import uk.ac.sanger.artemis.io.Feature;
import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.io.KeyVector;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierInfo;
import uk.ac.sanger.artemis.io.QualifierInfoException;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.sequence.NoSequenceException;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.InputStreamProgressEvent;
import uk.ac.sanger.artemis.util.InputStreamProgressListener;
import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.Options;

class ReadAndWriteEntry
{

  private static org.apache.log4j.Logger logger4j = 
    org.apache.log4j.Logger.getLogger(ReadAndWriteEntry.class);
 
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
    final DatabaseEntrySource entry_source = new DatabaseEntrySource();
    boolean promptUser = true;
    if(System.getProperty("read_only") != null)
      promptUser = false;

    if(!entry_source.setLocation(promptUser))
      return null;
    
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
   * Write entry to a file
   * @param entry
   * @param file
   * @param flatten Flatten the gene model and combine the qualifiers if true.
   *    If false it will write all features and qualifiers out.
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
                                              final boolean force,
                                              final int destination_type) 
         throws IOException, EntryInformationException
  {
    GeneUtils.lazyLoadAll(entry, null);

    EntryInformation artemis_entry_information = Options.getArtemisEntryInformation();
    if(!flatten)
    {
      final FeatureVector features = entry.getAllFeatures();
      for(int i=0; i<features.size(); i++)
        addAllKeysQualifiers(artemis_entry_information, features.elementAt(i).getEmblFeature());
    }
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
    final Key new_key = feature.getKey();
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
      Entry entry = ReadAndWriteEntry.readEntryFromDatabase("Pf3D7_03");
      ReadAndWriteEntry.writeDatabaseEntryToFile(
          entry, new File("Pf3D7_03.flatten"), true, false, 
          DocumentEntryFactory.EMBL_FORMAT);
      
      ReadAndWriteEntry.writeDatabaseEntryToFile(
          entry, new File("Pf3D7_03.not-flatten"), false, false,
          DocumentEntryFactory.EMBL_FORMAT);
    }
    catch(Exception e)
    {
      e.printStackTrace();
    }
  }
  
}