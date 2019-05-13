/* DocumentEntryFactory.java
 *
 * created: Sat Sep 11 1999
 *
 * This file is part of Artemis
 *
 * Copyright (C) 1999  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/DocumentEntryFactory.java,v 1.3 2007-04-11 13:47:09 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import uk.ac.sanger.artemis.util.Document;
import uk.ac.sanger.artemis.util.LinePushBackReader;

import java.io.File;
import java.io.IOException;

import javax.swing.JOptionPane;


/**
 *  This class contains the method makeDocumentEntry (), which creates a
 *  DocumentEntry object of the appropriate class from a Document object.
 *
 *  @author Kim Rutherford
 *  @version $Id: DocumentEntryFactory.java,v 1.3 2007-04-11 13:47:09 tjc Exp $
 **/

abstract public class DocumentEntryFactory 
{

  /** use if format of the entry is unknown   */
  final public static int UNKNOWN_FORMAT = 0;

  /** use if format of the entry is not important */
  final public static int ANY_FORMAT = UNKNOWN_FORMAT;

  /** use for an entry that is in EMBL format */
  final public static int EMBL_FORMAT = 1;

  /** use for an entry that is in GENBANK format */
  final public static int GENBANK_FORMAT = 2;

  /** use for an entry that is in GFF format  */
  final public static int GFF_FORMAT = 3;

  public static boolean REMOVE_PRODUCT_FROM_PSEUDOGENE = false;
  
  /** 
   * Boolean indicating whether or not the user should 
   * be asked if they want to use a GFF index. 
   */
  private static boolean displayIndexingQuestionForGffs = true;


  /**
   *  Read a DocumentEntry object from the given Document.
   *  @param entry_information The EntryInformation to use when reading.  This
   *    supplies the list of valid keys and qualifiers
   *  @param listener The object that will listen for ReadEvents.
   *  @exception EntryInformationException Thrown if an Entry using the given
   *    EntryInformation object cannot contain the Key, Qualifier or
   *    Key/Qualifier combination of one of the features in the Document.
   **/
  public static DocumentEntry makeDocumentEntry (final EntryInformation entry_information,
                                                 final Document document,
                                                 final ReadListener listener)
      throws IOException, EntryInformationException 
  {
    if(document.getInputStream() instanceof htsjdk.samtools.util.BlockCompressedInputStream)
    {
      if(IndexedGFFDocumentEntry.isIndexed( ((File)document.getLocation()) ))
      {
    	int sel = 0;
    	
    	if (getDisplayIndexingQuestionForGffs())
    	{
          Object[] possibleValues = { "Use index", "Ignore index" };
          sel = JOptionPane.showOptionDialog(null, 
              "Use the GFF index file (to increase the performance)\n"+ 
              "or concatenate the sequences together?",
              "Indexed GFF", 
              JOptionPane.DEFAULT_OPTION, JOptionPane.WARNING_MESSAGE,
              null, possibleValues, possibleValues[0]);
    	}
    	
        if(sel == 0)
          return new IndexedGFFDocumentEntry(document);
      }
    }
    
    final LinePushBackReader document_reader =
                        document.getLinePushBackReader();

    final String first_line = document_reader.readLine();

    if (first_line == null)  // empty file - create an empty EmblDocumentEntry
      return new EmblDocumentEntry(entry_information, document, listener);

    final int first_line_type = LineGroup.getLineType(first_line);
    document_reader.pushBack(first_line);

    switch (first_line_type) {
    case LineGroup.GFF_FEATURE:
    case LineGroup.GFF_MISC:
      {
        final SimpleDocumentEntry document_entry =
          new GFFDocumentEntry (document, listener);
        return document_entry;
      }
    case LineGroup.MSPCRUNCH_FEATURE:
      {
        final SimpleDocumentEntry document_entry =
          new MSPcrunchDocumentEntry (document, listener);
        return document_entry;
      }
    case LineGroup.BLAST_FEATURE:
      {
        final SimpleDocumentEntry document_entry =
          new BlastDocumentEntry (document, listener);
        return document_entry;
      }
    case LineGroup.GENBANK_MISC:
    case LineGroup.GENBANK_FEATURE:
      {
        final SimpleDocumentEntry document_entry =
          new GenbankDocumentEntry (entry_information, document, listener);
        return document_entry;
      }
//    case LineGroup.GAME_XML:
//      return new BioJavaEntry (document,
//                               new GAMEFormat ());
    default:
      {
        final SimpleDocumentEntry document_entry =
          new EmblDocumentEntry (entry_information, document, listener);
        return document_entry;
      }
    }
  }

  /**
   *  Make a new (nameless) DocumentEntry from the given Entry.
   *  @param entry_information The EntryInformation to use for the new object.
   *    This supplies the list of valid keys and qualifiers.  Note that this
   *    argument is ignored by some DocumentEntry types, because some formats
   *    (like GFF) have limitations on the possible keys and qualifiers.
   *  @param destination_type This parameter control the type of DocumentEntry
   *    that is created.  It should be a DocumentEntry type that can be
   *    constructed from any Entry eg. one of EMBL_FORMAT, GENBANK_FORMAT.
   *  @param force If true then invalid qualifiers and any features with
   *    invalid keys in the new Entry will be quietly thrown away.  "Invalid"
   *    means that the key/qualifier is not allowed to occur in an Entry of
   *    this type (probably determined by the EntryInformation object of this
   *    Entry).  If false an EntryInformationException will be thrown for
   *    invalid keys or qualifiers.
   *  @exception EntryInformationException Thrown if an Entry using the given
   *    EntryInformation object cannot contain the Key, Qualifier or
   *    Key/Qualifier combination of one of the features in the source Entry.
   **/
  public static DocumentEntry makeDocumentEntry (final EntryInformation
                                                   entry_information,
                                                 final Entry entry,
                                                 int destination_type,
                                                 final boolean force)
      throws EntryInformationException 
  {

    if(destination_type == ANY_FORMAT)  
    {
      if (entry instanceof EmblDocumentEntry)
        destination_type = EMBL_FORMAT;
      else
      {
        if (entry instanceof GenbankDocumentEntry) 
          destination_type = GENBANK_FORMAT;
        else
        {
          if (entry instanceof GFFDocumentEntry) 
            destination_type = GFF_FORMAT;
          else
            destination_type = EMBL_FORMAT;
        }
      }
    }

    switch (destination_type) {
    case EMBL_FORMAT:
    {
      EmblDocumentEntry ee = new EmblDocumentEntry (entry_information, entry, force);
      if(force)
        removeProductFromPseudogene(ee, destination_type);
      return ee;
    }
    case GENBANK_FORMAT:
      return new GenbankDocumentEntry (entry_information, entry, force);
    case GFF_FORMAT:
      return new GFFDocumentEntry (entry, force);
//     case BSML_FORMAT:
//     case GAME_FORMAT:
//     case AGAVE_FORMAT:
//     ...
    default:
      throw new Error ("internal error - unknown DocumentEntry type");
    }
  }
  
  /**
   * Remove product qualifier from CDS features with pseudogene qualifier
   * @param document_entry
   * @param destination_type
   * @throws EntryInformationException
   */
  private static void removeProductFromPseudogene(final DocumentEntry entry,
                                                  final int destination_type) 
      throws EntryInformationException
  {
    if( !REMOVE_PRODUCT_FROM_PSEUDOGENE )
      return;

    final FeatureVector features = entry.getAllFeatures();
    for(Feature f: features)
    {
      if(f.getKey().equals("CDS"))
      {
        Qualifier q = f.getQualifierByName("product");
        if(  q != null && 
            (f.getQualifierByName("pseudogene") != null || 
             f.getQualifierByName("pseudo") != null))
          f.getQualifiers().remove(q);
      }
    }
  }

  /**
   * Get the value of the displayIndexingQuestionForGffs property.
   * @return boolean
   */
  public static boolean getDisplayIndexingQuestionForGffs()
  {
	return displayIndexingQuestionForGffs;
  }

  /**
   * Set the value of the displayIndexingQuestionForGffs property.
   * @param displayIndexingQuestionForGffs boolean
   */
  public static void setDisplayIndexingQuestionForGffs(boolean displayIndexingQuestionForGffs)
  {
	DocumentEntryFactory.displayIndexingQuestionForGffs = displayIndexingQuestionForGffs;
  }
}
