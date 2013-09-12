/* 
 * This file is part of Artemis
 *
 * Copyright (C) 2013  Genome Research Limited
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
package uk.ac.sanger.artemis.io;

import java.io.IOException;
import java.net.URL;

import junit.framework.Assert;
import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureSegmentVector;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.Selection;
import uk.ac.sanger.artemis.SimpleEntryGroup;
import uk.ac.sanger.artemis.io.DocumentEntryFactory;
import uk.ac.sanger.artemis.io.Entry;
import uk.ac.sanger.artemis.io.EntryInformationException;
import uk.ac.sanger.artemis.sequence.AminoAcidSequence;
import uk.ac.sanger.artemis.sequence.Bases;
import uk.ac.sanger.artemis.sequence.NoSequenceException;
import uk.ac.sanger.artemis.util.Document;
import uk.ac.sanger.artemis.util.DocumentFactory;
import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.util.StringVector;

public class Utils
{
  public static Entry getEntry(final String fileName)
  {
    try
    {
      URL entryFile = Utils.class.getResource(fileName);
      final Document doc = DocumentFactory.makeDocument(entryFile.getFile());
      return DocumentEntryFactory.makeDocumentEntry(
          Options.getArtemisEntryInformation(),doc,null);
    }
    catch(EntryInformationException e) 
    {
      Assert.fail(e.getMessage());
    }
    catch(IOException e) 
    {
      Assert.fail(e.getMessage());
    }
    return null;
  }
  
  public static EntryGroup getEntryGroup(final String fileName)
  {
    try
    {
      final Entry new_embl_entry = getEntry(fileName);
      final uk.ac.sanger.artemis.Entry entry = new  uk.ac.sanger.artemis.Entry(new_embl_entry);
      final EntryGroup egrp = new SimpleEntryGroup(entry.getBases());
      egrp.add(entry);
      return egrp;
    }
    catch(OutOfRangeException e)
    {
      Assert.fail(e.getMessage());
    }
    catch(NoSequenceException e)
    {
      Assert.fail(e.getMessage());
    }
    return null;
  }
  
  /**
   * Method to change the translation table being used
   * @param n - genetic code table number
   */
  protected static void changeTranslationTable(String n)
  {
    StringVector options_file_table =
                 Options.getOptions().getOptionValues("translation_table_"+n);
    StringVector table = Options.getOptions().getOptionValues("translation_table_1");
    for(String cod_plus_aa: options_file_table)
    {
      final int codon_index = Bases.getIndexOfBase(cod_plus_aa.charAt(0)) * 16 +
                              Bases.getIndexOfBase(cod_plus_aa.charAt(1)) * 4 +
                              Bases.getIndexOfBase(cod_plus_aa.charAt(2));
      table.setElementAt(cod_plus_aa.substring(3), codon_index);
    }
        
    StringBuffer sbuff = new StringBuffer();
    for(int i = 0; i < 64; ++i)
      sbuff.append(table.elementAt(i)+" ");
    Options.getOptions().setGeneticCode(sbuff.toString());
      
    options_file_table =
          Options.getOptions().getOptionValues("start_codons_"+n);
    sbuff = new StringBuffer();
    for(String start: options_file_table)
      sbuff.append(start+" ");
        
    Options.getOptions().setProperty("start_codons",sbuff.toString());   
    AminoAcidSequence.setGeneCode();
  }
}