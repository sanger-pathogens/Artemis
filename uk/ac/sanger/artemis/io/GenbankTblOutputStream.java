/* GenbankTblOutputStream.java
 *
 * created: 2008
 *
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

import java.io.File;
import java.io.IOException;
import java.io.Writer;

import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JOptionPane;

import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.util.FileDocument;
import uk.ac.sanger.artemis.util.StringVector;


/**
 *  Handle writing tbl format:
 *  http://www.ncbi.nlm.nih.gov/Sequin/table.html
 */
public class GenbankTblOutputStream
{
  
  /**
   * Write out an entry as tbl format.
   * @param entry
   * @param f
   */
  public static void writeEntryAsTbl(final Entry entry, final JFrame f)
  {
    final JFileChooser chooser = new JFileChooser();
    chooser.setMultiSelectionEnabled(false);
    int returnVal = chooser.showSaveDialog(f);
    if(returnVal == JFileChooser.CANCEL_OPTION)
      return;

    final File fileSave = chooser.getSelectedFile();
    if(fileSave.exists())
    {
      returnVal = JOptionPane.showConfirmDialog(f, 
                          "File exists. Overwrite?", 
                          "File Exists",
                          JOptionPane.YES_NO_OPTION);
      if(returnVal == JOptionPane.CANCEL_OPTION ||
         returnVal == JOptionPane.NO_OPTION)
        return;
    }
    
    final FileDocument fileDocument = new FileDocument(fileSave);
    try
    {
      final Writer writer = fileDocument.getWriter();
      writer.write(">Feature "+entry.getName()+"\n");
      
      final FeatureVector features = entry.getAllFeatures();
      final EntryInformation entry_information = entry.getEntryInformation ();
      
      int count = 0;
      for(int i=0; i<features.size(); i++)
      {
        final Feature feature = features.elementAt(i);
        if(feature.getKey().getKeyString().equals("source"))
          continue;
        if(count > 0)
          writer.write("\n");
        
        count++;
        writeRanges(feature, writer);
        writeQualifiers(feature, entry_information, writer);
      }
      writer.close();
    }
    catch(IOException e1)
    {
      JOptionPane.showMessageDialog(f, e1.getMessage());
    }
  }
  
  /**
   * Write out ranges and feature key
   * @param feature
   * @param writer
   * @throws IOException
   */
  private static void writeRanges(final Feature feature, final Writer writer) 
          throws IOException
  {
    final RangeVector ranges = feature.getLocation().getRanges();
    if(!feature.isForwardFeature() && ranges.size() > 1)
      ranges.reverse();
    
    Range r = (Range)ranges.elementAt(0);
    if(feature.isForwardFeature())
      writer.write(r.getStart()+"\t"+r.getEnd());
    else
      writer.write(r.getEnd()+"\t"+r.getStart());
    
    writer.write("\t"+feature.getKey().getKeyString());
    
    for(int j=1; j<ranges.size(); j++)
    {
      writer.write("\n");
      r = (Range)ranges.elementAt(j);
      if(feature.isForwardFeature())
        writer.write(r.getStart()+"\t"+r.getEnd());
      else
        writer.write(r.getEnd()+"\t"+r.getStart());
    }
  }
  
  /**
   * Write out qualifiers
   * @param feature
   * @param entry_information
   * @param writer
   * @throws IOException 
   */
  private static void writeQualifiers(final Feature feature, 
                                      final EntryInformation entry_information,
                                      final Writer writer) throws IOException
  {
    QualifierVector qualifiers = feature.getQualifiers();
    for(int j=0; j<qualifiers.size(); j++)
    {
      Qualifier current_qualifier = (Qualifier) qualifiers.elementAt(j);

      if(feature.getKey().getKeyString().equals("CDS") &&
         current_qualifier.getName ().equals("translation"))
        continue;
      
      final QualifierInfo qualifier_info =
        entry_information.getQualifierInfo (current_qualifier.getName ());

      // this will contain one String for each /name=value pair
      final StringVector qualifier_strings =
        StreamQualifier.toStringVector (qualifier_info, current_qualifier);
      
      for(int k=0; k<qualifier_strings.size(); k++)
      {
        String qualifierStr = (String)qualifier_strings.get(k);
        int index = qualifierStr.indexOf("=");
        if(index < 1)
          continue;
        qualifierStr = qualifierStr.substring(index+1);
        if(qualifierStr.startsWith("\""))
          qualifierStr = qualifierStr.substring(1);
        if(qualifierStr.endsWith("\""))
          qualifierStr = qualifierStr.substring(0, qualifierStr.length()-1);
        
        writer.write("\n\t\t\t"+current_qualifier.getName()+"\t"+
            qualifierStr);
      }
    }
  }
}