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

    boolean isStartPartial = false;
    boolean isEndPartial   = false;
    if(feature.getEmblFeature() instanceof GFFStreamFeature)
    {
      try
      {
        if(feature.getQualifierByName("Start_range") != null)
          isStartPartial = true;
        if(feature.getQualifierByName("End_range") != null)
          isEndPartial = true;
      }
      catch (InvalidRelationException e){}
    }
    else
    {
      if(feature.getLocation().isPartial(true)) // 5prime
      {
        if(feature.isForwardFeature()) 
          isStartPartial = true;
        else
          isEndPartial = true;
      }
      if(feature.getLocation().isPartial(false)) // 3prime
      {
        if(feature.isForwardFeature())
          isEndPartial = true;
        else
          isStartPartial = true;
      }
    }

    writer.write( getPositionsStr(feature, ranges.elementAt(0), 
        isStartPartial, isEndPartial) );
    writer.write("\t"+feature.getKey().getKeyString());
    
    for(int j=1; j<ranges.size(); j++)
    {
      writer.write("\n");
      writer.write( getPositionsStr(feature, ranges.elementAt(j), 
          isStartPartial, isEndPartial) );
    }
  }
  
  /**
   * Get the start and end positions as a tab delimited string
   * @param feature
   * @param r
   * @param isStartRangePartial
   * @param isEndRangePartial
   * @return
   */
  private static String getPositionsStr(
      final Feature feature, 
      final Range r,
      final boolean isStartPartial,
      final boolean isEndPartial)
  {
    boolean firstBase;
    boolean lastBase;
    final int low_marker  = feature.getFirstBaseMarker().getPosition();
    final int high_marker = feature.getLastBaseMarker().getPosition();
    String low_pos;
    String high_pos;
    if(feature.isForwardFeature())
    {
      firstBase = (r.getStart() == low_marker);
      lastBase  = (r.getEnd()   == high_marker);
      low_pos   = Integer.toString(r.getStart());
      high_pos  = Integer.toString(r.getEnd());
    }
    else
    {
      firstBase = (r.getEnd()   == high_marker);
      lastBase  = (r.getStart() == low_marker);
      low_pos   = Integer.toString(r.getEnd());
      high_pos  = Integer.toString(r.getStart());
    }

    // set partials
    if(firstBase && isStartPartial)
      low_pos  = (feature.isForwardFeature() ? "<" : ">")+low_pos;
    if(lastBase && isEndPartial)
      high_pos = (feature.isForwardFeature() ? ">" : "<")+high_pos;
    return low_pos+"\t"+high_pos;
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
        {
          if(current_qualifier.getName().equals("pseudo"))
            writer.write("\n\t\t\tpseudo");
          continue;
        }
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