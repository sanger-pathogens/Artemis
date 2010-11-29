/* IOUtils
 *
 * created: 2010
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2010  Genome Research Limited
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
package uk.ac.sanger.artemis.components.variant;

import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.List;

import javax.swing.JOptionPane;

import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureEnumeration;
import uk.ac.sanger.artemis.FeatureKeyQualifierPredicate;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.components.MessageDialog;
import uk.ac.sanger.artemis.io.Key;

import net.sf.samtools.util.BlockCompressedInputStream;

class IOUtils
{

  private static void writeHeader(final String fileName, final Writer writer)
  {
    try
    {
      FileInputStream fileStream = new FileInputStream(fileName);
      BlockCompressedInputStream is = new BlockCompressedInputStream(fileStream);

      int c;
      int lc = -1;
      
      while((c = is.read()) >= 0)
      {
        if(lc == '\n' && c != '#')
          break;
        writer.write(c);
      }
      is.close();
    }
    catch (IOException e)
    {
      e.printStackTrace();
    }
  }
  
  /**
   * Write filtered uncompressed VCF. Uses the filter in VCFview to
   * determine if variants are written.
   * @param vcfFileName
   * @param writer
   * @param vcfView
   * @param features
   */
  protected static void writeVCF(final String vcfFileName, 
                                 final Writer writer, 
                                 final VCFview vcfView,
                                 final FeatureVector features)
  {
    writeHeader(vcfFileName, writer);
    try
    {
      TabixReader tr = new TabixReader(vcfFileName);
      String line;
      while ((line = tr.readLine()) != null)
      {
        if(line.startsWith("#"))
          continue;
        String parts[] = VCFview.tabPattern.split(line, 0);
        int basePosition = Integer.parseInt(parts[1]) + vcfView.getSequenceOffset(parts[0]);
        if( !vcfView.showVariant(parts[3], parts[4], features, basePosition, parts[5]) )
          continue;
        writer.write(line+'\n');
      }
      writer.close();
    }
    catch (IOException e)
    {
      e.printStackTrace();
    }
  }
  
  protected static void export(final EntryGroup entryGroup, 
                               final List<String> vcfFiles,
                               final VCFview vcfView)
  {
    // select all CDS features that do not have the /pseudo qualifier
    final FeatureKeyQualifierPredicate predicate =
      new FeatureKeyQualifierPredicate(Key.CDS, "pseudo", false);
    
    final FeatureVector features = new FeatureVector ();
    final FeatureEnumeration feature_enum = entryGroup.features ();

    while (feature_enum.hasMoreFeatures ())
    {
      final Feature current_feature = feature_enum.nextFeature ();
      if (predicate.testPredicate (current_feature)) 
        features.add (current_feature);
    }
    
    try
    {
      String filterFiles = "";
      for(int i=0; i<vcfFiles.size(); i++)
      {

        FileWriter writer = new FileWriter(vcfFiles.get(i)+".filter");
        IOUtils.writeVCF(vcfFiles.get(i), writer, vcfView, features);
        filterFiles += vcfFiles.get(i)+".filter\n";
      }

      new MessageDialog (null, "Saved Files", filterFiles, false);
    }
    catch (IOException e1)
    {
      JOptionPane.showMessageDialog(vcfView, e1.getMessage(),
          "Error", JOptionPane.ERROR_MESSAGE);
    }
  }
  
}