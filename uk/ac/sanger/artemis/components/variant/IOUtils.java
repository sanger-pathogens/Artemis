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
import java.util.List;

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
  /**
   * Write filtered uncompressed VCF. Uses the filter in VCFview to
   * determine if variants are written.
   * @param vcfFileName
   * @param vcfView
   * @param features
   */
  protected static void writeVCF(final String vcfFileName, 
                                 final VCFview vcfView,
                                 final FeatureVector features)
  {
    try
    {
      FileWriter writer = new FileWriter(vcfFileName+".filter");
      if(IOUtils.isBCF(vcfFileName))
      {
        BCFReader.writeVCF(writer, vcfFileName);
        return;
      }
      
      TabixReader tr = new TabixReader(vcfFileName);
      String line;
      while ((line = tr.readLine()) != null)
      {
        if(line.startsWith("#"))
        {
          writer.write(line+'\n');
          continue;
        }
        
        VCFRecord record = VCFRecord.parse(line);
        int basePosition = record.getPos() + vcfView.getSequenceOffset(record.getChrom());
        if( !vcfView.showVariant(record, features, basePosition) )
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
  
  /**
   * Export as a VCF based on the filtering applied in the VCFview.
   * @param entryGroup
   * @param vcfFiles
   * @param vcfView
   */
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
    
    String filterFiles = "";
    for(int i=0; i<vcfFiles.size(); i++)
    {
      IOUtils.writeVCF(vcfFiles.get(i), vcfView, features);
      filterFiles += vcfFiles.get(i)+".filter\n";
    }

    new MessageDialog (null, "Saved Files", filterFiles, false);
  }

  /**
   * Test if this is a BCF file.
   * @param fileName
   * @return
   * @throws IOException
   */
  protected static boolean isBCF(String fileName) throws IOException
  {
    FileInputStream fis = new FileInputStream(fileName);
    BlockCompressedInputStream is = new BlockCompressedInputStream(fis);
    byte[] magic = new byte[4];
    is.read(magic);
    fis.close();
    is.close();
    String line = new String(magic);
    if(line.equals("BCF\4"))
      return true;
    return false;
  }
  
}