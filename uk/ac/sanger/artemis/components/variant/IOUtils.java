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

import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.List;

import javax.swing.JFileChooser;

import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureEnumeration;
import uk.ac.sanger.artemis.FeatureKeyQualifierPredicate;
import uk.ac.sanger.artemis.FeaturePredicate;
import uk.ac.sanger.artemis.FeatureSegment;
import uk.ac.sanger.artemis.FeatureSegmentVector;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.components.MessageDialog;
import uk.ac.sanger.artemis.components.StickyFileChooser;
import uk.ac.sanger.artemis.components.variant.BCFReader.BCFReaderIterator;
import uk.ac.sanger.artemis.io.FastaStreamSequence;
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
  protected static File writeVCF(final String vcfFileName, 
                                 final VCFview vcfView,
                                 final FeatureVector features,
                                 final int nfiles)
  {
    try
    {
      File filterFile = getFile(vcfFileName, nfiles);
      FileWriter writer = new FileWriter(filterFile);
      if(IOUtils.isBCF(vcfFileName))
      {
        BCFReader.writeVCF(writer, vcfFileName);
        return filterFile;
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
      return filterFile;
    }
    catch (IOException e)
    {
      e.printStackTrace();
      return null;
    }
  }
  
  private static File getFile(final String vcfFileName, final int nfiles) throws IOException
  {
    if(nfiles > 1)
      return new File(vcfFileName+".filter");
    final StickyFileChooser file_dialog = new StickyFileChooser();

    file_dialog.setSelectedFile(new File(vcfFileName+".filter"));
    file_dialog.setDialogTitle("Choose save file ...");
    file_dialog.setDialogType(JFileChooser.SAVE_DIALOG);
    final int status = file_dialog.showSaveDialog(null);

    if(status != JFileChooser.APPROVE_OPTION ||
       file_dialog.getSelectedFile() == null) 
      return null;
    
    return file_dialog.getSelectedFile();
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
    // get all CDS features that do not have the /pseudo qualifier
    final FeatureVector features = getFeatures(
        new FeatureKeyQualifierPredicate(Key.CDS, "pseudo", false), entryGroup);
    
    String filterFiles = "";
    for(int i=0; i<vcfFiles.size(); i++)
    {
      File filterFile = IOUtils.writeVCF(vcfFiles.get(i), vcfView, features, vcfFiles.size());
      filterFiles += filterFile.getAbsolutePath()+"\n";
    }

    new MessageDialog (null, "Saved Files", filterFiles, false);
  }

  protected static void exportFasta(final EntryGroup entryGroup,
                                    final AbstractVCFReader vcfReaders[],
                                    final String chr,
                                    final VCFview vcfView,
                                    final boolean vcf_v4)
  {
    // get all CDS features that do not have the /pseudo qualifier
    final FeatureVector features = getFeatures(
        new FeatureKeyQualifierPredicate(Key.CDS, "pseudo", false), entryGroup);
    
    for (int i = 0; i < vcfReaders.length; i++)
    {
      String vcfFileName = vcfReaders[i].getFileName();
      try
      {
        File filterFile = getFile(vcfFileName, vcfReaders.length);
        FileWriter writer = new FileWriter(filterFile);
        for (int j = 0; j < features.size(); j++)
        {
          Feature f = features.elementAt(j);
          FeatureSegmentVector segs = f.getSegments();
          Hashtable<String, String> seqs = new Hashtable<String, String>();
          String bases = f.getBases();
          seqs.put("ref", bases);
          
          StringBuffer buff = new StringBuffer();
          
          for(int k=0; k<segs.size(); k++)
          {
            FeatureSegment seg = segs.elementAt(k);
            int sbeg = seg.getRawRange().getStart();
            int send = seg.getRawRange().getEnd();
            String segBases = seg.getBases();
            
            if (vcfReaders[i] instanceof BCFReader)
            {
              BCFReaderIterator it = ((BCFReader) vcfReaders[i]).query(chr, sbeg, send);
              VCFRecord bcfRecord = null;

              while ((bcfRecord = it.next()) != null)
                segBases = getSeqsVariations(bcfRecord, segBases, sbeg, f.isForwardFeature(), vcf_v4);
            }
            
            buff.append(segBases);
          }
 
          Enumeration<String> en = seqs.keys();
          while(en.hasMoreElements())
          {
            String key = en.nextElement();
            bases = seqs.get(key);
            FastaStreamSequence stream = new FastaStreamSequence(bases, f.getIDString()+":"+key);
            stream.writeToStream(writer);
          }
        }
        
        writer.close();
      }
      catch (IOException e)
      {
        e.printStackTrace();
      }
    }
  }
  
  private static int getNumberOfDeletions(VCFRecord vcfRecord, boolean vcf_v4)
  {
    if(vcf_v4)
      return vcfRecord.getRef().length()-vcfRecord.getAlt().length();
    
    int index = vcfRecord.getAlt().indexOf("D");
    int ndel = 0;
    try
    {
      ndel = Integer.parseInt( vcfRecord.getAlt().substring(index+1) );
    }
    catch(NumberFormatException e) { e.printStackTrace(); }
    return ndel;
  }
  
  protected static String getSeqsVariations(VCFRecord vcfRecord, 
      String bases, int sbeg, boolean isFwd, boolean vcf_v4)
  {   
    int position = vcfRecord.getPos()-sbeg;
    if(!isFwd)
      position = bases.length()-position;
    
    if(position > bases.length())
      return bases;

    if(vcfRecord.isDeletion(vcf_v4))
    {
      int ndel = getNumberOfDeletions(vcfRecord, vcf_v4);
      
      if(isFwd)
        return bases.substring(0,position)+bases.substring(position+ndel);
      else
        return bases.substring(0,position-ndel)+bases.substring(position);
    }
    else if(vcfRecord.isInsertion(vcf_v4))
    {
      
    }
    else
    {
      
    }
    return bases;
  }

  /**
   * Get all features in an entry group.
   * @param predicate
   * @param entryGroup
   * @return
   */
  private static FeatureVector getFeatures(FeaturePredicate predicate, EntryGroup entryGroup)
  {
    final FeatureVector features = new FeatureVector ();
    final FeatureEnumeration feature_enum = entryGroup.features ();
    while (feature_enum.hasMoreFeatures ())
    {
      final Feature current_feature = feature_enum.nextFeature ();
      if (predicate.testPredicate (current_feature)) 
        features.add (current_feature);
    }
    return features;
  }

  /**
   * Test if this is a BCF file.
   * @param fileName
   * @return
   * @throws IOException
   */
  protected static boolean isBCF(String fileName) throws IOException
  {
    InputStream ins;
    if(fileName.startsWith("http:"))
    {
      final URL urlBamIndexFile = new URL(fileName);
      ins = urlBamIndexFile.openStream();
    }
    else
      ins = new FileInputStream(fileName);
    BlockCompressedInputStream is = new BlockCompressedInputStream(ins);
    byte[] magic = new byte[4];
    is.read(magic);
    ins.close();
    is.close();
    String line = new String(magic);
    if(line.equals("BCF\4"))
      return true;
    return false;
  }
  
}