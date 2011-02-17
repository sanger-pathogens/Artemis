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
import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.sequence.Bases;

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
      File filterFile = getFile(vcfFileName, nfiles, "filter");
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
  
  private static File getFile(final String vcfFileName, final int nfiles, final String suffix) throws IOException
  {
    if(nfiles > 1)
      return new File(vcfFileName+suffix);
    final StickyFileChooser file_dialog = new StickyFileChooser();

    file_dialog.setSelectedFile(new File(vcfFileName+suffix));
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
    exportFasta(entryGroup, vcfReaders, chr, vcfView, vcf_v4, features); 
  }

  
  protected static void exportFasta(final EntryGroup entryGroup,
                                    final AbstractVCFReader vcfReaders[],
                                    final String chr,
                                    final VCFview vcfView,
                                    final boolean vcf_v4,
                                    final FeatureVector features)
  {
    String suffix = ".fasta";
    if(features.size() == 1)
      suffix = "."+features.elementAt(0).getIDString()+suffix;
    
    for (int i = 0; i < vcfReaders.length; i++)
    {
      String vcfFileName = vcfReaders[i].getFileName();
      try
      {
        File filterFile = getFile(vcfFileName, vcfReaders.length, suffix);
        FileWriter writer = new FileWriter(filterFile);
        for (int j = 0; j < features.size(); j++)
        {
          Feature f = features.elementAt(j);
          FeatureSegmentVector segs = f.getSegments();
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
              VCFRecord bcfRecord;
              while ((bcfRecord = it.next()) != null)
                segBases = getSeqsVariation(bcfRecord, segBases, sbeg, f.isForwardFeature(), vcf_v4);
            }
            buff.append(segBases);
          }
          
          StringBuffer header = new StringBuffer(f.getSystematicName());
          header.append(" "+f.getIDString()+" ");
          final String product = f.getProductString();
          header.append( (product == null ? "undefined product" : product) );
          header.append(" ").append(f.getWriteRange());

          writeSequence(writer, header.toString(), buff.toString()); 
        }
        
        writer.close();
      }
      catch (IOException e)
      {
        e.printStackTrace();
      }
    }
  }
  
  private static void writeSequence(FileWriter writer, String header, String bases) throws IOException
  {
    writer.write (">" + header + "\n");

    final int SEQUENCE_LINE_BASE_COUNT = 60;
    for(int k=0; k<bases.length(); k+=SEQUENCE_LINE_BASE_COUNT)
    {
      int end = k + SEQUENCE_LINE_BASE_COUNT;
      if(end > bases.length())
        end = bases.length();
      writer.write ( bases.substring(k,end) + "\n");
    }
  }
  
  private static String getSeqsVariation(VCFRecord vcfRecord, 
      String bases, int sbeg, boolean isFwd, boolean vcf_v4)
  {
    int position = vcfRecord.getPos()-sbeg;
    if(!isFwd)
      position = bases.length()-position-1;
    
    if(position > bases.length())
      return bases;
    else if(position < 0)
      return bases;

    if(position < bases.length()-1 && bases.charAt(position) == '-')
      return bases;
    
    StringBuffer buff = new StringBuffer();
    buff.append(bases.substring(0,position)); 
    
    if(vcfRecord.getAlt().isDeletion(vcf_v4))
    {
      int ndel = vcfRecord.getAlt().getNumberOfDeletions(vcf_v4);
      
      if(isFwd)
        position+=ndel;
      else
        buff.delete(position-ndel+1, position);
      
      for(int i=0; i<ndel; i++)
        buff.append("-");
    }
    else if(vcfRecord.getAlt().isInsertion(vcf_v4))
    {
      
    }
    else if(vcfRecord.getAlt().isMultiAllele())
    {
      
    }
    else if(vcfRecord.getAlt().isNonVariant())                   // non-variant
    {
      if(isFwd)
        buff.append(vcfRecord.getRef().toUpperCase());
      else
        buff.append(Bases.complement(vcfRecord.getRef()).toUpperCase());
    }
    else
    {
      String alt = vcfRecord.getAlt().toString().toLowerCase();  // SNP   
      if(isFwd)
        buff.append(alt);
      else
        buff.append(Bases.complement(alt));
    }
    
    if(isFwd && position < bases.length())
      buff.append(bases.substring(position+1));
    else if(!isFwd)
    {
      buff.append(bases.substring(position+1));
    }
    return buff.toString();
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