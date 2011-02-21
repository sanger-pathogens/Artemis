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

import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureEnumeration;
import uk.ac.sanger.artemis.FeatureKeyQualifierPredicate;
import uk.ac.sanger.artemis.FeaturePredicate;
import uk.ac.sanger.artemis.FeatureSegment;
import uk.ac.sanger.artemis.FeatureSegmentVector;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.Selection;
import uk.ac.sanger.artemis.components.FileViewer;
import uk.ac.sanger.artemis.components.MessageDialog;
import uk.ac.sanger.artemis.components.SequenceViewer;
import uk.ac.sanger.artemis.components.StickyFileChooser;
import uk.ac.sanger.artemis.components.variant.BCFReader.BCFReaderIterator;
import uk.ac.sanger.artemis.io.DocumentEntry;
import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.sequence.Bases;
import uk.ac.sanger.artemis.sequence.MarkerRange;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.FileDocument;
import uk.ac.sanger.artemis.util.RemoteFileDocument;

import net.sf.samtools.util.BlockCompressedInputStream;

class IOUtils
{ 
  
  private static final int MAXIMUM_SELECTED_FEATURES = 25;
  
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

  /**
   * Write out fasta for a selected base range
   * @param entryGroup
   * @param vcfReaders
   * @param chr
   * @param vcfView
   * @param vcf_v4
   * @param selection
   * @param view
   */
  protected static void exportFastaByRange(
                                    final EntryGroup entryGroup,
                                    final AbstractVCFReader vcfReaders[],
                                    final String chr,
                                    final VCFview vcfView,
                                    final boolean vcf_v4,
                                    final Selection selection,
                                    final boolean view)
  {
    MarkerRange marker = selection.getMarkerRange();
    Range range = marker.getRawRange();
    int direction = ( marker.isForwardMarker() ? Bases.FORWARD : Bases.REVERSE);
    FeatureVector features = entryGroup.getAllFeatures();
    FileWriter writer = null;
    String fastaFiles = "";
    
    String name = entryGroup.getActiveEntries().elementAt(0).getName();
    int sbeg = range.getStart();
    int send = range.getEnd();

    StringBuffer buffSeq = null;
    try
    {
      
      if(!view)
      {
        File newfile = new File(
           getBaseDirectoryFromEntry(entryGroup.getActiveEntries().elementAt(0)),
               name);
      
        File f = getFile(newfile.getAbsolutePath(), vcfReaders.length, ".fasta");
        writer = new FileWriter(f);
        fastaFiles += f.getAbsolutePath()+"\n";
      }
      else
        buffSeq = new StringBuffer();

      for (int i = 0; i < vcfReaders.length; i++)
      {
        String basesStr = entryGroup.getBases().getSubSequence(marker.getRange(), direction);
        if (vcfReaders[i] instanceof BCFReader)
        {
          BCFReaderIterator it = ((BCFReader) vcfReaders[i]).query(chr, sbeg, send);
          VCFRecord record;
          while ((record = it.next()) != null)
          {
            int basePosition = record.getPos() + vcfView.getSequenceOffset(record.getChrom());
            if(vcfView.showVariant(record, features, basePosition) )
              basesStr = getSeqsVariation(record, basesStr, sbeg, marker.isForwardMarker(), vcf_v4);
          }
        }
          
        StringBuffer header = new StringBuffer(name+" ");
        header.append(sbeg+":"+send+ (marker.isForwardMarker() ? "" : " reverse"));
        header.append(" (").append(vcfReaders[i].getName()).append(")");

        if(view) // sequence viewer
        {
          buffSeq.append(">");
          buffSeq.append(header.toString());
          buffSeq.append("\n");
          buffSeq.append(basesStr);
          buffSeq.append("\n");
        }
        else    // write to file
          writeSequence(writer, header.toString(), basesStr);
      }
      
      if(writer != null)
        writer.close();
    }
    catch(IOException e)
    {
      e.printStackTrace();
    }
    
    if(!view)
      new MessageDialog (null, "Saved Files", fastaFiles, false);
    else
    {
      FileViewer viewer = new FileViewer ("Feature base viewer for selected range: " + 
          sbeg+":"+send+(marker.isForwardMarker() ? "" : " reverse"), true);
      viewer.getTextPane().setText(buffSeq.toString());
    }
  }

  
  protected static void exportFasta(final AbstractVCFReader vcfReaders[],
                                    final String chr,
                                    final VCFview vcfView,
                                    final boolean vcf_v4,
                                    final FeatureVector features,
                                    final boolean view)
  {
    if(view && features.size () > MAXIMUM_SELECTED_FEATURES)
      new MessageDialog (null,
                        "warning: only viewing the sequences for " +
                        "the first " + MAXIMUM_SELECTED_FEATURES +
                        " selected features");
    
    String suffix = ".fasta";
    if(features.size() == 1)
      suffix = "."+features.elementAt(0).getIDString()+suffix;
    
    FileWriter writer = null;
    String fastaFiles = "";
    
    for (int i = 0; i < vcfReaders.length; i++)
    {
      try
      {
        if(!view)
        {
          File f = getFile(vcfReaders[i].getFileName(), vcfReaders.length, suffix);
          writer = new FileWriter(f);
          fastaFiles += f.getAbsolutePath()+"\n";
        }
        
        for (int j = 0; j < features.size() && (!view || j < MAXIMUM_SELECTED_FEATURES); j++)
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
              VCFRecord record;
              while ((record = it.next()) != null)
              {
                int basePosition = record.getPos() + vcfView.getSequenceOffset(record.getChrom());
                if(vcfView.showVariant(record, features, basePosition) )
                  segBases = getSeqsVariation(record, segBases, sbeg, f.isForwardFeature(), vcf_v4);
              }
            }
            buff.append(segBases);
          }

          StringBuffer header = new StringBuffer(f.getSystematicName());
          header.append(" "+f.getIDString()+" ");
          final String product = f.getProductString();
          header.append( (product == null ? "undefined product" : product) );
          header.append(" ").append(f.getWriteRange());
          header.append(" (").append(vcfReaders[i].getName()).append(")");

          if(view) // sequence viewer
          {
            SequenceViewer viewer =
              new SequenceViewer ("Feature base viewer for feature:" + 
                f.getIDString (), false);  
            viewer.setSequence(">"+header.toString(), buff.toString());
          }
          else    // write to file
            writeSequence(writer, header.toString(), buff.toString());
        }
        
        if(writer != null)
          writer.close();
      }
      catch (IOException e)
      {
        e.printStackTrace();
      }
    }
    
    if(!view )
      new MessageDialog (null, "Saved Files", fastaFiles, false);
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
      if(!vcfRecord.getAlt().toString().equals(".") && isFwd)
      {
        buff.append(getBase(vcfRecord.getAlt().toString(), isFwd));
        position+=vcfRecord.getAlt().toString().length();
      }
      
      if(isFwd)
        position+=ndel-1;
      else
      {
        if(position-ndel+1 < 0)
          buff.delete(0, position);
        else   
          buff.delete(position-ndel+1, position);
      }
      
      for(int i=0; i<ndel; i++)
        buff.append("-");
    }
    else if(vcfRecord.getAlt().isInsertion(vcf_v4))
      buff.append(getBase(vcfRecord.getAlt().toString(), isFwd));
    else if(vcfRecord.getAlt().isMultiAllele())
    {
      String base = MultipleAlleleVariant.getIUBCode(vcfRecord);
      if(base != null)
        buff.append(base);
      else
        buff.append(bases.charAt(position));
    }
    else if(vcfRecord.getAlt().isNonVariant())                   // non-variant
      buff.append(getBase(vcfRecord.getRef(), isFwd).toUpperCase());
    else
      buff.append(getBase(vcfRecord.getAlt().toString().toLowerCase(), isFwd));
    
    if(isFwd && position < bases.length())
      buff.append(bases.substring(position+1));
    else if(!isFwd && position < bases.length())
      buff.append(bases.substring(position+1));

    return buff.toString();
  }
  
  /**
   * Get the actual bases by reverse complementing if on the
   * reverse strand.
   * @param baseStr
   * @param isFwd
   * @return
   */
  private static String getBase(String baseStr, boolean isFwd)
  {
    if(isFwd)
      return baseStr;
    return Bases.reverseComplement(baseStr);
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
  
  /**
   *  Return the dirtectory that the given entry was read from.
   **/
  private static File getBaseDirectoryFromEntry(final Entry entry) 
  {
    final uk.ac.sanger.artemis.io.Entry embl_entry = entry.getEMBLEntry();

    if(embl_entry instanceof DocumentEntry) 
    {
      final DocumentEntry document_entry =(DocumentEntry) embl_entry;

      if(document_entry.getDocument() instanceof FileDocument) 
      {
        final FileDocument file_document =
         (FileDocument) document_entry.getDocument();

        if(file_document.getFile().getParent() != null) 
          return new File(file_document.getFile().getParent());
      }
    }
    if(((DocumentEntry)entry.getEMBLEntry()).getDocument()
       instanceof RemoteFileDocument ||
       ((DocumentEntry)entry.getEMBLEntry()).getDocument()
       instanceof DatabaseDocument)
      return new File(System.getProperty("user.dir"));

    return null;
  }
  
}