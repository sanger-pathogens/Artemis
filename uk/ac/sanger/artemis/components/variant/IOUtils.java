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

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.io.Writer;
import java.net.URL;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Vector;

import javax.swing.Box;
import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;

import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureEnumeration;
import uk.ac.sanger.artemis.FeatureKeyQualifierPredicate;
import uk.ac.sanger.artemis.FeaturePredicate;
import uk.ac.sanger.artemis.FeaturePredicateConjunction;
import uk.ac.sanger.artemis.FeatureSegment;
import uk.ac.sanger.artemis.FeatureSegmentVector;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.Selection;
import uk.ac.sanger.artemis.components.FileViewer;
import uk.ac.sanger.artemis.components.MessageDialog;
import uk.ac.sanger.artemis.components.SequenceViewer;
import uk.ac.sanger.artemis.components.StickyFileChooser;
import uk.ac.sanger.artemis.io.DocumentEntry;
import uk.ac.sanger.artemis.io.EntryInformationException;
import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.io.Location;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.sequence.Bases;
import uk.ac.sanger.artemis.sequence.MarkerRange;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.FileDocument;
import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.util.ReadOnlyException;
import uk.ac.sanger.artemis.util.RemoteFileDocument;

import net.sf.samtools.util.BlockCompressedInputStream;

class IOUtils
{ 
  
  private static int MAXIMUM_SELECTED_FEATURES = 25;
  private static int SEQUENCE_LINE_BASE_COUNT = 60;
  
  /**
   * Write filtered uncompressed VCF. Uses the filter in VCFview to
   * determine if variants are written.
   * @param manualHash
   * @param vcfFileName
   * @param vcfView
   * @param features
   * @param nfiles
   * @return
   */
  protected static File writeVCF(final Map<String, Boolean> manualHash,
                                 final String vcfFileName,
                                 final int vcfIndex,
                                 final VCFview vcfView,
                                 final FeatureVector features,
                                 final int nfiles)
  {
    try
    {
      File filterFile = getFile(vcfFileName, nfiles, ".filter", null);
      if(filterFile == null)
        return null;
      FileWriter writer = new FileWriter(filterFile);
      AbstractVCFReader.write(manualHash, vcfFileName, vcfIndex, writer, vcfView, features);

      return filterFile;
    }
    catch (IOException e)
    {
      e.printStackTrace();
      return null;
    }
  }
  
  private static File getFile(final String vcfFileName, final int nfiles,
      final String suffix, final JComponent comp) throws IOException
  {
    if(nfiles > 1)
    {
      if(vcfFileName.startsWith("http"))
      {
        int ind = vcfFileName.lastIndexOf('/')+1;
        return new File(vcfFileName.substring(ind)+suffix);
      }
      else
        return new File(vcfFileName+suffix);
    }

    final StickyFileChooser file_dialog = new StickyFileChooser();
    file_dialog.setSelectedFile(new File(vcfFileName+suffix));
    file_dialog.setDialogTitle("Choose save file ...");
    file_dialog.setDialogType(JFileChooser.SAVE_DIALOG);
    if(comp != null)
      file_dialog.setAccessory(comp);
    
    final int status = file_dialog.showSaveDialog(null);

    if(status != JFileChooser.APPROVE_OPTION ||
       file_dialog.getSelectedFile() == null) 
      return null;
    
    return file_dialog.getSelectedFile();
  }
  
  /**
   * Export as a VCF based on the filtering applied in the VCFview.
   * @param manualHash
   * @param vcfFiles
   * @param vcfView
   */
  protected static void export(final Map<String, Boolean> manualHash,
                               final List<String> vcfFiles,
                               final VCFview vcfView)
  {
    // get all CDS features that do not have the /pseudo or /pseudogene qualifier
    final FeatureVector features = getFeatures(
        new FeaturePredicateConjunction(
            new FeatureKeyQualifierPredicate(Key.CDS, "pseudo", false),
            new FeatureKeyQualifierPredicate(Key.CDS, "pseudogene", false),
            FeaturePredicateConjunction.AND), vcfView.getEntryGroup());
    
    String filterFiles = "";
    for(int i=0; i<vcfFiles.size(); i++)
    {
      File filterFile = IOUtils.writeVCF(manualHash, vcfFiles.get(i), i, vcfView, features, vcfFiles.size());
      if(filterFile == null)
        return;
      filterFiles += filterFile.getAbsolutePath()+"\n";
    }

    new MessageDialog (null, "Saved Files", filterFiles, false);
  }

  /**
   * Export all variant sites to a multiple fasta file.
   * @param vcfView
   */
  protected static void exportVariantFasta(final VCFview vcfView)
  {
    final EntryGroup entryGroup = vcfView.getEntryGroup();
    final String name = entryGroup.getActiveEntries().elementAt(0).getName();
    final File newfile = new File( 
        getBaseDirectoryFromEntry(entryGroup.getActiveEntries().elementAt(0)),
            name);

    try
    {
      final File f = getFile(newfile.getAbsolutePath(), 1, ".fasta", null);
      if(f == null)
        return;

      final FileOutputStream fos = new FileOutputStream(f);
      exportVariantFasta(vcfView, 
          new PrintWriter(fos), 
          entryGroup.getSequenceLength(), 
          entryGroup.getAllFeatures(), 
          entryGroup.getBases());
      fos.close();
    }
    catch(IOException ioe)
    {
      ioe.printStackTrace();
    }
  }
  
  /**
   * Export all variant sites to a multiple fasta file.
   * @param vcfView
   * @param pw
   * @param length
   * @param features
   * @param bases
   */
  protected static void exportVariantFasta(final VCFview vcfView, 
                                           final PrintWriter pw, 
                                           final int length, 
                                           final FeatureVector features, 
                                           final Bases bases)
  {
    try
    {
      int ntotalSamples = 0;
      for (int i = 0; i < vcfView.getVcfReaders().length; i++)
        ntotalSamples += vcfView.getVcfReaders()[i].getNumberOfSamples();
      
      final Writer[] writer = new Writer[ntotalSamples+1];
      final File[] tmpFiles = new File[ntotalSamples+1];
      
      for (int i = 0; i < vcfView.getVcfReaders().length; i++)
      {
        final String names[] = vcfView.getVcfReaders()[i].sampleNames;
        for(int j=0; j<names.length; j++)
        {
          final String fn = (names[j].equals("") ? (j+1)+"_sample" : names[j].replaceAll("[/\\:]", "_"));
          tmpFiles[i+j] = File.createTempFile(fn, "art");
          writer[i+j] = new FileWriter( tmpFiles[i+j] );
          writer[i+j].write(">"+fn);
        }
      }
      
      // include reference bases
      final String refName = vcfView.getEntryGroup().getActiveEntries().elementAt(0).getName();
      tmpFiles[ntotalSamples] = File.createTempFile("ref", "art");
      writer[ntotalSamples] = new FileWriter( tmpFiles[ntotalSamples] );
      writer[ntotalSamples].write(">"+refName);


      final int MAX_BASE_CHUNK = (10000/ntotalSamples)*SEQUENCE_LINE_BASE_COUNT;
      final HashMap<Integer, VCFRecord> records[] = new HashMap[MAX_BASE_CHUNK];
      int baseCount = 0;
      
      // write variant sites to tmp files
      for(int i=0; i<length; i+=MAX_BASE_CHUNK)
      {
        int start = i+1;
        int end = i+MAX_BASE_CHUNK;
        if(end > length)
          end = length;

        storeVCFRecords(vcfView, records, start, end);
        baseCount = writeVariants(vcfView, records, writer, features, 
                        ntotalSamples, start, end, bases, baseCount);
      }

      for(int i=0; i<writer.length; i++)
        writer[i].close();
      
      // concatenate the single fasta files into a multiple fasta file
      for (int i = 0; i < tmpFiles.length; i++) 
      {
        final BufferedReader br = new BufferedReader(new FileReader(tmpFiles[i].getPath()));
        String line;
        while ( (line = br.readLine()) != null) 
          pw.println(line);      
        br.close();
        tmpFiles[i].delete();
      }
      pw.close();
    }
    catch(IOException e)
    { 
      e.printStackTrace();
    }
    catch (OutOfRangeException e)
    {
      e.printStackTrace();
    }
  }
  
  /**
   * For a given range store the VFFRecord in a Hashtable.
   * @param vcfView
   * @param records
   * @param start
   * @param end
   * @throws IOException
   */
  private static void storeVCFRecords(final VCFview vcfView, 
                                      final HashMap<Integer, VCFRecord> records[], 
                                      final int start, 
                                      final int end) throws IOException
  {
    
    for (int i = 0; i < vcfView.getVcfReaders().length; i++)
      records[i] = new HashMap<Integer, VCFRecord> ();
    
    for (int i = 0; i < vcfView.getVcfReaders().length; i++)
    {
      AbstractVCFReader reader = vcfView.getVcfReaders()[i];
      if(vcfView.isConcatenate())
      {
        String[] contigs = reader.getSeqNames();
        for(int j=0; j<contigs.length; j++)
        {
          int offset = vcfView.getSequenceOffset(contigs[j]);
          int nextOffset;
          if(j<contigs.length-1)
            nextOffset = vcfView.getSequenceOffset(contigs[j+1]);
          else
            nextOffset = vcfView.seqLength;
          
          if( (start >= offset && start <= nextOffset) ||
              (end   >= offset && end <= nextOffset) )
          {
            int thisStart = start - offset;
            if(thisStart < 1)
              thisStart = 1;
            loadRecords(records[i], reader, contigs[j], thisStart, end - offset, offset);
          }
        }
      }
      else
        loadRecords(records[i], reader, vcfView.getChr(), start, end, 0);
    }
  }
  
  private static void loadRecords(HashMap<Integer, VCFRecord> records,
                                  final AbstractVCFReader reader, 
                                  final String contig, 
                                  final int start, 
                                  final int end,
                                  final int offset) throws IOException
  {
    VCFRecord record;
    while((record = reader.getNextRecord(contig, start, end)) != null)
    {
      if(records == null)
        records = new HashMap<Integer, VCFRecord> ();
      records.put(record.getPos()+offset, record);
    }
  }
  
  /**
   * Write variant positions.
   * @param vcfView
   * @param records
   * @param writer
   * @param features
   * @param ntotalSamples
   * @param start
   * @param end
   * @param bases
   * @param bc
   * @return
   * @throws IOException
   * @throws OutOfRangeException
   */
  private static int writeVariants(final VCFview vcfView,
                                   final HashMap<Integer, VCFRecord> records[],
                                   final Writer writer[],
                                   final FeatureVector features,
                                   final int ntotalSamples,
                                   final int start, 
                                   final int end,
                                   final Bases bases,
                                   int bc) throws IOException, OutOfRangeException
  {
    final char basesC[] = bases.getSubSequenceC(new Range(start, end), Bases.FORWARD);
    for (int i = start; i < end; i++)
    {
      int thisSample = 0;
      final String[] thisBase = new String[ntotalSamples+1];

      boolean seenSNP     = false;
      int insertionLength = 0;

      // loop over each VCF file
      for (int j = 0; j < vcfView.getVcfReaders().length; j++)
      {
        AbstractVCFReader reader = vcfView.getVcfReaders()[j];
        VCFRecord record = records[j].get(i);

        if(record == null)
          continue;

        boolean vcf_v4 = reader.isVcf_v4();
        int nsamples = reader.getNumberOfSamples();
        // loop over each sample
        for(int k=0; k<nsamples; k++)
        {
          // look at each type of variant
          if(vcfView.showVariant(record, features, i, reader, k, j) )
          {
            if(record.getAlt().isDeletion(vcf_v4))
            {
              // note: do not write out if just deletion
              thisBase[thisSample] = "-";
              
              /*if( thisBase[ntotalSamples] == null || 
                  thisBase[ntotalSamples].length() < record.getRef().length() )
              {
                thisBase[ntotalSamples] = record.getRef();
                if(!record.getAlt().toString().equals("."))
                  thisBase[thisSample] = record.getAlt().toString();
                else
                  thisBase[thisSample] = "";
                
                int padLength = thisBase[ntotalSamples].length() - thisBase[thisSample].length();
                for(int ipad=0; ipad<padLength; ipad++)
                  thisBase[thisSample] += "-";
              }*/
            }
            else if(record.getAlt().isInsertion(vcf_v4))
            {
              String in = record.getAlt().toString();
              if(in.startsWith("I"))
                in = in.substring(1);
              thisBase[thisSample] = in;
              if(in.length() > insertionLength)
                insertionLength = in.length();
              seenSNP = true;
              
              if( (thisBase[ntotalSamples] == null || 
                   thisBase[ntotalSamples].length() < record.getRef().length()) &&
                   in.toLowerCase().startsWith(record.getRef().toLowerCase()))
                thisBase[ntotalSamples] = record.getRef();
            }
            else if(record.getAlt().isMultiAllele(k))
            {
              String base = MultipleAlleleVariant.getIUBCode(record);
              if(base != null)
              {
                thisBase[thisSample] = base;
                seenSNP = true;
              }
            }
            else if(record.getAlt().isNonVariant()) 
            {
              thisBase[thisSample] = ".";
            }
            else
            {
              if(record.getAlt().toString().length() == record.getRef().length() )
              {
                thisBase[thisSample] = record.getAlt().toString();
                seenSNP = true;
                
                if(thisBase[ntotalSamples] == null || 
                   thisBase[ntotalSamples].length() < record.getRef().length())
                  thisBase[ntotalSamples]  = record.getRef();
              }
            }
          }
          else
            thisBase[thisSample] = "N"; // filtered out

          thisSample++;
        }
      }
      
      if(seenSNP)
      {
        // look-up reference base
        if(thisBase[ntotalSamples] == null)
          thisBase[ntotalSamples] = String.valueOf( basesC[i-start] );

        int remainder = 0;
        for(int j=0; j<thisBase.length; j++)
        {         
          if(thisBase[j] != null)
          {
            for(int k=0; k<thisBase[j].length(); k++)
            {
              remainder = (bc+k)%SEQUENCE_LINE_BASE_COUNT;
              if(remainder == 0)
                writer[j].write(System.getProperty("line.separator"));
              writer[j].write(thisBase[j].charAt(k));
            }
          }
          else
          {
            remainder = bc%SEQUENCE_LINE_BASE_COUNT;
            if(remainder == 0)
              writer[j].write(System.getProperty("line.separator"));
            writer[j].write("N");
          }
          
          if(insertionLength > 0)
          {
            int ins;
            if(thisBase[j] != null)
              ins = insertionLength-thisBase[j].length();
            else
              ins = insertionLength-1;
            
            int rem = remainder+1;
            for(int k=0; k<ins; k++)
            {
              remainder = (rem+k)%SEQUENCE_LINE_BASE_COUNT;
              if(remainder == 0)
                writer[j].write(System.getProperty("line.separator"));
              writer[j].write("-");
            }
          }
        }

        if(insertionLength > 0)
          bc+=insertionLength;
        else
          bc++;
      }
    }
    
    return bc;
  }
  
  
  /**
   * Write out FASTA for a selected base range
   * @param vcfView
   * @param selection
   * @param view
   */
  protected static void exportFastaByRange(
                                    final VCFview vcfView,
                                    final Selection selection,
                                    final boolean view,
                                    Writer writer)
  {
    if(selection.getMarkerRange() == null)
    {
      JOptionPane.showMessageDialog(null, 
          "No base range selected.", 
          "Warning", JOptionPane.WARNING_MESSAGE);
      return;  
    }

    AbstractVCFReader vcfReaders[] = vcfView.getVcfReaders();
    MarkerRange marker = selection.getMarkerRange();
    Range range = marker.getRawRange();
    String fastaFiles = "";

    EntryGroup entryGroup = vcfView.getEntryGroup();
    String name = entryGroup.getActiveEntries().elementAt(0).getName();
    int sbeg = range.getStart();
    int send = range.getEnd();

    StringBuffer buffSeq = null;
    try
    {
      final JCheckBox useNs = new JCheckBox("Use N for filtered out sites", true);
      useNs.setToolTipText("Mask filtered sites.");
      final JCheckBox useMask = new JCheckBox("Use N for sites without non-variant", true);
      useMask.setToolTipText("Mask sites that are not confirmed by a non-variant record.");
      Box yBox = Box.createVerticalBox();
      yBox.add(useNs);
      yBox.add(useMask);
      if(writer == null)
        JOptionPane.showMessageDialog(null, yBox, "Options", JOptionPane.INFORMATION_MESSAGE);
      else
        useMask.setSelected(false);

      if(!view && writer == null)
      {
        File newfile = new File(
           getBaseDirectoryFromEntry(entryGroup.getActiveEntries().elementAt(0)),
               name);
      
        File f = getFile(newfile.getAbsolutePath(), 1, ".fasta", null);
        if(f == null)
          return;
        writer = new FileWriter(f);
        fastaFiles += f.getAbsolutePath()+"\n";
      }
      else
        buffSeq = new StringBuffer();

      Bases bases = entryGroup.getSequenceEntry().getBases();
      // reference
      writeOrViewRange(null, -1, sbeg, send, writer, buffSeq, 
          marker, bases, name, vcfView, entryGroup, useNs.isSelected(), useMask.isSelected());

      // vcf sequences
      for (int i = 0; i < vcfReaders.length; i++)
        writeOrViewRange(vcfReaders[i], i, sbeg, send, writer, buffSeq,
            marker, bases, name, vcfView, entryGroup, useNs.isSelected(), useMask.isSelected());

      if(writer != null)
        writer.close();
    }
    catch(IOException e)
    {
      e.printStackTrace();
    }
    catch (OutOfRangeException e)
    {
      e.printStackTrace();
    }

    if(!view)
    {
      if(writer instanceof FileWriter)
        new MessageDialog (null, "Saved Files", fastaFiles, false);
    }
    else
    {
      FileViewer viewer = new FileViewer ("Feature base viewer for selected range: " + 
          sbeg+":"+send+(marker.isForwardMarker() ? "" : " reverse"), true, false, true);
      viewer.getTextPane().setText(buffSeq.toString());
    }
  }

  /**
   * Write the FASTA sequence out for the given features for each of the
   * VCF/BCF files.
   * @param vcfView
   * @param features
   * @param view
   */
  protected static void exportFasta(final VCFview vcfView,
                                    final FeatureVector features,
                                    final boolean view,
                                    Writer writer)
  {
    if(features.size () < 1)
    {
      JOptionPane.showMessageDialog(null, 
          "No features selected.", 
          "Warning", JOptionPane.WARNING_MESSAGE);
      return;  
    }

    if(view && features.size () > MAXIMUM_SELECTED_FEATURES)
      new MessageDialog (null,
                        "warning: only viewing the sequences for " +
                        "the first " + MAXIMUM_SELECTED_FEATURES +
                        " selected features");

    String suffix = ".fasta";
    if(features.size() == 1)
      suffix = "."+features.elementAt(0).getIDString()+suffix;

    String fastaFiles = "";
    final AbstractVCFReader vcfReaders[] = vcfView.getVcfReaders();

    final JCheckBox single = new JCheckBox("Single FASTA", true);
    final JCheckBox combineFeats = new JCheckBox("Combine feature sequences", true);
    final JCheckBox useNs = new JCheckBox("Use N for filtered out sites", true);
    useNs.setToolTipText("Mask filtered sites.");
    final JCheckBox useMask = new JCheckBox("Use N for sites without non-variant", true);
    useMask.setToolTipText("Mask sites that are not confirmed by a non-variant record.");
    
    if(writer != null)
      useMask.setSelected(false);
    
    Box yBox = Box.createVerticalBox();
    if(!view && vcfReaders.length > 1)
      yBox.add(single);
    yBox.add(combineFeats);
    yBox.add(useNs);
    yBox.add(useMask);
    
    final String name = vcfView.getEntryGroup().getActiveEntries().elementAt(0).getName();
    try
    {
      if(!view && writer == null)
      {
        File newfile = new File(
            getBaseDirectoryFromEntry(vcfView.getEntryGroup().getActiveEntries().elementAt(0)),
                name);
        File f = getFile(newfile.getAbsolutePath(), 1, suffix, yBox);
        if(f == null)
          return;
        writer = new FileWriter(f);
        fastaFiles += f.getAbsolutePath()+"\n";
      }
      else if(writer == null)
        JOptionPane.showMessageDialog(null, yBox, "View Option(s)", JOptionPane.INFORMATION_MESSAGE);

      // reference sequence
      StringBuffer buff = new StringBuffer();
      for (int j = 0; j < features.size() && (!view || j < MAXIMUM_SELECTED_FEATURES); j++)
      {
        Feature f = features.elementAt(j);
        buff.append( f.getBases() );
        if(!combineFeats.isSelected())
        {
          writeOrView(null, f, writer, buff, "");
          buff = new StringBuffer();
        }
      }
      if(combineFeats.isSelected())
        writeOrView(null, null, writer, buff, name);
      if(writer != null && !single.isSelected())
        writer.close();
      
      // 
      for (int i = 0; i < vcfReaders.length; i++)
      {
        if(!view && !single.isSelected())
        {
          File f = getFile(vcfReaders[i].getFileName(), vcfReaders.length, suffix, null);
          writer = new FileWriter(f);
          fastaFiles += f.getAbsolutePath()+"\n";
        }
        buff = new StringBuffer();
        
        for (int j = 0; j < features.size() && (!view || j < MAXIMUM_SELECTED_FEATURES); j++)
        {
          Feature f = features.elementAt(j);
          FeatureSegmentVector segs = f.getSegments();
          
          for(int k=0; k<segs.size(); k++)
          {
            FeatureSegment seg = segs.elementAt(k);
            int sbeg = seg.getRawRange().getStart();
            int send = seg.getRawRange().getEnd();
            buff.append( getAllBasesInRegion(vcfReaders[i], i, sbeg, send, seg.getBases(),
                  features, vcfView, f.isForwardFeature(), useNs.isSelected(), useMask.isSelected()) );
          }

          if(!combineFeats.isSelected())
          {
            writeOrView(vcfReaders[i], f, writer, buff, "");
            buff = new StringBuffer();
          }
        }
        
        if(combineFeats.isSelected())
          writeOrView(vcfReaders[i], null, writer, buff, "");

        if(writer != null && !single.isSelected())
          writer.close();
      }
      
      if(writer != null && single.isSelected())
        writer.close();
    } 
    catch(IOException e)
    {
      e.printStackTrace(); 
    }
    
    if(!view && writer instanceof FileWriter)
      new MessageDialog (null, "Saved Files", fastaFiles, false);
  }
  
  private static StringBuffer getHeader(AbstractVCFReader reader, 
      MarkerRange marker, String seqName, 
      int sbeg, int send)
  {
    StringBuffer header = new StringBuffer();
    if(reader != null)
      header.append(reader.getName()).append(" ");
    header.append(seqName).append(" ");
    header.append(sbeg).append(":").append(send);
    if(marker != null)
      header.append((marker.isForwardMarker() ? "" : " reverse"));
    return header;
  }
  
  private static void writeOrViewRange(AbstractVCFReader reader,
                                       final int vcfIndex,
                                       int sbeg, int send,
                                       Writer writer, StringBuffer buffSeq, 
                                       MarkerRange marker, Bases bases, 
                                       String name,
                                       VCFview vcfView,
                                       final EntryGroup entryGroup, 
                                       final boolean useNs,
                                       final boolean useMask) throws IOException, OutOfRangeException
  {
    int direction = ( marker.isForwardMarker() ? Bases.FORWARD : Bases.REVERSE);
    int length = send-sbeg+1;
    int MAX_BASE_CHUNK = 2000*SEQUENCE_LINE_BASE_COUNT;
    String basesStr;
    StringBuffer header = getHeader(reader, marker, name, sbeg, send);
    int linePos = 0;

    for(int i=0; i<length; i+=MAX_BASE_CHUNK)
    {
      int sbegc = sbeg+i;
      int sendc = sbeg+i+MAX_BASE_CHUNK-1;
      if(i+MAX_BASE_CHUNK-1 > length)
        sendc = send;

      int sbegc_raw = sbegc;
      int sendc_raw = sendc;
      if(direction == Bases.REVERSE)
      {
        sendc = bases.getLength () - sbegc_raw + 1;
        sbegc = bases.getLength () - sendc_raw + 1;
      }

      MarkerRange m = new MarkerRange(marker.getStrand(), sbegc, sendc);
      basesStr = bases.getSubSequence(m.getRange(), direction);
      FeatureVector features = entryGroup.getFeaturesInRange(m.getRange());
      //System.out.println((reader == null ? "" : reader.getName())+" "+sbegc+".."+sendc);
      if(reader != null)
        basesStr = getAllBasesInRegion(reader, vcfIndex, sbegc_raw, sendc_raw, basesStr,
                       features, vcfView, marker.isForwardMarker(), useNs, useMask);
      else
        basesStr = basesStr.toUpperCase();

      linePos = writeOrView(writer, header, basesStr, buffSeq, linePos);
      header = null;
    }
  }
  
  private static int writeOrView(Writer writer, 
                                 StringBuffer header, 
                                 String basesStr, 
                                 StringBuffer buff,
                                 int linePos) throws IOException
  {
    if(writer == null) // sequence viewer
    {
      if(header != null)
        buff.append(">").append(header.toString()).append("\n");
      wrapString(basesStr, buff);
    }
    else    // write to file
      return writeSequence(writer, header, basesStr, linePos);

    return 0;
  }
  
  /**
   * Construct a header and write or view the sequence.
   * @param reader
   * @param f
   * @param writer
   * @param buff
   * @throws IOException
   */
  private static void writeOrView(AbstractVCFReader reader, Feature f,
      Writer writer, StringBuffer buff, String hdr)
          throws IOException
  {
    StringBuffer header = new StringBuffer(hdr);
    final String basesStr;
    
    if(reader != null)
    {
      header.append(reader.getName()).append(" ");
      basesStr = buff.toString();
    }
    else
      basesStr = buff.toString().toUpperCase();
    
    if(f != null)
    {
      header.append(f.getSystematicName()).append(" ");
      header.append(f.getIDString()).append(" ");
      final String product = f.getProductString();
      header.append( (product == null ? "undefined product" : product) );
      header.append(" ").append(f.getWriteRange());
    }
    
    if(writer == null) // sequence viewer
    {
      SequenceViewer viewer =
          new SequenceViewer ("Feature base viewer for feature(s)", false);  
      viewer.setSequence(">"+header.toString(), basesStr);
    }
    else    // write to file
      writeSequence(writer, header, basesStr, 0);
  }
  
  /**
   * For a given VCF file change the sequence in a range and return the
   * base sequence as a string.
   * @param reader
   * @param vcfIndex
   * @param sbeg
   * @param send
   * @param basesStr
   * @param features
   * @param vcfView
   * @param isFwd
   * @param useNs
   * @param useMask
   * @return
   * @throws IOException
   */
  private static String getAllBasesInRegion(final AbstractVCFReader reader,
      final int vcfIndex,
      final int sbeg,
      final int send,
      String basesStr,
      final FeatureVector features,
      final VCFview vcfView,
      final boolean isFwd,
      final boolean useNs,
      final boolean useMask) throws IOException
  {
    if(vcfView.isConcatenate())
    {
      String[] contigs = reader.getSeqNames();
      for(int j=0; j<contigs.length; j++)
      {
        int offset = vcfView.getSequenceOffset(contigs[j]);
        int nextOffset;
        if(j<contigs.length-1)
          nextOffset = vcfView.getSequenceOffset(contigs[j+1]);
        else
          nextOffset = vcfView.seqLength;
        
        if( (offset >= sbeg && offset < send) ||
            (offset < sbeg && sbeg < nextOffset) )
        {
          int thisStart = sbeg - offset;
          if(thisStart < 1)
            thisStart = 1;
          int thisEnd   = send - offset;
          basesStr = getBasesInRegion(reader, vcfIndex, contigs[j], thisStart, thisEnd, 
              basesStr, features, vcfView, isFwd, useNs, useMask);
        }
      }
    }
    else
      basesStr = getBasesInRegion(reader, vcfIndex, vcfView.getChr(), sbeg, send, 
          basesStr, features, vcfView, isFwd, useNs, useMask);

    return basesStr;
  }
  
  /**
   * For a given VCF file change the sequence in a range and return the
   * base sequence as a string.
   * @param reader
   * @param vcfIndex
   * @param chr
   * @param sbeg
   * @param send
   * @param basesStr
   * @param features
   * @param vcfView
   * @param isFwd
   * @param useNs
   * @param useMask
   * @return
   * @throws IOException
   */
  private static String getBasesInRegion(final AbstractVCFReader reader,
                                         final int vcfIndex,
                                         final String chr,
                                         int sbeg,
                                         final int send,
                                         String basesStr,
                                         final FeatureVector features,
                                         final VCFview vcfView,
                                         final boolean isFwd,
                                         final boolean useNs, 
                                         final boolean useMask) throws IOException
  {
    boolean vcf_v4 = reader.isVcf_v4();
    int len = basesStr.length();
    int baseNum = sbeg;
    try
    {
      VCFRecord record;
      while ((record = reader.getNextRecord(chr, sbeg, send)) != null)
      {
        //
        // mask regions with N where there are no records
        if(useMask && baseNum < record.getPos())
          basesStr = maskSites(sbeg, record.getPos(), baseNum, isFwd, basesStr);
        baseNum = record.getPos()+1;

        int basePosition = record.getPos() + vcfView.getSequenceOffset(record.getChrom());
        if(vcfView.showVariant(record, features, basePosition, reader, -1, vcfIndex) )
          basesStr = getSeqsVariation(record, basesStr, sbeg, isFwd, vcf_v4);
        else if(useNs && isSNPorNonVariant(record))
        {
          int position = record.getPos()-sbeg;
          if(!isFwd)
            position = basesStr.length()-position-1;
          basesStr = basesStr.substring(0, position) + 'n' +
                     basesStr.substring(position+1);
        }

        // adjust for insertions
        if(basesStr.length() > len) 
        {
          sbeg -= (basesStr.length()-len);
          len = basesStr.length();
        }
      }
    }
    catch(NullPointerException e)
    {
      System.err.println(chr+":"+sbeg+"-"+send+"\n"+e.getMessage());
    }

    if(useMask && baseNum-sbeg < len)
      basesStr = maskSites(sbeg, len+sbeg, baseNum, isFwd, basesStr);

    return basesStr;
  }
  
  /**
   * Mask sequence sites
   * @param sbeg
   * @param endMask
   * @param baseNum
   * @param isFwd
   * @param basesStr
   * @return
   */
  private static String maskSites(final int sbeg, final int endMask, 
      final int baseNum, final boolean isFwd, String basesStr)
  {
    for(int i=baseNum; i<endMask; i++)
    {
      int position = i-sbeg;
      if(!isFwd)
        position = basesStr.length()-position-1;
      basesStr = basesStr.substring(0, position) + 'n' +
                 basesStr.substring(position+1);
    }
    return basesStr;
  }
  

  protected static void countVariants(final VCFview vcfView,
                                      final FeatureVector features) throws IOException
  {
    if(features.size () < 1)
    {
      JOptionPane.showMessageDialog(null, 
          "No features selected.", 
          "Warning", JOptionPane.WARNING_MESSAGE);
      return;  
    }
    
    String[] columnNames = { 
        "VCF", "Name", "Variant", "Non-variant", "Deletion", "Insertion", "Synonymous", "Non-synonymous"};
    Vector<String> columnData = new Vector<String>();
    for(String col: columnNames)
      columnData.add(col);
    Vector<Vector<Object>> rowData = new Vector<Vector<Object>>();
    
    AbstractVCFReader vcfReaders[] = vcfView.getVcfReaders();
    for(int vcfIndex=0; vcfIndex<vcfReaders.length; vcfIndex++)
    {
      AbstractVCFReader reader = vcfReaders[vcfIndex];
      for (int j = 0; j < features.size(); j++)
      {
        int count[] = new int[6];
        for(int c: count)
          c = 0;
        
        Feature f = features.elementAt(j);
        FeatureSegmentVector segs = f.getSegments();
        
        for(int k=0; k<segs.size(); k++)
        {
          FeatureSegment seg = segs.elementAt(k);
          int sbeg = seg.getRawRange().getStart();
          int send = seg.getRawRange().getEnd();

          if(vcfView.isConcatenate())
          {
            String[] contigs = reader.getSeqNames();
            for(int i=0; i<contigs.length; i++)
            {
              int offset = vcfView.getSequenceOffset(contigs[i]);
              int nextOffset;
              if(i<contigs.length-1)
                nextOffset = vcfView.getSequenceOffset(contigs[i+1]);
              else
                nextOffset = vcfView.seqLength;
              
              if( (offset >= sbeg && offset < send) ||
                  (offset < sbeg && sbeg < nextOffset) )
              {
                int thisStart = sbeg - offset;
                if(thisStart < 1)
                  thisStart = 1;
                int thisEnd   = send - offset;
              
                VCFRecord record;
                while ((record = reader.getNextRecord(vcfView.getChr(), thisStart, thisEnd)) != null)
                  count(record, count, features, reader, vcfIndex, vcfView);
              }
            }
          }
          else
          {
            VCFRecord record;
            while ((record = reader.getNextRecord(vcfView.getChr(), sbeg, send)) != null)
              count(record, count, features, reader, vcfIndex, vcfView);
          }
        }

        Object row[] = { 
            reader.getName(), f.getSystematicName(), count[0], count[1], count[2], count[3], count[4], count[5] };

        Vector<Object> thisRow = new Vector<Object>();
        for(Object obj: row)
          thisRow.add(obj);
        rowData.add(thisRow);
      }
    }
   
    TableViewer tab = new TableViewer(rowData, columnData, "Variant Overview");
    for(int i=2; i< columnData.size(); i++)
      tab.setIntegerRowSorter(i);
  }
  
  private static void count(VCFRecord record, int count[], FeatureVector features, AbstractVCFReader reader, int vcfIndex, VCFview vcfView)
  {
    int basePosition = record.getPos() + vcfView.getSequenceOffset(record.getChrom());
    if(!vcfView.showVariant(record, features, basePosition, reader, -1, vcfIndex) )
      return;
    
    if(record.getAlt().isNonVariant())
    {
      count[1]++;
      return;
    }
    else
      count[0]++;
   
    if(record.getAlt().isDeletion(reader.isVcf_v4()))
      count[2]++;
    else if(record.getAlt().isInsertion(reader.isVcf_v4()))
      count[3]++;
    
    if(record.getAlt().length() == 1 && record.getRef().length() == 1)
    {
      short synFlag = record.getSynFlag(features, record.getPos());
      switch(synFlag)
      {
        case 1:  count[4]++; break;  // synonymous
        default: count[5]++; break;  // non-synonymous
      }
    }  
  }
  
  private static boolean isSNPorNonVariant(VCFRecord record)
  {
    return (record.getRef().length() == 1 && record.getAlt().length() == 1) || record.getAlt().isNonVariant();
  }
  
  protected static void wrapString(String bases, StringBuffer buff)
  {
    final int SEQUENCE_LINE_BASE_COUNT = 60;
    for(int k=0; k<bases.length(); k+=SEQUENCE_LINE_BASE_COUNT)
    {
      int end = k + SEQUENCE_LINE_BASE_COUNT;
      if(end > bases.length())
        end = bases.length();
      buff.append ( bases.substring(k,end) ).append("\n");
    }
  }
  
  private static int writeSequence(Writer writer, 
                                   StringBuffer header, 
                                   String bases, 
                                   int startPos) throws IOException
  {
    if(header != null)
      writer.write (">" + header.toString() + "\n");
    int k = 0;
    for(k=0; k<bases.length(); k+=SEQUENCE_LINE_BASE_COUNT)
    {
      int end = k + SEQUENCE_LINE_BASE_COUNT - startPos;
      if(end > bases.length())
        end = bases.length();
      writer.write ( bases.substring(k,end) );
      
      if(k < bases.length() -1)
        writer.write("\n");
      
      startPos = 0;
    }
    
    return k % SEQUENCE_LINE_BASE_COUNT;
  }
  
  /**
   * Change the bases to reflect a variation record.
   * @param vcfRecord
   * @param bases
   * @param sbeg
   * @param isFwd
   * @param vcf_v4
   * @return
   */
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
      int ndel = vcfRecord.getAlt().getNumberOfIndels(vcf_v4);
      if(isFwd &&
         !vcfRecord.getAlt().toString().equals(".") && 
         !vcfRecord.getAlt().toString().startsWith("D"))
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
    {
      if(!isFwd)
        buff.delete(position-vcfRecord.getRef().length()+1, position);
      
      String in = vcfRecord.getAlt().toString();
      if(in.startsWith("I"))
        in = in.substring(1);
      buff.append(getBase(in, isFwd));

      if(isFwd)
        position+=(vcfRecord.getRef().toString().length()-1);
    }
    else if(vcfRecord.getAlt().isMultiAllele(-1))
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
   * Create features for each variant that has not been filtered out.
   * @param vcfView
   * @param entryGroup
   */
  protected static void createFeatures(final VCFview vcfView,
                                       final EntryGroup entryGroup)
  {
    final Entry newEntry = entryGroup.createEntry("VCF");

    int sbeg = 1;
    int send = entryGroup.getSequenceLength();
    int MAX_BASE_CHUNK = 1000 * SEQUENCE_LINE_BASE_COUNT;

    Bases bases = entryGroup.getSequenceEntry().getBases();

    for (int i = 0; i < send; i += MAX_BASE_CHUNK)
    {
      int sbegc = sbeg + i;
      int sendc = sbeg + i + MAX_BASE_CHUNK - 1;
      if (i + MAX_BASE_CHUNK - 1 > send)
        sendc = send;

      try
      {
        Range range = new Range(sbegc, sendc);
        FeatureVector features = entryGroup.getFeaturesInRange(range);
        String chr = vcfView.getChr();
        AbstractVCFReader vcfReaders[] = vcfView.getVcfReaders();
        for(int vcfIndex=0; vcfIndex<vcfReaders.length; vcfIndex++)
        {
          AbstractVCFReader reader = vcfReaders[vcfIndex];
          if(vcfView.isConcatenate())
          {
            for(String contig: reader.getSeqNames())
              makeFeatures(reader, vcfIndex, contig, sbegc, sendc, features, vcfView, bases, newEntry);
          }
          else
            makeFeatures(reader, vcfIndex, chr, sbegc, sendc, features, vcfView, bases, newEntry);
        }
      }
      catch (IOException ioe)
      {
        ioe.printStackTrace();
      }
      catch (OutOfRangeException e)
      {
        e.printStackTrace();
      }
    }
  }
  
  private static void makeFeatures(
      final AbstractVCFReader reader,
      final int vcfIndex,
      final String chr, 
      final int sbegc, 
      final int sendc, 
      final FeatureVector features, 
      final VCFview vcfView, 
      final Bases bases, 
      final Entry entry) throws IOException, OutOfRangeException
  {
    Key variantKey = new Key("misc_difference");
    try
    {
      VCFRecord record;
      while( (record = reader.getNextRecord(chr, sbegc, sendc)) != null)
      {
        makeFeature(record, reader.getName(), vcfView, features, bases, entry, variantKey, reader, vcfIndex);
      }
    }
    catch (NullPointerException e)
    {
      System.err.println(chr + ":" + sbegc + "-" + sendc + "\n"
          + e.getMessage());
    }
  }

  private static void makeFeature(
      final VCFRecord record,
      final String vcfFileName,
      final VCFview vcfView, 
      final FeatureVector features, 
      final Bases bases, 
      final Entry entry, 
      final Key variantKey, 
      final AbstractVCFReader vcfReader,
      final int vcfIndex) throws OutOfRangeException, ReadOnlyException
  {
    int basePosition = record.getPos() + vcfView.getSequenceOffset(record.getChrom());
    if (vcfView.showVariant(record, features, basePosition, vcfReader, -1, vcfIndex))
    {
      MarkerRange marker = new MarkerRange(bases.getForwardStrand(),
          basePosition, basePosition);
      Location location = marker.createLocation();
      QualifierVector qualifiers = new QualifierVector();
      String qualifierStr = record.getRef()+"->"+record.getAlt().toString()+
                            "; "+vcfFileName+"; score="+record.getQuality();
      if(record.getAlt().isMultiAllele(-1))
        qualifierStr += "; MULTI-ALLELE";
      else if(record.getAlt().isDeletion(vcfReader.isVcf_v4()))
        qualifierStr += "; DELETION";
      else if(record.getAlt().isInsertion(vcfReader.isVcf_v4()))
        qualifierStr += "; INSERTION";
      else if(record.getAlt().isNonVariant())
        return;
      
      try
      {
        FeatureVector fs = entry.getFeaturesInRange(marker.getRange());
        if(fs.size() > 0)
        {
          for(int i=0; i<fs.size(); i++)
          {
            Feature f = fs.elementAt(i);
            if(f.getKey().compareTo(variantKey) == 0)
            {
              f.getQualifiers().addQualifierValues(
                  new Qualifier("note", qualifierStr));
              return;
            }
          }
        }
        
        qualifiers.addQualifierValues(new Qualifier("note", qualifierStr));
        entry.createFeature(variantKey, location, qualifiers);
      }
      catch (EntryInformationException e)
      {
        e.printStackTrace();
      }
    }
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
    if(fileName.startsWith("http:") || fileName.startsWith("ftp:"))
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