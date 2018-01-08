/* BamUtils
 *
 * created: 2017
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2011 Genome Research Limited
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
package uk.ac.sanger.artemis.components.alignment;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.net.URL;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.List;
import java.util.Vector;

import javax.swing.JProgressBar;

import htsjdk.samtools.AlignmentBlock;
import htsjdk.samtools.BAMIndex;
import htsjdk.samtools.BAMIndexer;
import htsjdk.samtools.BamFileIoUtils;
import htsjdk.samtools.BamIndexValidator;
import htsjdk.samtools.CRAMCRAIIndexer;
import htsjdk.samtools.SAMException;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamFileValidator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.SamReaderFactory.Option;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.cram.CRAIIndex;
import htsjdk.samtools.cram.build.CramIO;
import htsjdk.samtools.seekablestream.SeekableStreamFactory;
import htsjdk.samtools.util.CloseableIterator;
import htsjdk.samtools.util.CloserUtil;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureSegmentVector;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.util.FTPSeekableStream;

/**
 * Utility methods for BamView.
 * 
 * @author kp11
 *
 */
class BamUtils
{
  private static org.apache.log4j.Logger logger = 
		    org.apache.log4j.Logger.getLogger(BamUtils.class);

  protected static float getFeatureLength(Feature f)
  {
    FeatureSegmentVector segs = f.getSegments();
    int len = 0;
    for(int i=0; i<segs.size(); i++)
    {
      Range r = segs.elementAt(i).getRawRange();
      len += r.getEnd()-r.getStart()+1;
    }
    return (float)len;
  }
  
  /**
   * Count the reads in a range.
   * @param start
   * @param end
   * @param bam
   * @param refName
   * @param samFileReaderHash
   * @param seqNames
   * @param offsetLengths
   * @param concatSequences
   * @param seqLengths
   * @param samRecordFlagPredicate
   * @param samRecordMapQPredicate
   * @param contained
   * @param useStrandTag - strand specific tag
   * @return
   */
  protected static float[] getCount(
      final BamView bamView,
      final int start,
      final int end,
      final String bam,
      final boolean contained,
      final boolean useStrandTag)
  {
    final Vector<String> seqNames = bamView.getSeqNames();
    final HashMap<String, Integer> offsetLengths = bamView.getOffsetLengths();
    final HashMap<String, Integer> seqLengths = bamView.getSeqLengths();

    int cnt[] = new int[2];
    cnt[0] = 0;
    cnt[1] = 0;
    if(bamView.isConcatSequences())
    {
      int len = 0;
      int lastLen = 1;
      for(String name : seqNames)
      {
        int thisLength = seqLengths.get(name);
        len += thisLength;

        if( (lastLen >= start && lastLen < end) ||
            (len >= start && len < end) ||
            (start >= lastLen && start < len) ||
            (end > lastLen && end < len) )
        {
          int offset = offsetLengths.get(name); 
          int thisStart = start - offset;
          if(thisStart < 1)
            thisStart = 1;
          int thisEnd   = end - offset;
          if(thisEnd > thisLength)
            thisEnd = thisLength;

          cnt = count(bamView, bam, thisStart, thisEnd, contained, true, useStrandTag);

        }
        lastLen = len;
      }
    }
    else
    {
      cnt = count(bamView, bam, start, end, contained, true, useStrandTag);
    }
    
    float cntf[] = new float[2];
    cntf[0] = cnt[0];
    cntf[1] = cnt[1];
    return cntf;
  }

  protected static int[] count(
          final BamView bamView,
          final String bam, 
          final int start,
          final int end,
          final boolean contained,
          final boolean byStrand,
          final boolean useStrandTag)
  {
    final String refName = (String) bamView.getCombo().getSelectedItem();
    final Hashtable<String, SamReader> samFileReaderHash = bamView.getSamFileReaderHash();
    final SAMRecordPredicate samRecordFlagPredicate = bamView.getSamRecordFlagPredicate();
    final SAMRecordPredicate samRecordMapQPredicate = bamView.getSamRecordMapQPredicate();

    int cnt[] = new int[2];
    cnt[0] = 0;
    cnt[1] = 0;
    
    SamReader inputSam = samFileReaderHash.get(bam);
    final CloseableIterator<SAMRecord> it = inputSam.query(refName, start, end, contained);

    while ( it.hasNext() )
    {
      SAMRecord samRecord = it.next();
      if( samRecordFlagPredicate == null ||
          !samRecordFlagPredicate.testPredicate(samRecord))
       {
         if(samRecordMapQPredicate == null ||
            samRecordMapQPredicate.testPredicate(samRecord))
         {
           if(byStrand && BamView.isNegativeStrand(samRecord, useStrandTag))
             cnt[1]++;
           else
             cnt[0]++;
         }
       }
    }
    it.close();
    return cnt;
  }
  
  protected static int[] calc(
      final BamView bamView, 
      final String refName, 
      final int sequenceLength,
      final boolean useStrandTag,
      final JProgressBar progressBar)
  {
    int mappedReads[] = new int[bamView.bamList.size()];
    int MAX_BASE_CHUNK = 2000 * 60;
    boolean contained = false;
    for (int i = 0; i < sequenceLength; i += MAX_BASE_CHUNK)
    {
      if(progressBar != null)
        progressBar.setValue(i);
      int sbegc = i;
      int sendc = i + MAX_BASE_CHUNK - 1;

      for (int j = 0; j < bamView.bamList.size(); j++)
      {
        String bam = bamView.bamList.get(j);
        if (bamView.isConcatSequences())
        {
          int len = 0;
          int lastLen = 1;
          for (String name : bamView.getSeqNames())
          {
            int thisLength = bamView.getSeqLengths().get(name);
            len += thisLength;

            if ((lastLen >= sbegc && lastLen < sendc)
                || (len >= sbegc && len < sendc)
                || (sbegc >= lastLen && sbegc < len)
                || (sendc >= lastLen && sendc < len))
            {
              int offset = bamView.getOffsetLengths().get(name);
              int thisStart = sbegc - offset;
              if (thisStart < 1)
                thisStart = 1;
              int thisEnd = sendc - offset;
              if (thisEnd > thisLength)
                thisEnd = thisLength;

              mappedReads[j] += BamUtils.count(bamView, bam, thisStart, thisEnd, 
                  contained, false, useStrandTag)[0];
            }
            lastLen = len;
          }
        }
        else
        {
          mappedReads[j] += BamUtils.count(bamView, bam, sbegc, sendc,
              contained, false, useStrandTag)[0];
        }
      }
    }
    return mappedReads;
  }

  /**
   * Return the coverage for each base in a range for the forward and
   * reverse strand.
   * @param bamView
   * @param bamFile
   * @param start
   * @param end
   * @param concatShift
   * @param cnt
   * @return
   */
  protected static int[][] countOverRange(
      final BamView bamView,
      final String bamFile, 
      final int start, 
      final int end, 
      final int concatShift, 
      final int cnt[][])
  {
    final String refName = (String) bamView.getCombo().getSelectedItem();
    final Hashtable<String, SamReader> samFileReaderHash = bamView.getSamFileReaderHash();
    final SAMRecordPredicate samRecordFlagPredicate = bamView.getSamRecordFlagPredicate();
    final SAMRecordPredicate samRecordMapQPredicate = bamView.getSamRecordMapQPredicate();

    SamReader inputSam = samFileReaderHash.get(bamFile);
    final CloseableIterator<SAMRecord> it = 
        inputSam.query(refName, start, end, false);

    while (it.hasNext())
    {
      SAMRecord samRecord = it.next();
      if (samRecordFlagPredicate == null
          || !samRecordFlagPredicate.testPredicate(samRecord))
      {
        if (samRecordMapQPredicate == null
            || samRecordMapQPredicate.testPredicate(samRecord))
        {
          List<AlignmentBlock> blocks = samRecord.getAlignmentBlocks();
          boolean isFwd = !samRecord.getReadNegativeStrandFlag();
          
          for(int j=0; j<blocks.size(); j++)
          {
            AlignmentBlock block = blocks.get(j);
            int refStart = block.getReferenceStart();
            for(int i=0; i<block.getLength(); i++)
            {
              int pos = refStart + i + concatShift;
              int bin = pos - start;
              if(bin < 0 || bin > cnt.length-1)
                continue;
              
              if(isFwd)
                cnt[bin][0]++;
              else
                cnt[bin][1]++;
            } 
          }
        }
      }
    }
    it.close();
    return cnt;
  }

  /**
   * For a list of features calculate the read count for each
   * @param bamView
   * @param features
   * @param contained
   * @param useIntrons
   * @param useStrandTag
   * @param mappedReads
   * @param progressBar
   * @return
   */
  protected static Hashtable<String, List<ReadCount>> calculateMappedReads(
      final BamView bamView,
      final FeatureVector features,
      final boolean contained, 
      final boolean useIntrons,
      final boolean useStrandTag,
      final int mappedReads[],
      final JProgressBar progressBar)
  {
    final Hashtable<String, List<ReadCount>> featureReadCount = 
        new Hashtable<String, List<ReadCount>>();
    for (int i = 0; i < features.size(); i++)
    {
      final Feature f = features.elementAt(i);
      if(progressBar != null)
        progressBar.setValue(i);

      int start = f.getRawFirstBase();
      int end = f.getRawLastBase();
      final float fLen = BamUtils.getFeatureLength(f);
      List<ReadCount> sampleCounts = new Vector<ReadCount>();

      for (int j = 0; j < bamView.bamList.size(); j++)
      {
        final String bam = bamView.bamList.get(j);
        float cnt[] = new float[2];

        cnt = BamUtils.getCount(bamView, start, end, bam, contained, useStrandTag);
        if (!useIntrons && f.getSegments().size() > 1)
        {
          // remove reads contained by intron
          for (int k = 0; k < f.getSegments().size()-1; k++)
          {
            int seg = k;
            int nextSeg = k+1;
            if(!f.isForwardFeature())
            {
              seg = f.getSegments().size()-k-1;
              nextSeg = seg-1;
            }

            start = f.getSegments().elementAt(seg).getRawRange().getEnd();
            end = f.getSegments().elementAt(nextSeg).getRawRange().getStart();

            float tmpcnt[] = new float[2];
            tmpcnt = BamUtils.getCount(bamView, start, end, bam, true, useStrandTag);
            cnt[0] -= tmpcnt[0];
            cnt[1] -= tmpcnt[1];
          }
        }
        
        if (mappedReads != null)
        {
          cnt[0] = (cnt[0] / (((float) mappedReads[j] / 1000000.f) * (fLen / 1000.f)));
          cnt[1] = (cnt[1] / (((float) mappedReads[j] / 1000000.f) * (fLen / 1000.f)));
        }

        sampleCounts.add( new ReadCount(cnt, f.isForwardFeature()) );
      }
      featureReadCount.put(ReadCountDialog.getFeatureName(f), sampleCounts);
    }
    return featureReadCount;
  }
  
  /**
   * Check whether two SAM records are the same.
   * @param bamRec1 SamViewRecord
   * @param bamRec2 SamViewRecord
   * @return boolean - true if equality holds.
   */
  public static boolean samRecordEqualityCheck(SAMRecord rec1, SAMRecord rec2)
  {
	  boolean result = false;
	  
	  if (rec1 == null && rec2 == null)
		  return true;
	  
	  result = (
		(rec1.getReadName().equals(rec2.getReadName())) &&
		(rec1.getAlignmentStart() == rec2.getAlignmentStart()) &&
		(rec1.getAlignmentEnd() == rec2.getAlignmentEnd()) &&
		(rec1.getFlags() == rec2.getFlags())
	  );
	  
	  return result;
  }
  
  /**
   * Determine if the given file name is expected to be a CRAM file
   * based on its extension.
   * @param filename String
   * @return boolean
   */
  public static boolean isCramFile(String filename) 
  {
	  boolean result = false;
	  
	  try 
	  {
		  result = ( filename.endsWith(CramIO.CRAM_FILE_EXTENSION) );
	  }
	  catch (Exception e)
	  {
		  result = false;
	  }
	  
	  return result;
  }
  
  /**
   * Determine if the given file name is expected to be a BAM file
   * based on its extension.
   * @param filename String
   * @return boolean
   */
  public static boolean isBamFile(String filename) 
  {
	  boolean result = false;
	  
	  try 
	  {
		  result = ( filename.endsWith(BamFileIoUtils.BAM_FILE_EXTENSION) );
	  }
	  catch (Exception e)
	  {
		  result = false;
	  }
	  
	  return result;
  }
  
  /**
   * Construct the index file name given the alignment file name.
   * @param filename String
   * @return String
   */
  public static String constructIndexFilename(String alignmentFilename) 
  {
	  String result = null;
	  
	  if (isBamFile(alignmentFilename))
		  result = alignmentFilename + BAMIndex.BAMIndexSuffix;
	  else if (isCramFile(alignmentFilename))
		  result = alignmentFilename + CRAIIndex.CRAI_INDEX_SUFFIX;
	  
	  return result;
  }
  
  /**
   * Create a temporary local index file.
   * @param alignmentFileName String
   * @return String
   * @throws IOException
   */
  public static File createTempIndexFile(URL alignmentFileName) throws IOException
  {
	  File file = null;
	  String ext = null;
	  
	  if (isBamFile(alignmentFileName.getFile()))
		  ext = BAMIndex.BAMIndexSuffix;
	  
	  else if (isCramFile(alignmentFileName.getFile()))
		  ext = CRAIIndex.CRAI_INDEX_SUFFIX;
	  
	  file = File.createTempFile(alignmentFileName.getFile().replaceAll("[\\/\\s]", "_"), ext);
	  file.deleteOnExit();
	      
	  return file;
  }
  
  /**
   * Create a local index file for the given alignment file.
   * @param alignmentFilename String
   * @param indexFile File to create
   * @throws IOException
   */
  public static void createIndexFileFromScratch(String alignmentFilename, File indexFile) throws IOException
  {
	  String indexFilePath = indexFile.getAbsolutePath();
	  
	  if (isBamFile(alignmentFilename))
	  {
		  // Only bother with local indexes
		  
		  logger.debug("Creating BAM index file from scratch: " + indexFilePath);
		  
		  final SamReaderFactory factory = SamReaderFactory.makeDefault();
		  factory.disable(Option.EAGERLY_DECODE);
		  factory.enable(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS);
		  factory.validationStringency(ValidationStringency.SILENT);
		  final SamReader bam = factory.open(new File(alignmentFilename));
		  
		  if (bam.type() != SamReader.Type.BAM_TYPE) 
		  {
			  throw new IOException("Input file must be a valid BAM file.");
	      }

		  if (!bam.getFileHeader().getSortOrder().equals(SAMFileHeader.SortOrder.coordinate)) 
		  {
			  throw new IOException("Input BAM file must be sorted by coordinate.");
		  }

		  BAMIndexer.createIndex(bam, indexFile);
		  CloserUtil.close(bam);
		 
	  }
	  else if (isCramFile(alignmentFilename))
	  {
		  logger.debug("Creating CRAM index file from scratch: " + indexFilePath);
		  
		  if (alignmentFilename.startsWith("ftp")) 
		  {
			  // special case for FTP
			  CRAMCRAIIndexer.writeIndex(
					  new FTPSeekableStream(new URL(alignmentFilename), getFtpSocketTimeout(), getFtpBufferSize()), 
					  new FileOutputStream(indexFile));
		  } 
		  else 
		  {
			  CRAMCRAIIndexer.writeIndex(
					  SeekableStreamFactory.getInstance().getStreamFor(alignmentFilename), 
					  new FileOutputStream(indexFile));
		  }
	  } 

  }
  
  /**
   * Validate a BAM/CRAM file,
   * @param reader SamReader
   * @param ref CRAMReferenceSequenceFile
   * @param throwOnError boolean
   * @throws SamException
   */
  public static void validateSAMFile(SamReader reader, CRAMReferenceSequenceFile ref, boolean throwOnError) throws SAMException
  { 
	  logger.debug("Performing alignment file validation...");
	  
	  try 
	  {
		  SamFileValidator validator = new SamFileValidator(new PrintWriter(System.out), 20);
		  validator.setIndexValidationStringency(BamIndexValidator.IndexValidationStringency.LESS_EXHAUSTIVE );
		  validator.validateSamFileSummary(reader, ref);
	  }
	  catch(SAMException e)
	  {
		 logger.warn("Validation errors: " + e.getMessage());
		 
		 if (throwOnError)
			 throw e;
	  }
  }
  
  /**
   * Read a SAM/BAM/CRAM index file and if one does not exist then
   * try to create it.
   * 
   * @param alignmentFile String
   * @return File index file
   * @throws IOException
   * @throws SAMException
   */
  static File getIndexFile(String alignmentFile) throws SAMException, IOException
  {
    File indexFile = null;
    boolean indexFileExists = false;
    boolean alignmentFileIsRemote = false;
    
    
    String indexFilename = constructIndexFilename(alignmentFile);
    if (indexFilename == null) {
    		// We cannot determine the alignment file type
  	  	throw new SAMException("Alignment file is not a known file type. The file suffix may be incorrect.");
    }
    
    if (alignmentFile.startsWith("http") || alignmentFile.startsWith("ftp"))
    {
      alignmentFileIsRemote = true;
    	  final URL urlIndexFile = new URL(indexFilename);
    	  
    	  /*
    	   * Create a local copy of the remote index file if possible.
    	   * [copy it to a temp file].
    	   */
      indexFile = createTempIndexFile(new URL(alignmentFile));

      InputStream is = null;
      FileOutputStream out = null;
      try 
      {
	      is = urlIndexFile.openStream();
	      out = new FileOutputStream(indexFile);
	      int c;
	      
	      while ((c = is.read()) != -1)
	        out.write(c);
	      
	      out.flush();
	      
	      indexFileExists = true;
	      logger.debug("Created local index file from remote... " + indexFile.getAbsolutePath());
      }
      catch (FileNotFoundException e)
      {
    	  	indexFileExists = false;
    	  	logger.debug("Cannot locate a remote index file: " + indexFile.getAbsolutePath());
      } 
      finally
      {
    	  	if (out != null)
    	  		out.close();
    	  	
    	  	if (is != null)
    	  		is.close();
      }
 
    }
    else
    {
    		// Use the local file index
    	
    		indexFile = new File(indexFilename);
    		indexFileExists = indexFile.exists();
    		
    		logger.debug("Using index file... " + indexFile.getAbsolutePath());
    }
    
    
    /*
     * We can't find an index file so try and create one locally
     * if the alignment file is local. If it is remote, then just flag
     * an error immediately.
     */
    if(!indexFileExists)
    {
	    	if (alignmentFileIsRemote) 
	    	{
	    		// Remote file
	    		
	    		String ls = System.getProperty("line.separator");
	    		String msg = ls+
		            "Failed to find an index file. Remote BAM/CRAM files must be indexed."+ls+
		            "Please download the file to your local file system."+ls+
		            "Indexing can then be done using samtools (http://www.htslib.org). For example: "+ls+ls+
		            "samtools sort <in.bam> -o <sorted.bam>"+ls+
		            "samtools index <sorted.bam>"+ls;
		        
		    logger.error("Failed to find an index file for remote alignment file: " + alignmentFile);
		    throw new SAMException(msg);
	    	}
	    	else 
	    	{
	    		// Local file
	    		
	    		try
	    		{
	    			logger.warn("Index file not found, so we will create one using the alignment file.");
	        
	    			createIndexFileFromScratch(alignmentFile, indexFile);
	    		}
	    		catch(Exception e)
	    		{
	    			String ls = System.getProperty("line.separator");
	    			String msg = ls+
		            "An index file could not be created."+ls+
		            "The BAM/CRAM file needs to be accessible, sorted and indexed."+ls+
		            "Indexing and sorting can be done using samtools (http://www.htslib.org). For example: "+ls+ls+
		            "samtools sort <in.bam> -o <sorted.bam>"+ls+
		            "samtools index <sorted.bam>"+ls;
		        
	    			logger.error("Failed to create index file " + indexFile.getAbsolutePath(), e);
		      
	    			throw new SAMException(msg);
	    		}
	    	}
    }

    return indexFile;
  }
  
  /**
   * Get an FTP timeout from options, if specified.
   * @return int
   */
  public static int getFtpSocketTimeout() 
  {
	  int timeout = 10000;
	  
	  if(Options.getOptions().getIntegerProperty("bamview_ftp_socket_timeout") != null)
	  { 
		  timeout = Options.getOptions().getIntegerProperty("bamview_ftp_socket_timeout");
    	  
		  logger.debug("BAM VIEW FTP SOCKET TIMEOUT=" + timeout);
	  }
	  
	  return timeout;
  }
  
  /**
   * Get an FTP buffer size from options, if specified.
   * @return int
   */
  public static int getFtpBufferSize() 
  {
	  int bufferSize = -1;  // automatically set
	  
	  if(Options.getOptions().getIntegerProperty("bamview_ftp_buffer_size") != null)
	  { 
		  bufferSize = Options.getOptions().getIntegerProperty("bamview_ftp_buffer_size");
    	  
		  logger.debug("BAM VIEW FTP BUFFER SIZE=" + bufferSize);
	  }
	  
	  return bufferSize;
  }
  
}

