/* BamUtils
 *
 * created: 2011
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

import java.util.HashMap;
import java.util.Hashtable;
import java.util.List;
import java.util.Vector;

import javax.swing.JProgressBar;

import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureSegmentVector;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.io.Range;

class BamUtils
{

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
    final Hashtable<String, SAMFileReader> samFileReaderHash = bamView.getSamFileReaderHash();
    final SAMRecordPredicate samRecordFlagPredicate = bamView.getSamRecordFlagPredicate();
    final SAMRecordPredicate samRecordMapQPredicate = bamView.getSamRecordMapQPredicate();

    int cnt[] = new int[2];
    cnt[0] = 0;
    cnt[1] = 0;
    
    SAMFileReader inputSam = samFileReaderHash.get(bam);
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
    final Hashtable<String, SAMFileReader> samFileReaderHash = bamView.getSamFileReaderHash();
    final SAMRecordPredicate samRecordFlagPredicate = bamView.getSamRecordFlagPredicate();
    final SAMRecordPredicate samRecordMapQPredicate = bamView.getSamRecordMapQPredicate();

    SAMFileReader inputSam = samFileReaderHash.get(bamFile);
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
}

