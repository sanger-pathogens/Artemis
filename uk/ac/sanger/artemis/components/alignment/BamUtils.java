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

import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureSegmentVector;
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
      final int start,
      final int end,
      final String bam,
      final String refName,
      final Hashtable<String, SAMFileReader> samFileReaderHash,
      final Vector<String> seqNames,
      final HashMap<String, Integer> offsetLengths,
      final boolean concatSequences, 
      final HashMap<String, Integer> seqLengths,
      final SAMRecordPredicate samRecordFlagPredicate,
      final SAMRecordMapQPredicate samRecordMapQPredicate,
      final boolean contained,
      final boolean useStrandTag)
  {
    int cnt[] = new int[2];
    cnt[0] = 0;
    cnt[1] = 0;
    if(concatSequences)
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

          cnt = count(bam, samFileReaderHash, name, thisStart, thisEnd, 
              samRecordFlagPredicate, samRecordMapQPredicate, contained, true, useStrandTag);

        }
        lastLen = len;
      }
    }
    else
    {
      cnt = count(bam, samFileReaderHash, refName, start, end, 
          samRecordFlagPredicate, samRecordMapQPredicate, contained, true, useStrandTag);
    }
    
    float cntf[] = new float[2];
    cntf[0] = cnt[0];
    cntf[1] = cnt[1];
    return cntf;
  }

  protected static int[] count(final String bam, 
                    final Hashtable<String, SAMFileReader> samFileReaderHash, 
                    final String refName, 
                    final int start, 
                    final int end,
                    final SAMRecordPredicate samRecordFlagPredicate,
                    final SAMRecordPredicate samRecordMapQPredicate,
                    final boolean contained,
                    final boolean byStrand,
                    final boolean useStrandTag)
  {
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

  /**
   * Return the coverage for each base in a range for the forward and
   * reverse strand.
   * @param bamFile
   * @param samFileReaderHash
   * @param refName
   * @param start
   * @param end
   * @param samRecordFlagPredicate
   * @param samRecordMapQPredicate
   * @return
   */
  protected static int[][] countOverRange(final String bamFile,
                                          final Hashtable<String, SAMFileReader> samFileReaderHash, 
                                          final String refName,
                                          final int start, 
                                          final int end,
                                          final int concatShift,
                                          final int cnt[][],
                                          final SAMRecordPredicate samRecordFlagPredicate,
                                          final SAMRecordPredicate samRecordMapQPredicate)
  {
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
}

