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

import java.util.Hashtable;
import java.util.Vector;

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
   * @return
   */
  protected static int getCount(
      final int start,
      final int end,
      final String bam,
      final String refName,
      final Hashtable<String, SAMFileReader> samFileReaderHash,
      final Vector<String> seqNames,
      final Hashtable<String, Integer> offsetLengths,
      final boolean concatSequences, 
      final Hashtable<String, Integer> seqLengths,
      final SAMRecordPredicate samRecordFlagPredicate,
      final SAMRecordMapQPredicate samRecordMapQPredicate,
      final boolean contained)
  {
    int cnt = 0;
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
            (end >= lastLen && end < len) )
        {
          int offset = offsetLengths.get(name); 
          int thisStart = start - offset;
          if(thisStart < 1)
            thisStart = 1;
          int thisEnd   = end - offset;
          if(thisEnd > thisLength)
            thisEnd = thisLength;

          cnt = count(bam, samFileReaderHash, name, thisStart, thisEnd, 
              samRecordFlagPredicate, samRecordMapQPredicate, contained);

        }
        lastLen = len;
      }
    }
    else
    {
      cnt = count(bam, samFileReaderHash, refName, start, end, 
          samRecordFlagPredicate, samRecordMapQPredicate, contained);
    }
    return cnt;
  }

  protected static int count(String bam, 
                    Hashtable<String, SAMFileReader> samFileReaderHash, 
                    String refName, 
                    int start, 
                    int end,
                    SAMRecordPredicate samRecordFlagPredicate,
                    SAMRecordPredicate samRecordMapQPredicate,
                    boolean contained)
  {
    int cnt = 0;
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
            cnt++;
         }
       }
    }
    it.close();
    return cnt;
  }

}

