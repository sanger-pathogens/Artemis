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

import java.awt.Container;
import java.awt.Frame;
import java.text.DecimalFormat;
import java.util.Hashtable;
import java.util.List;
import java.util.Vector;

import javax.swing.JFrame;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureSegmentVector;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.components.EntryEdit;
import uk.ac.sanger.artemis.components.FileViewer;
import uk.ac.sanger.artemis.components.MultiComparator;
import uk.ac.sanger.artemis.io.Range;

class BamUtils
{
  /**
   * Read count for selected reads.
   * @param features
   * @param refName
   * @param samFileReaderHash
   * @param bamList
   * @param seqNames
   * @param offsetLengths
   * @param concatSequences
   * @param seqLengths
   * @param samRecordFlagPredicate
   * @param samRecordMapQPredicate
   * @param contained
   * @param useIntrons
   * @param mappedReads
   */
  protected static void countReads(final FeatureVector features,
                                   final String refName,
                                   final Hashtable<String, SAMFileReader> samFileReaderHash,
                                   final List<String> bamList, 
                                   final Vector<String> seqNames,
                                   final Hashtable<String, Integer> offsetLengths,
                                   final boolean concatSequences, 
                                   final Hashtable<String, Integer> seqLengths,
                                   final SAMRecordFlagPredicate samRecordFlagPredicate,
                                   final SAMRecordMapQPredicate samRecordMapQPredicate,
                                   final boolean contained,
                                   final boolean useIntrons,
                                   final int mappedReads[])
  {
    Hashtable<String, List<Float>> featureReadCount = new Hashtable<String, List<Float>>();

    for(int i=0; i<features.size(); i++)
    {     
      Feature f = features.elementAt(i);
      
      int start  = f.getFirstBase();
      int end    = f.getLastBase();
      float fLenKb = getFeatureLength(f);
      List<Float> sampleCounts = new Vector<Float>();
        
      for(int j=0; j<bamList.size(); j++)
      {
        String bam = bamList.get(j);
        float cnt = 0;
        if(!useIntrons && f.getSegments().size() > 1)
        {
          for(int k=0; k<f.getSegments().size(); k++)
          {
            start = f.getSegments().elementAt(k).getStart().getPosition();
            end   = f.getSegments().elementAt(k).getEnd().getPosition();
            cnt += getCount(start, end, bam, refName, samFileReaderHash, 
              seqNames, offsetLengths, concatSequences, seqLengths, 
              samRecordFlagPredicate, samRecordMapQPredicate, contained); 
          }
        }
        else
          cnt = getCount(start, end, bam, refName, samFileReaderHash, 
              seqNames, offsetLengths, concatSequences, seqLengths, 
              samRecordFlagPredicate, samRecordMapQPredicate, contained); 
       
        if(mappedReads != null)
          cnt = cnt / ((((float)mappedReads[j]) / 1000000.f)*fLenKb);
        
        sampleCounts.add(cnt);
      }
      featureReadCount.put(f.getSystematicName(), sampleCounts);
    }

    DecimalFormat df = new DecimalFormat("0.00##");
    StringBuffer buff = new StringBuffer();
    for(int j=0; j<bamList.size(); j++)
    {
      String bam = bamList.get(j);
      buff.append("#BAM: "+bam);
      if(mappedReads != null)
        buff.append(" Mapped Reads/million: "+ df.format( ((float)mappedReads[j]) / 1000000.f) );
      buff.append("\n");
    }
    buff.append("\n");
    
    for (String fId : featureReadCount.keySet() ) {
      buff.append(fId+"\t");
      List<Float> cnts = featureReadCount.get(fId);
      for(int i=0; i<cnts.size(); i++)
        buff.append(df.format(cnts.get(i)) + (i<cnts.size()-1 ? "\t" : ""));
      buff.append("\n");
    }

    FileViewer viewer;
    if(mappedReads != null)
      viewer = new FileViewer ("RPKM", true, false, true);
    else
      viewer = new FileViewer ("Read Count", true, false, true);
    viewer.getTextPane().setText(buff.toString());
  }

  private static float getFeatureLength(Feature f)
  {
    FeatureSegmentVector segs = f.getSegments();
    int len = 0;
    for(int i=0; i<segs.size(); i++)
    {
      Range r = segs.elementAt(i).getRawRange();
      len += r.getEnd()-r.getStart();
    }
    return ((float)len) / 1000.f;
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
  private static int getCount(
      final int start,
      final int end,
      final String bam,
      final String refName,
      final Hashtable<String, SAMFileReader> samFileReaderHash,
      final Vector<String> seqNames,
      final Hashtable<String, Integer> offsetLengths,
      final boolean concatSequences, 
      final Hashtable<String, Integer> seqLengths,
      final SAMRecordFlagPredicate samRecordFlagPredicate,
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

  /**
   * Calculate the total number of mapped reads.
   * @param refName
   * @param samFileReaderHash
   * @param bamList
   * @param seqNames
   * @param offsetLengths
   * @param concatSequences
   * @param seqLengths
   * @param sequenceLength
   * @return
   */
  protected static int[] getTotalMappedReads(
      final String refName,
      final Hashtable<String, SAMFileReader> samFileReaderHash,
      final List<String> bamList, 
      final Vector<String> seqNames,
      final Hashtable<String, Integer> offsetLengths,
      final boolean concatSequences, 
      final Hashtable<String, Integer> seqLengths,
      final int sequenceLength)
  {
    int MAX_BASE_CHUNK = 2000*60;
    int mapped[] = new int[bamList.size()];
    boolean contained = false;
    SAMRecordFlagPredicate samRecordFlagPredicate = new SAMRecordFlagPredicate(SAMRecordFlagPredicate.READ_UNMAPPED_FLAG); 
    SAMRecordMapQPredicate samRecordMapQPredicate = new SAMRecordMapQPredicate(-1);

    for (int i = 0; i < sequenceLength; i += MAX_BASE_CHUNK)
    {
      int sbegc = i;
      int sendc = i + MAX_BASE_CHUNK - 1;

      for (int j=0; j<bamList.size(); j++)
      {
        String bam = bamList.get(j);
        if (concatSequences)
        {
          int len = 0;
          int lastLen = 1;
          for (String name : seqNames)
          {
            int thisLength = seqLengths.get(name);
            len += thisLength;

            if ((lastLen >= sbegc && lastLen < sendc)
                || (len >= sbegc && len < sendc)
                || (sbegc >= lastLen && sbegc < len)
                || (sendc >= lastLen && sendc < len))
            {
              int offset = offsetLengths.get(name);
              int thisStart = sbegc - offset;
              if (thisStart < 1)
                thisStart = 1;
              int thisEnd = sendc - offset;
              if (thisEnd > thisLength)
                thisEnd = thisLength;

              mapped[j] += count(bam, samFileReaderHash, name, thisStart, thisEnd,
                  samRecordFlagPredicate, samRecordMapQPredicate, contained);

            }
            lastLen = len;
          }
        }
        else
        {
          mapped[j] += count(bam, samFileReaderHash, refName, sbegc, sendc,
              samRecordFlagPredicate, samRecordMapQPredicate, contained);
        }
      }
    }
    return mapped;
  }
  
  private static int count(String bam, 
                    Hashtable<String, SAMFileReader> samFileReaderHash, 
                    String refName, 
                    int start, 
                    int end,
                    SAMRecordFlagPredicate samRecordFlagPredicate,
                    SAMRecordMapQPredicate samRecordMapQPredicate,
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
  
  protected static Container getBamContainer(BamView bamView)
  {
    Frame fs[] = JFrame.getFrames();
    for(Frame f: fs)
    {
      if( f instanceof JFrame && 
         ((JFrame)f) instanceof EntryEdit ||
         ((JFrame)f) instanceof MultiComparator)
        return ((JFrame)f).getContentPane();
    }
    return bamView;
  }
}

