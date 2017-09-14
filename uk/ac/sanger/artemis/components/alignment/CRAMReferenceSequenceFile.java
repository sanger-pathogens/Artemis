/* BamView
 *
 * created: 2012
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2012  Genome Research Limited
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

import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.samtools.SAMSequenceDictionary;
import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.sequence.Bases;
import uk.ac.sanger.artemis.util.OutOfRangeException;


  class CRAMReferenceSequenceFile implements ReferenceSequenceFile
  {
    private Entry sequence;
    private BamView bamView;
    
    CRAMReferenceSequenceFile(final Entry sequence, final BamView bamView)
    {
      this.sequence = sequence;
      this.bamView  = bamView;
    }
    
    /**
    * Retrieves the complete sequence described by this contig.
    * @param contig contig whose data should be returned.
    * @return The full sequence associated with this contig.
    * @throws UnsupportedOperationException if !sIndexed.
    */
    public ReferenceSequence getSequence(String contig)
    {
      return getSubsequenceAt(contig, 1L, sequence.getBases().getLength());
      
      //return new ReferenceSequence(sequence.getName(), 0, 
      //    sequence.getBases().getForwardStrand().getStrandBases().getBytes());
    }

    public SAMSequenceDictionary getSequenceDictionary()
    {
      return null;
    }

    /**
    * Gets the subsequence of the contig in the range [start,stop]
    * @param contig Contig whose subsequence to retrieve.
    * @param start inclusive, 1-based start of region.
    * @param stop inclusive, 1-based stop of region.
    * @return The partial reference sequence associated with this range.
    * @throws UnsupportedOperationException if !sIndexed.
    */
    public ReferenceSequence getSubsequenceAt(String contig, long start, long stop )
    {
      try
      {
        if(bamView.isConcatSequences())
        {
          int offset = bamView.getSequenceOffset(contig);
          start += offset;
          stop  += offset;
        }
        
        return new ReferenceSequence(sequence.getName(), 0, 
            sequence.getBases().getSubSequence(new Range((int)start, (int)stop), Bases.FORWARD).getBytes());
      }
      catch (OutOfRangeException e)
      {
        e.printStackTrace();
      }
      return null;
    }

    /**
    * @return true if getSequence and getSubsequenceAt methods are allowed.
    */
    public boolean isIndexed()
    {
      return true;
    }

    public ReferenceSequence nextSequence()
    {
      return null;
    }

    public void reset()
    {
    }

    public void close()
    {
    }
  }
