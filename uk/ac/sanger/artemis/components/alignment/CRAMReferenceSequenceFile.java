/* CRAMReferenceSequenceFile.java
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

import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.SAMSequenceDictionary;
import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.io.IndexFastaStream;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.sequence.Bases;
import uk.ac.sanger.artemis.util.OutOfRangeException;

/**
 * Holds the reference data required for CRAM file functionality.
 * 
 * @author kp11
 *
 */
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
    		ReferenceSequence result = null;
    	
    		if(!isReferenceIndexed())
        {
    			result = getSubsequenceAt(contig, 1L, sequence.getBases().getLength());
        }
    		else
    		{
    			// We can use the fasta index...
        	    //
    			IndexFastaStream is = (IndexFastaStream)sequence.getEMBLEntry().getSequence();
    			result = is.getIndexSeqFile().getSequence(contig);
    		}
    	
    		return result;
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
    	
    	  ReferenceSequence result = null;
    	  
      try
      {
    	  	if(!isReferenceIndexed())
        {
        		int offset = bamView.getSequenceOffset(contig);
        		start += offset;
        		stop  += offset;
          
        		result = new ReferenceSequence(sequence.getName(), 0, 
                  sequence.getBases().getSubSequence(new Range((int)start, (int)stop), Bases.FORWARD).getBytes());
        }
        else
        {
        		// We can use the fasta index...
        	    //
        	  	IndexFastaStream is = (IndexFastaStream)sequence.getEMBLEntry().getSequence();
        	  	result = is.getIndexSeqFile().getSubsequenceAt(contig, (int)start, (int)stop);
        }
      }
      catch (OutOfRangeException e)
      {
        e.printStackTrace();
      }
      
      return result;
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
    
    /**
     * Return true if the reference fasta has an .fai index file.
     * @return boolean
     */
    protected boolean isReferenceIndexed() 
    {
    		boolean result = true;
    		
    		try
    		{
    			result = ( sequence.getEMBLEntry().getSequence() instanceof IndexFastaStream );
    		} catch (Exception e) 
    		{
    			result = false;
    		}
    		
    		return result;
    }
}
