/* SAMRecordPositionComparator
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

package uk.ac.sanger.artemis.components.alignment;

import java.util.Comparator;

import net.sf.samtools.SAMRecord;
import uk.ac.sanger.artemis.components.alignment.BamViewRecord;

 class SAMRecordPositionComparator implements Comparator<Object>
  {
    public BamView bamView;
    public SAMRecordPositionComparator(BamView bamView)
    {
      this.bamView = bamView;
    }
    
    public int compare(Object o1, Object o2) 
    {
      SAMRecord pr1 = ((BamViewRecord) o1).sam;
      SAMRecord pr2 = ((BamViewRecord) o2).sam;
      
      int offset1 = bamView.getSequenceOffset(pr1.getReferenceName());
      int offset2 = bamView.getSequenceOffset(pr2.getReferenceName());
      
      if(pr1.getAlignmentStart()+offset1 < pr2.getAlignmentStart()+offset2)
        return -1;
      else if(pr1.getAlignmentStart()+offset1 > pr2.getAlignmentStart()+offset2)
        return 1;  
      return 0;
    }
  }