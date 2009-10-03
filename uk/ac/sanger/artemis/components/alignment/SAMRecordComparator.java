/* SAMRecordComparator
 *
 * created: 2009
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2009  Genome Research Limited
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

 class SAMRecordComparator implements Comparator<Object>
  {
    public int compare(Object o1, Object o2) 
    {
      SAMRecord pr1 = (SAMRecord) o1;
      SAMRecord pr2 = (SAMRecord) o2;
      
      int cmp = pr1.getReadName().compareTo(pr2.getReadName());
      if(cmp == 0)
      {
        if(pr1.getAlignmentStart() < pr2.getAlignmentStart())
          return -1;
        else
          return 1;
      }   
      return cmp;
    }
  }