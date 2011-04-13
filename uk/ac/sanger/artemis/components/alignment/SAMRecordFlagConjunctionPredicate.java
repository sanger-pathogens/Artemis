/* SAMRecordFlagConjunctionPredicate
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 2011  Genome Research Limited
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 **/

package uk.ac.sanger.artemis.components.alignment;

import net.sf.samtools.SAMRecord;

/**
 *  Test the SAMRecord flag.
 **/
public class SAMRecordFlagConjunctionPredicate implements SAMRecordPredicate
{
  private SAMRecordPredicate p1;
  private SAMRecordPredicate p2;
  private int type;  // AND or OR
  private boolean is1 = true;
  private boolean is2 = true;
  private static int AND = 0;
  
  public SAMRecordFlagConjunctionPredicate(SAMRecordPredicate p1, SAMRecordPredicate p2, int type, boolean is1, boolean is2)
  {
    this.p1 = p1;
    this.p2 = p2;
    this.type = type;
    this.is1 = is1;
    this.is2 = is2;
  }
  
  public boolean testPredicate(SAMRecord samRecord)
  {
    if(type == AND)
    {
      if(!is1 && !is2)
      {
        if(!p1.testPredicate(samRecord) && !p2.testPredicate(samRecord))
          return true;
      }
      else if(!is1)
      {
        if(!p1.testPredicate(samRecord) && p2.testPredicate(samRecord))
          return true;
      }
      else if(!is2)
      {
        if(p1.testPredicate(samRecord) && !p2.testPredicate(samRecord))
          return true;
      }
      else if(is1 && is2)
      {
        if(p1.testPredicate(samRecord) && p2.testPredicate(samRecord))
          return true;
      }
    }
    else
    {
      if(!is1 && !is2)
      {
        if(!p1.testPredicate(samRecord) || !p2.testPredicate(samRecord))
          return true;
      }
      else if(!is1)
      {
        if(!p1.testPredicate(samRecord) || p2.testPredicate(samRecord))
          return true;
      }
      else if(!is2)
      {
        if(p1.testPredicate(samRecord) || !p2.testPredicate(samRecord))
          return true;
      }
      else if(is1 && is2)
      {
        if(p1.testPredicate(samRecord) || p2.testPredicate(samRecord))
          return true;
      }
    }
    return false;
  }
  
}