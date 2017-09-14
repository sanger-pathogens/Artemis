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

import java.util.Vector;

import net.sf.samtools.SAMRecord;

/**
 *  Test the SAMRecord flag.
 **/
public class SAMRecordFlagConjunctionPredicate implements SAMRecordPredicate
{
  private Vector<SAMRecordPredicate> predicates;
  private int type;  // AND or OR
  protected static int AND = 0;
  protected static int OR = 1;
  
  public SAMRecordFlagConjunctionPredicate(SAMRecordPredicate p1, SAMRecordPredicate p2, int type)
  {
    predicates = new Vector<SAMRecordPredicate>();
    predicates.add(p1);
    predicates.add(p2);
    this.type = type;
  }
  
  public SAMRecordFlagConjunctionPredicate(Vector<SAMRecordPredicate> predicates, int type)
  {
    this.predicates = predicates;
    this.type = type;
  }

  public boolean testPredicate(SAMRecord samRecord)
  {
    for(SAMRecordPredicate predicate: predicates)
    {
      if(predicate.testPredicate(samRecord))
      {
        if(type == OR)
          return true;
      }
      else
      {
        if(type == AND)
          return false;
      }
    }

    if(type == AND)
      return true;
    return false;
  }

}