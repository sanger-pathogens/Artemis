/* SAMRecordPredicate
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 2009  Genome Research Limited
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

import htsjdk.samtools.SAMRecord;

/**
 *  Test the SAMRecord flag.
 **/
public class SAMRecordFlagPredicate implements SAMRecordPredicate
{
  private int flag;
  private boolean isSet;
  
  protected static final int READ_PAIRED_FLAG = 0x1;
  protected static final int PROPER_PAIR_FLAG = 0x2;
  protected static final int READ_UNMAPPED_FLAG = 0x4;
  protected static final int MATE_UNMAPPED_FLAG = 0x8;
  protected static final int READ_STRAND_FLAG = 0x10;
  protected static final int MATE_STRAND_FLAG = 0x20;
  protected static final int FIRST_OF_PAIR_FLAG = 0x40;
  protected static final int SECOND_OF_PAIR_FLAG = 0x80;
  protected static final int NOT_PRIMARY_ALIGNMENT_FLAG = 0x100;
  protected static final int READ_FAILS_VENDOR_QUALITY_CHECK_FLAG = 0x200;
  protected static final int DUPLICATE_READ_FLAG = 0x400;
  
  protected static final String[] FLAGS_DESCRIPTION =
  {
    "Read Paired",
    "Proper Pair",
    "Read Unmapped",
    "Mate Unmapped",
    "Read on Negative Strand",
    "Mate on Negative Strand",
    "First of Pair",
    "Second of Pair",
    "Not Primary Alignment",
    "Read Fails Vendor Quality Check",
    "Duplicate Read"
  };
  
  protected static int[] FLAGS =
  {
    READ_PAIRED_FLAG,
    PROPER_PAIR_FLAG,
    READ_UNMAPPED_FLAG,
    MATE_UNMAPPED_FLAG,
    READ_STRAND_FLAG,
    MATE_STRAND_FLAG,
    FIRST_OF_PAIR_FLAG,
    SECOND_OF_PAIR_FLAG,
    NOT_PRIMARY_ALIGNMENT_FLAG,
    READ_FAILS_VENDOR_QUALITY_CHECK_FLAG,
    DUPLICATE_READ_FLAG
  };

  
  public SAMRecordFlagPredicate(int flag)
  {
    this(flag, true);
  }
  
  public SAMRecordFlagPredicate(int flag, boolean isSet)
  {
    this.flag = flag;
    this.isSet = isSet;
  }
  
  /**
   *  Test the given SAMRecord against this predicate.
   *  @param feature The SAMRecord to test the predicate against.
   *  @return Return true if and only if this predicate is true for the given
   *    SAMRecord.
   **/
  public boolean testPredicate (final SAMRecord samRecord)
  {
    if(!isSet)
      return !isFlagSet(samRecord.getFlags());
    return isFlagSet(samRecord.getFlags());
  }
  
  private boolean isFlagSet(int thisFlag)
  {
    for(int i=0; i<FLAGS.length; i++)
    {
      if((flag & FLAGS[i]) == FLAGS[i])
      {
        if((thisFlag & FLAGS[i]) == FLAGS[i])
          return true;
      }
    }
    return false;
  }
}


