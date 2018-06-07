/* SequenceChangeEvent.java
 *
 * created: Wed Oct 28 1998
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 1998,1999,2000  Genome Research Limited
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
 *
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/sequence/SequenceChangeEvent.java,v 1.3 2005-11-28 16:46:38 tjc Exp $
 */

package uk.ac.sanger.artemis.sequence;

import uk.ac.sanger.artemis.io.Range;

/**
 *  This event is sent when the sequence of bases in a strand changes.
 *
 *  @author Kim Rutherford
 *  @version $Id: SequenceChangeEvent.java,v 1.3 2005-11-28 16:46:38 tjc Exp $
 *
 **/

public class SequenceChangeEvent extends uk.ac.sanger.artemis.ChangeEvent 
{
  /**
   *  The type that was passed to the constructor.
   **/
  /* final */ private int type;

  /**
   *  The position that was passed to the constructor.
   **/
  /* final */ private int position;

  /**
   *  The sub sequence of bases that was passed to the constructor.
   **/
  /* final */ private String sub_sequence;

  private int length;

  private Range range;

  final public static int DELETION = 1;

  final public static int INSERTION = 2;

  final public static int REVERSE_COMPLEMENT = 3;

  final public static int CONTIG_REVERSE_COMPLEMENT = 4;

  final public static int CONTIG_REORDER = 5;

  /**
   *  Create a new SequenceChangeEvent object.
   *  @param bases The Bases object that generated the event.
   *  @param type The type of the event (INSERTION or DELETION).
   *  @param position The position of the first base that was deleted or the
   *    position base immediately before the insertion.  If the insertion is
   *    at the start of the sequence the position will be 0.
   *  @param sub_sequence The bases that were inserted or deleted.
   **/
  public SequenceChangeEvent(final Bases bases,
                             final int type,
                             final int position,
                             final String sub_sequence) 
  {
    super(bases);
    this.type = type;
    this.position = position;
    this.sub_sequence = sub_sequence;
  }

  /**
   *  Create a new SequenceChangeEvent object (CONTIG_REORDER).
   **/
  public SequenceChangeEvent(final int type,
                             final int position,
                             final Range range)
  {
    super(range);
    this.type = type;
    this.position = position;
    this.range  = range;
  }

  public SequenceChangeEvent(final Bases bases,
                             final int type,
                             final Range range,
                             final int length) 
  {
    super(bases);
    this.type = type;
    this.range  = range;
    this.length = length;
  }


  /**
   *  Create a new SequenceChangeEvent object.
   *  @param bases The Bases object that generated the event.
   *  @param type The type of the event (should be REVERSE_COMPLEMENT).
   **/
  public SequenceChangeEvent(final Bases bases,
                             final int type) 
  {
    super (bases);
    this.type = type;
    this.position = 0;
    this.sub_sequence = null;
  }

  /**
   *  Return the Bases reference that was passed to the constructor.
   **/
  public Bases getBases() 
  {
    return (Bases)getSource();
  }

  /**
   *  Return the type of this event (INSERTION or DELETION).
   **/
  public int getType() 
  {
    return type;
  }

  /**
   *  Return the position of the deletion or insertion (as passed to the
   *  constructor).
   **/
  public int getPosition()
  {
    return position;
  }

  /**
   *  Return a String containing the bases that were inserted or deleted.
   **/
  public String getSubSequence()
  {
    return sub_sequence;
  }

  public int getLength()
  {
    return length;
  }

  public Range getRange()
  {
    return range;
  }

}


