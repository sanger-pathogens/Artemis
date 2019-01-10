/* MarkerChangeEvent.java
 *
 * created: Sun Mar 21 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/sequence/MarkerChangeEvent.java,v 1.1 2004-06-09 09:52:18 tjc Exp $
 */

package uk.ac.sanger.artemis.sequence;

/**
 *  This event is sent when the position or Strand of a Marker changes.
 *
 *  @author Kim Rutherford
 *  @version $Id: MarkerChangeEvent.java,v 1.1 2004-06-09 09:52:18 tjc Exp $
 **/

public class MarkerChangeEvent extends java.util.EventObject {
  /**
   *  Create a new MarkerChange event
   *  @param marker The Marker that has changed.
   *  @param strand The Strand that the Marker was associated with before
   *    the MarkerChange event.
   *  @param position The position on the Strand of the Marker before the
   *    event.
   **/
  public MarkerChangeEvent (final Marker marker,
                            final Strand strand, final int position) {
    super (marker);
    this.marker = marker;
    this.strand = strand;
    this.position = position;
  }

  /**
   *  Return the Marker that was passed to the constructor.
   **/
  public Marker getMarker () {
    return marker;
  }

  /**
   *  Return the Strand that was passed to the constructor.
   **/
  public Strand getStrand () {
    return strand;
  }

  /**
   *  The Marker that was passed to the constructor.
   **/
  final private Marker marker;

  /**
   *  The Strand that was passed to the constructor.
   **/
  final private Strand strand;

  /**
   *  The position that was passed to the constructor.
   **/
  final private int position;
}


