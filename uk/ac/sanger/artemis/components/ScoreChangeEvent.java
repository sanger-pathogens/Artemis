/* ScoreChangeEvent.java
 *
 * created: Thu Oct 21 1999
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 1999  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/ScoreChangeEvent.java,v 1.1 2004-06-09 09:47:22 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

/**
 *  The adjustment event emitted when a score changes.
 *
 *  @author Kim Rutherford
 *  @version $Id: ScoreChangeEvent.java,v 1.1 2004-06-09 09:47:22 tjc Exp $
 **/

public class ScoreChangeEvent extends java.util.EventObject {
  /**
   *  Create new ScoreChangeEvent.
   **/
  public ScoreChangeEvent (final Object source, final int value) {
    super (source);
    this.value = value;
  }

  /**
   *  Return the new score value that caused this event (as passed to the
   *  constructor).
   **/
  public int getValue () {
    return value;
  }
  

  /**
   *  The new score value that caused this event (as passed to the
   *  constructor).
   **/
  final int value;
}
