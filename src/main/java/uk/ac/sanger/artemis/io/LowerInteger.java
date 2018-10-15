/* LowerInteger.java
 *
 * created: Mon May  3 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/LowerInteger.java,v 1.1 2004-06-09 09:49:55 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

/**
 *  This class is used to represent a lower end of a range like this <100..200
 *
 *  @author Kim Rutherford
 *  @version $Id: LowerInteger.java,v 1.1 2004-06-09 09:49:55 tjc Exp $
 **/

public class LowerInteger {
  /**
   *  Create a new LowerInteger from the given Integer.
   **/
  public LowerInteger (final Integer position) {
    this.position = position.intValue ();
  }

  /**
   *  Create a new LowerInteger from the given int.
   **/
  public LowerInteger (final int position) {
    this.position = position;
  }

  /**
   *  Return the position that was passed to the constructor.
   **/
  public int getPosition () {
    return position;
  }

  /**
   *  Return a String representing this object.  It will look something like
   *  this: "<100"
   **/
  public String toString () {
    return "<" + getPosition ();
  }

  /**
   *  The position that was passed to the constructor
   **/
  private int position;
}


