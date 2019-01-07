/* UpperInteger.java
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/UpperInteger.java,v 1.1 2004-06-09 09:50:39 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

/**
 *  This class is used to represent a upper end of a range like this 100..>200
 *  @author Kim Rutherford
 *  @version $Id: UpperInteger.java,v 1.1 2004-06-09 09:50:39 tjc Exp $
 **/

public class UpperInteger {
  /**
   *  Create a new UpperInteger from the given int.
   **/
  protected UpperInteger (final Integer position) {
    this.position = position.intValue ();
  }

  /**
   *  Create a new UpperInteger from the given Integer.
   **/
  protected UpperInteger (final int position) {
    this.position = position;
  }

  /**
   *  Return the position that was passed to the constructor.
   **/
  protected int getPosition () {
    return position;
  }

  /**
   *  Return a String representing this.  It will look something like
   *  this: ">100"
   **/
  public String toString () {
    return ">" + getPosition ();
  }

  /** position that was passed to the constructor */
  private int position;
}
