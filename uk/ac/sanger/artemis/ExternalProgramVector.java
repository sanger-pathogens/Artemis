/* ExternalProgramVector.java
 *
 * created: Tue Jan 26 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/ExternalProgramVector.java,v 1.1 2004-06-09 09:44:35 tjc Exp $
 **/

package uk.ac.sanger.artemis;

import java.util.Vector;

/**
 *  This class implements a Vector of ExternalProgram objects.
 *
 *  @author Kim Rutherford
 *  @version $Id: ExternalProgramVector.java,v 1.1 2004-06-09 09:44:35 tjc Exp $
 **/

public class ExternalProgramVector {
  /**
   *  Create a new (empty) ExternalProgramVector.
   **/
  public ExternalProgramVector () {

  }

  /**
   *  Performs the same function as Vector.addElement ()
   */
  public void add (final ExternalProgram program) {
    vector.addElement (program);
  }

  /**
   *  Performs the same function as Vector.elementAt ()
   */
  public ExternalProgram elementAt (final int index) {
    return (ExternalProgram) vector.elementAt (index);
  }

  /**
   *  Return true if this object contains the given ExternalProgram.
   **/
  public boolean contains (final ExternalProgram program) {
    if (indexOf (program) == -1) {
      return false;
    } else {
      return true;
    }
  }

  /**
   *  Performs the same function as Vector.removeElement ()
   **/
  public int indexOf (final ExternalProgram program) {
    return vector.indexOf (program);
  }

  /**
   *  Performs the same function as Vector.size ()
   */
  public int size () {
    return vector.size ();
  }

  /**
   *  Storage for ExternalProgram objects.
   */
  private Vector vector = new Vector ();
}


