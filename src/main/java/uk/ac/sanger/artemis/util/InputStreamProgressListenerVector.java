/* InputStreamProgressListenerVector.java
 *
 * created: Wed Aug  7 2002
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 2001  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/util/InputStreamProgressListenerVector.java,v 1.1 2004-06-09 09:53:04 tjc Exp $
 */

package uk.ac.sanger.artemis.util;

import java.util.Vector;

/**
 *  A Vector of InputStreamProgressListener objects.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: InputStreamProgressListenerVector.java,v 1.1 2004-06-09 09:53:04 tjc Exp $
 **/

public class InputStreamProgressListenerVector {
  /**
   *  Create a new (empty) vector of InputStreamProgressListener objects.
   **/
  public InputStreamProgressListenerVector () {
    vector = new Vector ();
  }

  /**
   *  Performs the same function as Vector.addElement ()
   **/
  public void add (final InputStreamProgressListener node) {
    vector.addElement (node);
  }

  /**
   *  Performs the same function as Vector.elementAt ()
   **/
  public InputStreamProgressListener elementAt (final int index) {
    return (InputStreamProgressListener) vector.elementAt (index);
  }

  /**
   *  Performs the same function as Vector.size ()
   **/
  public int size () {
    return vector.size ();
  }

  /**
   *  Storage for InputStreamProgressListener objects.
   **/
  private Vector vector;
}
