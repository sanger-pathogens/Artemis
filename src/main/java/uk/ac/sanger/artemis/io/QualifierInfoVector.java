/* QualifierInfoVector.java
 *
 * created: Sun Feb 21 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/QualifierInfoVector.java,v 1.1 2004-06-09 09:50:11 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import java.util.Vector;

/**
 *  A Vector of QualifierInfo objects.
 *
 *  @author Kim Rutherford
 *  @version $Id: QualifierInfoVector.java,v 1.1 2004-06-09 09:50:11 tjc Exp $
 **/

public class QualifierInfoVector {
  /**
   *  Create a new (empty) QualifierInfoVector object.
   **/
  public QualifierInfoVector () {

  }

  /**
   *  Performs the same function as Vector.addElement ()
   */
  public void addElement (final QualifierInfo info) {
    vector.addElement (info);
  }
  
  /**
   *  Performs the same function as addElement ()
   */
  public void add (final QualifierInfo info) {
    vector.addElement (info);
  }
  
  /**
   *  Performs the same function as Vector.setElementAt ()
   **/
  public void setElementAt (final QualifierInfo qualifier_info,
                            final int index) {
    vector.setElementAt (qualifier_info, index);
  }

  /**
   *  Performs the same function as Vector.elementAt ()
   */
  public QualifierInfo elementAt (final int index) {
    return (QualifierInfo) vector.elementAt (index);
  }

  /**
   *  Return a new copy of this object.
   **/
  public QualifierInfoVector copy () { 
    final QualifierInfoVector new_vector = new QualifierInfoVector ();

    new_vector.vector = (Vector) vector.clone ();

    return new_vector;
  }

  /**
   *  Performs the same function as Vector.size ()
   */
  public int size () {
    return vector.size ();
  }
  
  /**
   *  Storage for QualifierInfo objects.
   */
  private Vector vector = new Vector ();
}


