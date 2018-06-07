/* EntryVector.java
 *
 * created: Sat Oct 17 1998
 *
 * This file is part of Artemis
 * 
 * Copyright(C) 1998,1999,2000  Genome Research Limited
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or(at your option) any later version.
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/EntryVector.java,v 1.1 2004-06-09 09:44:27 tjc Exp $
 */

package uk.ac.sanger.artemis;

import java.util.Vector;

/**
 *  This class is a Vector of Entry objects.
 *
 *  @author Kim Rutherford
 *  @version $Id: EntryVector.java,v 1.1 2004-06-09 09:44:27 tjc Exp $
 *
 **/

public class EntryVector 
{

  /** Storage for Entry objects. */
  private Vector vector = new Vector();

  /**
   *  Create a new(empty) EntryVector object.
   **/
  public EntryVector() 
  {
  }

  /**
   *  Appends the given Entry object to the vector if and only if it isn't
   *  already in the vector.
   **/
  protected void addElement(Entry entry) 
  {
    if(indexOf(entry) == -1) 
      vector.addElement(entry);
  }
  
  /**
   *  Appends the given Entry object to the vector if and only if it isn't
   *  already in the vector. (same as addElement()).
   **/
  protected void add(Entry entry) 
  {
    addElement(entry);
  }
  
  /**
   *  Performs the same function as Vector.elementAt()
   **/
  public Entry elementAt(int index) 
  {
    return(Entry) vector.elementAt(index);
  }

  /**
   *  Performs the same function as Vector.removeElement()
   **/
  protected boolean removeElement(Entry entry) 
  {
    return vector.removeElement(entry);
  }

  /**
   *  Return true if this object contains the given Entry.
   **/
  public boolean contains(Entry entry) 
  {
    if(indexOf(entry) == -1) 
      return false;
    else 
      return true;
  }

  /**
   *  Performs the same function as Vector.removeAllElements()
   **/
  public void removeAllElements() 
  {
    vector.removeAllElements();
  }

  /**
   *  Performs the same function as Vector.removeElement()
   **/
  public int indexOf(Entry entry) 
  {
    return vector.indexOf(entry);
  }

  /**
   *  Performs the same function as Vector.size()
   **/
  public int size() 
  {
    return vector.size();
  }
  
  /**
   *  Create a new EntryVector with the same contents as this one.
   **/
  public Object clone() 
  {
    final EntryVector return_vector = new EntryVector();
    return_vector.vector = (Vector)vector.clone();
    return return_vector;
  }

}
