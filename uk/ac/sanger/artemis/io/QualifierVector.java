/* QualifierVector.java
 *
 * created: Tue Oct 13 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/QualifierVector.java,v 1.7 2006-12-06 17:27:12 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import java.util.Vector;

/**
 *  This class implements a vector of Qualifier objects.  It behaves
 *  differently to the Vector class (see addElement() and replaceElement()).
 *
 *  @author Kim Rutherford
 *  @version $Id: QualifierVector.java,v 1.7 2006-12-06 17:27:12 tjc Exp $
 *
 */

public class QualifierVector extends Vector<Qualifier>
{
  private static final long serialVersionUID = 1L;

  /**
   *  Create a new (empty) vector of Qualifier objects.
   */
  public QualifierVector() 
  {
    super(7);
  }


  /**
   *  Add the values from the given qualifier to the Qualifier object with the
   *  same name in this QualifierVector or otherwise add a copy of the
   *  argument.
   *  @param qualifier This object contians name and values to add.
   *  @return The Qualifier that was changed or created.
   **/
  public Qualifier addQualifierValues(final Qualifier qualifier) 
  {
    if(qualifier.getName() == null) 
      throw new Error("");

    final int index_of_qualifier =
      indexOfQualifierWithName(qualifier.getName());

    if(index_of_qualifier == -1) 
    {
      addElement(qualifier.copy());
      return null;
    } 
    else
    {
      final Qualifier current_qualifier = elementAt(index_of_qualifier);
      current_qualifier.addValues(qualifier.getValues());
      return current_qualifier;
    }
  }

  /**
   *  Add the given Qualifier to this QualifierVector, replacing any exisiting
   *  Qualifier with the same name.
   *  @param qualifier The Qualifier to add.
   **/
  public void setQualifier(final Qualifier qualifier) 
  {
    final int index = indexOfQualifierWithName(qualifier.getName());

    if(index == -1) 
      addElement(qualifier);
    else 
    {
      removeQualifierByName(qualifier.getName());
      addElement(qualifier.copy());
    }
  }

  /**
   *  Remove the Qualifier with the given name.  If there is no Qualifier with
   *  that name return immediately.
   *  @param name The Qualifier name to look for.
   **/
  public void removeQualifierByName(final String name) 
  {
    final int index = indexOfQualifierWithName(name);

    if(index != -1)
      removeElementAt(index);
  }


  /**
   *  Returns true if and only if this QualifierVector contains a qualifier
   *  with given name.
   **/
  public boolean contains(final String name) 
  {
    if(indexOfQualifierWithName (name) == -1) 
      return false;
    else 
      return true;
  }
  
  /**
   *  Returns the index of the Qualifier in this QualifierVector with the
   *  given name of -1 if no such Qualifier exists.
   *  @param name The Qualifier name to look for.
   **/
  public int indexOfQualifierWithName(String name) 
  {
    final int vsize = size();
    for(int i = 0; i < vsize; ++i) 
    {
      if(elementAt(i).getName().equals(name)) 
        return i;
    }
    return -1;
  }
  
  /**
   *  Returns the Qualifier in this QualifierVector with the given name or
   *  null if no such Qualifier exists.
   *  @param name The Qualifier name to look for.
   **/
  public Qualifier getQualifierByName(String name) 
  {
    final int index_of_named_qualifier = indexOfQualifierWithName(name);
    if(index_of_named_qualifier == -1)
      return null; 
    else
      return elementAt(index_of_named_qualifier);
  }
  
  /**
   *  Return the reference of a new copy of this QualifierVector.  All of the
   *  Qualifier objects in the vector will be copied too.
   **/
  public QualifierVector copy() 
  {
    return (QualifierVector)super.clone();
//  final QualifierVector return_vector = new QualifierVector();
//  final int vsize = size();  
//  for(int i = 0 ; i < vsize; ++i) 
//    return_vector.addElement(((Qualifier)elementAt(i)).copy());

//  return return_vector;
  }

}
