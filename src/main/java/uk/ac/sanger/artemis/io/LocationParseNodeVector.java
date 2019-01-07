/* LocationParseNodeVector.java
 *
 * created: Wed Oct  7 1998
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
 */

package uk.ac.sanger.artemis.io;

import java.util.Vector;

/**
 *  LocationParseNodeVector class is a Vector of objects of type
 *  LocationParseNode.
 *  @author Kim Rutherford
 * */

public class LocationParseNodeVector extends Vector<LocationParseNode> 
{
  private static final long serialVersionUID = 1L;

  /**
   *  Add the given node to the end of the vector.
   **/
  protected void addElementAtEnd (LocationParseNode node) 
  {
    insertElementAt (node, size ());
  }
}


