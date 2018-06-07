/* FeatureVector.java
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
 */

package uk.ac.sanger.artemis.io;

import java.util.Vector;

/**
 *  This class implements a Vector of Feature objects.
 *  @author Kim Rutherford
 *  @version $Id: FeatureVector.java,v 1.2 2004-11-24 11:55:52 tjc Exp $
 */

public class FeatureVector extends Vector<Feature>
{
  private static final long serialVersionUID = 1L;

  /**
   *  Create a new vector of Feature objects with an initial capacity of 100.
   **/
  public FeatureVector() 
  {
    super(100);
  }

  /**
   *  Performs the same function as Vector.elementAt()
   **/
  public Feature featureAt(int index) 
  {
    return elementAt(index);
  }
}

