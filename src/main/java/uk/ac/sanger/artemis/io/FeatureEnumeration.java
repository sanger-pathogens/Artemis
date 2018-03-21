/* FeatureEnumeration.java
 *
 * created: Sun Jan 31 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/io/FeatureEnumeration.java,v 1.1 2004-06-09 09:49:23 tjc Exp $
 */

package uk.ac.sanger.artemis.io;

import java.util.NoSuchElementException;

/**
 *  An object that implements the FeatureEnumeration interface generates a
 *  series of Feature objects, one at a time. Successive calls to the
 *  nextFeature method return successive elements of the series.
 *
 *  @author Kim Rutherford
 *  @version $Id: FeatureEnumeration.java,v 1.1 2004-06-09 09:49:23 tjc Exp $
 *
 **/

public interface FeatureEnumeration {
  /**
   *  Tests if this enumeration contains more features.
   *  @return true if and only if this enumeration object contains at least
   *    one more feature to provide; false otherwise.
   **/
  boolean hasMoreFeatures ();
  
  /**
   *  Returns the next Feature object of this enumeration if this enumeration
   *  object has at least one more element to provide.
   *  @return The next Feature of this enumeration.
   *  @exception NoSuchElementException If no more elements exist.
   **/
  Feature nextFeature ()
     throws NoSuchElementException;
}
