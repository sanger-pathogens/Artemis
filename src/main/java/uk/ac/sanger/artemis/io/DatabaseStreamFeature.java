/* DatbaseStreamFeature.java
 *
 * created: Mar 2005
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2005  Genome Research Limited
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
 */

package uk.ac.sanger.artemis.io;

/**
 *  This is an implementation of Feature that can read and write itself to a
 *  CHADO stream.
 *
 *  @version $Id: DatabaseStreamFeature.java,v 1.7 2007-06-05 09:33:41 tjc Exp $
 **/

public class DatabaseStreamFeature
    extends GFFStreamFeature
    implements DocumentFeature, StreamFeature, ComparableFeature
{

  public DatabaseStreamFeature(final Key key, final Location location,
                          final QualifierVector qualifiers)
  {
    super(key, location, qualifiers);
  }

  /**
   *  Create a new DatabaseStreamFeature with the same key, location and
   *  qualifiers as the given feature.  The feature should be added to an
   *  Entry (with Entry.add ()).
   *  @param feature The feature to copy.
   **/
  public DatabaseStreamFeature (final Feature feature) 
  {
    super(feature);
  }

  /**
   *  Return the reference of a new copy of this Feature.
   **/
  public Feature copy() 
  {
    final Feature return_value = new DatabaseStreamFeature(this);
    return return_value;
  }

}