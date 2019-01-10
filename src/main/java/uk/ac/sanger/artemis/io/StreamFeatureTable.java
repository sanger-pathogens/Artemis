/* StreamFeatureTable.java
 *
 * created: Wed Sep 15 1999
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
 */

package uk.ac.sanger.artemis.io;

import java.io.IOException;
import java.io.Writer;


/**
 *  This object contains all the features in an entry.  This object can write
 *  itself to a stream.
 *
 *  @author Kim Rutherford
 */

class StreamFeatureTable extends FeatureTable {
  /**
   *  Create a new empty StreamFeatureTable object.
   **/
  StreamFeatureTable () {

  }

  /**
   *  Create a new StreamFeatureTable object from the given FeatureTable.
   **/
  StreamFeatureTable (final FeatureTable feature_table) {
    final FeatureEnumeration feature_enum = feature_table.features ();

    while (feature_enum.hasMoreFeatures ()) {
      final Feature new_feature = feature_enum.nextFeature ();

      add (new_feature);
    }
  }

  /**
   *  Write this StreamFeatureTable to the given stream.
   *  @param writer The stream to write to.
   **/
  public void writeToStream (final Writer writer)
      throws IOException {

    final FeatureEnumeration feature_enum = features ();

    while (feature_enum.hasMoreFeatures ()) {
      final StreamFeature this_feature =
        (StreamFeature) feature_enum.nextFeature ();

      this_feature.writeToStream (writer);

      Thread.yield ();
    }
  }
}
