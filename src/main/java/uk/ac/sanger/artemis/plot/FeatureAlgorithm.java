/* FeatureAlgorithm.java
 *
 * created: Wed Dec 16 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/plot/FeatureAlgorithm.java,v 1.1 2004-06-09 09:51:25 tjc Exp $
 */

package uk.ac.sanger.artemis.plot;

import uk.ac.sanger.artemis.Feature;

/**
 *  The FeatureAlgorithm class is the base class for algorithms that work
 *  on features.  A FeatureAlgorithm has a name and is specific to one
 *  Feature, meaning the algorithm can't change feature.
 *
 *  @author Kim Rutherford
 *  @version $Id: FeatureAlgorithm.java,v 1.1 2004-06-09 09:51:25 tjc Exp $
 **/

public abstract class FeatureAlgorithm extends Algorithm {
  /**
   *  Create a new FeatureAlgorithm object.
   *  @param feature The feature to do the calculation on.
   *  @param algorithm_name A String used to identify this algorithm to the
   *    user.
   *  @param algorithm_short_name A String used to identify this algorithm
   *    internally.  See the Algorithm constructor for more details.
   **/
  public FeatureAlgorithm (final Feature feature, final String algorithm_name,
                           final String algorithm_short_name) {
    super (algorithm_name, algorithm_short_name);
    this.feature = feature;
  }

  /**
   *  Returns the strand we will do the calculation on.
   **/
  public Feature getFeature () {
    return feature;
  }

  /**
   *  Return the value of the function between a pair of amino acids.
   *  @param start The start amino acid (included in the range).
   *  @param end The end amino acid (included in the range).
   *  @param values The results are returned in this array, hence it should be
   *    allocated at the size given by getValueCount ().
   **/
  public abstract void getValues (int start, int end, float [] values);

  /**
   *  Return the number of values a call to getValues () will return.
   **/
  public abstract int getValueCount ();

  /**
   *  The Feature we will do the calculation on.
   **/
  private Feature feature;
}


