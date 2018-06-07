/* HydrophilicityAlgorithm.java
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/plot/HydrophilicityAlgorithm.java,v 1.1 2004-06-09 09:51:33 tjc Exp $
 */

package uk.ac.sanger.artemis.plot;

import uk.ac.sanger.artemis.*;

/**
 *  HydrophilicityAlgorithm class
 *
 *  @author Kim Rutherford
 *  @version $Id: HydrophilicityAlgorithm.java,v 1.1 2004-06-09 09:51:33 tjc Exp $
 **/

public class HydrophilicityAlgorithm extends HydroAlgorithm {
  /**
   *  Create a new HydrophilicityAlgorithm object.
   *  @param feature The feature to do future calculation on.
   **/
  public HydrophilicityAlgorithm (Feature feature) {
    super (feature, "Hopp-Woods Hydrophilicity", "hydrophilicity", data_array);
  }

  /**
   *  Return the number of values a call to getValues () will return - three
   *  in this case.
   **/
  public int getValueCount () {
    return 1;
  }

  /**
   *  Return the maximum value of this algorithm.
   *  @return The maximum is 4.
   **/
  protected Float getMaximumInternal () {
    return new Float (3);
  } 

  /**
   *  Return the minimum value of this algorithm.
   *  @return The minimum is -2.
   **/
  protected Float getMinimumInternal () {
    return new Float (-2);
  }

  /**
   *  Hydrophobicity value for each amino acid.  This array is indexed by
   *  comparing to amino_acid_array.  
   **/
  private final static float [] data_array = {
    +3.00F, +3.00F, +3.00F, +0.20F, +0.20F,
    +3.00F, -0.50F, +0.30F, -0.40F, +0.00F,
    -2.30F, -1.00F, +0.00F, -0.50F, -1.30F,
    -3.40F, -1.80F, -1.50F, -2.50F, -1.80F
  };
}
