/* HydrophobicityAlgorithm.java
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/plot/HydrophobicityAlgorithm.java,v 1.1 2004-06-09 09:51:34 tjc Exp $
 */

package uk.ac.sanger.artemis.plot;

import uk.ac.sanger.artemis.*;

/**
 *  Objects of this class have one useful method - getValues (), which takes a
 *  range of bases and returns a single floating point number.  The feature to
 *  use is set in the constructor.
 *
 *  @author Kim Rutherford
 *  @version $Id: HydrophobicityAlgorithm.java,v 1.1 2004-06-09 09:51:34 tjc Exp $
 **/

public class HydrophobicityAlgorithm extends HydroAlgorithm {
  /**
   *  Create a new HydrophobicityAlgorithm object.
   *  @param feature The feature to do future calculation on.
   **/
  public HydrophobicityAlgorithm (Feature feature) {
    super (feature, "Kyte-Doolittle Hydrophobicity", "hydrophobicity",
           data_array);
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
    return new Float (4);
  } 

  /**
   *  Return the minimum value of this algorithm.
   *  @return The minimum is -4..
   **/
  protected Float getMinimumInternal () {
    return new Float (-3);
  }

  /**
   *  Hydrophobicity value for each amino acid.  This array is indexed by
   *  comparing to amino_acid_array.  
   **/
  private final static float [] data_array = {
    -4.50F, -3.90F, -3.50F, -3.50F, -3.50F,
    -3.50F, -3.20F, -0.80F, -0.70F, -1.60F,
    -1.30F, +2.50F, -0.40F, +1.80F, +1.90F,
    -0.90F, +3.80F, +4.20F, +2.80F, +4.50F
  };

}
