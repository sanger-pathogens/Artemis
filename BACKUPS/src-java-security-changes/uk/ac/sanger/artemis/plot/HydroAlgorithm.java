/* HydroAlgorithm.java
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/plot/HydroAlgorithm.java,v 1.2 2005-09-06 07:23:15 tjc Exp $
 */

package uk.ac.sanger.artemis.plot;

import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.sequence.Bases;

/**
 *  This is the base class for HydrophobicityAlgorithm and
 *  HydrophilicityAlgorithm.
 *
 *  @author Kim Rutherford
 *  @version $Id: HydroAlgorithm.java,v 1.2 2005-09-06 07:23:15 tjc Exp $
 **/

public abstract class HydroAlgorithm extends FeatureAlgorithm {
  /**
   *  Create a new HydroAlgorithm object.
   *  @param feature The feature to do future calculation on.
   *  @param algorithm_name A String used to identify this algorithm to the
   *    user.
   *  @param data_array An array of data values, one value for each amino acid.
   *    The data should correspond to the amino acids in this order:
   *    RKDQNEHSTPYCGAMWLVFI.
   **/
  public HydroAlgorithm (final Feature feature, final String algorithm_name,
                         final String algorithm_short_name,
                         final float [] data_array) {
    super (feature, algorithm_name, algorithm_short_name);
    this.data_array = data_array;
  }
  
  /**
   *  Return the hydrophobicity between a pair of amino acids.
   *  @param start The start amino acid (included in the range), indexed from
   *    zero.
   *  @param end The end amino acid (included in the range), indexed from 0.
   *  @param values The one return value for this algorithm is returned in
   *    this array.
   **/
  public void getValues (int start, int end, float [] values) {
    final String translation =
      getFeature().getTranslation().toString().substring(start, end).toLowerCase();

//  System.out.println(translation+"   "+start + " " + end + " " +
//                            getFeature ().getTranslation ().length ());

    float total = 0;

    for (int i = 0 ; i < translation.length () ; ++i) {
      total += getCodonValue (data_array, translation.charAt (i));
    }

    values[0] = total / translation.length ();

//  System.out.println("foo: " + values[0] + " " + total + " " +
//                     translation.length ());
  }

  /**
   *  Return the default or optimal window size.
   *  @return null is returned if this algorithm doesn't have optimal window
   *    size.
   **/
  public Integer getDefaultWindowSize () {
    final Integer super_window_size = super.getDefaultWindowSize ();
    if (super_window_size != null) {
      // the superclass version of getDefaultWindowSize () returns non-null
      // iff the user has set the window size in the options file
      return super_window_size;
    }
    return new Integer (7);
  } 

  /**
   *  Return the default maximum window size for this algorithm.
   *  @return null is returned if this algorithm doesn't have maximum window
   *    size.
   **/
  public Integer getDefaultMaxWindowSize () {
    final Integer super_max_window_size = super.getDefaultMaxWindowSize ();
    if (super_max_window_size != null) {
      // the superclass version of getDefaultMaxWindowSize () returns non-null
      // iff the user has set the max window size in the options file
      return super_max_window_size;
    }
    return new Integer (100);
  } 

  /**
   *  Return the default minimum window size for this algorithm.
   *  @return null is returned if this algorithm doesn't have minimum window
   *    size.
   **/
  public Integer getDefaultMinWindowSize () {
    final Integer super_min_window_size = super.getDefaultMinWindowSize ();
    if (super_min_window_size != null) {
      // the superclass version of getDefaultMinWindowSize () returns non-null
      // iff the user has set the minimum window size in the options file
      return super_min_window_size;
    }
    return new Integer (7);
  } 

  /**
   *  Return the default or optimal step size.
   *  @return null is returned if this algorithm doesn't have optimal step
   *    size.
   **/
  public Integer getDefaultStepSize (int window_size) {
    return new Integer (1);
  } 

  /**
   *  Return the average value of function over the whole strand.
   *  @return null is returned if this algorithm doesn't have an average or if
   *    the average can't be calculated.
   **/
  public Float getAverage () {
    return new Float (0);
  } 

  /**
   *  Return the value in data_array that corresponds to the given amino acid
   *  symbol (uses amino_acid_array to do this).
   **/
  private float getCodonValue (float [] data_array, char amino_acid_symbol) {
    final int index = amino_acid_array.indexOf (amino_acid_symbol);

//  System.out.println("getCodonValue (): " + index + " " +
//                      amino_acid_symbol) ; // + " " + data_array[index]);
    
    if (index == -1) {
      // return the average value in this case.
      return 0;
    } else {
      return data_array[index];
    }
  }
  
  /**
   *  Data values for each amino acid.  This array is indexed by comparing to
   *  amino_acid_array.
   **/
  /* final */ private float [] data_array;

  /**
   *  The index of the amino acid symbol in this String gives the index to
   *  look up in data_array to get the hydrophobicity value.  eg. the symbol
   *  "r" has index 0 hence has value +3.00
   **/
  private final static String amino_acid_array = "rkdqnehstpycgamwlvfi";
}


