/* BaseAlgorithm.java
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/plot/BaseAlgorithm.java,v 1.1 2004-06-09 09:51:16 tjc Exp $
 */

package uk.ac.sanger.artemis.plot;

import uk.ac.sanger.artemis.sequence.*;

/**
 *  The BaseAlgorithm class is the base class for algorithms that work
 *  directly on bases.  A BaseAlgorithm has a name and is specific to one
 *  Strand of DNA, meaning the algorithm can't change strand part way along.
 *
 *  @author Kim Rutherford
 *  @version $Id: BaseAlgorithm.java,v 1.1 2004-06-09 09:51:16 tjc Exp $
 **/

public abstract class BaseAlgorithm extends Algorithm {
  /**
   *  Create a new BaseAlgorithm object.
   *  @param strand The strand to do the calculation on.
   *  @param algorithm_name A String used to identify this algorithm to the
   *    user.
   *  @param algorithm_short_name A String used to identify this algorithm
   *    internally.  See the Algorithm constructor for more details.
   **/
  public BaseAlgorithm (final Strand strand, final String algorithm_name,
                        final String algorithm_short_name) {
    super (algorithm_name, algorithm_short_name);
    this.bases = strand.getBases ();

    if (strand.isForwardStrand ()) {
      forward_flag = true;
    } else {
      forward_flag = false;
    }
  }

  /**
   *  Return the Bases object of the Strand that was passed to the
   *  constructor.
   **/
  public Bases getBases () {
    return bases;
  }

  /**
   *  Returns the strand we will do the calculation on.
   **/
  public Strand getStrand () {
    if (forward_flag ^ rev_comp_display) {
      return getBases ().getForwardStrand ();
    } else {
      return getBases ().getReverseStrand ();
    }
  }
  
  /**
   *  If rev_comp_display is true all calculations will be performs on the
   *  opposite Strand to the strand that was passed to the constructor.
   **/
  public void setRevCompDisplay (final boolean rev_comp_display) {
    this.rev_comp_display = rev_comp_display;
  }

  /**
   *  Returns true if the FeatureDisplay is reverse complemented.  All
   *  calculations should be performed on the opposite Strand to the strand
   *  that was passed to the constructor.
   **/
  public boolean isRevCompDisplay () {
    return rev_comp_display;
  }

  /**
   *  Return the value of the function between a pair of bases.
   *  @param start The start base (included in the range).
   *  @param end The end base (included in the range).
   *  @param values The results are returned in this array, hence it should be
   *    allocated at the size given by getValueCount ().
   **/
  public abstract void getValues (int start, int end, final float [] values);

  /**
   *  Return the number of values a call to getValues () will return.
   **/
  public abstract int getValueCount ();
  
  /**
   *  The Bases we will do the calculation on.
   **/
  private Bases bases;

  /**
   *  If rev_comp_display is true all calculations will be performed on the
   *  opposite Strand to the strand that was passed to the constructor.
   **/
  private boolean rev_comp_display = false;

  /**
   *  true if and only if the calculations should be done on the forward
   *  Strand.
   **/
  private boolean forward_flag;
}
