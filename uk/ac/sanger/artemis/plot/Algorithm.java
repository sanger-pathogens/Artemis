/* Algorithm.java
 *
 * created: Tue Dec 15 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/plot/Algorithm.java,v 1.2 2009-04-07 09:11:29 tjc Exp $
 */

package uk.ac.sanger.artemis.plot;

import uk.ac.sanger.artemis.Options;

/**
 *  This class represents an algorithm that can be plotted.
 *
 *  @author Kim Rutherford
 *  @version $Id: Algorithm.java,v 1.2 2009-04-07 09:11:29 tjc Exp $
 **/

public abstract class Algorithm {
  
  
  /**
   *  The name of this algorithm, as passed to the constructor.
   **/
  private String algorithm_name;

  /**
   *  The short (one word) name of this algorithm, as passed to the
   *  constructor.
   **/
  private String algorithm_short_name;

  /**
   *  Set by disableMaxAndMin () and enableMaxAndMin ().
   **/
  protected boolean max_min_disabled = false;

  /**
   *  Set by the constructor by looking at the options with this name:
   *  getAlgorithmShortName () + "_default_window_size"
   **/
  private Integer options_window_size = null;

  /**
   *  Set by the constructor by looking at the options with this name:
   *  getAlgorithmShortName () + "_default_min_window_size"
   **/
  private Integer options_max_window_size = null;

  /**
   *  Set by the constructor by looking at the options with this name:
   *  getAlgorithmShortName () + "_default_min_window_size"
   **/
  private Integer options_min_window_size = null;
  
  private boolean userMaxMin = false;
  private float userMin = Float.MIN_VALUE;
  private float userMax = Float.MAX_VALUE;
  
  /**
   *  Create a new Algorithm object.
   *  @param algorithm_name A String used to identify this algorithm to the
   *    user.
   *  @param algorithm_short_name A String used to identify this algorithm
   *    internally.  This is used to make the names of the keys for looking up
   *    the window size options.  eg. long name: "GC Content (%)", short name:
   *    gc_content, which gives these options names:
   *    gc_content_default_window_size, gc_content_default_max_window_size and
   *    gc_content_default_min_window_size.
   **/
  public Algorithm (final String algorithm_name,
                    final String algorithm_short_name) 
  {
    this.algorithm_name = algorithm_name;
    this.algorithm_short_name = algorithm_short_name;

    options_min_window_size =
      Options.getOptions ().getIntegerProperty (getAlgorithmShortName () +
                                                "_default_min_window");

    if (options_min_window_size != null &&
        options_min_window_size.intValue () < 1) {
      options_min_window_size = null;
    }

    options_max_window_size =
      Options.getOptions ().getIntegerProperty (getAlgorithmShortName () +
                                                "_default_max_window");

    if (options_max_window_size != null) {
      if (options_min_window_size == null) {
        if (options_max_window_size.intValue () < 1) {
          options_max_window_size = null;
        }
      } else {
        if (options_max_window_size.intValue () <
            options_min_window_size.intValue ()) {
          options_max_window_size = options_min_window_size;
        }
      }
    }

    options_window_size =
      Options.getOptions ().getIntegerProperty (getAlgorithmShortName () +
                                                "_default_window_size");

    if (options_window_size != null) {
      if (options_min_window_size == null) {
        if (options_window_size.intValue () < 1) {
          options_window_size = null;
        }
      } else {
        if (options_window_size.intValue () <
            options_min_window_size.intValue ()) {
          options_window_size = options_min_window_size;
        }
      }
    }
  }

  /**
   *  Return the name of this algorithm.
   **/
  public String getAlgorithmName () 
  {
     return algorithm_name;
  }
  
  /**
   *  Return the name of this algorithm.
   **/
  public void setAlgorithmName (String algorithm_name) 
  {
     this.algorithm_name = algorithm_name;
  }

  /**
   *  Return the short (one word) name of this algorithm, as passed tp the
   *  constructor.
   **/
  public String getAlgorithmShortName () 
  {
     return algorithm_short_name;
  }

  /**
   *  Return the default or optimal window size.
   *  @return null is returned if this algorithm doesn't have optimal window
   *    size.
   **/
  public Integer getDefaultWindowSize () 
  {
    return options_window_size;
  } 

  /**
   *  Return the default maximum window size for this algorithm.
   *  @return null is returned if this algorithm doesn't have maximum window
   *    size.
   **/
  public Integer getDefaultMaxWindowSize ()
  {
    return options_max_window_size;
  } 

  /**
   *  Return the default minimum window size for this algorithm.
   *  @return null is returned if this algorithm doesn't have minimum window
   *    size.
   **/
  public Integer getDefaultMinWindowSize () 
  {
    return options_min_window_size;
  } 

  /**
   *  Return the default or optimal step size.
   *  @return null is returned if this algorithm doesn't have optimal step
   *    size.
   **/
  public Integer getDefaultStepSize (int window_size) 
  {
    return null;
  } 

  /**
   *  Return the maximum value of this algorithm.    Subclasses of this class
   *  should override getMaximumInternal () rather than this method.
   *  @return null is returned if this algorithm doesn't have a fixed maximum
   *    or if maxAndMinDisabled () is true.
   **/
  final public Float getMaximum () 
  {
    if (scalingFlag ())
      return null;
    else if(isUserMaxMin())
      return getUserMax();
    else
      return getMaximumInternal ();
  } 

  /**
   *  Return the maximum value of this algorithm.
   *  @return null is returned if this algorithm doesn't have a fixed maximum.
   **/
  protected Float getMaximumInternal () 
  {
    return null;
  } 

  /**
   *  Return the minimum value of this algorithm.  Subclasses of this class
   *  should override getMinimumInternal () rather than this method.
   *  @return null is returned if this algorithm doesn't have a fixed minimum
   *    or if maxAndMinDisabled () is true.
   **/
  final public Float getMinimum () 
  {
    if (scalingFlag ()) 
      return null;
    else if(isUserMaxMin())
      return getUserMin();
    else
      return getMinimumInternal ();
  } 

  /**
   *  Return the minimum value of this algorithm.
   *  @return null is returned if this algorithm doesn't have a fixed minimum.
   **/
  protected Float getMinimumInternal () 
  {
    return null;
  } 

  /**
   *  Return the average value of function over the whole strand.
   *  @return null is returned if this algorithm doesn't have an average or if
   *    the average can't be calculated.
   **/
  public Float getAverage () 
  {
    return null;
  }

  /**
   *  Force getMinimum () and getMaximum () to return null if and only if
   *  enable_scaling is true.
   **/
  public void setScalingFlag (final boolean enable_scaling) 
  {
    max_min_disabled = enable_scaling;
  }

  /**
   *  Return true if and only if getMaximum () and getMinimum () should return
   *  null.
   **/
  public boolean scalingFlag ()
  {
    return max_min_disabled;
  }
  
  public float getUserMin()
  {
    return userMin;
  }

  public void setUserMin(float userMin)
  {
    this.userMin = userMin;
  }

  public float getUserMax()
  {
    return userMax;
  }

  public void setUserMax(float userMax)
  {
    this.userMax = userMax;
  }
  
  public boolean isUserMaxMin()
  {
    return userMaxMin;
  }

  public void setUserMaxMin(boolean userMaxMin)
  {
    this.userMaxMin = userMaxMin;
  }
}
