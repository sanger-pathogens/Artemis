/* CoilFeatureAlgorithm.java
 *
 * created: Mon Dec 21 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/plot/CoilFeatureAlgorithm.java,v 1.1 2004-06-09 09:51:24 tjc Exp $
 **/

package uk.ac.sanger.artemis.plot;

import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.sequence.*;
import uk.ac.sanger.artemis.io.Range;

/**
 *  Objects of this class have one useful method - getValues (), which takes a
 *  range of bases and returns a single floating point number.  See
 *  "Predicting Coiled Coils from Protein Sequences", Science Vol. 252 page
 *  1162 for details of the alogrithm used.
 *
 *  @author Kim Rutherford
 *  @version $Id: CoilFeatureAlgorithm.java,v 1.1 2004-06-09 09:51:24 tjc Exp $
 **/

public class CoilFeatureAlgorithm extends FeatureAlgorithm {
  /**
   *  Create a new CoilFeatureAlgorithm object for the given Feature.
   **/
  public CoilFeatureAlgorithm (Feature feature) {
    super (feature, "Coiled Coils", "coiled_coil");
  }

  private static final int WINDOW_SIZE = 28;

  /**
   *  Return the coiled coil score of the (28) residues in the given range.
   *  The start and end values should be 28*3 bases apart (inclusive).
   *  @param start The start base (included in the range).
   *  @param end The end base (included in the range).
   *  @param values The one return value for this algorithm is returned in
   *    this array.
   **/
  public void getValues (int start, int end, float [] values) {
    
    final int number_of_weights = weight_array[0].length;

    final AminoAcidSequence translation = getFeature ().getTranslation ();

    final String translation_string =
      translation.toString ().substring (start, end + 1);

    final float [] scores = new float [WINDOW_SIZE];

    float max_score = -1.0F;

    int k = 0;

    for (int frame = 0 ; frame < 7 ; ++frame) {
      k = frame - start % 7;
      if (k < 0) {
        k = k + 7;
      }
      for (int j = 0 ; j < WINDOW_SIZE ; ++j) {
//        System.out.println (k + " " + frame + " " + j);

        final char this_char = translation_string.charAt (j);
        final int index =
          AminoAcidSequence.getSymbolIndex (this_char);

        final float this_score = weight_array [index][k];

        scores [j] = this_score;
          
        ++k;

        if (k >= 7) {
          k = 0;
        }
      }

      final float score = geometricMean (scores);

      if (score > max_score) {
        max_score = score;
      }
    }

    values[0] = probCoil (max_score);
  }

  /**
   *  Return the geometric mean of the floating point values in the argument.
   **/
  private float geometricMean (final float [] scores) {
    float total = 1;

    for (int i = 0 ; i < scores.length ; ++i) {
      total *= scores [i];
    }

    return (float) Math.pow (total, 1.0 / scores.length);
  }


  /**
   *  Returns probability of coiled-coil for the given score from the
   *  statistics in the original paper.
   **/
  private float probCoil (final float score) {
    final float gcc_mean = 1.63F;
    final float gg_mean = 0.77F;
    final float gcc_sd = 0.24F;
    final float gg_sd =0.20F;
    final float gcc = gauss (gcc_mean, gcc_sd, score);
    final float gg = gauss (gg_mean,  gg_sd,  score);

    return (float) (gcc/(30.*gg+gcc));
  }
  
  /**
   *  Calculates probability based on a Gaussian distribution.
   **/
  private float gauss (final float mean, final float sd, final float score) {
    return
      (float) (Math.pow (Math.E, -0.5 * Math.pow ((score - mean) / sd, 2)) /
      (sd * 2.0 * Math.PI));
  }
  
  /**
   *  Return the number of values a call to getValues () will return - one
   *  in this case.
   **/
  public int getValueCount () {
    return 1;
  }

  /**
   *  Return the default or optimal window size.
   *  @return null is returned if this algorithm doesn't have optimal window
   *    size.
   **/
  public Integer getDefaultWindowSize () {
    return new Integer (WINDOW_SIZE);
  }

  /**
   *  Return the default maximum window size for this algorithm.
   *  @return null is returned if this algorithm doesn't have maximum window
   *    size.
   **/
  public Integer getDefaultMaxWindowSize () {
    return new Integer (WINDOW_SIZE);
  }

  /**
   *  Return the default minimum window size for this algorithm.
   *  @return null is returned if this algorithm doesn't have minimum window
   *    size.
   **/
  public Integer getDefaultMinWindowSize () {
    return new Integer (WINDOW_SIZE);
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
   *  Return the maximum value of this algorithm.
   *  @return The maximum is 4.
   **/
  protected Float getMaximumInternal () {
    return new Float (1.01);
  } 

  /**
   *  Return the minimum value of this algorithm.
   *  @return The minimum is -4..
   **/
  protected Float getMinimumInternal () {
    return new Float (-0.01);
  }

  /**
   *  Return the average value of function over the whole strand.
   *  @return null is returned if this algorithm doesn't have an average or if
   *    the average can't be calculated.
   **/
  public Float getAverage () {
    return null;
  }

  /**
   *  The weighting matrix from the Science paper.
   **/
  private final static float [] [] weight_array = {
    { 1.297F,1.551F,1.084F,2.612F,0.377F,1.248F,0.877F, },  // Ala A    7.59
    { 0.659F,1.163F,1.210F,0.031F,1.358F,1.937F,1.798F, },  // arg R    5.39  
    { 0.835F,1.475F,1.534F,0.039F,1.722F,2.456F,2.280F, },  // Asn N    4.25
    { 0.030F,2.352F,2.268F,0.237F,0.663F,1.620F,1.448F, },  // Asp D    5.03
    { 0.824F,0.022F,0.308F,0.152F,0.180F,0.156F,0.044F, },  // Cys C    1.86
    { 0.179F,2.114F,1.778F,0.631F,2.550F,1.578F,2.526F, },  // Gln Q    4.27
    { 0.262F,3.496F,3.108F,0.998F,5.685F,2.494F,3.048F, },  // Glu E    6.10
    { 0.045F,0.275F,0.578F,0.216F,0.211F,0.426F,0.156F, },  // Gly G    7.10
    { 0.347F,0.275F,0.679F,0.395F,0.294F,0.579F,0.213F, },  // His H    2.25
    { 2.597F,0.098F,0.345F,0.894F,0.514F,0.471F,0.431F, },  // Ile I    5.35
    { 3.167F,0.297F,0.398F,3.902F,0.585F,0.501F,0.483F, },  // Leu L    9.33
    { 1.375F,2.639F,1.763F,0.191F,1.815F,1.961F,2.795F, },  // Lys K    5.72
    { 2.240F,0.370F,0.480F,1.409F,0.541F,0.772F,0.663F, },  // Met M    2.34
    { 0.531F,0.076F,0.403F,0.662F,0.189F,0.106F,0.013F, },  // Phe F    3.88
    { 0.0F,  0.008F,0.0F,  0.013F,0.0F,  0.0F,  0.0F,   },  // Pro P    5.28  
    { 0.382F,0.583F,1.052F,0.419F,0.525F,0.916F,0.628F, },  // Ser S    7.28
    { 0.169F,0.702F,0.955F,0.654F,0.791F,0.843F,0.647F, },  // Thr T    5.97
    { 0.240F,0.0F,  0.0F,  0.456F,0.019F,0.0F,  0.0F,   },  // Trp W    1.41  
    { 1.417F,0.090F,0.122F,1.659F,0.190F,0.130F,0.155F, },  // Tyr Y    3.16
    { 1.665F,0.403F,0.386F,0.949F,0.211F,0.342F,0.360F, },  // Val V    6.42
    { 0.000F,0.000F,0.000F,0.000F,0.000F,0.000F,0.000F, },  // Opl *    0.00
    { 0.000F,0.000F,0.000F,0.000F,0.000F,0.000F,0.000F, },  // Ocr #    0.00
    { 0.000F,0.000F,0.000F,0.000F,0.000F,0.000F,0.000F, },  // Amb +    0.00
    { 0.000F,0.000F,0.000F,0.000F,0.000F,0.000F,0.000F, },  // --- .    0.00
  };
}
