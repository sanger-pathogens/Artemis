/* FeaturePatternPredicate.java
 *
 * created: Tue Aug  6 2002
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2001  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/FeaturePatternPredicate.java,v 1.1 2004-06-09 09:44:43 tjc Exp $
 */

package uk.ac.sanger.artemis;

import uk.ac.sanger.artemis.sequence.AminoAcidSequence;

/**
 *  Each object of this class can be used to test Feature objects to see if
 *  they contain a given sequence pattern.  See FeaturePredicate.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: FeaturePatternPredicate.java,v 1.1 2004-06-09 09:44:43 tjc Exp $
 **/

public class FeaturePatternPredicate implements FeaturePredicate {
  /**
   *  Creata a new FeaturePatternPredicate object which will test for the
   *  presence of the given AminoAcidSequence in the translation of the Feature
   **/
  public FeaturePatternPredicate (final AminoAcidSequence motif_pattern) {
    this.motif_pattern = motif_pattern;
  }

  /**
   *  Test the given Feature against this FeaturePredicate
   *  @param feature The Feature to test the predicate against.
   *  @return Return true if and only if the given Feature contains the given
   *    sequence pattern.
   **/
  public boolean testPredicate (final Feature feature) {
    return motif_pattern.checkForMatch (feature.getTranslation ());
  }

  /**
   *  The AminoAcidSequence that was passed to the constructor.
   **/
  private AminoAcidSequence motif_pattern;
}
