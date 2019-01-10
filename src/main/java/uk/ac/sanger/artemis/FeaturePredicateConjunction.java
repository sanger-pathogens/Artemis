/* FeaturePredicateConjunction.java
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/FeaturePredicateConjunction.java,v 1.1 2004-06-09 09:44:45 tjc Exp $
 */

package uk.ac.sanger.artemis;

/**
 *  A class for combining two predicates.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: FeaturePredicateConjunction.java,v 1.1 2004-06-09 09:44:45 tjc Exp $
 **/

public class FeaturePredicateConjunction implements FeaturePredicate {
  /**
   *  Create a new FeaturePredicateConjunction object.
   *  @param type the type of conjunction to use when testPredicate() is
   *    called.
   **/
  public FeaturePredicateConjunction (final FeaturePredicate predicate1,
                                      final FeaturePredicate predicate2,
                                      final int type) {
    this.type = type;

    predicates.add (predicate1);
    predicates.add (predicate2);

    if (!(type == OR || type == AND)) {
      throw new Error ("internal error - illegal type given to " +
                       "FeaturePredicateConjunction constructor");
    }
  }

  /**
   *  Create a new FeaturePredicateConjunction object.
   *  @param type the type of conjunction to use when testPredicate() is
   *    called.
   **/
  public FeaturePredicateConjunction (final FeaturePredicateVector predicates,
                                      final int type) {
    this.predicates = predicates.copy ();
    this.type = type;

    if (!(type == OR || type == AND)) {
      throw new Error ("internal error - illegal type given to " +
                       "FeaturePredicateConjunction constructor");
    }

    if (predicates.size () == 0) {
      throw new Error ("internal error - no predicates given to " +
                       "FeaturePredicateConjunction constructor");
    }
  }

  /**
   *  The && predicate type.
   **/
  final public static int AND = 1;

  /**
   *  The || predicate type.
   **/
  final public static int OR = 0;

  /**
   *  Test the given Feature against this predicate.
   *  @param feature The Feature to test the predicate against.
   *  @return predicate1 && predicate2 if AND was passed as a type to the
   *    constructor.  Returns predicate1 || predicate2 if OR was passed to the
   *    constructor.
   **/
  public boolean testPredicate (final Feature feature) {
    for (int i = 0 ; i < predicates.size () ; ++i) {
      if (type == AND) {
        if (!predicates.elementAt (i).testPredicate (feature)) {
          return false;
        }
      } else {
        if (predicates.elementAt (i).testPredicate (feature)) {
          return true;
        }
      }
    }

    if (type == AND) {
      return true;
    } else {
      return false;
    }
  }

  /**
   *  The FeaturePredicates passed to the constructor.
   **/
  private FeaturePredicateVector predicates = new FeaturePredicateVector ();

  /**
   *  The type of conjunction passed to the constructor.
   **/
  private int type;
}
