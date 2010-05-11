/* FeatureKeyQualifierPredicate.java
 *
 * created: Tue Mar 30 1999
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
 *
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/FeatureKeyQualifierPredicate.java,v 1.1 2004-06-09 09:44:42 tjc Exp $
 */

package uk.ac.sanger.artemis;

import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.InvalidRelationException;
import uk.ac.sanger.artemis.util.StringVector;

/**
 *  Each object of this class can be used to test Feature objects to see if
 *  they have a given Key and a given Qualifier.
 *
 *  @author Kim Rutherford
 *  @version $Id: FeatureKeyQualifierPredicate.java,v 1.1 2004-06-09 09:44:42 tjc Exp $
 **/
public class FeatureKeyQualifierPredicate
    implements FeaturePredicate {
  /**
   *  Create a new FeatureKeyPredicate object.
   *  @param key testPredicate () will return false if this Key isn't the
   *    same as the Key of the test Feature.
   *  @param qualifier_name testPredicate () will return false if the test
   *    Feature does not contain a Qualifier with this name.  If null then
   *    match any qualifier name.
   *  @param qualifier_value testPredicate () will return false if the test
   *    Feature does not contain a Qualifier named qualifier_name with this
   *    value.
   **/
  public FeatureKeyQualifierPredicate (final Key key,
                                       final String qualifier_name,
                                       final String qualifier_value) {
    this.key                  = key;
    this.qualifier_name       = qualifier_name;
    this.qualifier_value      = qualifier_value;
    this.qualifier_must_exist = false;
  }

  /**
   *  Create a new FeatureKeyPredicate object.
   *  @param key testPredicate () will return false if this Key isn't the
   *    same as the Key of the test Feature.
   *  @param qualifier_name testPredicate () will return false if the test
   *    Feature does not contain a Qualifier with this name.  If null then
   *    match any qualifier.
   *  @param qualifier_value testPredicate () will return false if the test
   *    Feature does not contain a Qualifier named qualifier_name with this
   *    value.
   *  @param sub_string_match if true then qualifier_value need only match a
   *    substring for the predicate to be true.  If false then qualifier_value
   *    must match the full length of the target for the predicate to be true.
   *  @param ignore_case If true then case will be ignored when searching for
   *    qualifier_value.
   **/
  public FeatureKeyQualifierPredicate (final Key key,
                                       final String qualifier_name,
                                       final String qualifier_value,
                                       final boolean sub_string_match,
                                       final boolean ignore_case) {
    this(key, qualifier_name, qualifier_value, sub_string_match, ignore_case, false);
  }
  
  /**
   *  Create a new FeatureKeyPredicate object.
   *  @param key testPredicate () will return false if this Key isn't the
   *    same as the Key of the test Feature.
   *  @param qualifier_name testPredicate () will return false if the test
   *    Feature does not contain a Qualifier with this name.  If null then
   *    match any qualifier.
   *  @param qualifier_value testPredicate () will return false if the test
   *    Feature does not contain a Qualifier named qualifier_name with this
   *    value.
   *  @param sub_string_match if true then qualifier_value need only match a
   *    substring for the predicate to be true.  If false then qualifier_value
   *    must match the full length of the target for the predicate to be true.
   *  @param ignore_case If true then case will be ignored when searching for
   *    qualifier_value.
   **/
  public FeatureKeyQualifierPredicate (final Key key,
                                       final String qualifier_name,
                                       final String qualifier_value,
                                       final boolean sub_string_match,
                                       final boolean ignore_case,
                                       final boolean deleteQualifier) {
    this.key                  = key;
    this.qualifier_name       = qualifier_name;
    this.qualifier_must_exist = false;
    this.sub_string_match     = sub_string_match;
    this.ignore_case          = ignore_case;
    this.deleteQualifier      = deleteQualifier;

    if (ignore_case) {
      this.qualifier_value = qualifier_value.toLowerCase ();
    } else {
      this.qualifier_value = qualifier_value;
    }
  }

  /**
   *  Create a new FeatureKeyPredicate object.
   *  @param key testPredicate () will return false if this Key isn't the
   *    same as the Key of the test Feature.
   *  @param qualifier_name testPredicate () will return false if the test
   *    Feature does not contain a Qualifier with this name.  If null then
   *    match any qualifier.
   **/
  public FeatureKeyQualifierPredicate (final Key key,
                                       final String qualifier_name) {
    this (key, qualifier_name, true);
  }

  /**
   *  Create a new FeatureKeyPredicate object.
   *  @param key testPredicate () will return false if this Key isn't the
   *    same as the Key of the test Feature.
   *  @param qualifier_name testPredicate () will return false if the test
   *    Feature does not contain a Qualifier with this name.
   *  @param qualifier_must_exist If false then the predicate is true only if
   *    there is NO Qualifier with the given value in the Feature.  If true
   *    then the predicate is true only if there IS a Qualifier with the
   *    given value in the Feature.
   **/
  public FeatureKeyQualifierPredicate (final Key key,
                                       final String qualifier_name,
                                       final boolean qualifier_must_exist) {
    this.key                  = key;
    this.qualifier_name       = qualifier_name;
    this.qualifier_value      = null;
    this.qualifier_must_exist = qualifier_must_exist;
  }

  /**
   *  Test the given Feature against this FeatureKeyPredicate.
   *  @param feature The Feature to test the predicate against.
   *  @return Return true if and only if the given Feature has the same key
   *    as the one the was passed to the constructor and contains a
   *    Qualifier with the name and value that was passed to the
   *    constructor.
   **/
  public boolean testPredicate (final Feature feature) {
    if (key != null && !feature.getKey ().equals (key)) {
      return false;
    }

    if (qualifier_value == null) {
      // we are only testing for existence in this case
      Qualifier qualifier = null;
      try {
        qualifier = feature.getQualifierByName (qualifier_name);
      } catch (InvalidRelationException e) {
        // ignore - qualifier is null
      }

      final boolean has_qualifier = qualifier != null;

      if (qualifier_must_exist && has_qualifier ||
          !qualifier_must_exist && !has_qualifier) {
        return true;
      } else {
        return false;
      }
    } else {
      final StringVector qualifier_names_to_search;

      if (qualifier_name == null) {
        qualifier_names_to_search = null;
      } else {
        qualifier_names_to_search = new StringVector (qualifier_name);
      }

      return feature.findOrReplaceText (qualifier_value, ignore_case,
                                   sub_string_match, deleteQualifier,
                                   qualifier_names_to_search, null);
    }
  }

  /**
   *  The Key that was passed to the constructor.
   **/
  private Key key = null;

  /**
   *  The Qualifier name that was passed to the constructor.
   **/
  private String qualifier_name = null;

  /**
   *  The Qualifier value that was passed to the constructor.
   **/
  private String qualifier_value = null;

  /**
   *  If false then the predicate is true only if there is NO Qualifier with
   *  the given value in the Feature.  If true then the predicate is true
   *  only if there IS a Qualifier with the given value in the Feature.
   **/
  private boolean qualifier_must_exist = true;

  /**
   *  If true then qualifier_value need only match a substring for the
   *  predicate to be true.  If false then qualifier_value must match the full
   *  length of the target for the predicate to be true.  (set by the
   *  constructor).
   **/
  private boolean sub_string_match = false;

  /**
   *  If true then case will be ignored when searching for qualifier_value.
   **/
  private boolean ignore_case = false;
  
  /** If true delete the qualifier value */
  private boolean deleteQualifier = false;
}
