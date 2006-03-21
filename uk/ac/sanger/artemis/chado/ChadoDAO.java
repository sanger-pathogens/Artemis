/* 
 *
 * created: 2006
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2006  Genome Research Limited
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
 */

package uk.ac.sanger.artemis.chado;

import java.sql.*;
import java.io.*;
import java.util.List;

public interface ChadoDAO
{

  /**
   *
   * Get the residues of a feature.
   * @param feature_id  id of feature to query
   * @param schema      schema/organism name or null
   *
   */
  public ChadoFeature getSequence(final int feature_id,
                                  final String schema)
                        throws SQLException;

  /**
   *
   * Get feature name given the feature_id and schema.
   * @param feature_id  id of feature to query
   * @param schema      schema/organism name or null
   *
   */
  public String getFeatureName(final int feature_id,
                               final String schema)
                       throws SQLException;

  /**
   *
   * Get child feature properties for a given parent
   * feature to be able to construct a GFF like feature.
   *
   * @param parentFeatureID  id of parent feature to query
   * @param schema           schema/organism name or null
   *
   */
  public List getGff(final int parentFeatureID,
                     final String schema)
                     throws SQLException;

  /**
   *
   * Given a list of distict cvterm_id/type_id's of feature types
   * that have residues in the given schema and the schema name
   * return a list of features in the schema with residues.
   * @param cvterm_ids list of cvterm_id/type_id's
   * @param schema      schema/organism name or null
   *
   */
  public List getResidueFeatures(List cvterm_ids,
                                 final String schema)
                     throws SQLException;

  /**
   *
   * For a schema return the type_id's with residues.
   * @param schema      schema/organism name or null
   * @return list of type_id's
   * 
   */
  public List getResidueType(final String schema)
                     throws SQLException;

  /**
   * 
   * Get available schemas (as a List of ChadoFeature objects).
   *
   */
  public List getSchema()
              throws SQLException;

  /**
   * 
   * Get a lst of the cvterm's.
   * @param name cvterm name
   * @param cv_name ontology name (e.g. gene, sequence)
   *
   */
  public List getCvterm()
              throws SQLException;

//
// WRITE BACK
//
  /**
   *
   * @param schema schema to update.
   *
   */
  public void updateAttributes
                    (final String schema, final ChadoTransaction tsn)
                     throws SQLException;

  /**
   *
   * @param schema schema to update.
   *
   */
  public void insertAttributes
                    (final String schema, final ChadoTransaction tsn)
                     throws SQLException;

  /**
   *
   * @param schema schema to update.
   *
   */
  public void deleteAttributes
                    (final String schema, final ChadoTransaction tsn)
                     throws SQLException;

  /**
   *
   * Insert a feature into the database.
   * @param schema schema to update.
   *
   */
  public void insertFeature
                    (final String schema, final ChadoTransaction tsn,
                     final String srcfeature_id)
                     throws SQLException;

}
