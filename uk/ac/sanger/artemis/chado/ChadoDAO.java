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
import java.util.List;

/**
 *
 * Data access object
 * @see		uk.ac.sanger.artemis.chado.JdbcDAO
 * @see		uk.ac.sanger.artemis.chado.IBatisDAO
 *
 */
public interface ChadoDAO
{

  /**
   * Return the feature corresponding to this feature_id 
   * 
   * @param id the systematic id
   * @return the Feature, or null
   */
  public Feature getFeatureById(int id)
                      throws SQLException;

  /**
   * Return a features with this systematic id
   *  
   * @param name the systematic id
   * @return the Feature, or null
   */
  public Feature getFeatureByUniqueName(String name)
                      throws SQLException;
  
  /**
   * This can be used to get individual features or children.
   * If Feature.featureloc.srcfeature_id is set this is used
   * to return the children of that srcfeature_id.
   * @param feature     the feature used to query
   * @return	the <code>List</code> of child <code>Feature</code> objects
   * @throws SQLException
   */
   public List getFeaturesByLocatedOnFeature(Feature parent)
                         throws SQLException;
  
   /**
    * Return a list of features with any current (ie non-obsolete) name or synonym
    *  
    * @param name the lookup name
    * @return a (possibly empty) List<Feature> of children with this current name
    */
   public List getFeaturesByAnyCurrentName(String name)
                 throws SQLException;
   
   /**
    * Return a list of features with this name or synonym (including obsolete names)
    *  
    * @param name the lookup name
    * @return a (possibly empty) List<Feature> of children with this name
    */
   public List getFeaturesByAnyName(String name, String featureType);
   
  /**
   * Given a list of distict cvterm_id/type_id's of feature types
   * that have residues (from getResidueType()) in the given schema
   * and the schema name return a list of features in the schema
   * with residues.
   * @param cvterm_ids list of cvterm_id/type_id's
   * @param schema      schema/organism name or null
   * @return	the <code>List</code> of <code>Feature</code> objects
   * @throws SQLException
   */
  public List getResidueFeatures(List cvterm_ids,
                                 final String schema)
                     throws SQLException;

  /**
   * For a schema return the type_id's with residues.
   * @param schema      schema/organism name or null
   * @return 	the <code>List</code> of type_id's as <code>String</code>
   *            objects
   * @throws SQLException
   */
  public List getResidueType(String schema)
                     throws SQLException;

  /**
   * Get available schemas (as a <code>List</code> of <code>String</code>       
   * objects).
   * @return    the available schemas
   * @throws SQLException
   */
  public List getSchema()
              throws SQLException;

  /**
   * Get the full list of cvterm_id and name as a <code>List</code> of 
   * <code>Cvterm</code> objects.
   * @return	the full list of cvterm_id and name
   * @throws SQLException
   */
  public List getCvterm()
              throws SQLException;
  

  /**
   * Get dbxref for a feature.
   * @param uniquename  the unique name for the feature. If set to NULL
   *                    all <code>FeatureDbxref</code> are returned.
   * @return a <code>List</code> of feature_dbxrefs.
   * @throws SQLException
   */
  public List getFeatureDbxrefByUniquename(final String uniquename)
              throws SQLException;
  
  /**
   * Return a list of FeatureSynonyms for a uniquename
   * @param uniquename  the unique name for the feature. If set to NULL
   *                    all <code>FeatureSynonym</code> are returned.
   * @return
   * @throws SQLException
   */
  public List getFeatureSynonymsByUniquename(final String uniquename)
         throws SQLException;
  
  /**
   * Return a synonym of the given name and type if it exists
   * 
   * @param name the name to lookup
   * @param type the type of the Synonym
   * @return a Synonym, or null  
   */
  public Synonym getSynonymByNameAndCvTerm(String name, Cvterm type)
         throws SQLException;
  
  
  /**
   * Return a list of FeatureSynonyms which link a given Feature and Synonym
   * 
   * @param feature the test Feature
   * @param synonym the test Synonym
   * @return a (possibly empty) List<FeatureSynonym>
   */
  public List getFeatureSynonymsByFeatureAndSynonym(
         Feature feature, Synonym synonym)
         throws SQLException;
  
  
//
// WRITE BACK
//
  /**
   * Update attributes defined by the <code>ChadoTransaction</code>.
   * @param tsn		the <code>ChadoTransaction</code>
   * @return    number of rows changed
   * @throws SQLException
   */
  public int updateAttributes
                    (final ChadoTransaction tsn)
                     throws SQLException;

  /**
   * Insert attributes defined by the <code>ChadoTransaction</code>.
   * @param tsn         the <code>ChadoTransaction</code>
   * @throws SQLException
   */
  public void insertAttributes
                    (final ChadoTransaction tsn)
                     throws SQLException;

  /**
   * Delete attributes defined by the <code>ChadoTransaction</code>.
   * @param tsn         the <code>ChadoTransaction</code>
   * @throws SQLException
   */
  public void deleteAttributes
                    (final ChadoTransaction tsn)
                     throws SQLException;

  /**
   * Insert a feature into the database defined by the <code>ChadoTransaction</code>.
   * @param tsn         	the <code>ChadoTransaction</code>
   * @parma srcfeature_id	the parent feature identifier
   * @throws SQLException
   */
  public void insertFeature
                    (final ChadoTransaction tsn,
                     final String srcfeature_id)
                     throws SQLException;

  /**
   * Delete a feature from the database defined by the <code>ChadoTransaction</code>.
   * @param tsn         the <code>ChadoTransaction</code>
   * @return    number of rows deleted
   * @throws SQLException
   */
  public int deleteFeature
                    (final ChadoTransaction tsn)
                     throws SQLException;
  
  /**
   * Insert a dbxref for a feature.
   * @param tsn           the <code>ChadoTransaction</code>
   * @return    number of rows changed
   * @throws SQLException
   */
  public int insertFeatureDbxref(final ChadoTransaction tsn)
                     throws SQLException;
  
  /**
   * Delete a dbxref for a feature.
   * @param tsn           the <code>ChadoTransaction</code>
   * @return    number of rows changed
   * @throws SQLException
   */
  public int deleteFeatureDbxref(final ChadoTransaction tsn)
                     throws SQLException;
  
  /**
   * Insert a synonym for a feature.
   * @param tsn           the <code>ChadoTransaction</code>
   * @return    number of rows changed
   * @throws SQLException
   */
  public int insertFeatureAlias(final ChadoTransaction tsn)
                     throws SQLException;

  /**
   * Delete a synonym for a feature.
   * @param tsn           the <code>ChadoTransaction</code>
   * @return    number of rows changed
   * @throws SQLException
   */
  public int deleteFeatureAlias(final ChadoTransaction tsn)
                     throws SQLException;

  
  /**
   * Update feature_relationship for a feature.
   * @param tsn           the <code>ChadoTransaction</code>
   * @return    number of rows changed
   * @throws SQLException
   */
  public void updateFeatureRelationshipsForSubjectId(
      final ChadoTransaction tsn)
                     throws SQLException;
  /**
   * 
   * Write the time a feature was last modified
   * @param uniquename	the unique name of the feature
   * @param timestamp   the time stamp to use, 
   *                    if NULL use CURRENT_TIMESTAMP
   * @return  number of rows changed
   * @throws SQLException
   */
  public int writeTimeLastModified
                    (final String uniquename,
                     final Timestamp timestamp)
                     throws SQLException;

  /**
   *
   * Write the time a feature was last accessed
   * @param uniquename  the unique name of the feature
   * @param timestamp   the time stamp to use, 
   *                    if NULL use CURRENT_TIMESTAMP
   * @return  number of rows changed
   * @throws SQLException
   */
  public int writeTimeAccessioned
                    (final String uniquename,
                     final Timestamp timestamp)
                     throws SQLException;
}
