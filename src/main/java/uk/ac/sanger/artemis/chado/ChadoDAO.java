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

import org.gmod.schema.sequence.Synonym;
import org.gmod.schema.sequence.FeatureCvTerm;
import org.gmod.schema.cv.CvTerm;

/**
 * Data access object
 * @see		uk.ac.sanger.artemis.chado.JdbcDAO
 * @see		uk.ac.sanger.artemis.chado.IBatisDAO
 */
public interface ChadoDAO
{

  /**
   * Return the feature corresponding to this feature_id 
   * 
   * @param id the systematic id
   * @return the Feature, or null
   */
  public org.gmod.schema.sequence.Feature getFeatureById(int id);

  /**
   * Return a features with this systematic id
   *  
   * @param name the systematic id
   * @return the Feature, or null
   */
  public org.gmod.schema.sequence.Feature getFeatureByUniqueName(String name);
  
  /**
   * This can be used to get individual features or children.
   * If Feature.featureloc.srcfeature_id is set this is used
   * to return the children of that srcfeature_id.
   * @param feature     the feature used to query
   * @return	the <code>List</code> of child <code>Feature</code> objects
   * @throws SQLException
   */
   public List getFeaturesByLocatedOnFeature(org.gmod.schema.sequence.Feature parent);
  
   /**
    * Return a list of features with any current (ie non-obsolete) name or synonym
    *  
    * @param name the lookup name
    * @return a (possibly empty) List<Feature> of children with this current name
    */
   public List getFeaturesByAnyCurrentName(String name);
   
   /**
    * Return a list of features with this name or synonym (including obsolete names)
    *  
    * @param name the lookup name
    * @return a (possibly empty) List<Feature> of children with this name
    */
   public List getFeaturesByAnyName(String name, String featureType);
   
   /**
    * Return all the FeatureDbXRefs for a given feature, <b>specified by name</b>, or all if 
    * <code>null</code> is passed
    * 
    * @param uniqueName the uniquename of a Feature, or null for all FeatureDbXRefs
    * @return a (possibly empty) List<FeatureDbXRefI> 
    */
   public List getFeatureDbXRefsByFeatureUniquename(final String uniqueName);
   
   /**
    * Return the list of FeatureSynonyms for a given Feature, <b>specified by name</b>, or all if 
    * <code>null</code> is passed
    * 
    * @param uniqueName the uniquename of a Feature, or null for all
    * @return a (possibly empty) List<FeatureSynonymI> of matching synonyms
    */
   public List getFeatureSynonymsByFeatureUniquename(final String uniqueName);

  /**
   * Given a list of distict cvterm_id/type_id's of feature types
   * that have residues (from getResidueType()) in the given schema
   * and the schema name return a list of features in the schema
   * with residues.
   * @param cvterm_ids list of cvterm_id/type_id's
   * @param schema      schema/organism name or null
   * @return	the <code>List</code> of <code>Feature</code> objects
   */
  public List getResidueFeatures(List cvterm_ids,
                                 final String schema);

  /**
   * For a schema return the type_id's with residues.
   * @param schema      schema/organism name or null
   * @return 	the <code>List</code> of type_id's as <code>String</code>
   *            objects
   */
  public List getResidueType(String schema);

  /**
   * Get available schemas (as a <code>List</code> of <code>String</code>       
   * objects).
   * @return    the available schemas
   */
  public List getSchema();

  /**
   * Get the all the <code>CvTerm</code> objects as a<code>List</code>.
   * @return	the full list of cvterm_id and name
   */
  public List getCvTerms();
  
  /**
   * Return a synonym of the given name and type if it exists
   * 
   * @param name the name to lookup
   * @param type the type of the Synonym
   * @return a Synonym, or null  
   */
  public Synonym getSynonymByNameAndCvTerm(String name, CvTerm type);
  
  
  /**
   * Return a list of FeatureSynonyms which link a given Feature and Synonym
   * 
   * @param feature the test Feature
   * @param synonym the test Synonym
   * @return a (possibly empty) List<FeatureSynonym>
   */
  public List getFeatureSynonymsByFeatureAndSynonym(
         org.gmod.schema.sequence.Feature feature, Synonym synonym);
  
  
  /**
   * Return the FeatureCvTerm that links a given Feature and CvTerm, 
   * with a given value of 'not'
   * 
   * @param feature the Feature to test the link for
   * @param cvTerm the CvTerm to test the link for
   * @param not test for the not flag in the FeatureCvTerm 
   * @return the Feature, or null
   */
  public FeatureCvTerm getFeatureCvTermByFeatureAndCvTerm(
          org.gmod.schema.sequence.Feature feature,
          CvTerm cvTerm, boolean not);
  
  public List getOrganisms();

//
// WRITE BACK
//
  
  /**
   * Merge (update) an already persistent object back to the database 
   * (at the end of the current transaction, or depending upon flush mode). 
   * This method is defined in all the DAOs. It's recommended to call it 
   * through an appropriate one eg SequenceDaoI for FeatureI 
   * @param o The object to merge
   */
  public void merge(Object o);
  
  /**
   * Save the object to the database (at the end of the current transaction, 
   * or depending upon flush mode). This method is defined in all the DAOs. 
   * It's recommended to call it through an appropriate one eg SequenceDaoI 
   * for FeatureI 
   * @param o The object to store
   */
  public void persist(Object o);
  
  /**
   * Remove the object from the database (at the end of the current transaction, 
   * or depending upon flush mode). This method is defined in all the DAOs. 
   * It's recommended to call it through an appropriate one eg SequenceDaoI for 
   * FeatureI 
   * @param o The object to delete
   */
  public void delete(Object o);
  
}
