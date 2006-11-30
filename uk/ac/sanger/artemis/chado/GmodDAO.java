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

import java.util.List;

import org.gmod.schema.cv.Cv;
import org.gmod.schema.cv.CvTerm;
import org.gmod.schema.dao.*;
import org.gmod.schema.general.DbXRef;
import org.gmod.schema.organism.Organism;
import org.gmod.schema.sequence.Feature;
import org.gmod.schema.sequence.FeatureDbXRef;

public abstract class GmodDAO implements SequenceDaoI, SchemaDaoI, OrganismDaoI, CvDaoI
{

  public abstract void merge(Object obj);
  public abstract void persist(Object obj);
  public abstract void delete(Object obj);
  
  public abstract List getFeatureCvTermsByFeature(Feature feature);

  public abstract List getPubDbXRef();
  
  public abstract List getFeatureCvTermDbXRefByFeature(Feature feature);
  
  public abstract List getFeatureCvTermPubByFeature(Feature feature);
  
  /**
   * Return the list of all feature_synonyms as Feature.featureSynonyms 
   * 
   * @return a (possibly empty) List<Features> of matching synonyms
   */
  public abstract List getAllFeatureSynonymsAsFeature();
  
  /**
   * Return a list of features that have this particular cvterm 
   * @param cvTermName the CvTerm name
   * @return a (possibly empty) List<Feature> of children
   */
  public List getFeaturesByCvTermNameAndCvName(String cvTermName, String cvName)
  {
    return null;
  }
  
  //////
  ////// OrganismDaoI
  //////
  //////
  /**
   * Get the organism corresponding to this id
   * 
   * @param id the organism id (primary key) to lookup by
   * @return the corresponding organism, or null
   */
  public Organism getOrganismById(int id)
  {
    return null;
  }

  /**
   * Get the organism corresponding to this common name 
   * 
   * @param commonName the short name to look up
   * @return the corresponding organism, or null
   */
  public Organism getOrganismByCommonName(String commonName)
  {
    return null;
  }

  /**
   * Get a list of the common name of all the organisms.
   * 
   * @return a (possibly empty) List<String> of all the organisms' common names
   */
  public List findAllOrganismCommonNames()
  {
    return null;
  }
 
  
  //////
  ////// CvDaoI
  //////
  //////
  
  /**
   * Get a CV by id
   * 
   * @param id the cv id (primary key)
   * @return the corresponding Cv, or null
   */
  public Cv getCvById(int id)
  {
    return null;
  }

  // TODO Should this return a list or just one?
  /**
   * Retrieve a controlled vocabulary by its name
   * 
   * @param name the name to lookup
   * @return the List<Cv> of matches, or null
   */
  public List getCvByName(String name)
  {
    return null;
  }

  /**
   * Retrieve a CvTerm by id
   * 
   * @param id then cvterm id (primary key)
   * @return the corresponding CvTerm, or null
   */
  public CvTerm getCvTermById(int id)
  {
    return null;
  }


  // TODO Should this return a list or just one?
  /**
   * Retrieve a named CvTerm from a given Cv
   * 
   * @param cvTermName the name of the cvterm
   * @param cv the controlled vocabulary this cvterm is part of
   * @return a (possibly empty) list of matching cvterms
   */
  public List getCvTermByNameInCv(String cvTermName, Cv cv)
  {
    return null;
  }


  /**
   * Retrieve a CvTerm from the Gene Ontology
   * 
   * @param value the 
   * @return the corresponding CvTerm, or null
   */
  public CvTerm getGoCvTermByAcc(String value)
  {
    return null;
  }


  /**
   * Retrieve a CvTerm from the Gene Ontology via it's database entry
   * 
   * @param id the database name eg GO:123456
   * @return the corresponding CvTerm, or null
   */
  public CvTerm getGoCvTermByAccViaDb(final String id)
  {
    return null;
  }
  

  public List getFeatureCvTermsByFeatureAndCvTermAndNot(Feature arg0, CvTerm arg1, boolean arg2)
  {
    // TODO Auto-generated method stub
    return null;
  }


  public FeatureDbXRef getFeatureDbXRefByFeatureAndDbXRef(Feature arg0, DbXRef arg1)
  {
    // TODO Auto-generated method stub
    return null;
  }


  public boolean existsNameInOntology(String arg0, Cv arg1)
  {
    // TODO Auto-generated method stub
    return false;
  }
  
}