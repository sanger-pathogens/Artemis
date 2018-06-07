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

import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import org.gmod.schema.analysis.AnalysisFeature;
import org.gmod.schema.cv.Cv;
import org.gmod.schema.cv.CvTerm;
import org.gmod.schema.dao.*;
import org.gmod.schema.general.Db;
import org.gmod.schema.general.DbXRef;
import org.gmod.schema.organism.Organism;
import org.gmod.schema.pub.Pub;
import org.gmod.schema.pub.PubDbXRef;
import org.gmod.schema.sequence.Feature;
import org.gmod.schema.sequence.FeatureCvTerm;
import org.gmod.schema.sequence.FeatureCvTermDbXRef;
import org.gmod.schema.sequence.FeatureCvTermProp;
import org.gmod.schema.sequence.FeatureCvTermPub;
import org.gmod.schema.sequence.FeatureDbXRef;
import org.gmod.schema.sequence.FeatureLoc;

import uk.ac.sanger.artemis.util.DatabaseDocument;

public abstract class GmodDAO 
       implements SequenceDaoI, SchemaDaoI, OrganismDaoI, CvDaoI, PubDaoI, GeneralDaoI
{
  private static org.apache.log4j.Logger logger4j = 
    org.apache.log4j.Logger.getLogger(GmodDAO.class);

  public abstract Graph getGraph(final Integer graphId);
  public abstract List getGraphs(final Integer featureId);
  public abstract List getTableColumns(final String tableName);
  
  public abstract List getOrganismsContainingSrcFeatures();
  public abstract List getSimilarityMatchesByFeatureIds(final List featureIds);
  public abstract List getSimilarityMatches(final Integer srcFeatureId);
  public abstract List getClustersByFeatureIds(final List featureIds);
  public abstract List getParentFeaturesByChildFeatureIds(final List subjectIds);
  public abstract Feature getLazyFeatureNoResiduesById(final Integer featureId);
  
  
  public abstract List getFeatureLocsByFeatureId(int featureId);
  public abstract List getFeatureLocsBySrcFeatureId(int srcFeatureId);
  
  /**
   * Return a <code>List</code> of featureLoc's corresponding for a
   * <code>List</code> of feature_id's.
   * @param featureIds the list of featureIds to search
   * @return a (possibly empty) List<Feature>
   */
  public abstract List getFeatureLocsByListOfIds(final Collection featureIds);
  
  /**
   * Return all the Feature.featureDbXRefs for a <code>List</code> of feature_id's. 
   * These are grouped by their feature_id and returned in a <code>List</code>
   * of Feature's.
   * @param featureIds the list of featureIds to search
   * @return a (possibly empty) List<Feature> 
   */
  public abstract List getFeatureDbXRefsByFeatureId(final List featureIds);
  
  /**
   * Return a <code>List</code> of feature's corresponding for a
   * <code>List</code> of feature_id's.
   * @param featureIds the list of featureIds to search
   * @return a (possibly empty) List<Feature>
   */
  public abstract List getFeaturesByListOfIds(final List featureIds);
  
  /**
   * Return a list of chado features with residues
   * with residues.
   * @return    the <code>List</code> of <code>Feature</code> objects
   */
  public abstract List getResidueFeatures(final Integer organismId);
  
  /**
   * Return a list of chado features with residues
   * with residues.
   * @param commonName - organism common name
   * @return    the <code>List</code> of <code>Feature</code> objects
   */
  public abstract List getResidueFeaturesByOrganismCommonName(final String commonName);
  
  /**
   * Get the residues (sub-sequence) for a feature given it's 
   * uniquename
   * @param uniqueName
   * @return
   */
  public abstract List getResiduesByUniqueName(final String uniqueName);
  
  /**
   * Return all the Feature.featureProps for a <code>List</code> of feature_id's. 
   * These are grouped by their feature_id and returned in a <code>List</code>
   * of Feature's.
   * @param featureIds
   * @return a (possibly empty) List<Feature> 
   */
  public abstract List getFeaturePropByFeatureIds(final List featureIds);
  
  /**
   * Return a <code>List</code> of FeatureCvTerm's for all Feature's on
   * the given srcFeatureId
   * @param srcFeatureId  srcfeature_id
   * @return
   */
  public abstract List getFeatureCvTermsBySrcFeature(Feature srcFeature);
  
  /**
   * Return the FeatureCvTremDbXRef's for all Feature's given their srcfeature_id
   * @param srcFeatureId
   * @return
   */
  public abstract List getFeatureCvTermDbXRefBySrcFeature(Feature srcFeature);
  
  /**
   * Return the FeatureCvTermPub's for all Feature's given their srcfeature_id
   * @param srcfeature_id
   * @return
   */
  public abstract List getFeatureCvTermPubBySrcFeature(Feature srcFeature);
  
  /**
   * Return the FeaturePub's for all Feature's given their srcfeature_id
   * @param srcFeatureId
   * @return
   */
  public abstract List getFeaturePubsBySrcFeature(Feature srcFeature);
  
  /**
   * Return the FeaturePub's for a Feature
   * @param feature
   * @return
   */
  public abstract List getFeaturePubsByFeature(Feature feature);
  
  /**
   * Return the FeatureSynonym's for all Feature's given their srcfeature_id
   * @param srcFeatureId
   * @return
   */
  public abstract List getFeatureSynonymsBySrcFeature(Feature srcFeature);
  
  /**
   * Return the FeatureSynonym's for all Feature's given their feature_id's
   * @param 
   * @return
   */
  public abstract List getFeatureSynonymsByFeatureIds(List featureIds);
  
  /**
   * Return a list of features that have this particular cvterm 
   * @param cvTermName the CvTerm name
   * @return a (possibly empty) List<Feature> of children
   */
  public List getFeaturesByCvTermNameAndCvName(String cvTermName, String cvName)
  {
    return null;
  }
 
  
  public List getFeaturesByOrganism(Organism organism)
  {
    return null;
  }
  
  public List getFeaturesByUniqueNames(List name)
  {
    return null;  
  }
  
  //////
  ////// SchemaDaoI
  //////
  //////
  
  
  /**
   * Return the FeatureDbXRef's for all Feature's given their srcfeature_id
   * @param srcFeatureId
   * @return
   */
  public abstract List getFeatureDbXRefsBySrcFeature(Feature srcFeature);
  
  public List getFeaturesByAnyNameAndOrganism(String arg0, String arg1, String arg2)
  {
    // TODO Auto-generated method stub
    return null;
  }

  public List getPossibleMatches(String arg0, CvTerm arg1, int arg2)
  {
    // TODO Auto-generated method stub
    return null;
  }

  public AnalysisFeature getAnalysisFeatureFromFeature(Feature arg0)
  {
    // TODO Auto-generated method stub
    return null;
  }
  
  //////
  ////// GeneralDaoI
  //////
  //////
  /**
   * Retrieve a database by name
   * 
   * @param name the name to lookup
   * @return the corresponding db, or null
   */
  public Db getDbByName(String name)
  {
    return null;
  }

  /**
   * Retrieve the db xref corresponding to a given DB and accession number
   * 
   * @param db the db the dbxref refers to
   * @param accession the accession "number" the dbxref refers to
   * @return the dbxref, or null
   */
  public DbXRef getDbXRefByDbAndAcc(Db db, String accession)
  {
    return null;
  }
  
  public abstract List getDbs();
  
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
  public abstract Organism getOrganismByCommonName(String commonName);

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
  
  public abstract List getAllCvs();
  
  public List getAllTermsInCvWithCount(Cv arg0)
  {
    // TODO Auto-generated method stub
    return null;
  }


  public CvTerm getCvTermByDbXRef(DbXRef arg0)
  {
    // TODO Auto-generated method stub
    return null;
  }
  
  public List getPossibleMatches(String arg0, Cv arg1, int arg2)
  {
    // TODO Auto-generated method stub
    return null;
  }
  
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
  
  
  
  //
  //
  // PubDaoI
  //
  
  public Pub getPubByDbXRef(DbXRef arg0)
  {
    // TODO Auto-generated method stub
    return null;
  }

  public Pub getPubById(int arg0)
  {
    // TODO Auto-generated method stub
    return null;
  }

  public Pub getPubByUniqueName(String arg0)
  {
    // TODO Auto-generated method stub
    return null;
  }

  public List getPubPropByPubAndCvTerm(Pub arg0, CvTerm arg1)
  {
    // TODO Auto-generated method stub
    return null;
  }
  
  
  
  //
  //
  //  Common functions
  //
  
  
  /**
   * Find a dbxref in the database and retrieve the associated 
   * dbxref_id and db_id
   * @param dbXRef
   * @return
   */
  protected DbXRef loadDbXRef(DbXRef dbXRef)
  {
    Integer db_id = getDbId(dbXRef.getDb());

    if(db_id == null)
      throw new RuntimeException("No database called " + 
          dbXRef.getDb().getName() +
          " found -check the spelling!");
    
    dbXRef.getDb().setDbId(db_id.intValue());
     
    Integer dbxref_id = getDbXRefId(dbXRef);
    if(dbxref_id == null)
    {
      dbXRef.setVersion("1");
      // create a new accession entry in dbxref
      insertDbXRef(dbXRef);
      // now get the new dbxref_id
      dbxref_id = getDbXRefId(dbXRef);
    }

    dbXRef.setDbXRefId(dbxref_id.intValue());
    return dbXRef;
  }
  
  /**
   * Find a Pub if it exists. If the Pub does not exist then create
   * one.
   * @param pub
   * @return
   */
  protected Pub loadPub(Pub pub)
  {
    Pub pubResult = getPubByUniqueName(pub);
    
    if(pubResult == null)
    {
      // define the pub.type_id 
      //
      CvTerm cvTerm;
      
      if(pub.getUniqueName().startsWith("PMID:") ||
         pub.getUniqueName().startsWith("PubMed:"))
        cvTerm = DatabaseDocument.getCvTermByCvAndCvTerm("unfetched", "genedb_literature");
      else
        cvTerm = DatabaseDocument.getCvTermByCvAndCvTerm("unknown", "genedb_literature");
      pub.setCvTerm(cvTerm);
      
      insertPub(pub);
      pubResult = getPubByUniqueName(pub);
      
      // add PubDbXRef
      loadPubDbXRef(pubResult);
    }
    
    return pubResult;
  }
  
  /**
   * Create PubDbXRef for new Pub's and to handle links to
   * links to eg, pubmed.
   * @param pub
   */
  private void loadPubDbXRef(final Pub pub)
  {
    try
    {
      int index = pub.getUniqueName().indexOf(':');
      if (index > -1)
      {
        DbXRef dbXRef = new DbXRef();
        dbXRef.setAccession(pub.getUniqueName().substring(index + 1));
        Db db = new Db();
        db.setName(pub.getUniqueName().substring(0, index));
        dbXRef.setDb(db);
        dbXRef = loadDbXRef(dbXRef);

        PubDbXRef pubDbXRef = new PubDbXRef();
        pubDbXRef.setDbXRef(dbXRef);
        pubDbXRef.setPub(pub);
        insertPubDbXRef(pubDbXRef);
      }
    }
    catch (Exception e)
    {
      logger4j.warn("GmodDAO.loadPubDbXRef() :: "+e.getMessage());
    }
  }
  
  /**
   * Insert a feature_cvterm and associated feature_cvtermprop's,
   * feature_cvterm_dbxref's and feature_cvterm_pub.
   * @param feature_cvterm
   */
  protected void insertAllFeatureCvTerm(final FeatureCvTerm feature_cvterm)
  {
    // get the pub_id and create a new Pub if necessary
    if(feature_cvterm.getPub() != null)
      feature_cvterm.setPub( loadPub(feature_cvterm.getPub()) );
    
    insertFeatureCvTerm(feature_cvterm); 
    
    //
    // get the current feature_id sequence value
    int feature_cvterm_id = getCurrval("feature_cvterm_feature_cvterm_id_seq");
    feature_cvterm.setFeatureCvTermId(feature_cvterm_id);
    
    if(feature_cvterm.getFeatureCvTermProps() != null)
    {
      Collection featureCvTermProps = feature_cvterm.getFeatureCvTermProps();
      Iterator it = featureCvTermProps.iterator();
      while(it.hasNext())
      {
        FeatureCvTermProp featureCvTermProp = (FeatureCvTermProp)it.next();
        featureCvTermProp.setFeatureCvTerm(feature_cvterm);

        insertFeatureCvTermProp(featureCvTermProp);
      }
    }
    
    // feature_cvterm_pub's
    if(feature_cvterm.getFeatureCvTermPubs() != null)
    {
      Collection featureCvTermPubs = feature_cvterm.getFeatureCvTermPubs();
      Iterator it = featureCvTermPubs.iterator();
      while(it.hasNext())
      {
        FeatureCvTermPub featureCvTermPub = (FeatureCvTermPub)it.next();
        featureCvTermPub.setFeatureCvTerm(feature_cvterm);

        // get the pub_id and create a new Pub if necessary
        featureCvTermPub.setPub( loadPub(featureCvTermPub.getPub()) );
        
        insertFeatureCvTermPub(featureCvTermPub);
      }
    
    }
    // feature_cvterm_dbxref's
    if(feature_cvterm.getFeatureCvTermDbXRefs() != null)
    {
      Collection featureCvTermDbXRefs = feature_cvterm.getFeatureCvTermDbXRefs();
      Iterator it = featureCvTermDbXRefs.iterator();
      while(it.hasNext())
      {
        FeatureCvTermDbXRef featureCvTermDbXRef = (FeatureCvTermDbXRef)it.next();
        featureCvTermDbXRef.setFeatureCvTerm(feature_cvterm);
        
        // look for dbxref in the database
        DbXRef dbxref = loadDbXRef(featureCvTermDbXRef.getDbXRef());
        insertFeatureCvTermDbXRef(featureCvTermDbXRef);
      }
    }
  }
  
  protected abstract Integer getDbId(Db db);
  protected abstract Integer getDbXRefId(DbXRef dbXRef);
  protected abstract void insertDbXRef(DbXRef dbXRef);
  protected abstract Pub getPubByUniqueName(Pub pub);
  protected abstract void insertPub(Pub pub);
  protected abstract void insertPubDbXRef(PubDbXRef pubDbXRef);
  protected abstract void insertFeatureCvTerm(final FeatureCvTerm feature_cvterm);
  protected abstract int getCurrval(String seq_id);
  protected abstract void insertFeatureCvTermProp(FeatureCvTermProp featureCvTermProp);
  protected abstract void insertFeatureCvTermPub(FeatureCvTermPub featureCvTermPub);
  protected abstract void insertFeatureCvTermDbXRef(FeatureCvTermDbXRef featureCvTermDbXRef);
}