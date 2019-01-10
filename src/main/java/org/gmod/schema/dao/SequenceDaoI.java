/*
 * Copyright (c) 2006 Genome Research Limited.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Library General Public License as published
 * by  the Free Software Foundation; either version 2 of the License or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public License
 * along with this program; see the file COPYING.LIB.  If not, write to
 * the Free Software Foundation Inc., 59 Temple Place - Suite 330,
 * Boston, MA  02111-1307 USA
 */

package org.gmod.schema.dao;


import org.gmod.schema.cv.CvTerm;
import org.gmod.schema.general.DbXRef;
import org.gmod.schema.organism.Organism;
import org.gmod.schema.sequence.Feature;
import org.gmod.schema.sequence.FeatureCvTerm;
import org.gmod.schema.sequence.FeatureCvTermDbXRef;
import org.gmod.schema.sequence.FeatureCvTermPub;
import org.gmod.schema.sequence.FeatureDbXRef;
import org.gmod.schema.sequence.FeatureSynonym;
import org.gmod.schema.sequence.Synonym;
import org.gmod.schema.utils.CountedName;

import java.util.List;

public interface SequenceDaoI extends BaseDaoI {

    /**
     * Return the feature corresponding to this feature_id 
     * 
     * @param id the systematic id
     * @return the Feature, or null
     */
    public Feature getFeatureById(int id);

    /**
     * Return a list of features with this uniquename
     *  
     * @param name the uniquename
     * @return the Feature, or null
     */
    public List<Feature> getFeaturesByUniqueName(String name);
    
    /**
     * 
     * @param name the uniquename
     * @param featureType the type of feature to return eg "gene". <b>NB</> String, not a type argument
     * @return
     */
    public Feature getFeatureByUniqueName(String name, String featureType);
    
    /**
     * Return a list of features with any current (ie non-obsolete) name or synonym
     *  
     * @param name the lookup name
     * @return a (possibly empty) List<Feature> of children with this current name
     */
    public List<Feature> getFeaturesByAnyCurrentName(String name);
    
    /**
     * Return a list of features with this name or synonym (including obsolete names). The 
     * name can contain an SQL wildcard (%) 
     *  
     * @param name the lookup name
     * @param featureType the type of feature to return eg "gene"
     * @return a (possibly empty) List<Feature> of children with this name
     */
    public List<Feature> getFeaturesByAnyName(String name, String featureType);

    
    // TODO Document overlap behaviour
    /**
     * Return a list of features located on a source Feature, within a given range
     *  
     * @param min the minimum (interbase) coordinate
     * @param max the maximum (interbase) coordinate
     * @param strand 
     * @param parent the source feature
     * @param type 
     * @return a List<Feature> which ??? this range
     */
    public List<Feature> getFeaturesByRange(int min, int max, int strand,
            Feature parent, String type);

    /**
     * Return a list of features located on a source Feature 
     *  
     * @param parent the parent feature
     * @return a (possibly empty) List<Feature> of children located on this parent
     */
    public List<Feature> getFeaturesByLocatedOnFeature(Feature parent);

    /**
     * Return the FeatureCvTerm that links a given Feature and CvTerm, with a given value of 'not'
     * 
     * @param feature the Feature to test the link for
     * @param cvTerm the CvTerm to test the link for
     * @param not test for the not flag in the FeatureCvTerm 
     * @return the Feature, or null
     */
    public List<FeatureCvTerm> getFeatureCvTermsByFeatureAndCvTermAndNot(Feature feature,
            CvTerm cvTerm, boolean not);

    /**
     * Return a list of FeatureCvterm's for a Feature, or a list
     * of all FeatureCvTerm's if Feature is null.
     * @param feature the Feature to retrieve associated FeatureCvTerm's
     * @return the FeatureCvTerm's
     */
    public List<FeatureCvTerm> getFeatureCvTermsByFeature(Feature feature);
    
    /**
     * Get a list of all FeatureCvTermDbXRef's for a Feature, or a list
     * of all FeatureCvTermDbXRef's if Feature is null.
     * @param feature the Feature to retrieve associated FeatureCvTermDbXRef's
     * @return the FeatureCvTermDbXRef's
     */
    public List<FeatureCvTermDbXRef> getFeatureCvTermDbXRefByFeature(Feature feature);
    
    /**
     * Get a list of all FeatureCvTermPub's for a Feature, or a list
     * of all FeatureCvTermPub's if Feature is null.
     * @param feature the Feature to retrieve associated FeatureCvTermPub's
     * @return the FeatureCvTermPub's
     */
    public List<FeatureCvTermPub> getFeatureCvTermPubByFeature(Feature feature);
    
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
    public List<FeatureSynonym> getFeatureSynonymsByFeatureAndSynonym(
            Feature feature, Synonym synonym);
    
    public FeatureDbXRef getFeatureDbXRefByFeatureAndDbXRef(final Feature feature, final DbXRef dbXRef);
    
    /**
     * Return all the FeatureDbXRefs for a given feature, <b>specified by name</b>, or all if 
     * <code>null</code> is passed
     * 
     * @param uniqueName the uniquename of a Feature, or null for all FeatureDbXRefs
     * @return a (possibly empty) List<FeatureDbXRefI> 
     */
    public List<FeatureDbXRef> getFeatureDbXRefsByFeatureUniquename(final String uniqueName);
    
    /**
     * Return the list of FeatureSynonyms for a given Feature, <b>specified by name</b>, or all if 
     * <code>null</code> is passed
     * 
     * @param uniqueName the uniquename of a Feature, or null for all
     * @return a (possibly empty) List<FeatureSynonymI> of matching synonyms
     */
    public List<FeatureSynonym> getFeatureSynonymsByFeatureUniquename(final String uniqueName);
    
    /**
     * Return the list of all feature_synonyms as Feature.featureSynonyms 
     * 
     * @return a (possibly empty) List<Features> of matching synonyms
     */
    public List<Feature> getAllFeatureSynonymsAsFeature();
    
    /**
     * Return the list of Features for a given GO number 
     * 
     * 
     * @param go the GO number
     * @return a (possibly empty) List<Feature> of matching genes
     */
    public List<List> getFeatureByGO(final String go);
    
     /**
     * Return a list of features that have this particular cvterm 
     * 
     *  
     * @param cvTermName the CvTerm name
     * @return a (possibly empty) List<Feature> of children
     */
    public List<Feature> getFeaturesByCvTermName(String cvTermName);
    
    /**
     * Return a list of features that have this particular cvterm 
     * 
     *  
     * @param cvTermName the CvTerm name
     * @return a (possibly empty) List<Feature> of children
     */
    public List<Feature> getFeaturesByCvTermNameAndCvName(String cvTermName, String cvName);
    
    /**
     * Return a list of top-level features 
     * 
     *  
     * @return a (possibly empty) List<Feature> of children
     */
    public List<Feature> getTopLevelFeatures();
    
    public List<CountedName> getProducts();
    
    /**
     * Return a list of feature uniquename based on cvterm for auto-completion 
     * 
     * @param name the Feature uniquename
     * @param cvTerm the CvTerm
     * @param limit the number of maximum results to return
     * @return a (possibly empty) List<String> of feature uniquename
     */
    public List<String> getPossibleMatches(String name,CvTerm cvTerm, int limit);
    
    /**
     * Return a list of feature uniquename based on cvterm for auto-completion 
     * 
     * @param name the Feature uniquename
     * @param orgNames the comma seperated organism common names
     * @param featureType the type of Features to return e.g gene
     * @param limit the number of maximum results to return
     * @return a (possibly empty) List<Feature> of Feature
     */
    public List<Feature> getFeaturesByAnyNameAndOrganism(String name,String orgNames,String featureType);
    
    /**
     * Return a list of feature based on organism 
     * 
     * @param organism the Organism
     * @return a (possibly empty) List<String> of feature
     */
    public List<Feature> getFeaturesByOrganism(Organism org);
    
    /**
     * Return the features corresponding to uniquenames in the list
     * 
     * @param names the list of uniquenames
     * @return the list of Features, or null
     */
    public List<Feature> getFeaturesByUniqueNames(List<String> names);

}
