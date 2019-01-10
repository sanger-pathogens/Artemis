package org.gmod.schema.cv;


import org.gmod.schema.analysis.AnalysisProp;
import org.gmod.schema.general.DbXRef;
import org.gmod.schema.organism.OrganismProp;
import org.gmod.schema.phylogeny.Phylonode;
import org.gmod.schema.phylogeny.PhylonodeProp;
import org.gmod.schema.phylogeny.PhylonodeRelationship;
import org.gmod.schema.phylogeny.Phylotree;
import org.gmod.schema.pub.Pub;
import org.gmod.schema.pub.PubProp;
import org.gmod.schema.pub.PubRelationship;
import org.gmod.schema.sequence.Feature;
import org.gmod.schema.sequence.FeatureCvTerm;
import org.gmod.schema.sequence.FeatureCvTermProp;
import org.gmod.schema.sequence.FeatureProp;
import org.gmod.schema.sequence.FeatureRelationship;
import org.gmod.schema.sequence.FeatureRelationshipProp;
import org.gmod.schema.sequence.Synonym;

import java.io.Serializable;
import java.util.Collection;













public class CvTerm implements Serializable {

    private Collection<Phylotree> phylotrees;
    private Collection<PhylonodeProp> phylonodeProps;
    private Collection<PhylonodeRelationship> phylonodeRelationships;
    private Collection<Phylonode> phylonodes;
   
    


    public Collection<Phylotree> getPhylotrees() {
        return this.phylotrees;
    }
    
    public void setPhylotrees(Collection<Phylotree> phylotrees) {
        this.phylotrees = phylotrees;
    }


    public Collection<PhylonodeProp> getPhylonodeProps() {
        return this.phylonodeProps;
    }
    
    public void setPhylonodeProps(Collection<PhylonodeProp> phylonodeProps) {
        this.phylonodeProps = phylonodeProps;
    }




    public Collection<PhylonodeRelationship> getPhylonodeRelationships() {
        return this.phylonodeRelationships;
    }
    
    public void setPhylonodeRelationships(Collection<PhylonodeRelationship> phylonodeRelationships) {
        this.phylonodeRelationships = phylonodeRelationships;
    }


    public Collection<Phylonode> getPhylonodes() {
        return this.phylonodes;
    }
    
    public void setPhylonodes(Collection<Phylonode> phylonodes) {
        this.phylonodes = phylonodes;
    }
    
    
    
    
    
    
    
    
    
    
    
    // Fields    


     private int cvTermId;
    


     private DbXRef dbXRef;
     


     private Cv cv;
     

     private String name;
     

     private String definition;
     

     private int isObsolete;
     

     private int isRelationshipType;
     

     private Collection<AnalysisProp> analysisProps;
     

     private Collection<CvTermProp> cvTermPropsForTypeId;
     

     private Collection<CvTermProp> cvTermPropsForCvTermId;
     

     private Collection<DbXRefProp> dbXRefProps;
     

     private Collection<Synonym> synonyms;
     

     private Collection<CvTermDbXRef> cvTermDbXRefs;
     

     private Collection<CvTermPath> cvTermPathsForTypeId;
     

     private Collection<FeatureCvTermProp> featureCvTermProps;
     

     private Collection<FeatureCvTerm> featureCvTerms;
     

     private Collection<CvTermRelationship> cvTermRelationshipsForTypeId;
     

     private Collection<CvTermRelationship> cvTermRelationshipsForObjectId;
     

     private Collection<PubProp> pubProps;
     

     private Collection<OrganismProp> organismProps;
     

     private Collection<CvTermRelationship> cvTermRelationshipsForSubjectId;
     

     private Collection<CvTermSynonym> cvTermSynonymsForCvTermId;
     

     private Collection<FeatureProp> featureProps;
     

     private Collection<CvTermPath> cvTermPathsForSubjectId;
     

     private Collection<CvTermPath> cvTermPathsForObjectId;
     

     private Collection<CvTermSynonym> cvTermSynonymsForTypeId;
     

     private Collection<Pub> pubs;
     

     private Collection<FeatureRelationshipProp> featureRelationshipProps;
     

     private Collection<Feature> features;
     

     private Collection<PubRelationship> pubRelationships;
     

     private Collection<FeatureRelationship> featureRelationships;

     // Constructors

    /** default constructor */
    public CvTerm() {
    	// Deliberately empty default constructor
    }
    
    /** useful constructor! */
    public CvTerm(Cv cv, DbXRef dbXRef, String name, String definition) {
       this.dbXRef = dbXRef;
       this.cv = cv;
       this.name = name;
       this.definition = definition;
    }
    
   
    // Property accessors
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#getCvTermId()
     */
    public int getCvTermId() {
        return this.cvTermId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#setCvTermId(int)
     */
    public void setCvTermId(int cvTermId) {
        this.cvTermId = cvTermId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#getDbXRef()
     */
    public DbXRef getDbXRef() {
        return this.dbXRef;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#setDbXRef(org.genedb.db.jpa.DbXRef)
     */
    public void setDbXRef(DbXRef dbXRef) {
        this.dbXRef = dbXRef;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#getCv()
     */
    public Cv getCv() {
        return this.cv;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#setCv(org.gmod.schema.cv.CvI)
     */
    public void setCv(Cv cv) {
        this.cv = cv;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#getName()
     */
    public String getName() {
        return this.name;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#setName(java.lang.String)
     */
    public void setName(String name) {
        this.name = name;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#getDefinition()
     */
    public String getDefinition() {
        return this.definition;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#setDefinition(java.lang.String)
     */
    public void setDefinition(String definition) {
        this.definition = definition;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#getIsObsolete()
     */
    public int getIsObsolete() {
        return this.isObsolete;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#setIsObsolete(int)
     */
    public void setIsObsolete(int isObsolete) {
        this.isObsolete = isObsolete;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#getIsRelationshipType()
     */
    public int getIsRelationshipType() {
        return this.isRelationshipType;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#setIsRelationshipType(int)
     */
    public void setIsRelationshipType(int isRelationshipType) {
        this.isRelationshipType = isRelationshipType;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#getAnalsisProps()
     */
    private Collection<AnalysisProp> getAnalysisProps() {
        return this.analysisProps;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#setAnalsisProps(java.util.Set)
     */
    private void setAnalysisProps(Collection<AnalysisProp> analysisProps) {
        this.analysisProps = analysisProps;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#getCvTermPropsForTypeId()
     */
    private Collection<CvTermProp> getCvTermPropsForTypeId() {
        return this.cvTermPropsForTypeId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#setCvTermPropsForTypeId(java.util.Set)
     */
    private void setCvTermPropsForTypeId(Collection<CvTermProp> cvTermPropsForTypeId) {
        this.cvTermPropsForTypeId = cvTermPropsForTypeId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#getCvTermPropsForCvTermId()
     */
    private Collection<CvTermProp> getCvTermPropsForCvTermId() {
        return this.cvTermPropsForCvTermId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#setCvTermPropsForCvTermId(java.util.Set)
     */
    private void setCvTermPropsForCvTermId(Collection<CvTermProp> cvTermPropsForCvTermId) {
        this.cvTermPropsForCvTermId = cvTermPropsForCvTermId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#getDbXRefProps()
     */
    private Collection<DbXRefProp> getDbXRefProps() {
        return this.dbXRefProps;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#setDbXRefProps(java.util.Set)
     */
    private void setDbXRefProps(Collection<DbXRefProp> dbXRefProps) {
        this.dbXRefProps = dbXRefProps;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#getSynonyms()
     */
    private Collection<Synonym> getSynonyms() {
        return this.synonyms;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#setSynonyms(java.util.Set)
     */
    private void setSynonyms(Collection<Synonym> synonyms) {
        this.synonyms = synonyms;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#getCvTermDbXRefs()
     */
    private Collection<CvTermDbXRef> getCvTermDbXRefs() {
        return this.cvTermDbXRefs;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#setCvTermDbXRefs(java.util.Set)
     */
    private void setCvTermDbXRefs(Collection<CvTermDbXRef> cvTermDbXRefs) {
        this.cvTermDbXRefs = cvTermDbXRefs;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#getCvTermPathsForTypeId()
     */
    private Collection<CvTermPath> getCvTermPathsForTypeId() {
        return this.cvTermPathsForTypeId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#setCvTermPathsForTypeId(java.util.Set)
     */
    private void setCvTermPathsForTypeId(Collection<CvTermPath> cvTermPathsForTypeId) {
        this.cvTermPathsForTypeId = cvTermPathsForTypeId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#getFeatureCvTermProps()
     */
    private Collection<FeatureCvTermProp> getFeatureCvTermProps() {
        return this.featureCvTermProps;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#setFeatureCvTermProps(java.util.Set)
     */
    private void setFeatureCvTermProps(Collection<FeatureCvTermProp> featureCvTermProps) {
        this.featureCvTermProps = featureCvTermProps;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#getFeatureCvTerms()
     */
    private Collection<FeatureCvTerm> getFeatureCvTerms() {
        return this.featureCvTerms;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#setFeatureCvTerms(java.util.Set)
     */
    private void setFeatureCvTerms(Collection<FeatureCvTerm> featureCvTerms) {
        this.featureCvTerms = featureCvTerms;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#getCvTermRelationshipsForTypeId()
     */
    private Collection<CvTermRelationship> getCvTermRelationshipsForTypeId() {
        return this.cvTermRelationshipsForTypeId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#setCvTermRelationshipsForTypeId(java.util.Set)
     */
    private void setCvTermRelationshipsForTypeId(Collection<CvTermRelationship> cvTermRelationshipsForTypeId) {
        this.cvTermRelationshipsForTypeId = cvTermRelationshipsForTypeId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#getCvTermRelationshipsForObjectId()
     */
    public Collection<CvTermRelationship> getCvTermRelationshipsForObjectId() {
        return this.cvTermRelationshipsForObjectId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#setCvTermRelationshipsForObjectId(java.util.Set)
     */
    public void setCvTermRelationshipsForObjectId(Collection<CvTermRelationship> cvTermRelationshipsForObjectId) {
        this.cvTermRelationshipsForObjectId = cvTermRelationshipsForObjectId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#getPubProps()
     */
    private Collection<PubProp> getPubProps() {
        return this.pubProps;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#setPubProps(java.util.Set)
     */
    private void setPubProps(Collection<PubProp> pubProps) {
        this.pubProps = pubProps;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#getOrganismProps()
     */
    private Collection<OrganismProp> getOrganismProps() {
        return this.organismProps;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#setOrganismProps(java.util.Set)
     */
    private void setOrganismProps(Collection<OrganismProp> organismProps) {
        this.organismProps = organismProps;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#getCvTermRelationshipsForSubjectId()
     */
    public Collection<CvTermRelationship> getCvTermRelationshipsForSubjectId() {
        return this.cvTermRelationshipsForSubjectId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#setCvTermRelationshipsForSubjectId(java.util.Set)
     */
    public void setCvTermRelationshipsForSubjectId(Collection<CvTermRelationship> cvTermRelationshipsForSubjectId) {
        this.cvTermRelationshipsForSubjectId = cvTermRelationshipsForSubjectId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#getCvTermSynonymsForCvTermId()
     */
    private Collection<CvTermSynonym> getCvTermSynonymsForCvTermId() {
        return this.cvTermSynonymsForCvTermId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#setCvTermSynonymsForCvTermId(java.util.Set)
     */
    private void setCvTermSynonymsForCvTermId(Collection<CvTermSynonym> cvTermSynonymsForCvTermId) {
        this.cvTermSynonymsForCvTermId = cvTermSynonymsForCvTermId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#getFeatureProps()
     */
    private Collection<FeatureProp> getFeatureProps() {
        return this.featureProps;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#setFeatureProps(java.util.Set)
     */
    private void setFeatureProps(Collection<FeatureProp> featureProps) {
        this.featureProps = featureProps;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#getCvTermPathsForSubjectId()
     */
    private Collection<CvTermPath> getCvTermPathsForSubjectId() {
        return this.cvTermPathsForSubjectId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#setCvTermPathsForSubjectId(java.util.Set)
     */
    private void setCvTermPathsForSubjectId(Collection<CvTermPath> cvTermPathsForSubjectId) {
        this.cvTermPathsForSubjectId = cvTermPathsForSubjectId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#getCvTermPathsForObjectId()
     */
    private Collection<CvTermPath> getCvTermPathsForObjectId() {
        return this.cvTermPathsForObjectId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#setCvTermPathsForObjectId(java.util.Set)
     */
    private void setCvTermPathsForObjectId(Collection<CvTermPath> cvTermPathsForObjectId) {
        this.cvTermPathsForObjectId = cvTermPathsForObjectId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#getCvTermSynonymsForTypeId()
     */
    private Collection<CvTermSynonym> getCvTermSynonymsForTypeId() {
        return this.cvTermSynonymsForTypeId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#setCvTermSynonymsForTypeId(java.util.Set)
     */
    private void setCvTermSynonymsForTypeId(Collection<CvTermSynonym> cvTermSynonymsForTypeId) {
        this.cvTermSynonymsForTypeId = cvTermSynonymsForTypeId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#getPubs()
     */
    private Collection<Pub> getPubs() {
        return this.pubs;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#setPubs(java.util.Set)
     */
    private void setPubs(Collection<Pub> pubs) {
        this.pubs = pubs;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#getFeatureRelationshipProps()
     */
    private Collection<FeatureRelationshipProp> getFeatureRelationshipProps() {
        return this.featureRelationshipProps;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#setFeatureRelationshipProps(java.util.Set)
     */
    private void setFeatureRelationshipProps(Collection<FeatureRelationshipProp> featureRelationshipProps) {
        this.featureRelationshipProps = featureRelationshipProps;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#getFeatures()
     */
    private Collection<Feature> getFeatures() {
        return this.features;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#setFeatures(java.util.Set)
     */
    private void setFeatures(Collection<Feature> features) {
        this.features = features;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#getPubRelationships()
     */
    private Collection<PubRelationship> getPubRelationships() {
        return this.pubRelationships;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#setPubRelationships(java.util.Set)
     */
    private void setPubRelationships(Collection<PubRelationship> pubRelationships) {
        this.pubRelationships = pubRelationships;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#getFeatureRelationships()
     */
    private Collection<FeatureRelationship> getFeatureRelationships() {
        return this.featureRelationships;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermI#setFeatureRelationships(java.util.Set)
     */
    private void setFeatureRelationships(Collection<FeatureRelationship> featureRelationships) {
        this.featureRelationships = featureRelationships;
    }


}

