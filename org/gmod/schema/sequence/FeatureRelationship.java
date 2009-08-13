package org.gmod.schema.sequence;

import org.gmod.schema.cv.CvTerm;

import java.io.Serializable;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

public class FeatureRelationship implements Serializable {

    // Fields    
    private int featureRelationshipId;
    private Feature featureBySubjectId;
    private Feature featureByObjectId;
    private CvTerm cvTerm;
    private String value;
    private int rank;
    private Set<FeatureRelationshipProp> featureRelationshipProps = new HashSet<FeatureRelationshipProp>(0);
    private Set<FeatureRelationshipPub> featureRelationshipPubs = new HashSet<FeatureRelationshipPub>(0);

     // Constructors

    /** default constructor */
    public FeatureRelationship() {
    	// Deliberately empty default constructor
    }

	/** minimal constructor */
    public FeatureRelationship(Feature featureBySubjectId, Feature featureByObjectId, CvTerm cvTerm, int rank) {
        this.featureBySubjectId = featureBySubjectId;
        this.featureByObjectId = featureByObjectId;
        this.cvTerm = cvTerm;
        this.rank = rank;
    }
    
   
    // Property accessors


    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureRelationshipI#getFeatureRelationshipId()
     */
    public int getFeatureRelationshipId() {
        return this.featureRelationshipId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureRelationshipI#setFeatureRelationshipId(int)
     */
    public void setFeatureRelationshipId(int featureRelationshipId) {
        this.featureRelationshipId = featureRelationshipId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureRelationshipI#getFeatureBySubjectId()
     */
    public Feature getFeatureBySubjectId() {
        return this.featureBySubjectId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureRelationshipI#setFeatureBySubjectId(org.genedb.db.jpa.Feature)
     */
    public void setFeatureBySubjectId(Feature featureBySubjectId) {
        this.featureBySubjectId = featureBySubjectId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureRelationshipI#getFeatureByObjectId()
     */
    public Feature getFeatureByObjectId() {
        return this.featureByObjectId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureRelationshipI#setFeatureByObjectId(org.genedb.db.jpa.Feature)
     */
    public void setFeatureByObjectId(Feature featureByObjectId) {
        this.featureByObjectId = featureByObjectId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureRelationshipI#getCvterm()
     */
    public CvTerm getCvTerm() {
        return this.cvTerm;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureRelationshipI#setCvterm(org.gmod.schema.cv.CvTermI)
     */
    public void setCvTerm(CvTerm cvTerm) {
        this.cvTerm = cvTerm;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureRelationshipI#getValue()
     */
    public String getValue() {
        return this.value;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureRelationshipI#setValue(java.lang.String)
     */
    private void setValue(String value) {
        this.value = value;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureRelationshipI#getRank()
     */
    public int getRank() {
        return this.rank;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureRelationshipI#setRank(int)
     */
    public void setRank(int rank) {
        this.rank = rank;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureRelationshipI#getFeatureRelationshipprops()
     */
    private Collection<FeatureRelationshipProp> getFeatureRelationshipProps() {
        return this.featureRelationshipProps;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureRelationshipI#setFeatureRelationshipprops(java.util.Set)
     */
    private void setFeatureRelationshipProps(Set<FeatureRelationshipProp> featureRelationshipProps) {
        this.featureRelationshipProps = featureRelationshipProps;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureRelationshipI#getFeatureRelationshipPubs()
     */
    private Collection<FeatureRelationshipPub> getFeatureRelationshipPubs() {
        return this.featureRelationshipPubs;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureRelationshipI#setFeatureRelationshipPubs(java.util.Set)
     */
    private void setFeatureRelationshipPubs(Set<FeatureRelationshipPub> featureRelationshipPubs) {
        this.featureRelationshipPubs = featureRelationshipPubs;
    }




}


