package org.gmod.schema.sequence;


import org.gmod.schema.cv.CvTerm;
import org.gmod.schema.utils.propinterface.PropertyI;




import java.io.Serializable;
import java.util.HashSet;
import java.util.Set;

public class FeatureProp implements Serializable, PropertyI {

    // Fields    
    private int featurePropId;
    public CvTerm cvTerm;
    private Feature feature;    
    private String value;
    private int rank;  
    private Set<FeaturePropPub> featurePropPubs = new HashSet<FeaturePropPub>(0);

     // Constructors
    /** default constructor */
    public FeatureProp() {
    	// Deliberately empty default constructor
    }
    
    /** useful constructor ! */
    public FeatureProp(Feature feature, CvTerm cvTerm, String value, int rank) {
       this.cvTerm = cvTerm;
       this.feature = feature;
       this.value = value;
       this.rank = rank;
    }

    // Property accessors

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeaturePropI#getCvterm()
     */
    public CvTerm getCvTerm() {
        return this.cvTerm;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeaturePropI#setCvterm(org.gmod.schema.cv.CvTermI)
     */
    public void setCvTerm(CvTerm cvTerm) {
        this.cvTerm = cvTerm;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeaturePropI#getFeature()
     */
    public Feature getFeature() {
        return this.feature;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeaturePropI#setFeature(org.genedb.db.jpa.Feature)
     */
    public void setFeature(Feature feature) {
        this.feature = feature;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeaturePropI#getValue()
     */
    public String getValue() {
        return this.value;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeaturePropI#setValue(java.lang.String)
     */
    public void setValue(String value) {
        this.value = value;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeaturePropI#getRank()
     */
    public int getRank() {
        return this.rank;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeaturePropI#setRank(int)
     */
    public void setRank(int rank) {
        this.rank = rank;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeaturePropI#getFeaturepropPubs()
     */
    private Set<FeaturePropPub> getFeaturePropPubs() {
        return this.featurePropPubs;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeaturePropI#setFeaturepropPubs(java.util.Set)
     */
    public void setFeaturePropPubs(Set<FeaturePropPub> featurePropPubs) {
        this.featurePropPubs = featurePropPubs;
    }

    public int getFeaturePropId() {
        return this.featurePropId;
    }

    public void setFeaturePropId(final int featurePropId) {
        this.featurePropId = featurePropId;
    }

}


