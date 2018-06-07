package org.gmod.schema.sequence;

import org.gmod.schema.pub.Pub;

import java.io.Serializable;

public class FeaturePub implements Serializable {

    // Fields    
    private int featurePubId;
    private Feature feature;
    private Pub pub;

     // Constructors

    /** default constructor */
    public FeaturePub() {
    	// Deliberately empty default constructor
    }

    /** full constructor */
    public FeaturePub(Feature feature, Pub pub) {
       this.feature = feature;
       this.pub = pub;
    }
    
   
    // Property accessors

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeaturePubI#getFeaturePubId()
     */
    private int getFeaturePubId() {
        return this.featurePubId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeaturePubI#setFeaturePubId(int)
     */
    private void setFeaturePubId(int featurePubId) {
        this.featurePubId = featurePubId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeaturePubI#getFeature()
     */
    public Feature getFeature() {
        return this.feature;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeaturePubI#setFeature(org.genedb.db.jpa.Feature)
     */
    public void setFeature(Feature feature) {
        this.feature = feature;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeaturePubI#getPub()
     */
    public Pub getPub() {
        return this.pub;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeaturePubI#setPub(org.gmod.schema.pub.PubI)
     */
    public void setPub(Pub pub) {
        this.pub = pub;
    }




}


