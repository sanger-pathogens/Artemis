package org.gmod.schema.sequence;

import org.gmod.schema.general.DbXRef;

import java.io.Serializable;

public class FeatureDbXRef implements Serializable {

    // Fields    
    private int featureDbXRefId;
    private DbXRef dbXRef;
    private Feature feature;
    private boolean current;

     // Constructors

    /** default constructor */
    public FeatureDbXRef() {
    	// Deliberately empty default constructor
    }

    /** full constructor */
    public FeatureDbXRef(DbXRef dbXRef, Feature feature, boolean current) {
       this.dbXRef = dbXRef;
       this.feature = feature;
       this.current = current;
    }
    
   
    // Property accessors
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureDbXRefI#getFeatureDbXRefId()
     */
    public int getFeatureDbXRefId() {
        return this.featureDbXRefId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureDbXRefI#setFeatureDbXRefId(int)
     */
    public void setFeatureDbXRefId(int featureDbXRefId) {
        this.featureDbXRefId = featureDbXRefId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureDbXRefI#getDbxref()
     */
    public DbXRef getDbXRef() {
        return this.dbXRef;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureDbXRefI#setDbxref(org.gmod.schema.general.DbXRefI)
     */
    public void setDbXRef(DbXRef dbXRef) {
        this.dbXRef = dbXRef;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureDbXRefI#getFeature()
     */
    public Feature getFeature() {
        return this.feature;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureDbXRefI#setFeature(org.genedb.db.jpa.Feature)
     */
    public void setFeature(Feature feature) {
        this.feature = feature;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureDbXRefI#isCurrent()
     */
    public boolean isCurrent() {
        return this.current;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureDbXRefI#setCurrent(boolean)
     */
    public void setCurrent(boolean current) {
        this.current = current;
    }




}


