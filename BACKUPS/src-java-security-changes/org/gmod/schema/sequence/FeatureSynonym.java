package org.gmod.schema.sequence;

import org.gmod.schema.pub.Pub;

import java.io.Serializable;

public class FeatureSynonym implements Serializable {

    // Fields    
    private int featureSynonymId;
    private Synonym synonym;
    private Feature feature;
    private Pub pub;
    private boolean current;
    private boolean internal;

     // Constructors

    /** default constructor */
    public FeatureSynonym() {
    	// Deliberately empty default constructor
    }

    /** full constructor */
    public FeatureSynonym(Synonym synonym, Feature feature, Pub pub, boolean current, boolean internal) {
       this.synonym = synonym;
       this.feature = feature;
       this.pub = pub;
       this.current = current;
       this.internal = internal;
    }
    
   
    // Property accessors

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureSynonymI#getFeatureSynonymId()
     */
    public int getFeatureSynonymId() {
        return this.featureSynonymId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureSynonymI#setFeatureSynonymId(int)
     */
    public void setFeatureSynonymId(int featureSynonymId) {
        this.featureSynonymId = featureSynonymId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureSynonymI#getSynonym()
     */
    public Synonym getSynonym() {
        return this.synonym;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureSynonymI#setSynonym(org.gmod.schema.sequence.SynonymI)
     */
    public void setSynonym(Synonym synonym) {
        this.synonym = synonym;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureSynonymI#getFeature()
     */
    public Feature getFeature() {
        return this.feature;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureSynonymI#setFeature(org.genedb.db.jpa.Feature)
     */
    public void setFeature(Feature feature) {
        this.feature = feature;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureSynonymI#getPub()
     */
    public Pub getPub() {
        return this.pub;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureSynonymI#setPub(org.gmod.schema.pub.PubI)
     */
    public void setPub(Pub pub) {
        this.pub = pub;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureSynonymI#isCurrent()
     */
    public boolean isCurrent() {
        return this.current;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureSynonymI#setCurrent(boolean)
     */
    public void setCurrent(boolean current) {
        this.current = current;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureSynonymI#isInternal()
     */
    public boolean isInternal() {
        return this.internal;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureSynonymI#setInternal(boolean)
     */
    public void setInternal(boolean internal) {
        this.internal = internal;
    }




}


