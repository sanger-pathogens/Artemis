package org.gmod.schema.sequence;



import org.gmod.schema.pub.Pub;

import java.io.Serializable;

public class FeatureCvTermPub implements Serializable {

    // Fields    
    private int featureCvTermPubId;
    private Pub pub;
    private FeatureCvTerm featureCvTerm;

     // Constructors

    /** default constructor */
    public FeatureCvTermPub() {
    	// Deliberately empty default constructor
    }

    /** full constructor */
    public FeatureCvTermPub(Pub pub, FeatureCvTerm featureCvTerm) {
       this.pub = pub;
       this.featureCvTerm = featureCvTerm;
    }
    
   
    // Property accessors
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureCvTermPubI#getFeatureCvTermPubId()
     */
    public int getFeatureCvTermPubId() {
        return this.featureCvTermPubId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureCvTermPubI#setFeatureCvTermPubId(int)
     */
    public void setFeatureCvTermPubId(int featureCvTermPubId) {
        this.featureCvTermPubId = featureCvTermPubId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureCvTermPubI#getPub()
     */
    public Pub getPub() {
        return this.pub;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureCvTermPubI#setPub(org.gmod.schema.pub.PubI)
     */
    public void setPub(Pub pub) {
        this.pub = pub;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureCvTermPubI#getFeatureCvterm()
     */
    public FeatureCvTerm getFeatureCvTerm() {
        return this.featureCvTerm;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureCvTermPubI#setFeatureCvterm(org.gmod.schema.sequence.FeatureCvTermI)
     */
    public void setFeatureCvTerm(FeatureCvTerm featureCvTerm) {
        this.featureCvTerm = featureCvTerm;
    }




}


