package org.gmod.schema.sequence;



import org.gmod.schema.pub.Pub;

import java.io.Serializable;

public class FeaturePropPub implements Serializable {

    // Fields    
    private int featurePropPubId;
    private FeatureProp featureProp;
    private Pub pub;
 
    // Property accessors

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeaturePropPubI#getFeaturePropPubId()
     */
    private int getFeaturePropPubId() {
        return this.featurePropPubId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeaturePropPubI#setFeaturePropPubId(int)
     */
    private void setFeaturePropPubId(int featurePropPubId) {
        this.featurePropPubId = featurePropPubId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeaturePropPubI#getFeatureprop()
     */
    private FeatureProp getFeatureProp() {
        return this.featureProp;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeaturePropPubI#setFeatureprop(org.genedb.db.jpa.FeatureProp)
     */
    private void setFeatureProp(FeatureProp featureProp) {
        this.featureProp = featureProp;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeaturePropPubI#getPub()
     */
    private Pub getPub() {
        return this.pub;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeaturePropPubI#setPub(org.gmod.schema.pub.PubI)
     */
    private void setPub(Pub pub) {
        this.pub = pub;
    }




}


