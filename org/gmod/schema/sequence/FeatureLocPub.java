package org.gmod.schema.sequence;

import org.gmod.schema.pub.Pub;

import java.io.Serializable;

public class FeatureLocPub implements Serializable {

    // Fields    
    private int featureLocPubId;
    private FeatureLoc featureLoc;
    private Pub pub;

    // Property accessors
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureLocPubI#getFeatureLocPubId()
     */
    private int getFeatureLocPubId() {
        return this.featureLocPubId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureLocPubI#setFeatureLocPubId(int)
     */
    private void setFeatureLocPubId(int featureLocPubId) {
        this.featureLocPubId = featureLocPubId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureLocPubI#getFeatureloc()
     */
    private FeatureLoc getFeatureloc() {
        return this.featureLoc;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureLocPubI#setFeatureloc(org.gmod.schema.sequence.FeatureLocI)
     */
    private void setFeatureloc(FeatureLoc featureloc) {
        this.featureLoc = featureloc;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureLocPubI#getPub()
     */
    private Pub getPub() {
        return this.pub;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureLocPubI#setPub(org.gmod.schema.pub.PubI)
     */
    private void setPub(Pub pub) {
        this.pub = pub;
    }




}


