package org.gmod.schema.sequence;

import org.gmod.schema.pub.Pub;

import java.io.Serializable;

public class FeatureRelationshipPub implements Serializable {

    // Fields
    private int featureRelationshipPubId;
    private Pub pub;
    private FeatureRelationship featureRelationship;
   
    // Property accessors

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureRelationshipPubI#getFeatureRelationshipPubId()
     */
    private int getFeatureRelationshipPubId() {
        return this.featureRelationshipPubId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureRelationshipPubI#setFeatureRelationshipPubId(int)
     */
    private void setFeatureRelationshipPubId(int featureRelationshipPubId) {
        this.featureRelationshipPubId = featureRelationshipPubId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureRelationshipPubI#getPub()
     */
    private Pub getPub() {
        return this.pub;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureRelationshipPubI#setPub(org.gmod.schema.pub.PubI)
     */
    private void setPub(Pub pub) {
        this.pub = pub;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureRelationshipPubI#getFeatureRelationship()
     */
    private FeatureRelationship getFeatureRelationship() {
        return this.featureRelationship;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureRelationshipPubI#setFeatureRelationship(org.genedb.db.jpa.FeatureRelationship)
     */
    private void setFeatureRelationship(FeatureRelationship featureRelationship) {
        this.featureRelationship = featureRelationship;
    }




}


