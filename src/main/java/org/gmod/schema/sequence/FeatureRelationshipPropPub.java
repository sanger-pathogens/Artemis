package org.gmod.schema.sequence;

import org.gmod.schema.pub.Pub;

import java.io.Serializable;

public class FeatureRelationshipPropPub implements Serializable {

    // Fields    
    private int featureRelationshipPropPubId;
    private FeatureRelationshipProp featureRelationshipProp;
    private Pub pub;  
   
    // Property accessors

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureRelationshipPropPubI#getFeatureRelationshipPropPubId()
     */
    private int getFeatureRelationshipPropPubId() {
        return this.featureRelationshipPropPubId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureRelationshipPropPubI#setFeatureRelationshipPropPubId(int)
     */
    private void setFeatureRelationshipPropPubId(int featureRelationshipPropPubId) {
        this.featureRelationshipPropPubId = featureRelationshipPropPubId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureRelationshipPropPubI#getFeatureRelationshipprop()
     */
    private FeatureRelationshipProp getFeatureRelationshipProp() {
        return this.featureRelationshipProp;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureRelationshipPropPubI#setFeatureRelationshipprop(org.genedb.db.jpa.FeatureRelationshipProp)
     */
    private void setFeatureRelationshipProp(FeatureRelationshipProp featureRelationshipProp) {
        this.featureRelationshipProp = featureRelationshipProp;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureRelationshipPropPubI#getPub()
     */
    private Pub getPub() {
        return this.pub;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureRelationshipPropPubI#setPub(org.gmod.schema.pub.PubI)
     */
    private void setPub(Pub pub) {
        this.pub = pub;
    }




}


