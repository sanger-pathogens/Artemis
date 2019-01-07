package org.gmod.schema.sequence;

import org.gmod.schema.cv.CvTerm;
import org.gmod.schema.utils.propinterface.PropertyI;

import java.io.Serializable;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

public class FeatureRelationshipProp implements Serializable, PropertyI {

    // Fields    
    private int featureRelationshipPropId;
    private CvTerm cvTerm;
    private FeatureRelationship featureRelationship;
    private String value;
    private int rank;
    private Set<FeatureRelationshipPropPub> featureRelationshipPropPubs = new HashSet<FeatureRelationshipPropPub>(0);
   
    // Property accessors

    private int getFeatureRelationshipPropId() {
        return this.featureRelationshipPropId;
    }
    
    private void setFeatureRelationshipPropId(int featureRelationshipPropId) {
        this.featureRelationshipPropId = featureRelationshipPropId;
    }

    public CvTerm getCvTerm() {
        return this.cvTerm;
    }
    
    private void setCvTerm(CvTerm cvTerm) {
        this.cvTerm = cvTerm;
    }

    private FeatureRelationship getFeatureRelationship() {
        return this.featureRelationship;
    }
    
    private void setFeatureRelationship(FeatureRelationship featureRelationship) {
        this.featureRelationship = featureRelationship;
    }
    

    private String getValue() {
        return this.value;
    }
    
    private void setValue(String value) {
        this.value = value;
    }
    

    private int getRank() {
        return this.rank;
    }
    
    private void setRank(int rank) {
        this.rank = rank;
    }

    private Collection<FeatureRelationshipPropPub> getFeatureRelationshipPropPubs() {
        return this.featureRelationshipPropPubs;
    }
    
    private void setFeatureRelationshipPropPubs(Set<FeatureRelationshipPropPub> featureRelationshipPropPubs) {
        this.featureRelationshipPropPubs = featureRelationshipPropPubs;
    }




}


