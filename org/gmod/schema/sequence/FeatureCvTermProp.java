package org.gmod.schema.sequence;


import org.gmod.schema.cv.CvTerm;
import org.gmod.schema.utils.propinterface.PropertyI;
import org.gmod.schema.utils.Rankable;

import java.io.Serializable;

/**
 * This represents a key/value pair attached to a FeatureCvTerm. The key is itself a CvTerm, with 
 * a free-text value
 * 
 * Database constraints: feature_cvtermprop_c1 unique (feature_cvterm_id, type_id, rank)
 * 
 * @author art
 */
/**
 * @author art
 *
 */


public class FeatureCvTermProp implements Serializable, PropertyI, Rankable {
 
    // Fields    

    
     /**
     * Database unique primary key 
     */
     private int featureCvTermPropId;
     

    /**
     * The CvTerm that acts as the key in this map of properties
     */
     private CvTerm cvTerm;
     
     /**
     * The FeatureCvTerm to which this property is attached
     */
     private FeatureCvTerm featureCvTerm;
     
     /**
     * The value of this property
     */
     private String value;
     
     /**
     * The rank is used to distinguish multiple 
     * values for the same key eg /foo="value1", /foo="value2" in an EMBL file could be stored as two FeatureCvTerm with 
     * different ranks. The default is 0;
     */
     private int rank;

     // Constructors

    /** default constructor */
    public FeatureCvTermProp() {
    	// Deliberately empty default constructor
    }

	/** minimal constructor */
    public FeatureCvTermProp(CvTerm cvTerm, FeatureCvTerm featureCvTerm, int rank) {
        this.cvTerm = cvTerm;
        this.featureCvTerm = featureCvTerm;
        this.rank = rank;
    }
    /** full constructor */
    public FeatureCvTermProp(CvTerm cvTerm, FeatureCvTerm featureCvTerm, String value, int rank) {
       this.cvTerm = cvTerm;
       this.featureCvTerm = featureCvTerm;
       this.value = value;
       this.rank = rank;
    }
    
   
    // Property accessors

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureCvTermPropI#getFeatureCvTermpropId()
     */
    private int getFeatureCvTermPropId() {
        return this.featureCvTermPropId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureCvTermPropI#setFeatureCvTermpropId(int)
     */
    private void setFeatureCvTermPropId(int featureCvTermPropId) {
        this.featureCvTermPropId = featureCvTermPropId;
    }

    /**
     * Accessor for featureCvTerm
     * 
     * @see featureCvTerm
     */
    public CvTerm getCvTerm() {
        return this.cvTerm;
    }
    
    /**
     * Accessor for cvTerm
     * 
     * @see cvTerm
     */
    public void setCvTerm(CvTerm cvTerm) {
        this.cvTerm = cvTerm;
    }

    /**
     * Accessor for featureCvTerm
     * 
     * @see featureCvTerm
     */
    public FeatureCvTerm getFeatureCvTerm() {
        return this.featureCvTerm;
    }
    

    /**
     * Accessor for featureCvTerm
     * 
     * @see featureCvTerm
     */
    public void setFeatureCvTerm(FeatureCvTerm featureCvTerm) {
        this.featureCvTerm = featureCvTerm;
    }
    

    /**
     * Accessor for value
     * 
     * @see value
     */
    public String getValue() {
        return this.value;
    }
    
    /**
     * Accessor for value
     * 
     * @see value
     */
    public void setValue(String value) {
        this.value = value;
    }
    

    /**
     * Accessor for rank
     * 
     * @see rank
     */
    public int getRank() {
        return this.rank;
    }
    
    /**
     * Accessor for rank
     * 
     * @see rank
     */
    public void setRank(int rank) {
        this.rank = rank;
    }




}


