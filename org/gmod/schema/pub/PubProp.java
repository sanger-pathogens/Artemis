package org.gmod.schema.pub;


import org.gmod.schema.cv.CvTerm;
import org.gmod.schema.utils.propinterface.PropertyI;

import java.io.Serializable;












public class PubProp implements Serializable, PropertyI {

    // Fields    

    

     private int pubPropId;
     

         

     private CvTerm cvTerm;
     

         

     private Pub pub;
     

     private String value;
     

     private Integer rank;

     // Constructors

    /** default constructor */
    public PubProp() {
    	// Deliberately empty default constructor
    }

	/** minimal constructor */
    public PubProp(CvTerm cvTerm, Pub pub, String value) {
        this.cvTerm = cvTerm;
        this.pub = pub;
        this.value = value;
    }
    /** full constructor */
    public PubProp(CvTerm cvTerm, Pub pub, String value, Integer rank) {
       this.cvTerm = cvTerm;
       this.pub = pub;
       this.value = value;
       this.rank = rank;
    }
    
   
    // Property accessors

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubPropI#getPubPropId()
     */
    private int getPubPropId() {
        return this.pubPropId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubPropI#setPubPropId(int)
     */
    private void setPubPropId(int pubPropId) {
        this.pubPropId = pubPropId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubPropI#getCvTerm()
     */
    public CvTerm getCvTerm() {
        return this.cvTerm;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubPropI#setCvTerm(org.gmod.schema.cv.CvTermI)
     */
    private void setCvTerm(CvTerm cvTerm) {
        this.cvTerm = cvTerm;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubPropI#getPub()
     */
    private Pub getPub() {
        return this.pub;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubPropI#setPub(org.gmod.schema.pub.PubI)
     */
    private void setPub(Pub pub) {
        this.pub = pub;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubPropI#getValue()
     */
    private String getValue() {
        return this.value;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubPropI#setValue(java.lang.String)
     */
    private void setValue(String value) {
        this.value = value;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubPropI#getRank()
     */
    public Integer getRank() {
        return this.rank;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubPropI#setRank(java.lang.Integer)
     */
    private void setRank(Integer rank) {
        this.rank = rank;
    }




}


