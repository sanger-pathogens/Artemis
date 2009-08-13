package org.gmod.schema.cv;



import java.io.Serializable;












public class CvTermProp implements Serializable {

    // Fields    


     private int cvTermPropId;
     


     private CvTerm cvTermByCvTermId;
     


     private CvTerm cvTermByTypeId;
     

     private String value;
     

     private int rank;
    
   
    // Property accessors
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermPropI#getCvTermpropId()
     */
    private int getCvTermPropId() {
        return this.cvTermPropId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermPropI#setCvTermpropId(int)
     */
    private void setCvTermPropId(int cvTermPropId) {
        this.cvTermPropId = cvTermPropId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermPropI#getCvTermByCvTermId()
     */
    private CvTerm getCvTermByCvTermId() {
       return this.cvTermByCvTermId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermPropI#setCvTermByCvTermId(org.gmod.schema.cv.CvTermI)
     */
    private void setCvTermByCvTermId(CvTerm cvTermByCvTermId) {
        this.cvTermByCvTermId = cvTermByCvTermId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermPropI#getCvTermByTypeId()
     */
    private CvTerm getCvTermByTypeId() {
        return this.cvTermByTypeId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermPropI#setCvTermByTypeId(org.gmod.schema.cv.CvTermI)
     */
    private void setCvTermByTypeId(CvTerm cvTermByTypeId) {
        this.cvTermByTypeId = cvTermByTypeId;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermPropI#getValue()
     */
    private String getValue() {
        return this.value;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermPropI#setValue(java.lang.String)
     */
    private void setValue(String value) {
        this.value = value;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermPropI#getRank()
     */
    private int getRank() {
        return this.rank;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermPropI#setRank(int)
     */
    private void setRank(int rank) {
        this.rank = rank;
    }




}


