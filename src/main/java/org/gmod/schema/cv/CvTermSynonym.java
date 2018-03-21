package org.gmod.schema.cv;




import java.io.Serializable;











public class CvTermSynonym implements Serializable {

    // Fields    


     private int cvTermSynonymId;
     


     private CvTerm cvTermByCvTermId;
     


     private CvTerm cvTermByTypeId;
     

     private String synonym;
    
   
    // Property accessors

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermSynonymI#getCvTermSynonymId()
     */
    private int getCvTermSynonymId() {
        return this.cvTermSynonymId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermSynonymI#setCvTermSynonymId(int)
     */
    private void setCvTermSynonymId(int cvTermSynonymId) {
        this.cvTermSynonymId = cvTermSynonymId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermSynonymI#getCvTermByCvTermId()
     */
    private CvTerm getCvTermByCvTermId() {
        return this.cvTermByCvTermId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermSynonymI#setCvTermByCvTermId(org.gmod.schema.cv.CvTermI)
     */
    private void setCvTermByCvTermId(CvTerm cvTermByCvTermId) {
        this.cvTermByCvTermId = cvTermByCvTermId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermSynonymI#getCvTermByTypeId()
     */
    private CvTerm getCvTermByTypeId() {
        return this.cvTermByTypeId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermSynonymI#setCvTermByTypeId(org.gmod.schema.cv.CvTermI)
     */
    private void setCvTermByTypeId(CvTerm cvTermByTypeId) {
        this.cvTermByTypeId = cvTermByTypeId;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermSynonymI#getSynonym()
     */
    private String getSynonym() {
        return this.synonym;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermSynonymI#setSynonym(java.lang.String)
     */
    private void setSynonym(String synonym) {
        this.synonym = synonym;
    }




}


