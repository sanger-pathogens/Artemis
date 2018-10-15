package org.gmod.schema.cv;



import java.io.Serializable;











public class CvTermPath implements Serializable {

    // Fields    


     private int cvTermPathId;
    


     private CvTerm cvTermBySubjectId;
     

        

     private CvTerm cvTermByObjectId;
     

        

     private CvTerm cvTermByTypeId;
     

        

     private Cv cv;
     

     private Integer pathDistance;
    
   
    // Property accessors

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermPathI#getCvTermPathId()
     */
    private int getCvTermPathId() {
        return this.cvTermPathId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermPathI#setCvTermPathId(int)
     */
    private void setCvTermPathId(int cvTermPathId) {
        this.cvTermPathId = cvTermPathId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermPathI#getCvTermBySubjectId()
     */
    private CvTerm getCvTermBySubjectId() {
        return this.cvTermBySubjectId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermPathI#setCvTermBySubjectId(org.gmod.schema.cv.CvTermI)
     */
    private void setCvTermBySubjectId(CvTerm cvTermBySubjectId) {
        this.cvTermBySubjectId = cvTermBySubjectId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermPathI#getCvTermByObjectId()
     */
    private CvTerm getCvTermByObjectId() {
        return this.cvTermByObjectId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermPathI#setCvTermByObjectId(org.gmod.schema.cv.CvTermI)
     */
    private void setCvTermByObjectId(CvTerm cvTermByObjectId) {
        this.cvTermByObjectId = cvTermByObjectId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermPathI#getCvTermByTypeId()
     */
    private CvTerm getCvTermByTypeId() {
        return this.cvTermByTypeId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermPathI#setCvTermByTypeId(org.gmod.schema.cv.CvTermI)
     */
    private void setCvTermByTypeId(CvTerm cvTermByTypeId) {
        this.cvTermByTypeId = cvTermByTypeId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermPathI#getCv()
     */
    private Cv getCv() {
        return this.cv;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermPathI#setCv(org.gmod.schema.cv.CvI)
     */
    private void setCv(Cv cv) {
        this.cv = cv;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermPathI#getPathDistance()
     */
    private Integer getPathDistance() {
        return this.pathDistance;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermPathI#setPathDistance(java.lang.Integer)
     */
    private void setPathDistance(Integer pathDistance) {
        this.pathDistance = pathDistance;
    }




}


