package org.gmod.schema.cv;



import java.io.Serializable;











public class CvTermRelationship implements Serializable {

    // Fields    


     private int cvTermRelationshipId;
     


     private CvTerm cvTermBySubjectId;
     


     private CvTerm cvTermByObjectId;
     


     private CvTerm cvTermByTypeId;

     // Constructors

    /** default constructor */
    public CvTermRelationship() {
    	// Deliberately empty default constructor
    }

    /** full constructor */
    public CvTermRelationship(CvTerm cvTermBySubjectId, CvTerm cvTermByObjectId, CvTerm cvTermByTypeId) {
       this.cvTermBySubjectId = cvTermBySubjectId;
       this.cvTermByObjectId = cvTermByObjectId;
       this.cvTermByTypeId = cvTermByTypeId;
    }
    
   
    // Property accessors
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermRelationshipI#getCvTermRelationshipId()
     */
    private int getCvTermRelationshipId() {
        return this.cvTermRelationshipId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermRelationshipI#setCvTermRelationshipId(int)
     */
    private void setCvTermRelationshipId(int cvTermRelationshipId) {
        this.cvTermRelationshipId = cvTermRelationshipId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermRelationshipI#getCvTermBySubjectId()
     */
    private CvTerm getCvTermBySubjectId() {
        return this.cvTermBySubjectId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermRelationshipI#setCvTermBySubjectId(org.gmod.schema.cv.CvTermI)
     */
    private void setCvTermBySubjectId(CvTerm cvTermBySubjectId) {
        this.cvTermBySubjectId = cvTermBySubjectId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermRelationshipI#getCvTermByObjectId()
     */
    private CvTerm getCvTermByObjectId() {
        return this.cvTermByObjectId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermRelationshipI#setCvTermByObjectId(org.gmod.schema.cv.CvTermI)
     */
    private void setCvTermByObjectId(CvTerm cvTermByObjectId) {
        this.cvTermByObjectId = cvTermByObjectId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermRelationshipI#getCvTermByTypeId()
     */
    private CvTerm getCvTermByTypeId() {
        return this.cvTermByTypeId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermRelationshipI#setCvTermByTypeId(org.gmod.schema.cv.CvTermI)
     */
    private void setCvTermByTypeId(CvTerm cvTermByTypeId) {
        this.cvTermByTypeId = cvTermByTypeId;
    }




}


