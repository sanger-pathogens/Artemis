package org.gmod.schema.pub;

import org.gmod.schema.cv.CvTerm;

import java.io.Serializable;











public class PubRelationship implements Serializable {

    // Fields    

    

     private int pubRelationshipId;
     

         

     private Pub pubBySubjectId;
     

         

     private Pub pubByObjectId;
     

         

     private CvTerm cvTerm;
    
   
    // Property accessors

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubRelationshipI#getPubRelationshipId()
     */
    private int getPubRelationshipId() {
        return this.pubRelationshipId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubRelationshipI#setPubRelationshipId(int)
     */
    private void setPubRelationshipId(int pubRelationshipId) {
        this.pubRelationshipId = pubRelationshipId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubRelationshipI#getPubBySubjectId()
     */
    private Pub getPubBySubjectId() {
        return this.pubBySubjectId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubRelationshipI#setPubBySubjectId(org.gmod.schema.pub.PubI)
     */
    private void setPubBySubjectId(Pub pubBySubjectId) {
        this.pubBySubjectId = pubBySubjectId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubRelationshipI#getPubByObjectId()
     */
    private Pub getPubByObjectId() {
        return this.pubByObjectId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubRelationshipI#setPubByObjectId(org.gmod.schema.pub.PubI)
     */
    private void setPubByObjectId(Pub pubByObjectId) {
        this.pubByObjectId = pubByObjectId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubRelationshipI#getCvTerm()
     */
    private CvTerm getCvTerm() {
        return this.cvTerm;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubRelationshipI#setCvTerm(org.gmod.schema.cv.CvTermI)
     */
    private void setCvTerm(CvTerm cvTerm) {
        this.cvTerm = cvTerm;
    }




}


