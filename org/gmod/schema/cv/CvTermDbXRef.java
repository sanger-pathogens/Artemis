package org.gmod.schema.cv;



import org.gmod.schema.general.DbXRef;

import java.io.Serializable;












public class CvTermDbXRef implements Serializable {

    // Fields    


     private int cvTermDbXRefId;
    


     private CvTerm cvTerm;
     


     private DbXRef dbXRef;
     

     private int isForDefinition;
    
   
    // Property accessors

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermDbXRefI#getCvTermDbXRefId()
     */
    private int getCvTermDbXRefId() {
        return this.cvTermDbXRefId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermDbXRefI#setCvTermDbXRefId(int)
     */
    private void setCvTermDbXRefId(int cvTermDbXRefId) {
        this.cvTermDbXRefId = cvTermDbXRefId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermDbXRefI#getCvTerm()
     */
    private CvTerm getCvTerm() {
        return this.cvTerm;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermDbXRefI#setCvTerm(org.gmod.schema.cv.CvTermI)
     */
    private void setCvTerm(CvTerm cvTerm) {
        this.cvTerm = cvTerm;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermDbXRefI#getDbXRef()
     */
    private DbXRef getDbXRef() {
        return this.dbXRef;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermDbXRefI#setDbXRef(org.gmod.schema.general.DbXRefI)
     */
    private void setDbXRef(DbXRef dbXRef) {
        this.dbXRef = dbXRef;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermDbXRefI#getIsForDefinition()
     */
    private int getIsForDefinition() {
        return this.isForDefinition;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvTermDbXRefI#setIsForDefinition(int)
     */
    private void setIsForDefinition(int isForDefinition) {
        this.isForDefinition = isForDefinition;
    }




}


