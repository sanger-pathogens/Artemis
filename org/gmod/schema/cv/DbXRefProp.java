package org.gmod.schema.cv;




import org.gmod.schema.general.DbXRef;
import org.gmod.schema.utils.propinterface.PropertyI;

import java.io.Serializable;











public class DbXRefProp implements Serializable, PropertyI {

    // Fields    


     private int dbXRefPropId;
     


     private CvTerm cvTerm;
     


     private DbXRef dbXRef;
     

     private String value;
     

     private int rank;
   
    // Property accessors
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.DbXRefPropI#getDbXRefpropId()
     */
    private int getDbXRefPropId() {
        return this.dbXRefPropId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.DbXRefPropI#setDbXRefpropId(int)
     */
    private void setDbXRefPropId(int dbXRefPropId) {
        this.dbXRefPropId = dbXRefPropId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.DbXRefPropI#getCvTerm()
     */
    public CvTerm getCvTerm() {
        return this.cvTerm;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.DbXRefPropI#setCvTerm(org.gmod.schema.cv.CvTermI)
     */
    private void setCvTerm(CvTerm cvTerm) {
        this.cvTerm = cvTerm;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.DbXRefPropI#getDbXRef()
     */
    private DbXRef getDbXRef() {
        return this.dbXRef;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.DbXRefPropI#setDbXRef(org.gmod.schema.general.DbXRefI)
     */
    private void setDbXRef(DbXRef dbXRef) {
        this.dbXRef = dbXRef;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.DbXRefPropI#getValue()
     */
    private String getValue() {
        return this.value;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.DbXRefPropI#setValue(java.lang.String)
     */
    private void setValue(String value) {
        this.value = value;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.DbXRefPropI#getRank()
     */
    private int getRank() {
        return this.rank;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.DbXRefPropI#setRank(int)
     */
    private void setRank(int rank) {
        this.rank = rank;
    }




}


