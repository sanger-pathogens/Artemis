package org.gmod.schema.sequence;

import org.gmod.schema.general.DbXRef;

import java.io.Serializable;

public class FeatureCvTermDbXRef implements Serializable {

    // Fields    
    private int featureCvTermDbXRefId;
    private DbXRef dbXRef;
    private FeatureCvTerm featureCvTerm;

     // Constructors

    /** default constructor */
    public FeatureCvTermDbXRef() {
    	// Deliberately empty default constructor
    }

    /** full constructor */
    public FeatureCvTermDbXRef(DbXRef dbXRef, FeatureCvTerm featureCvTerm) {
       this.dbXRef = dbXRef;
       this.featureCvTerm = featureCvTerm;
    }
    
   
    // Property accessors
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureCvTermDbXRefI#getFeatureCvTermDbXrefId()
     */
    private int getFeatureCvTermDbXRefId() {
        return this.featureCvTermDbXRefId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureCvTermDbXRefI#setFeatureCvTermDbXrefId(int)
     */
    private void setFeatureCvTermDbXRefId(int featureCvTermDbXRefId) {
        this.featureCvTermDbXRefId = featureCvTermDbXRefId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureCvTermDbXRefI#getDbxref()
     */
    public DbXRef getDbXRef() {
        return this.dbXRef;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureCvTermDbXRefI#setDbxref(org.gmod.schema.general.DbXRefI)
     */
    public void setDbXRef(DbXRef dbXRef) {
        this.dbXRef = dbXRef;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureCvTermDbXRefI#getFeatureCvterm()
     */
    public FeatureCvTerm getFeatureCvTerm() {
        return this.featureCvTerm;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureCvTermDbXRefI#setFeatureCvterm(org.gmod.schema.sequence.FeatureCvTermI)
     */
    public void setFeatureCvTerm(FeatureCvTerm featureCvTerm) {
        this.featureCvTerm = featureCvTerm;
    }




}


