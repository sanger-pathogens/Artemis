package org.gmod.schema.pub;

import org.gmod.schema.general.DbXRef;

import java.io.Serializable;











public class PubDbXRef implements Serializable {

    // Fields    

    

     private int pubDbXRefId;
     

         

     private DbXRef dbXRef;
     

         

     private Pub pub;
     

     private boolean current;

     // Constructors

    /** default constructor */
    public PubDbXRef() {
    	// Deliberately empty default constructor
    }

    /** full constructor */
    public PubDbXRef(Pub pub, DbXRef dbXRef, boolean current) {
       this.dbXRef = dbXRef;
       this.pub = pub;
       this.current = current;
    }
    
   
    // Property accessors

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubDbXRefI#getPubDbXRefId()
     */
    private int getPubDbXRefId() {
        return this.pubDbXRefId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubDbXRefI#setPubDbXRefId(int)
     */
    private void setPubDbXRefId(int pubDbXRefId) {
        this.pubDbXRefId = pubDbXRefId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubDbXRefI#getDbXRef()
     */
    public DbXRef getDbXRef() {
        return this.dbXRef;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubDbXRefI#setDbXRef(org.gmod.schema.general.DbXRefI)
     */
    public void setDbXRef(DbXRef dbXRef) {
        this.dbXRef = dbXRef;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubDbXRefI#getPub()
     */
    public Pub getPub() {
        return this.pub;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubDbXRefI#setPub(org.gmod.schema.pub.PubI)
     */
    public void setPub(Pub pub) {
        this.pub = pub;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubDbXRefI#isCurrent()
     */
    public boolean isCurrent() {
        return this.current;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubDbXRefI#setCurrent(boolean)
     */
    public void setCurrent(boolean current) {
        this.current = current;
    }




}


