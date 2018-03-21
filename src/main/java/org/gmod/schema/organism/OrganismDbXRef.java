package org.gmod.schema.organism;

import org.gmod.schema.general.DbXRef;

import java.io.Serializable;











public class OrganismDbXRef implements Serializable {

    // Fields    

    

     private int organismDbXRefId;
     

         

     private Organism organism;
     

         

     private DbXRef dbXRef;
    
   
    // Property accessors

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.OrganismDbXRefI#getOrganismDbXRefId()
     */
    private int getOrganismDbXRefId() {
        return this.organismDbXRefId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.OrganismDbXRefI#setOrganismDbXRefId(int)
     */
    private void setOrganismDbXRefId(int organismDbXRefId) {
        this.organismDbXRefId = organismDbXRefId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.OrganismDbXRefI#getOrganism()
     */
    private Organism getOrganism() {
        return this.organism;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.OrganismDbXRefI#setOrganism(org.gmod.schema.organism.OrganismI)
     */
    private void setOrganism(Organism organism) {
        this.organism = organism;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.OrganismDbXRefI#getDbXRef()
     */
    private DbXRef getDbXRef() {
        return this.dbXRef;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.OrganismDbXRefI#setDbXRef(org.gmod.schema.general.DbXRefI)
     */
    private void setDbXRef(DbXRef dbXRef) {
        this.dbXRef = dbXRef;
    }




}


