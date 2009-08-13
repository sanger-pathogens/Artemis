package org.gmod.schema.organism;



import org.gmod.schema.cv.CvTerm;
import org.gmod.schema.utils.propinterface.PropertyI;

import java.io.Serializable;











public class OrganismProp implements Serializable, PropertyI {

    // Fields    

    

     private int organismPropId;
     

         

     private Organism organism;
     

         

     private CvTerm cvTerm;
     

     private String value;
     

     private int rank;
    
   
    // Property accessors

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.OrganismPropI#getOrganismPropId()
     */
    private int getOrganismPropId() {
        return this.organismPropId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.OrganismPropI#setOrganismPropId(int)
     */
    private void setOrganismPropId(int organismPropId) {
        this.organismPropId = organismPropId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.OrganismPropI#getOrganism()
     */
    private Organism getOrganism() {
        return this.organism;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.OrganismPropI#setOrganism(org.gmod.schema.organism.OrganismI)
     */
    private void setOrganism(Organism organism) {
        this.organism = organism;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.OrganismPropI#getCvTerm()
     */
    public CvTerm getCvTerm() {
        return this.cvTerm;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.OrganismPropI#setCvTerm(org.gmod.schema.cv.CvTermI)
     */
    private void setCvTerm(CvTerm cvTerm) {
        this.cvTerm = cvTerm;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.OrganismPropI#getValue()
     */
    public String getValue() {
        return this.value;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.OrganismPropI#setValue(java.lang.String)
     */
    private void setValue(String value) {
        this.value = value;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.OrganismPropI#getRank()
     */
    private int getRank() {
        return this.rank;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.OrganismPropI#setRank(int)
     */
    private void setRank(int rank) {
        this.rank = rank;
    }




}


