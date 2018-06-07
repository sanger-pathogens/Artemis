package org.gmod.schema.cv;



import java.io.Serializable;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;










public class Cv implements Serializable {

    // Fields    


     private int cvId;
    

     private String name;
    

     private String definition;
    

     private Set<CvTermPath> cvTermPaths = new HashSet<CvTermPath>(0);
    

     private Set<CvTerm> cvTerms = new HashSet<CvTerm>(0);
    
   
    // Property accessors

    public int getCvId() {
        return this.cvId;
    }
    
    public void setCvId(int cvId) {
        this.cvId = cvId;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvI#getName()
     */
    public String getName() {
        return this.name;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvI#setName(java.lang.String)
     */
    public void setName(String name) {
        this.name = name;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvI#getDefinition()
     */
    public String getDefinition() {
        return this.definition;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvI#setDefinition(java.lang.String)
     */
    public void setDefinition(String definition) {
        this.definition = definition;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvI#getCvTermPaths()
     */
    public Collection<CvTermPath> getCvTermPaths() {
        return this.cvTermPaths;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvI#setCvTermPaths(java.util.Set)
     */
    public void setCvTermPaths(Set<CvTermPath> cvTermPaths) {
        this.cvTermPaths = cvTermPaths;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvI#getCvTerms()
     */
    public Collection<CvTerm> getCvTerms() {
        return this.cvTerms;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.CvI#setCvTerms(java.util.Set)
     */
    public void setCvTerms(Set<CvTerm> cvTerms) {
        this.cvTerms = cvTerms;
    }

}


