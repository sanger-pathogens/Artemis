package org.gmod.schema.sequence;

import org.gmod.schema.cv.CvTerm;

import java.io.Serializable;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;

public class Synonym implements Serializable {

    // Fields    
    private int synonymId;
    private CvTerm cvTerm;
    private String name;
    private String synonymSgml;
    private Set<FeatureSynonym> featureSynonyms = new HashSet<FeatureSynonym>(0);

     // Constructors

    /** default constructor */
    public Synonym() {
    	// Deliberately empty default constructor
    }

	/** minimal constructor */
    public Synonym(CvTerm cvTerm, String name, String synonymSgml) {
        this.cvTerm = cvTerm;
        this.name = name;
        this.synonymSgml = synonymSgml;
    }
    
   
    // Property accessors

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.SynonymI#getSynonymId()
     */
    public int getSynonymId() {
        return this.synonymId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.SynonymI#setSynonymId(int)
     */
    public void setSynonymId(int synonymId) {
        this.synonymId = synonymId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.SynonymI#getCvterm()
     */
    public CvTerm getCvTerm() {
        return this.cvTerm;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.SynonymI#setCvterm(org.gmod.schema.cv.CvTermI)
     */
    public void setCvTerm(CvTerm cvterm) {
        this.cvTerm = cvterm;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.SynonymI#getName()
     */
    public String getName() {
        return this.name;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.SynonymI#setName(java.lang.String)
     */
    public void setName(String name) {
        this.name = name;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.SynonymI#getSynonymSgml()
     */
    private String getSynonymSgml() {
        return this.synonymSgml;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.SynonymI#setSynonymSgml(java.lang.String)
     */
    public void setSynonymSgml(String synonymSgml) {
        this.synonymSgml = synonymSgml;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.SynonymI#getFeatureSynonyms()
     */
    private Collection<FeatureSynonym> getFeatureSynonyms() {
        return this.featureSynonyms;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.SynonymI#setFeatureSynonyms(java.util.Set)
     */
    private void setFeatureSynonyms(Set<FeatureSynonym> featureSynonyms) {
        this.featureSynonyms = featureSynonyms;
    }




}


