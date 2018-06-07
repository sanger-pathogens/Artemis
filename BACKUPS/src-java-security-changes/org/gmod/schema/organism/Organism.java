package org.gmod.schema.organism;



import org.gmod.schema.phylogeny.PhylonodeOrganism;
import org.gmod.schema.sequence.Feature;

import java.io.Serializable;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;











public class Organism implements Serializable {
    
    private Set<PhylonodeOrganism> phylonodeOrganisms = new HashSet<PhylonodeOrganism>(0);

    public Set<PhylonodeOrganism> getPhylonodeOrganisms() {
        return this.phylonodeOrganisms;
    }
    
    public void setPhylonodeOrganisms(Set<PhylonodeOrganism> phylonodeOrganisms) {
        this.phylonodeOrganisms = phylonodeOrganisms;
    }

    // Fields    

    

     private int organismId;
     

     private String abbreviation;
     

     private String genus;
     

     private String species;
     

     private String commonName;
     

     private String comment;
     

     private Set<OrganismProp> organismProps = new HashSet<OrganismProp>(0);
     

     private Set<Feature> features = new HashSet<Feature>(0);
     

     private Set<OrganismDbXRef> organismDbXRefs = new HashSet<OrganismDbXRef>(0);
    
   
    // Property accessors

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.OrganismI#getOrganismId()
     */
    public int getOrganismId() {
        return this.organismId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.OrganismI#setOrganismId(int)
     */
    public void setOrganismId(int organismId) {
        this.organismId = organismId;
    }
    
    public String getFullName() {
        return getGenus()+' '+getSpecies();
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.OrganismI#getAbbreviation()
     */
    public String getAbbreviation() {
        return this.abbreviation;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.OrganismI#setAbbreviation(java.lang.String)
     */
    public void setAbbreviation(String abbreviation) {
        this.abbreviation = abbreviation;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.OrganismI#getGenus()
     */
    public String getGenus() {
        return this.genus;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.OrganismI#setGenus(java.lang.String)
     */
    public void setGenus(String genus) {
        this.genus = genus;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.OrganismI#getSpecies()
     */
    public String getSpecies() {
        return this.species;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.OrganismI#setSpecies(java.lang.String)
     */
    public void setSpecies(String species) {
        this.species = species;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.OrganismI#getCommonName()
     */
    public String getCommonName() {
        return this.commonName;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.OrganismI#setCommonName(java.lang.String)
     */
    public void setCommonName(String commonName) {
        this.commonName = commonName;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.OrganismI#getComment()
     */
    public String getComment() {
        return this.comment;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.OrganismI#setComment(java.lang.String)
     */
    public void setComment(String comment) {
        this.comment = comment;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.OrganismI#getOrganismProps()
     */
    public Set<OrganismProp> getOrganismProps() {
        return this.organismProps;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.OrganismI#setOrganismProps(java.util.Set)
     */
    private void setOrganismProps(Set<OrganismProp> organismProps) {
        this.organismProps = organismProps;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.OrganismI#getFeatures()
     */
    private Collection<Feature> getFeatures() {
        return this.features;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.OrganismI#setFeatures(java.util.Set)
     */
    private void setFeatures(Set<Feature> features) {
        this.features = features;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.OrganismI#getOrganismDbXRefs()
     */
    private Collection<OrganismDbXRef> getOrganismDbXRefs() {
        return this.organismDbXRefs;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.OrganismI#setOrganismDbXRefs(java.util.Set)
     */
    private void setOrganismDbXRefs(Set<OrganismDbXRef> organismDbXRefs) {
        this.organismDbXRefs = organismDbXRefs;
    }




}


