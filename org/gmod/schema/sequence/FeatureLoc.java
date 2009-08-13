package org.gmod.schema.sequence;

import java.io.Serializable;
import java.util.HashSet;
import java.util.Set;

public class FeatureLoc implements Serializable {

    // Fields 
    private int featureLocId;
    private Feature featureBySrcFeatureId;
    private Feature featureByFeatureId;
    private Integer fmin;
    private boolean fminPartial;
    private Integer fmax;
    private boolean fmaxPartial;
    private Short strand;
    private Integer phase;
    private String residueInfo;
    private int locGroup;
    private int rank;
    private Set<FeatureLocPub> featureLocPubs = new HashSet<FeatureLocPub>(0);

    // used by artemis
    private int srcFeatureId;
    
     // Constructors

    /** default constructor */
    public FeatureLoc() {
    	// Deliberately empty default constructor
    }

	/** minimal constructor */

	private FeatureLoc(Feature featureBySrcFeatureId, boolean fminPartial, boolean fmaxPartial, int locGroup, int rank) {
        this.featureBySrcFeatureId = featureBySrcFeatureId;
        this.fminPartial = fminPartial;
        this.fmaxPartial = fmaxPartial;
        this.locGroup = locGroup;
        this.rank = rank;
    }
    /** large constructor */
    public FeatureLoc(Feature featureBySrcfeatureId, Feature featureByFeatureId, Integer fmin, boolean fminPartial, Integer fmax, boolean fmaxPartial, Short strand, Integer phase, int locGroup, int rank) {
       this.featureBySrcFeatureId = featureBySrcfeatureId;
       this.featureByFeatureId = featureByFeatureId;
       this.fmin = fmin;
       this.fminPartial = fminPartial;
       this.fmax = fmax;
       this.fmaxPartial = fmaxPartial;
       this.strand = strand;
       this.phase = phase;
       this.locGroup = locGroup;
       this.rank = rank;
    }
    
   
    // Property accessors


    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureLocI#getFeatureBySrcfeatureId()
     */
    public Feature getFeatureBySrcFeatureId() {
        return this.featureBySrcFeatureId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureLocI#setFeatureBySrcfeatureId(org.genedb.db.jpa.Feature)
     */
    public void setFeatureBySrcFeatureId(Feature featureBySrcFeatureId) {
        this.featureBySrcFeatureId = featureBySrcFeatureId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureLocI#getFeatureByFeatureId()
     */
    public Feature getFeatureByFeatureId() {
        return this.featureByFeatureId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureLocI#setFeatureByFeatureId(org.genedb.db.jpa.Feature)
     */
    public void setFeatureByFeatureId(Feature featureByFeatureId) {
        this.featureByFeatureId = featureByFeatureId;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureLocI#getFmin()
     */
    public Integer getFmin() {
        return this.fmin;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureLocI#setFmin(java.lang.Integer)
     */
    public void setFmin(Integer fmin) {
        this.fmin = fmin;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureLocI#isFminPartial()
     */
    public boolean isFminPartial() {
        return this.fminPartial;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureLocI#setFminPartial(boolean)
     */
    public void setFminPartial(boolean fminPartial) {
        this.fminPartial = fminPartial;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureLocI#getFmax()
     */
    public Integer getFmax() {
        return this.fmax;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureLocI#setFmax(java.lang.Integer)
     */
    public void setFmax(Integer fmax) {
        this.fmax = fmax;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureLocI#isFmaxPartial()
     */
    public boolean isFmaxPartial() {
        return this.fmaxPartial;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureLocI#setFmaxPartial(boolean)
     */
    public void setFmaxPartial(boolean fmaxPartial) {
        this.fmaxPartial = fmaxPartial;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureLocI#getStrand()
     */
    public Short getStrand() {
        return this.strand;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureLocI#setStrand(java.lang.Short)
     */
    public void setStrand(Short strand) {
        this.strand = strand;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureLocI#getPhase()
     */
    public Integer getPhase() {
        return this.phase;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureLocI#setPhase(java.lang.Integer)
     */
    public void setPhase(Integer phase) {
        this.phase = phase;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureLocI#getResidueInfo()
     */
    public String getResidueInfo() {
        return this.residueInfo;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureLocI#setResidueInfo(java.lang.String)
     */
    public void setResidueInfo(String residueInfo) {
        this.residueInfo = residueInfo;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureLocI#getLocgroup()
     */
    public int getLocGroup() {
        return this.locGroup;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureLocI#setLocgroup(int)
     */
    public void setLocGroup(int locGroup) {
        this.locGroup = locGroup;
    }
    
 
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureLocI#getRank()
     */
    public int getRank() {
        return this.rank;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureLocI#setRank(int)
     */
    public void setRank(int rank) {
        this.rank = rank;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureLocI#getFeaturelocPubs()
     */
    public Set<FeatureLocPub> getFeatureLocPubs() {
        return this.featureLocPubs;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureLocI#setFeaturelocPubs(java.util.Set)
     */
    public void setFeaturelocPubs(Set<FeatureLocPub> featureLocPubs) {
        this.featureLocPubs = featureLocPubs;
    }


    public int getFeatureLocId() {
        return this.featureLocId;
    }


    public void setFeatureLocId(int featureLocId) {
        this.featureLocId = featureLocId;
    }

    public int getSrcFeatureId() {
      return srcFeatureId;
    }

    public void setSrcFeatureId(int srcFeatureId) {
      this.srcFeatureId = srcFeatureId;
    }

}
