package org.gmod.schema.sequence;

import org.gmod.schema.cv.CvTerm;
import org.gmod.schema.pub.Pub;
import org.gmod.schema.utils.Rankable;
import org.gmod.schema.utils.propinterface.PropertyI;

import java.io.Serializable;
import java.util.Collection;
import java.util.HashSet;

public class FeatureCvTerm implements Serializable, Rankable, PropertyI {

    // Fields 
    private int featureCvTermId;
    private CvTerm cvTerm;
    private Feature feature;
    private Pub pub;
    private boolean not;
    private int rank;
    private Collection<FeatureCvTermProp> featureCvTermProps = new HashSet<FeatureCvTermProp>(0);
    private Collection<FeatureCvTermPub> featureCvTermPubs = new HashSet<FeatureCvTermPub>(0);
    private Collection<FeatureCvTermDbXRef> featureCvTermDbXRefs = new HashSet<FeatureCvTermDbXRef>(0);

     // Constructors

    /** default constructor */
    public FeatureCvTerm() {
    	// Deliberately empty default constructor
    }

	/** minimal constructor */
    public FeatureCvTerm(CvTerm cvTerm, Feature feature, Pub pub, boolean not,int rank) {
        this.cvTerm = cvTerm;
        this.feature = feature;
        this.pub = pub;
        this.not = not;
        this.rank = rank;
    }
    
   
    // Property accessors
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureCvTermI#getFeatureCvTermId()
     */
    public int getFeatureCvTermId() {
        return this.featureCvTermId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureCvTermI#setFeatureCvTermId(int)
     */
    public void setFeatureCvTermId(int featureCvTermId) {
        this.featureCvTermId = featureCvTermId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureCvTermI#getCvterm()
     */
    public CvTerm getCvTerm() {
        return this.cvTerm;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureCvTermI#setCvterm(org.gmod.schema.cv.CvTermI)
     */
    public void setCvTerm(CvTerm cvTerm) {
        this.cvTerm = cvTerm;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureCvTermI#getFeature()
     */
    public Feature getFeature() {
        return this.feature;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureCvTermI#setFeature(org.genedb.db.jpa.Feature)
     */
    public void setFeature(Feature feature) {
        this.feature = feature;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureCvTermI#getPub()
     */
    public Pub getPub() {
        return this.pub;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureCvTermI#setPub(org.gmod.schema.pub.PubI)
     */
    public void setPub(Pub pub) {
        this.pub = pub;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureCvTermI#isNot()
     */
    public boolean isNot() {
        return this.not;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureCvTermI#setNot(boolean)
     */
    public void setNot(boolean not) {
        this.not = not;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureCvTermI#getFeatureCvtermprops()
     */
    public Collection<FeatureCvTermProp> getFeatureCvTermProps() {
        return this.featureCvTermProps;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureCvTermI#setFeatureCvtermprops(java.util.Set)
     */
    public void setFeatureCvTermProps(Collection<FeatureCvTermProp> featureCvTermProps) {
        this.featureCvTermProps = featureCvTermProps;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureCvTermI#getFeatureCvtermPubs()
     */
    public Collection<FeatureCvTermPub> getFeatureCvTermPubs() {
        return this.featureCvTermPubs;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureCvTermI#setFeatureCvtermPubs(java.util.Set)
     */
    public void setFeatureCvTermPubs(Collection<FeatureCvTermPub> featureCvTermPubs) {
        this.featureCvTermPubs = featureCvTermPubs;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureCvTermI#getFeatureCvtermDbxrefs()
     */
    public Collection<FeatureCvTermDbXRef> getFeatureCvTermDbXRefs() {
        return this.featureCvTermDbXRefs;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureCvTermI#setFeatureCvtermDbxrefs(java.util.Set)
     */
    public void setFeatureCvTermDbXRefs(Collection<FeatureCvTermDbXRef> featureCvTermDbXRefs) {
        this.featureCvTermDbXRefs = featureCvTermDbXRefs;
    }

	public int getRank() {
		return rank;
	}

	public void setRank(int rank) {
		this.rank = rank;
	}




}


