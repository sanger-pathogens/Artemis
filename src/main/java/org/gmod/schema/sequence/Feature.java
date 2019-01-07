package org.gmod.schema.sequence;


import org.gmod.schema.analysis.AnalysisFeature;
import org.gmod.schema.cv.CvTerm;
import org.gmod.schema.general.DbXRef;
import org.gmod.schema.organism.Organism;
import org.gmod.schema.phylogeny.Phylonode;
import org.gmod.schema.utils.CollectionUtils;




import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.Collection;
import java.util.Date;
import java.sql.Timestamp;

public class Feature implements java.io.Serializable {

    private int featureId;
    private Organism organism;
    private CvTerm cvTerm;
    private String name;
    private String uniqueName;
    private Integer seqLen = -1;
    private String md5Checksum;
    private boolean analysis;
    private boolean obsolete;
    private Timestamp timeAccessioned;
    private Timestamp timeLastModified;

    // -------------------------------------------------------------------------------
    // Unsorted properties below here


    private Collection<Phylonode> phylonodes;
    private DbXRef dbXRef;
    private byte residues[];
    private Collection<FeatureLoc> featureLocsForSrcFeatureId;
    private Collection<FeatureRelationship> featureRelationshipsForObjectId;
    private Collection<FeatureRelationship> featureRelationshipsForSubjectId;
    private Collection<FeatureDbXRef> featureDbXRefs;
    private Collection<FeatureLoc> featureLocsForFeatureId;
    private Collection<FeatureCvTerm> featureCvTerms;

    //@Cascade( {CascadeType.ALL, CascadeType.DELETE_ORPHAN} )
    private Collection<FeatureProp> featureProps;
    private Collection<FeaturePub> featurePubs;
    private Collection<AnalysisFeature> analysisFeatures;
    private Collection<FeatureSynonym> featureSynonyms; 
    private FeatureLoc featureLoc;
    
     // Constructors

    /** default constructor */
    public Feature() {
    	// Deliberately empty default constructor
    }

	/** minimal constructor */
    public Feature(Organism organism, CvTerm cvTerm, String uniqueName, boolean analysis, boolean obsolete, Timestamp timeAccessioned, Timestamp timeLastModified) {
        this.organism = organism;
        this.cvTerm = cvTerm;
        this.uniqueName = uniqueName;
        this.analysis = analysis;
        this.obsolete = obsolete;
        this.timeAccessioned = timeAccessioned;
        this.timeLastModified = timeLastModified;
    }
    
   
    // Property accessors

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#getFeatureId()
     */
    public int getFeatureId() {
        return this.featureId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#setFeatureId(int)
     */
    public void setFeatureId(int featureId) {
        this.featureId = featureId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#getOrganism()
     */
    public Organism getOrganism() {
        return this.organism;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#setOrganism(org.gmod.schema.organism.OrganismI)
     */
    public void setOrganism(Organism organism) {
        this.organism = organism;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#getCvTerm()
     */
    public CvTerm getCvTerm() {
        return this.cvTerm;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#setCvTerm(org.gmod.schema.cv.CvTermI)
     */
    public void setCvTerm(CvTerm cvTerm) {
        this.cvTerm = cvTerm;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#getDbxref()
     */
    public DbXRef getDbXRef() {
        return this.dbXRef;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#setDbxref(org.gmod.schema.general.DbXRefI)
     */
    public void setDbXRef(DbXRef dbXRef) {
        this.dbXRef = dbXRef;
    }
    

    /**
     * Get the human-readable form of the feature eg the gene name
     * 
     * @return the name, may be null
     */
    public String getName() {
        return this.name;
    }
    
    
    /**
     * Set the human-readable form of the feature eg the gene name
     * 
     * @param name the human-readable name
     */
    public void setName(String name) {
        this.name = name;
    }
    

    /**
     * Fetch the unique name (systematic id) for the feature 
     * 
     * @return the unique name, not null
     */
    /**
     * @return
     */
    public String getUniqueName() {
        return this.uniqueName;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#setUniquename(java.lang.String)
     */
    public void setUniqueName(String uniqueName) {
        this.uniqueName = uniqueName;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#getResidues()
     */
    public byte[] getResidues() {
        return this.residues;
    }
    
    /**
     * Fetch a subset of the sequence (may be lazy) 
     * 
     * @param min the lower bound, in interbase coordinates
     * @param max the upper bound, in interbase coordinates
     * @return
     */
    public byte[] getResidues(int min, int max) {
        byte[] results = new byte[max - min];
        System.arraycopy(getResidues(), 0, results, 0, max);
        return results;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#setResidues(java.lang.String)
     */
    public void setResidues(byte[] residues) {
        this.residues = residues;
        if (residues == null) {
            seqLen = 0;
            md5Checksum = "";
            return;
        }
        seqLen = residues.length;
        this.md5Checksum = calcMD5(this.residues);
    }
    

    /**
     * Fetch the length of the sequence. Find it from the parent if necessary 
     * 
     * @return the length
     */
    public int getSeqLen() {
        if (this.seqLen.intValue() == -1 && residues != null) {
            return getResidues().length;
        }
        return this.seqLen.intValue();
    }
    
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#setSeqlen(java.lang.Integer)
     */
    public void setSeqLen(Integer seqLen) {
    	if (seqLen == null) {
    		throw new IllegalArgumentException("Length of '"+uniqueName+"' attempted to be set to null");
    	}
        this.seqLen = seqLen;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#getMd5Checksum()
     */
    public String getMd5Checksum() {
        return this.md5Checksum;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#setMd5Checksum(java.lang.String)
     */
    public void setMd5Checksum(String md5Checksum) {
        this.md5Checksum = md5Checksum;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#isAnalysis()
     */
    public boolean isAnalysis() {
        return this.analysis;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#setAnalysis(boolean)
     */
    public void setAnalysis(boolean analysis) {
        this.analysis = analysis;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#isObsolete()
     */
    public boolean isObsolete() {
        return this.obsolete;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#setObsolete(boolean)
     */
    public void setObsolete(boolean obsolete) {
        this.obsolete = obsolete;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#getTimeAccessioned()
     */
    public Date getTimeAccessioned() {
        return this.timeAccessioned;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#setTimeAccessioned(java.util.Date)
     */
    public void setTimeAccessioned(Timestamp timeAccessioned) {
        this.timeAccessioned = timeAccessioned;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#getTimeLastModified()
     */
    public Timestamp getTimeLastModified() {
        return this.timeLastModified;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#setTimeLastModified(java.util.Date)
     */
    public void setTimeLastModified(Timestamp timeLastModified) {
        this.timeLastModified = timeLastModified;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#getFeaturelocsForSrcfeatureId()
     */
    public Collection<FeatureLoc> getFeatureLocsForSrcFeatureId() {
        return (featureLocsForSrcFeatureId = CollectionUtils.safeGetter(featureLocsForSrcFeatureId));
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#setFeaturelocsForSrcfeatureId(java.util.Set)
     */
    public void setFeatureLocsForSrcFeatureId(Collection<FeatureLoc> featureLocsForSrcFeatureId) {
        this.featureLocsForSrcFeatureId = featureLocsForSrcFeatureId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#getFeatureRelationshipsForObjectId()
     */
    public Collection<FeatureRelationship> getFeatureRelationshipsForObjectId() {
        return (featureRelationshipsForObjectId = CollectionUtils.safeGetter(featureRelationshipsForObjectId));
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#setFeatureRelationshipsForObjectId(java.util.Set)
     */
    public void setFeatureRelationshipsForObjectId(Collection<FeatureRelationship> featureRelationshipsForObjectId) {
        this.featureRelationshipsForObjectId = featureRelationshipsForObjectId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#getFeatureRelationshipsForSubjectId()
     */
    public Collection<FeatureRelationship> getFeatureRelationshipsForSubjectId() {
        return (featureRelationshipsForSubjectId = CollectionUtils.safeGetter(featureRelationshipsForSubjectId));
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#setFeatureRelationshipsForSubjectId(java.util.Set)
     */
    public void setFeatureRelationshipsForSubjectId(Collection<FeatureRelationship> featureRelationshipsForSubjectId) {
        this.featureRelationshipsForSubjectId = featureRelationshipsForSubjectId;
    }
    
    public void addFeatureRelationshipsForSubjectId(FeatureRelationship featureRelationshipForSubjectId) {
      featureRelationshipForSubjectId.setFeatureBySubjectId(this);
      this.featureRelationshipsForSubjectId.add(featureRelationshipForSubjectId);
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#getFeatureDbxrefs()
     */
    public Collection<FeatureDbXRef> getFeatureDbXRefs() {
        return this.featureDbXRefs;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#setFeatureDbxrefs(java.util.Set)
     */
    public void setFeatureDbXRefs(Collection<FeatureDbXRef> featureDbXRefs) {
        this.featureDbXRefs = featureDbXRefs;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#getFeaturelocsForFeatureId()
     */
    public Collection<FeatureLoc> getFeatureLocsForFeatureId() {
        return (featureLocsForFeatureId = CollectionUtils.safeGetter(featureLocsForFeatureId));
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#setFeaturelocsForFeatureId(java.util.Set)
     */
    public void addFeatureLocsForFeatureId(FeatureLoc featureLocForFeatureId) {
        featureLocForFeatureId.setFeatureByFeatureId(this);
        getFeatureLocsForFeatureId().add(featureLocForFeatureId);
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#getFeatureCvterms()
     */
    public Collection<FeatureCvTerm> getFeatureCvTerms() {
        return this.featureCvTerms;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#setFeatureCvterms(java.util.Set)
     */
    public void setFeatureCvTerms(Collection<FeatureCvTerm> featureCvTerms) {
        this.featureCvTerms = featureCvTerms;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#getFeatureProps()
     */
    public Collection<FeatureProp> getFeatureProps() {
        return (featureProps = CollectionUtils.safeGetter(featureProps));
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#setFeatureProps(java.util.Set)
     */
    public void addFeatureProp(FeatureProp featureProp) {
        featureProp.setFeature(this);
        getFeatureProps().add(featureProp);
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#getFeaturePubs()
     */
    public Collection<FeaturePub> getFeaturePubs() {
        return this.featurePubs;
    }
    
    public void setFeaturePubs(Collection<FeaturePub> featurePubs) {
        this.featurePubs = featurePubs;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#getAnalysisfeatures()
     */
    public Collection<AnalysisFeature> getAnalysisFeatures() {
        return this.analysisFeatures;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#setAnalysisfeatures(java.util.Set)
     */
    public void setAnalysisFeatures(Collection<AnalysisFeature> analysisFeatures) {
        this.analysisFeatures = analysisFeatures;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#getFeatureSynonyms()
     */
    public Collection<FeatureSynonym> getFeatureSynonyms() {
        return (featureSynonyms = CollectionUtils.safeGetter(featureSynonyms));
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.FeatureI#setFeatureSynonyms(java.util.Set)
     */
    public void setFeatureSynonyms(Collection<FeatureSynonym> featureSynonyms) {
        this.featureSynonyms = featureSynonyms;
    }

    /**
     * Get the display name for the gene, preferrably the name, 
     * otherwise the display name   
     * 
     * @return the preferred display name, never null 
     */
    public String getDisplayName() {
        return (getName() != null) ? getName() : getUniqueName(); 
    }


    public void setFeatureLocsForFeatureId(Collection<FeatureLoc> featureLocsForFeatureId) {
        this.featureLocsForFeatureId = featureLocsForFeatureId;
    }


    public void setFeatureProps(Collection<FeatureProp> featureProps) {
        this.featureProps = featureProps;
    }
    

    public Collection<Phylonode> getPhylonodes() {
        return this.phylonodes;
    }

    public void setPhylonodes(Collection<Phylonode> phylonodes) {
        this.phylonodes = phylonodes;
    }

    private String calcMD5(byte[] residue) {
        try {
            MessageDigest md5 = MessageDigest.getInstance("MD5");
//          MessageDigest tc1 = md.clone();
            byte[] md5Bytes = md5.digest(residue);

            StringBuilder hexValue = new StringBuilder();
            for (int i=0 ; i<md5Bytes.length ; i++) {
                int val = md5Bytes[i] & 0xff; 
                if (val < 16) hexValue.append("0");
                hexValue.append(Integer.toHexString(val));
            }
            return hexValue.toString();

        }
        catch (NoSuchAlgorithmException exp) {
            exp.printStackTrace(); // Shouldn't happen - MD5 is supported algorithm
        }
        return null;
    }
    
    public FeatureLoc getFeatureLoc()
    {
      return featureLoc;
    }

    public void setFeatureLoc(FeatureLoc featureLoc)
    {
      this.featureLoc = featureLoc;
    }

}
