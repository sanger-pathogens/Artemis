package org.gmod.schema.analysis;


import java.io.Serializable;
import java.util.Collection;
import java.util.Date;
import java.util.HashSet;
import java.util.Set;










public class Analysis implements Serializable {

    // Fields    


     private int analysisId;
     

     private String name;
     

     private String description;
     

     private String program;
     

     private String programVersion;
     

     private String algorithm;
     

     private String sourceName;
     

     private String sourceVersion;
     

     private String sourceUri;
     

     private Date timeExecuted;
     

     private Set<AnalysisFeature> analysisFeatures = new HashSet<AnalysisFeature>(0);
     

     private Set<AnalysisProp> analysisProps = new HashSet<AnalysisProp>(0);
    
   
    // Property accessors
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.AnalysisI#getAnalysisId()
     */
     public int getAnalysisId() {
        return this.analysisId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.AnalysisI#setAnalysisId(int)
     */
    public void setAnalysisId(int analysisId) {
        this.analysisId = analysisId;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.AnalysisI#getName()
     */
    public String getName() {
        return this.name;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.AnalysisI#setName(java.lang.String)
     */
    public void setName(String name) {
        this.name = name;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.AnalysisI#getDescription()
     */
    public String getDescription() {
        return this.description;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.AnalysisI#setDescription(java.lang.String)
     */
    public void setDescription(String description) {
        this.description = description;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.AnalysisI#getProgram()
     */
    public String getProgram() {
        return this.program;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.AnalysisI#setProgram(java.lang.String)
     */
    public void setProgram(String program) {
        this.program = program;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.AnalysisI#getProgramversion()
     */
    public String getProgramVersion() {
        return this.programVersion;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.AnalysisI#setProgramversion(java.lang.String)
     */
    public void setProgramVersion(String programVersion) {
        this.programVersion = programVersion;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.AnalysisI#getAlgorithm()
     */
    public String getAlgorithm() {
        return this.algorithm;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.AnalysisI#setAlgorithm(java.lang.String)
     */
    public void setAlgorithm(String algorithm) {
        this.algorithm = algorithm;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.AnalysisI#getSourcename()
     */
    public String getSourceName() {
        return this.sourceName;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.AnalysisI#setSourcename(java.lang.String)
     */
    public void setSourceName(String sourceName) {
        this.sourceName = sourceName;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.AnalysisI#getSourceversion()
     */
    public String getSourceVersion() {
        return this.sourceVersion;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.AnalysisI#setSourceversion(java.lang.String)
     */
    public void setSourceVersion(String sourceVersion) {
        this.sourceVersion = sourceVersion;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.AnalysisI#getSourceuri()
     */
    public String getSourceUri() {
        return this.sourceUri;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.AnalysisI#setSourceuri(java.lang.String)
     */
    public void setSourceUri(String sourceUri) {
        this.sourceUri = sourceUri;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.AnalysisI#getTimeexecuted()
     */
    public Date getTimeExecuted() {
        return this.timeExecuted;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.AnalysisI#setTimeexecuted(java.util.Date)
     */
    public void setTimeExecuted(Date timeExecuted) {
        this.timeExecuted = timeExecuted;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.AnalysisI#getAnalysisfeatures()
     */
    public Collection<AnalysisFeature> getAnalysisFeatures() {
        return this.analysisFeatures;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.AnalysisI#setAnalysisfeatures(java.util.Set)
     */
    public void setAnalysisFeatures(Set<AnalysisFeature> analysisFeatures) {
        this.analysisFeatures = analysisFeatures;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.AnalysisI#getAnalysisprops()
     */
    public Collection<AnalysisProp> getAnalysisProps() {
        return this.analysisProps;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.AnalysisI#setAnalysisprops(java.util.Set)
     */
    private void setAnalysisProps(Set<AnalysisProp> analysisProps) {
        this.analysisProps = analysisProps;
    }




}


