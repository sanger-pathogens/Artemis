package org.gmod.schema.analysis;




import org.gmod.schema.cv.CvTerm;
import org.gmod.schema.utils.propinterface.PropertyI;

import java.io.Serializable;











public class AnalysisProp implements Serializable, PropertyI {

    // Fields    



     private int analysisPropId;
    


     private Analysis analysis;
    


     private CvTerm cvTerm;
    

     private String value;
    
   
    // Property accessors

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.AnalysisPropI#getAnalysispropId()
     */
    private int getAnalysisPropId() {
        return this.analysisPropId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.AnalysisPropI#setAnalysispropId(int)
     */
    private void setAnalysisPropId(int analysisPropId) {
        this.analysisPropId = analysisPropId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.AnalysisPropI#getAnalysis()
     */
    private Analysis getAnalysis() {
        return this.analysis;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.AnalysisPropI#setAnalysis(org.genedb.db.jpa.Analysis)
     */
    private void setAnalysis(Analysis analysis) {
        this.analysis = analysis;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.AnalysisPropI#getCvterm()
     */
    public CvTerm getCvTerm() {
        return this.cvTerm;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.AnalysisPropI#setCvterm(org.gmod.schema.cv.CvTermI)
     */
    private void setCvTerm(CvTerm cvTerm) {
        this.cvTerm = cvTerm;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.AnalysisPropI#getValue()
     */
    private String getValue() {
        return this.value;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.AnalysisPropI#setValue(java.lang.String)
     */
    private void setValue(String value) {
        this.value = value;
    }




}


