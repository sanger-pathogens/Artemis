package org.gmod.schema.pub;



import java.io.Serializable;











public class PubAuthor implements Serializable {

    // Fields    

    

     private int pubAuthorId;
     


     private Pub pub;
     

     private int rank;
     

     private Boolean editor;
     

     private String surname;
     

     private String givenNames;
     

     private String suffix;

   
    // Property accessors

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubAuthorI#getPubAuthorId()
     */
    private int getPubAuthorId() {
        return this.pubAuthorId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubAuthorI#setPubAuthorId(int)
     */
    private void setPubAuthorId(int pubAuthorId) {
        this.pubAuthorId = pubAuthorId;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubAuthorI#getPub()
     */
    private Pub getPub() {
        return this.pub;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubAuthorI#setPub(org.gmod.schema.pub.PubI)
     */
    private void setPub(Pub pub) {
        this.pub = pub;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubAuthorI#getRank()
     */
    private int getRank() {
        return this.rank;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubAuthorI#setRank(int)
     */
    private void setRank(int rank) {
        this.rank = rank;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubAuthorI#getEditor()
     */
    private Boolean getEditor() {
        return this.editor;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubAuthorI#setEditor(java.lang.Boolean)
     */
    private void setEditor(Boolean editor) {
        this.editor = editor;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubAuthorI#getSurname()
     */
    private String getSurname() {
        return this.surname;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubAuthorI#setSurname(java.lang.String)
     */
    private void setSurname(String surname) {
        this.surname = surname;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubAuthorI#getGivenNames()
     */
    private String getGivenNames() {
        return this.givenNames;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubAuthorI#setGivenNames(java.lang.String)
     */
    private void setGivenNames(String givenNames) {
        this.givenNames = givenNames;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubAuthorI#getSuffix()
     */
    private String getSuffix() {
        return this.suffix;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.PubAuthorI#setSuffix(java.lang.String)
     */
    private void setSuffix(String suffix) {
        this.suffix = suffix;
    }




}


