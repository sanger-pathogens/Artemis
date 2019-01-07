package org.gmod.schema.general;


import java.io.Serializable;
import java.util.Collection;
import java.util.HashSet;
import java.util.Set;










public class Db implements Serializable {

    // Fields    

    

     private int dbId;
     

     private String name;
     

     private String description;
     

     private String urlPrefix;


     private String url;
     

     private Set<DbXRef> dbXRefs = new HashSet<DbXRef>(0);
    
   
    // Property accessors
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.DbI#getDbId()
     */
    public int getDbId() {
        return this.dbId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.DbI#setDbId(int)
     */
    public void setDbId(int dbId) {
        this.dbId = dbId;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.DbI#getName()
     */
    public String getName() {
        return this.name;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.DbI#setName(java.lang.String)
     */
    public void setName(String name) {
        this.name = name;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.DbI#getDescription()
     */
    public String getDescription() {
        return this.description;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.DbI#setDescription(java.lang.String)
     */
    public void setDescription(String description) {
        this.description = description;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.DbI#getUrlPrefix()
     */
    public String getUrlPrefix() {
        return this.urlPrefix;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.DbI#setUrlPrefix(java.lang.String)
     */
    public void setUrlPrefix(String urlPrefix) {
        this.urlPrefix = urlPrefix;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.DbI#getUrl()
     */
    public String getUrl() {
        return this.url;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.DbI#setUrl(java.lang.String)
     */
    public void setUrl(String url) {
        this.url = url;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.DbI#getDbXRefs()
     */
    public Collection<DbXRef> getDbXRefs() {
        return this.dbXRefs;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.DbI#setDbXRefs(java.util.Set)
     */
    public void setDbXRefs(Set<DbXRef> dbXRefs) {
        this.dbXRefs = dbXRefs;
    }




}


