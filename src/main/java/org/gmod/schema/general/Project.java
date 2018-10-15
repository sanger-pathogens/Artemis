package org.gmod.schema.general;


import java.io.Serializable;

public class Project implements Serializable {

    // Fields    
     private int projectId;
     private String name;
     private String description;
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.ProjectI#getProjectId()
     */
    private int getProjectId() {
        return this.projectId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.ProjectI#setProjectId(int)
     */
    private void setProjectId(int projectId) {
        this.projectId = projectId;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.ProjectI#getName()
     */
    private String getName() {
        return this.name;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.ProjectI#setName(java.lang.String)
     */
    private void setName(String name) {
        this.name = name;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.ProjectI#getDescription()
     */
    private String getDescription() {
        return this.description;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.ProjectI#setDescription(java.lang.String)
     */
    private void setDescription(String description) {
        this.description = description;
    }

}


