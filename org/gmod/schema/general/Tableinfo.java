package org.gmod.schema.general;


import java.io.Serializable;
import java.util.Date;

public class Tableinfo implements Serializable {

    // Fields    

     private int tableinfoId;
     private String name;
     private String primaryKeyColumn;
     private int isView;
     private Integer viewOnTableId;
     private Integer superclassTableId;
     private int isUpdateable;
     private Date modificationDate;
    
   
    // Property accessors
     /* (non-Javadoc)
     * @see org.genedb.db.jpa.TableInfoI#getTableinfoId()
     */

    

    private int getTableinfoId() {
        return this.tableinfoId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.TableInfoI#setTableinfoId(int)
     */
    private void setTableinfoId(int tableinfoId) {
        this.tableinfoId = tableinfoId;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.TableInfoI#getName()
     */
    private String getName() {
        return this.name;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.TableInfoI#setName(java.lang.String)
     */
    private void setName(String name) {
        this.name = name;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.TableInfoI#getPrimaryKeyColumn()
     */
    private String getPrimaryKeyColumn() {
        return this.primaryKeyColumn;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.TableInfoI#setPrimaryKeyColumn(java.lang.String)
     */
    private void setPrimaryKeyColumn(String primaryKeyColumn) {
        this.primaryKeyColumn = primaryKeyColumn;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.TableInfoI#getIsView()
     */
    private int getIsView() {
        return this.isView;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.TableInfoI#setIsView(int)
     */
    private void setIsView(int isView) {
        this.isView = isView;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.TableInfoI#getViewOnTableId()
     */
    private Integer getViewOnTableId() {
        return this.viewOnTableId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.TableInfoI#setViewOnTableId(java.lang.Integer)
     */
    private void setViewOnTableId(Integer viewOnTableId) {
        this.viewOnTableId = viewOnTableId;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.TableInfoI#getSuperclassTableId()
     */
    private Integer getSuperclassTableId() {
        return this.superclassTableId;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.TableInfoI#setSuperclassTableId(java.lang.Integer)
     */
    private void setSuperclassTableId(Integer superclassTableId) {
        this.superclassTableId = superclassTableId;
    }
    

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.TableInfoI#getIsUpdateable()
     */
    private int getIsUpdateable() {
        return this.isUpdateable;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.TableInfoI#setIsUpdateable(int)
     */
    private void setIsUpdateable(int isUpdateable) {
        this.isUpdateable = isUpdateable;
    }

    /* (non-Javadoc)
     * @see org.genedb.db.jpa.TableInfoI#getModificationDate()
     */
    private Date getModificationDate() {
        return this.modificationDate;
    }
    
    /* (non-Javadoc)
     * @see org.genedb.db.jpa.TableInfoI#setModificationDate(java.util.Date)
     */
    private void setModificationDate(Date modificationDate) {
        this.modificationDate = modificationDate;
    }




}


