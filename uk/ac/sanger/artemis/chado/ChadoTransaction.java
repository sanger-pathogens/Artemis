/* ChadoTransaction
 *
 * created: July 2006
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2006  Genome Research Limited
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 */

package uk.ac.sanger.artemis.chado;

import java.sql.Timestamp;

import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.util.DatabaseDocument;


/**
 * Store information about a SQL transaction <i>e.g</i> UPDATE, INSERT, DELETE
 * a feature, featureloc, featureprop, feature_dbxref, feature_synonym. 
 **/
public class ChadoTransaction
{
  //
  // FEATURE TRANSACTIONS
  /** update statement */
  public static final int UPDATE = 1;
  /** insert statement */
  public static final int INSERT = 2;
  /** delete statement */
  public static final int DELETE = 3;
  
  /** type of statement <i>e.g.</i> UPDATE, INSERT, DELETE, ... */
  protected int type;
  /** feature unique name */
  protected String old_uniquename;
  private String uniquename;

  /** last time feature was modified */
  private Timestamp lastmodified;
  /** the feature object */
  private Object feature_obj;
  /** the feature key */
  private String featureKey;
  
  private GFFStreamFeature gff_feature;
  
  private String logComment;
  
  /**
   * Create a transaction to represent a single insert, delete or
   * update.
   * @param type          the type of transcation, e.g. insert, delete, update.
   * @param feature_obj   the <code>Feature</code>, <code>FeatureLoc</code>,
   * <code>FeatureProp</code>, <code>FeatureRelationship</code>, <code>FeatureDbxref</code>,
   * or <code>FeatureSynonym</code>.
   * @param lastmodified  the last modified timestamp  
   * @param gff_feature   the artemis GFF feature effected
   */
  public ChadoTransaction(final int type, 
                          final Object feature_obj,
                          final Timestamp lastmodified,
                          final GFFStreamFeature gff_feature,
                          final String featureKey, 
                          final String logComment)
  {
    this.type = type;
    this.lastmodified = lastmodified;
    this.feature_obj  = feature_obj;
    this.gff_feature  = gff_feature;
    this.logComment   = logComment;
    
    if(featureKey != null &&
       featureKey.equals(DatabaseDocument.EXONMODEL))
      this.featureKey = "exon";
    else
      this.featureKey   = featureKey;
  }

  /**
   * Copy this transaction
   * @return
   */
  public ChadoTransaction copy()
  {
    final ChadoTransaction tsn = new ChadoTransaction(getType(), getFeatureObject(), 
        getLastModified(), getGff_feature(), 
        getFeatureKey(), logComment);
    
    if(uniquename != null)
      tsn.setUniquename(uniquename);
    
    tsn.setOldUniquename(old_uniquename);
    return tsn;
  }

  /**
   * The type of SQL transaction
   * @return 	the transaction type
   */
  public int getType()
  {
    return type;
  } 

  public String getTypeAsString()
  {
    if(type == UPDATE)
      return "UPDATE";
    else if(type == INSERT)
      return "INSERT";
    else if(type == DELETE)
      return "DELETE";
    return "";
  }

  /**
   * Set the old uniquename, used when updating the uniquename
   * @param uniquename
   */
  public void setOldUniquename(final String old_uniquename)
  {
    this.old_uniquename = old_uniquename;
  }

  /**
   * Get the unique names of features to change.
   */
  public String getOldUniquename()
  {
    return old_uniquename;
  }
  
  /**
   * Set a uniquename, e.g. to be used when changing a featureprop
   * @param uniquename
   */
  public void setUniquename(final String uniquename)
  {
    this.uniquename = uniquename;
  }
  
  public String getUniquename()
  {
    if(uniquename != null)
      return uniquename;
    
    if(getGff_feature() == null)
      return null;
    
    return (String)
      getGff_feature().getQualifierByName("ID").getValues().get(0);
  }

  /**
   * Get the last time modified time stamp.
   * @return  the <code>Timestamp</code>
   */
  public Timestamp getLastModified()
  {
    return this.lastmodified;
  }
  
  public Object getFeatureObject()
  {
    return feature_obj;
  }

  public GFFStreamFeature getGff_feature()
  {
    return gff_feature;
  }


  public String getFeatureKey()
  {
    return featureKey;
  }


  public String getLogComment()
  {
    String key = "";
    if(getFeatureKey() != null)
      key = " KEY=" + getFeatureKey();
    return "["+getTypeAsString()+"] "+logComment+key;
  }

}
