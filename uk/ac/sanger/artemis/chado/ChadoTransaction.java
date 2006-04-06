/* ChadoTransaction
 *
 * created: July 2005
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2005  Genome Research Limited
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

import java.util.List;
import java.util.Vector;
import java.util.StringTokenizer;

/**
*
* Store information about a SQL transaction <i>e.g</i> UPDATE, INSERT, DELETE
* a feature or a feature attribute. 
*
**/
public class ChadoTransaction
{
  //
  // ATTRIBUTE TRANSACTIONS
  /** update statement */
  public static final int UPDATE = 1;
  /** insert statement */
  public static final int INSERT = 2;
  /** delete statement */
  public static final int DELETE = 3;

  // FEATURE TRANSACTIONS
  /** insert feature statement */
  public static final int INSERT_FEATURE = 4;
  /** delete feature statement */
  public static final int DELETE_FEATURE = 5;
  
  // DBXREF TRANSACTIONS
  public static final int UPDATE_DBXREF = 6;
  public static final int INSERT_DBXREF = 7;
  public static final int DELETE_DBXREF = 8;
  
  /** properties store <i>e.g.</i> value=<quote>product=hypothetical protein</quote> */
  private List properties;
  /** constriant store <i>e.g.</i> type_id=21078 */
  protected List constraint;
  /** type of statement <i>e.g.</i> UPDATE, INSERT, DELETE, ... */
  protected int type;
  /** feature unique name */
  protected String uniquename;
  /** chado table for this transaction */
  protected String chadoTable;
  /** postgres schema */
  private String schema;
  /** feature id */
  private int feature_id;
  /** chado feature */
  private ChadoFeature chadoFeature;

  private Dbxref dbxref;
 

  /**
   * Used to construct a <code>ChadoTransaction</code> that can
   * be used to describe a SQL transaction.
   *
   * <i>e.g.</i> to DELETE a feature and where uniquename is a <code>String</code>:
   * <blockquote><pre>
   * ChadoTransaction tsn = new ChadoTransaction(ChadoTransaction.DELETE_FEATURE,
   *                                             uniquename, "feature");
   * </pre></blockquote>
   *
   * Or to INSERT a feature property:
   * <blockquote><pre>
   * ChadoTransaction tsn = new ChadoTransaction(ChadoTransaction.INSERT,
   *                                             feature_id, "featureprop");
   * tsn.addProperty("value", "'"+ qualifier_string +"'");
   * tsn.addProperty("type_id", "'"+cvterm_id+"'");
   * tsn.addProperty("rank", Integer.toString(value_index));
   * </pre></blockquote>
   * 
   * @param type	the transaction type
   * @param uniquename	the unique name of feature
   * @param chadoTable	the chado table used in transaction <i>e.g.</i> 
   *                 	feature, featureloc....
   */
  public ChadoTransaction(int type, String uniquename, 
                          String chadoTable)
  {
    this.type = type;
    this.uniquename = uniquename;
    this.chadoTable = chadoTable;
  }
 
  /**
   * Used to construct a <code>ChadoTransaction</code> that can
   * be used to describe a SQL transaction.
   *
   * <i>e.g.</i> to INSERT a feature and wher chado_feature is a 
   *        <code>ChadoFeature</code>
   * <blockquote><pre>
   * ChadoTransaction tsn = new ChadoTransaction(ChadoTransaction.INSERT_FEATURE,
   *                                             chado_feature);
   * </pre></blockquote>
   *
   * @param type          the transaction type
   * @param chado_feature the <code>ChadoFeature</code> to be inserted
   *
   */
  public ChadoTransaction(int type, ChadoFeature chadoFeature)
  {
    this.chadoFeature = chadoFeature;
    this.type = type;
  }

  /**
   * Used to construct a <code>ChadoTransaction</code> that can
   * be used to describe a SQL transaction.
   * @param type    the transaction type
   * @param dbxref  a feature dbxref
   */
  public ChadoTransaction(final int type, final String uniquename,
                          final Dbxref dbxref)
  {
    this.type = type;
    this.dbxref = dbxref;
    this.uniquename = uniquename;
  }

  /**
   * The <code>Dbxref</code> used in a transaction.
   * @return  the feature dbxref
   */
  public Dbxref getFeatureDbxref()
  {
    return dbxref;
  }
  
  /**
   * The <code>ChadoFeature</code> feature used in a transaction.
   * @return chado_feature the <code>ChadoFeature</code> 
   */
  public ChadoFeature getChadoFeature()
  {
    return chadoFeature;
  }

  /**
   * Set the postgres schema.
   * @param     the postgres schema
   */
  public void setSchema(final String schema)
  {
    this.schema = schema;
  }

  /**
   * Get the postgres schema.
   * @return    the postgres schema
   */
  public String getSchema()
  {
    return schema;
  }

  /**
   * The type of SQL transaction
   * @return 	the transaction type
   */
  public int getType()
  {
    return type;
  } 

  /**
   * Get the unique name for a feature.
   * @return    the unique name for a feature.
   */
  public String getUniqueName()
  {
    return uniquename;
  }


  /**
   * Get the chado table that is involved in the SQL command
   * @return 	the chado table
   */
  public String getChadoTable()
  {
    return chadoTable;
  }


  /**
   * Add a property to this transaction that will be updated
   * or inserted.
   * @param name	the property name
   * @param value	the property value
   */
  public void addProperty(String name, String value) 
  {
    if(properties == null)
      properties = new Vector();

    properties.add(name+"="+value);
  }


  /**
   * Set a constraint on this transaction.
   * @param name	the constraint name
   * @param value	the constraint value
   */
  public void setConstraint(String name, String value)
  {
    if(constraint == null)
      constraint = new Vector();
    constraint.add(name+"="+value);
  }


  /**
   * Get the properties for this transaction that will be changed.
   * @return	the transaction properties used to INSERT/UPDATE 
   */
  public List getProperties()
  {
    return properties;
  }
 
  /**
   * Set the properties for this transaction that will be changed.
   * @param properties	the transaction properties used to INSERT/UPDATE
   */
  public void setProperties(List properties)
  {
    this.properties = properties;
  }

  /**
   * Get the names of the properties for this transaction that will 
   * be changed.
   * @return	the names of the properties to INSERT/UPDATE
   */
  public List getPropertiesName()
  {
    List propertiesName = new Vector();
    String property;
    int index;
    for(int i=0; i<properties.size();i++)
    {
      property = (String)properties.get(i);
      index = property.indexOf("=");
      propertiesName.add(property.substring(0,index));
    }

    return propertiesName;
  }

  /**
   * Get the values of the properties for this transaction that will  
   * be changed..
   * @return	the values of the properties to INSERT/UPDATE
   */
  public List getPropertiesValue()
  {
    List propertiesValue = new Vector();
    String property;
    int index;
    for(int i=0; i<properties.size();i++)
    {
      property = (String)properties.get(i);
      index = property.indexOf("=");
      propertiesValue.add(property.substring(index+1));
    }

    return propertiesValue;
  }

  /**
   * Get the constraints on this transaction.
   * @return 	the transaction constraints	
   */
  public List getConstraint()
  {
    return constraint;
  }

  /**
   * Get the unique names of features to change.
   */
  public List getUniquename()
  {
    Vector names = new Vector();
    StringTokenizer tok = new StringTokenizer(uniquename,",");
    while(tok.hasMoreTokens())
    {
      names.add(tok.nextToken());
    }

    return names;
  }

  /**
   * Set the feature id
   * @param feature_id	the feature identifier
   */
  public void setFeature_id(final int feature_id)
  {
    this.feature_id = feature_id;
  }

  /**
   * Get the feature id
   * @return	the feature identifier
   */
  public int getFeature_id()
  {
    return feature_id;
  }

}
