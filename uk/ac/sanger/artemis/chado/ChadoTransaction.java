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
* Store information about a transaction.
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

  /** properties store */
  private List properties;
  /** old properties store */
  protected List constraint;
  /** type of statement */
  protected int type;
  /** feature uniquename */
  protected String uniquename;
  /** chado table for this transaction */
  protected String chadoTable;
  /** schema */
  private String schema;
  /** feature id */
  private int feature_id;
  /** chado feature */
  private ChadoFeature chadoFeature;

  /**
  *
  * @param type		transaction type
  * @param uniquename	uniquename of feature
  * @param chadoTable	chado table used in transaction
  *
  */
  public ChadoTransaction(int type, String uniquename, 
                          String chadoTable)
  {
    this.type = type;
    this.uniquename = uniquename;
    this.chadoTable = chadoTable;
  }
 
  /**
  *
  * @param type          transaction type
  * @param chado_feature ChadoFeature to be inserted
  *
  */
  public ChadoTransaction(int type, ChadoFeature chadoFeature)
  {
    this.chadoFeature = chadoFeature;
    this.type = type;
  }


  /**
  *
  *
  */
  public ChadoFeature getChadoFeature()
  {
    return chadoFeature;
  }

  /**
  *
  *
  */
  public void setSchema(final String schema)
  {
    this.schema = schema;
  }

  /**
  *
  *
  */
  public String getSchema()
  {
    return schema;
  }

  /**
  *
  * @return transaction type
  *
  */
  public int getType()
  {
    return type;
  } 

  /**
  *
  * @return uniquename of feature
  *
  */
  public String getUniqueName()
  {
    return uniquename;
  }


  /**
  *
  * @return chado table
  *
  */
  public String getChadoTable()
  {
    return chadoTable;
  }


  /**
  * 
  *  Add a property to this transaction that will be changed.
  *
  */
  public void addProperty(String name, String value) 
  {
    if(properties == null)
      properties = new Vector();

    properties.add(name+"="+value);
  }


  /**
  *
  *  Set a constraint on this transaction.
  *
  */
  public void setConstraint(String name, String value)
  {
    if(constraint == null)
      constraint = new Vector();
    constraint.add(name+"="+value);
  }


  /**
  *
  *  Get properties to this transaction that will be changed.
  *
  */
  public List getProperties()
  {
    return properties;
  }
 
  public void setProperties(List properties)
  {
    this.properties = properties;
  }

  /**
  *
  *  Get properties to this transaction that will be changed.
  *
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
  *
  *  Get properties to this transaction that will be changed.
  *
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
  *
  *  Get constraints on this transaction.
  *
  */
  public List getConstraint()
  {
    return constraint;
  }

  /**
  *
  *  Get uniquenames of features.
  *
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
  *
  * Set the feature id
  *
  */
  public void setFeature_id(final int feature_id)
  {
    this.feature_id = feature_id;
  }

  /**
  *
  * Get the feature id
  *
  */
  public int getFeature_id()
  {
    return feature_id;
  }

}
