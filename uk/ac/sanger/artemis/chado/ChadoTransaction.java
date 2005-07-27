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

import java.util.Properties;
import java.util.Hashtable;
import java.util.Vector;
import java.util.Enumeration;

/**
*
* Store information about a transaction.
*
**/
public class ChadoTransaction
{
  /** update statement */
  public static final int UPDATE = 1;
  /** insert statement */
  public static final int INSERT = 2;
  /** delete statement */
  public static final int DELETE = 3;

  /** properties store */
  protected Hashtable properties; 
  /** old properties store */
  protected Hashtable constraint;

  /** type of statement */
  protected int type;
  /** feature uniquename */
  protected String uniquename;
  /** chado table for this transaction */
  protected String chadoTable;

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
      properties = new Hashtable();
    properties.put(name, value);
  }

  public void setConstraint(String name, String value)
  {
    if(constraint == null)
      constraint = new Hashtable();
    constraint.put(name, value);
  }

  public String getSqlQuery()
  {
    StringBuffer sqlBuff = new StringBuffer();

    if(type == UPDATE)
    {
      sqlBuff.append("UPDATE "+chadoTable);
      sqlBuff.append(" SET ");

      String name; 
      String value;
      Enumeration enum_prop = properties.keys();
      while(enum_prop.hasMoreElements())
      {
        name  = (String)enum_prop.nextElement();
        value = (String)properties.get(name);
        sqlBuff.append(name+"="+value);
        if(enum_prop.hasMoreElements())
          sqlBuff.append(" , ");
      }
      sqlBuff.append(" WHERE feature.feature_id="+chadoTable+".feature_id AND feature.uniquename='");
      sqlBuff.append(uniquename+"'");

      if(constraint != null)
      {
        Enumeration enum_constraint = constraint.keys();
        while(enum_constraint.hasMoreElements())
        {
          name  = (String)enum_constraint.nextElement();
          value = (String)constraint.get(name);
          sqlBuff.append(" AND "+name+"="+value);
        }
      }
    }
//  else if(type == INSERT)
//  else
           
    return sqlBuff.toString();
  }
}
