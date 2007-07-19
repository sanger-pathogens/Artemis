/* Qualifier.java
 *
 * created: Feb 2007
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2007  Genome Research Limited
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

package uk.ac.sanger.artemis.io;

import java.util.Iterator;
import java.util.List;
import java.util.Vector;

import uk.ac.sanger.artemis.util.*;

/**
 *  Each object of this class represents a group of embl qualifiers
 *  (name/value pairs) with the same name, so in effect the object contains
 *  one name and many values.
 *
 *  @author Kim Rutherford
 *  @version $Id: QualifierLazyLoading.java,v 1.5 2007-07-19 09:49:32 tjc Exp $
 * */

public class QualifierLazyLoading extends Qualifier
{

  /** The name that was passed to the constructor. */
  private String name;

  /** The values that were passed to the constructor. */
  private List values;
  private boolean forceLoad = false;
  
  /**
   *  Create a new Qualifier object.  This object consists of a name and
   *  possibly some values.  In the raw embl file we have /name=value.
   *  @param name The name of this qualifier (ie. the text immediately after
   *    the / in the qualifier)
   * @param values The values of this qualifier (ie the text immediately after
   *    the = in the qualifiers).  This argument should be null if the
   *    qualifier has no values.  The values String objects should not include
   *    the quote characters (),"" or [].  For example if the original
   *    qualifier was /citation=[3] then the value String should be: 3.  This
   *    argument is copied by the constructor.
   **/
  public QualifierLazyLoading (final String name, final List set) 
  {
    initialiseLazy (name, set);
  }

  /**
   *  Create a new Qualifier object.  This object consists of a name and
   *  possibly some values.  In the raw embl file we have /name=value.
   *  @param name The name of this qualifier (ie. the text immediately after
   *    the / in the qualifier)
   *  @param value The value of this qualifier (see the other constructor for
   *    details).  Other values may be added later.
   **/
  public QualifierLazyLoading (final String name, final LazyQualifierValue value) 
  {
    List set = new Vector();
    set.add(value);
    initialiseLazy (name, set);
  }

  /**
   *  Called by the constructor to initialise the object.
   *  @param name The name of this qualifier (ie. the text immediately after
   *    the / in the qualifier)
   *  @param values The values of this qualifier (see the other constructor for
   *    details).  Other values may be added later.
   **/
  private void initialiseLazy(final String name, final List values) 
  {
    this.name = name;
    if(values == null) 
      this.values = null;
    else
    {
      if(values.size() == 0)
        throw new Error ("internal error - zero length values vector");

      this.values = values;
    }
  }
  
  /**
   *  Return the values of this qualifier.  The StringVector that is returned
   *  is a copy (use addValues (), removeValue () etc. to change the
   *  Qualifier), but the String objects inside the vector are not copied so
   *  the same String references will occur in all copies.
   **/
  public StringVector getValues() 
  {
    if(values == null)
      return null;
    else
    {
      StringVector v = new StringVector();
      
      for(int i=0; i<values.size(); i++)
      {
        LazyQualifierValue lazy = (LazyQualifierValue)values.get(i);
        if(forceLoad)
          lazy.setForceLoad(true);
        v.add(i, lazy.getString());
      }
      
      return v;
    }
  }
  
  /**
   * 
   * @param index
   * @return
   */
  public LazyQualifierValue getValue(int index)
  {
    return (LazyQualifierValue)values.get(index);
  }
  
  /**
   * Test if all lazy values are loaded
   * @return
   */
  public boolean isAllLazyValuesLoaded()
  {
    for(int i=0; i<values.size(); i++)
      if(!((LazyQualifierValue)values.get(i)).isLazyLoaded())
        return false;
    return true;
  }

  /**
   * 
   * @param index
   * @return
   */
  public List getLazyValues()
  {
    return values;
  }
  
  /**
   *  Add the given values to this object.
   **/
  public void addValues(final List new_values) 
  {
    if(values == null) 
    {
      values = new Vector();
      values.add ((String)null);
    }
    if(new_values == null) 
    {
      // if new_values is null then we have a qualifier with no value (like
      // /pseudo).  we add one null value for each occurrence of the
      // qualifier.
      values.add ((String)null);
    }
    else
    {
      Iterator it = new_values.iterator();
      while(it.hasNext()) 
        values.add(it.next());
    }
  }

  /**
   *  Add the given value to this object.
   **/
  public void addValue(final Object new_value)
  {
    if(values == null) 
      values = new Vector();
    
    values.add(new_value);
  }
  
  public void removeValue(final Object value)
  {
    values.remove(value);

    if(values.size() == 0) 
      values = null;
  }
  
  public String getName()
  {
    return name;
  }
  
  /**
   *  Return the reference of a new copy of this Qualifier.
   **/
  public Qualifier copy() 
  {
    return new QualifierLazyLoading(getName(), values);
  }

  public void setForceLoad(boolean forceLoad)
  {
    this.forceLoad = forceLoad;
  }
}