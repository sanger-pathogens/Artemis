/* CacheHashMap.java
 *
 * created: 2012
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2012  Genome Research Limited
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or(at your option) any later version.
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
 **/
package uk.ac.sanger.artemis.util;

import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;

/**
 * Cache using a HashMap.
 */
public class CacheHashMap extends HashMap<Object, Object> 
{
  private static final long serialVersionUID = 1L;
  private final int maxSize;
  private final int shrinkSize;
  private final LinkedList<Object> accessList = new LinkedList<Object>();

  /**
   * When the maximum size is exceeded key/value entries are removed
   * from the map.
   * @param maxSize     maximum size of the map
   * @param shrinkSize  number to shrink when maxSize is reached
   */
  public CacheHashMap(int maxSize, int shrinkSize) 
  {
    super(maxSize);
    this.maxSize = maxSize;
    this.shrinkSize = shrinkSize;
  }

  public Object put(Object key, Object val)
  {
    if(size() >= maxSize)
      shrink();
    updateAccessList(key);
    return super.put(key, val);
  }

  public Object get(Object key) 
  {
    updateAccessList(key);
    return super.get(key);
  }

  public Object getLastKey()
  {
    return accessList.getLast();
  }
  
  private void shrink() 
  {
    final Iterator<Object> it = accessList.iterator();
    for(int i = 0; i < shrinkSize; i++) 
    {
      if(!it.hasNext())
        return;
      Object key = it.next();
      this.remove(key);
      it.remove();
    }
  }
  
  private void updateAccessList(Object key) 
  {
    int idx = accessList.indexOf(key);
    if (idx >= 0)
      accessList.remove(idx);
    accessList.add(key);
  }
}