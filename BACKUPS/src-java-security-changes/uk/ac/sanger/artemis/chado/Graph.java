/* Graph.java
 *
 * created: 2009
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2009  Genome Research Limited
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

import java.io.Serializable;

/**
 * Database graph.
 */
public class Graph implements Serializable 
{
  private static final long serialVersionUID = 1L;
  // table columns
  private int graphId;
  private int featureId;
  private String name;
  private String description;
  private byte[] data;

  public Graph() 
  {
  }

  public int getGraphId()
  {
    return graphId;
  }

  public void setGraphId(int graphId)
  {
    this.graphId = graphId;
  }

  public int getFeatureId()
  {
    return featureId;
  }

  public void setFeatureId(int featureId)
  {
    this.featureId = featureId;
  }

  public String getName()
  {
    return name;
  }

  public void setName(String name)
  {
    this.name = name;
  }

  public String getDescription()
  {
    return description;
  }

  public void setDescription(String description)
  {
    this.description = description;
  }

  public byte[] getData()
  {
    return data;
  }

  public void setData(byte[] data)
  {
    this.data = data;
  }
}