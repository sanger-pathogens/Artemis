/* FeatureForUpdatingResidues
 *
 * created: July 2007
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

package uk.ac.sanger.artemis.chado;

import org.gmod.schema.sequence.Feature;

public class FeatureForUpdatingResidues extends Feature
{
  private int startBase;
  private int endBase;
  private int basesToEnd;
  private int length;
  private String newSubSequence;
  private boolean isResidueUpdate = false;
  
  public int getEndBase()
  {
    return endBase;
  }
  public void setEndBase(int endBase)
  {
    this.endBase = endBase;
  }
  public String getNewSubSequence()
  {
    return newSubSequence;
  }
  public void setNewSubSequence(String newSubSequence)
  {
    this.newSubSequence = newSubSequence;
  }
  public int getStartBase()
  {
    return startBase;
  }
  public void setStartBase(int startBase)
  {
    this.startBase = startBase;
  }
  public int getBasesToEnd()
  {
    return basesToEnd;
  }
  public void setBasesToEnd(int basesToEnd)
  {
    this.basesToEnd = basesToEnd;
  }
  public int getLength()
  {
    return length;
  }
  public void setLength(int length)
  {
    this.length = length;
  }
  public boolean isResidueUpdate()
  {
    return isResidueUpdate;
  }
  public void setResidueUpdate(boolean isResidueUpdate)
  {
    this.isResidueUpdate = isResidueUpdate;
  }
  
}