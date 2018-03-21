/*
 * Copyright (C) 2008  Genome Research Limited
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
 *  @author: Tim Carver
 */

package uk.ac.sanger.artemis.circular;

import uk.ac.sanger.artemis.io.InvalidRelationException;
import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.io.Location;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.util.OutOfRangeException;

public class Feature
{
  private int colour;
  private int start;
  private int end;
  private String emblKey;
  
  public Feature(final String emblKey, 
                 final int start, 
                 final int end, 
                 final int colour)
  {
    this.emblKey = emblKey;
    this.start = start;
    this.end = end;
    this.colour = colour;
  }

  public int getEnd()
  {
    return end;
  }

  public void setEnd(int end)
  {
    this.end = end;
  }

  public int getStart()
  {
    return start;
  }

  public void setStart(int start)
  {
    this.start = start;
  }
  
  public String getEmblKey()
  {
    return emblKey;
  }

  public void setEmblKey(String emblKey)
  {
    this.emblKey = emblKey;
  }
  
  public int getColour()
  {
    return colour;
  }

  public void setColour(int colour)
  {
    this.colour = colour;
  }
  
  public uk.ac.sanger.artemis.Feature getArtemisFeature()
  {
    Key key = new Key(emblKey);

    QualifierVector qualifiers = new QualifierVector();
    Qualifier qualifier = new Qualifier("colour", Integer.toString(getColour())); 
    qualifiers.add(qualifier);
    
    uk.ac.sanger.artemis.io.Feature f = null;
    try
    {
      Range r = new Range(getStart(), getEnd());
      Location location = new Location(r);
      
      f = new uk.ac.sanger.artemis.io.EmblStreamFeature(key, location, qualifiers);
      return new uk.ac.sanger.artemis.Feature(f);
    }
    catch(InvalidRelationException e)
    {
      e.printStackTrace();
    }
    catch(OutOfRangeException e)
    {
      e.printStackTrace();
    }
    return null;
  }

}