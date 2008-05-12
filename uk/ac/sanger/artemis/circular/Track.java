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

import java.awt.Color;

import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureKeyPredicate;
import uk.ac.sanger.artemis.FeatureKeyQualifierPredicate;
import uk.ac.sanger.artemis.FeaturePredicate;
import uk.ac.sanger.artemis.io.Key;

public class Track
{
  private double position = .9d;
  private float size = 10.f;
  private FeaturePredicate featurePredicate;
  private boolean showForward = true;
  private boolean showReverse = true;
  private boolean notQualifier = false;
  private boolean any = false;
  private String keyStr;
  private String qualifier;
  private String qualifierValue;
  private Entry entry;
  private Color colour;
  
  public Track(double position, Entry entry)
  {
    this.position = position;
    this.entry = entry;
  }
  
  public Track(double position, 
               String keyStr,
               String qualifier,
               boolean notQualifier,
               boolean showForward,
               boolean showReverse,
               Entry entry)
  {
    this.position = position;
    this.showForward = showForward;
    this.showReverse = showReverse;
    this.keyStr = keyStr;
    this.qualifier = qualifier;
    this.notQualifier = notQualifier;
    this.entry = entry;
    
    if(keyStr != null)
    {
      final Key key = new Key(keyStr);
      if(qualifier != null)
        featurePredicate = 
          new FeatureKeyQualifierPredicate(key, qualifier, isNotQualifier());
      else
        featurePredicate = new FeatureKeyPredicate(key);
    }
    
    if(featurePredicate == null)
      any = true;
  }
  
  
  public Track(double position, 
               String keyStr,
               boolean showForward, 
               boolean showReverse,
               Entry entry)
  {
    this(position, keyStr, null, true, showForward, showReverse, entry);
  }
  
  
  public double getPosition()
  {
    return position;
  }

  public void setPosition(double position)
  {
    this.position = position;
  }

  /**
   * Test if this feature is drawn on this track
   * @param this_feature
   * @return
   */
  public boolean isOnTrack(Feature this_feature)
  {
    if(getEntry() != null && !getEntry().contains(this_feature))
      return false;
    
    if(isAny())
      return true;
    
    if(featurePredicate == null)
      return false;
    
    if(featurePredicate.testPredicate(this_feature))
    {
      if( (this_feature.isForwardFeature()  && isShowForward()) ||
          (!this_feature.isForwardFeature() && isShowReverse()) )
        return true;
    }
    return false;
  }

  public boolean isShowForward()
  {
    return showForward;
  }

  public void setShowForward(boolean showForward)
  {
    this.showForward = showForward;
  }

  public boolean isShowReverse()
  {
    return showReverse;
  }

  public void setShowReverse(boolean showReverse)
  {
    this.showReverse = showReverse;
  }

  public boolean isAny()
  {
    return any;
  }

  public void setAny(boolean any)
  {
    this.any = any;
  }

  public FeaturePredicate getFeaturePredicate()
  {
    return featurePredicate;
  }

  public void setFeaturePredicate(FeaturePredicate featurePredicate)
  {
    this.featurePredicate = featurePredicate;
  }

  public String getKeyStr()
  {
    return keyStr;
  }

  public void setKeyStr(String keyStr)
  {
    this.keyStr = keyStr;
  }

  public String getQualifier()
  {
    return qualifier;
  }

  public void setQualifier(String qualifier)
  {
    this.qualifier = qualifier;
  }

  public boolean isNotQualifier()
  {
    return notQualifier;
  }

  public void setNotQualifier(boolean notQualifier)
  {
    this.notQualifier = notQualifier;
  }

  public String getQualifierValue()
  {
    return qualifierValue;
  }

  public void setQualifierValue(String qualifierValue)
  {
    this.qualifierValue = qualifierValue;
  }

  public float getSize()
  {
    return size;
  }

  public void setSize(float size)
  {
    this.size = size;
  }

  public Entry getEntry()
  {
    return entry;
  }

  public void setEntry(Entry entry)
  {
    this.entry = entry;
  }

  public Color getColour()
  {
    return colour;
  }

  public void setColour(Color colour)
  {
    this.colour = colour;
  }
}