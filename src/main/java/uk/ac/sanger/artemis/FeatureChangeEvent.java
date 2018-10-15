/* FeatureChangeEvent.java
 *
 * created: Sat Oct 17 1998
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 1998,1999,2000  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/FeatureChangeEvent.java,v 1.3 2006-07-19 16:10:31 tjc Exp $
 */

package uk.ac.sanger.artemis;

import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.io.Location;
import uk.ac.sanger.artemis.io.QualifierVector;

import java.io.*;

/**
 *  This event is sent when a change occurs within a feature.  eg the location
 *  or qualifier changes.
 *
 *  @author Kim Rutherford
 *  @version $Id: FeatureChangeEvent.java,v 1.3 2006-07-19 16:10:31 tjc Exp $
  */

public class FeatureChangeEvent extends ChangeEvent 
{
  /** Event type - feature location changed. */
  public static final int LOCATION_CHANGED = 1;

  /** Event type - feature qualifier changed. */
  public static final int QUALIFIER_CHANGED = 2;

  /** Event type - feature key changed. */
  public static final int KEY_CHANGED = 3;

  /**  Event type - feature has been re-read. */
  public static final int ALL_CHANGED = 4;
  
  /** feature segment added or deleted */
  public static final int SEGMENT_CHANGED = 5;
  
  /**
   *  The feature that this event was generated for.
   **/
  private Feature feature;

  /**
   *  This is the type of this event (eg LOCATION_CHANGED, QUALIFIER_CHANGED,
   *  etc), as passed to the constructor
   **/
  private int type;

  /**
   *  See getOldKey().
   **/
  private Key old_key;

  /**
   *  See getOldLocation().
   **/
  private Location old_location;

  /**
   *  See getOldQualifiers().
   **/
  private QualifierVector old_qualifiers;

  /**
   *  See getNewKey().
   **/
  private Key new_key;

  /**
   *  See getNewLocation().
   **/
  private Location new_location;

  /**
   *  See getNewQualifiers().
   **/
  private QualifierVector new_qualifiers;

  /**
   *  Create a new FeatureChangeEvent object.
   *  @param source The object that generated the event.
   *  @param feature This Feature object that this event refers to.
   *  @param key The key before the Feature changed or null if the key hasn't
   *    changed.
   *  @param location The Location before the Feature changed or null if the
   *    location hasn't changed.
   *  @param key The Qualifiers before the Feature changed or null if the
   *    qualifiers have't changed.
   *  @param type This type of the event.
   **/
  public FeatureChangeEvent (Object source,
                             Feature feature,
                             Key old_key,
                             Location old_location,
                             QualifierVector old_qualifiers,
                             int type) 
  {
    super (source);

    this.feature = feature;
    this.type = type;

    this.old_key = old_key;
    this.old_location = old_location;
    this.old_qualifiers = old_qualifiers;

    this.new_key = feature.getKey ();
    this.new_location = feature.getLocation ();
    this.new_qualifiers = feature.getQualifiers ();
  }

  /**
   *  Return the Feature of this event.
   **/
  public Feature getFeature () 
  {
    return feature;
  }
   
  /**
   *  Return the type of this event as passed to the constructor.
   **/
  public int getType () 
  {
    return type;
  }

  /**
   *  Return the old Key from before the event.  Could return null.
   **/
  public Key getOldKey ()
  {
    return old_key;
  }

  /**
   *  Return the old Location from before the event.  Could return null.
   **/
  public Location getOldLocation ()
  {
    return old_location;
  }

  /**
   *  Return the old qualifiers from before the event.  Could return null.
   **/
  public QualifierVector getOldQualifiers ()
  {
    return old_qualifiers;
  }

  /**
   *  Return the Key after the event (extracted from the Feature was is
   *  passed to the constructor).
   **/
  public Key getNewKey ()
  {
    return new_key;
  }

  /**
   *  Return the Location after the event (extracted from the Feature was is
   *  passed to the constructor).
   **/
  public Location getNewLocation ()
  {
    return new_location;
  }

  /**
   *  Return the qualifiers after the event (extracted from the Feature was is
   *  passed to the constructor).
   **/
  public QualifierVector getNewQualifiers ()
  {
    return new_qualifiers;
  }

  /**
   *  Return true if and only if this event caused an actual change in the
   *  Feature (rather than someone just pressing OK in a FeatureEdit for
   *  example).
   **/
  public boolean featureHasChanged ()
  {
    if(getOldKey() != null)
    {
      if(!getOldKey().equals(getNewKey()))
      {
        return true;
      }
    }

    if(getOldLocation() != null)
    {
      if(!getOldLocation().equals(getNewLocation()))
      {
        return true;
      }
    }

    if(getOldQualifiers() != null)
    {
      if(!getOldQualifiers().equals(getNewQualifiers()))
      {
        return true;
      }
    }

    return false;
  }
 
}
