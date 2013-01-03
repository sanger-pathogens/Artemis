/* FeatureBaseViewer.java
 *
 * created: Sat Dec 19 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/FeatureBaseViewer.java,v 1.2 2008-12-16 11:46:15 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import uk.ac.sanger.artemis.*;

/**
 *  A component for viewing the bases of a feature. Once created this
 *  component listens for FeatureChange events to keep the the sequence up to
 *  date.
 *  @author Kim Rutherford
 **/

public class FeatureBaseViewer
    implements EntryChangeListener, FeatureChangeListener 
{
  /** The Feature that this component is showing information about.   */
  private Feature feature = null;

  /** The SequenceViewer object that is displaying the feature bases. */
  private SequenceViewer sequence_viewer;

  /** The Entry that contains the Feature this object is displaying.  */
  private Entry entry;
  
  private FeatureSegmentVector segments;
  
  /**
   *  Create a new FeatureBaseViewer component to display the bases of the
   *  given Feature.
   *  @param feature The feature to view.
   *  @param include_numbers If true then the sequence will be numbered
   *    (every second line of the display will be numbers rather than
   *    sequence).
   **/
  public FeatureBaseViewer (final Feature feature,
                            final boolean include_numbers,
                            final FeatureSegmentVector segments) 
  {
    this.feature = feature;
    this.entry = feature.getEntry ();
    this.segments = segments;

    sequence_viewer =
      new SequenceViewer ("Feature base viewer for feature:" +
                          getFeature ().getIDString (), include_numbers);  
    redisplay ();
    getFeature ().getEntry ().addEntryChangeListener (this);
    getFeature ().addFeatureChangeListener (this);
    sequence_viewer.addWindowListener (new WindowAdapter () 
    {
      public void windowClosed (WindowEvent event)
      {
        stopListening ();
      }
    });
  }
  
  /**
   *  Remove this object as a entry and feature change listener.
   **/
  private void stopListening () 
  {
    getEntry ().removeEntryChangeListener (this);
    getFeature ().removeFeatureChangeListener (this);
  }
  
  /**
   *  Implementation of the EntryChangeListener interface.  We listen to
   *  EntryChange events so we can delete this component if the feature gets
   *  deleted.
   **/
  public void entryChanged (EntryChangeEvent event) 
  {
    switch (event.getType ()) 
    {
    case EntryChangeEvent.FEATURE_DELETED:
      if (event.getFeature () == getFeature ()) 
      {
        stopListening ();
        sequence_viewer.dispose ();
      }
      break;
    default:
      // do nothing;
      break;
    }
  }

  /**
   *  Implementation of the FeatureChangeListener interface.  We need to
   *  listen to feature change events from the Features in this object so that
   *  we can keep the display up to date.
   *  @param event The change event.
   **/
  public void featureChanged (FeatureChangeEvent event) 
  {
    // re-read the information from the feature
    redisplay ();
  }

  /**
   *  Redisplay the bases.
   **/
  private void redisplay () 
  {
    final String product = getFeature().getProductString();
    final StringBuilder hdr = new StringBuilder();
    hdr.append(getFeature().getSystematicName()).append(" ");
    hdr.append(getFeature().getIDString()).append(" ");
    hdr.append(product == null ? "undefined product" : product);
    
    final String bases;
    if(segments != null)  // display just for selected segments
    {
      final StringBuilder buffer = new StringBuilder();
      segments.sortByPosition();
      for(int i = 0; i < segments.size(); ++i) 
      {
        final FeatureSegment segment = segments.elementAt(i);
        buffer.append(segment.getBases());
        hdr.append(i == 0 ? " " : ",").append(segment.getRawRange().toString());
      }
      hdr.append(getFeature ().isForwardFeature() ? " forward" : " reverse");
      bases = buffer.toString();
    }
    else
    {
      hdr.append(" ").append(getFeature().getWriteRange());
      bases = getFeature ().getBases ();
    }

    sequence_viewer.setSequence (">" + hdr.toString(), bases.toUpperCase ());
  }

  /**
   *  Return the feature this component is showing information about.
   **/
  private Feature getFeature () 
  {
    return feature;
  }

  /**
   *  Return the Entry that contains the Feature this object is displaying.
   **/
  private Entry getEntry () 
  {
    return entry;
  }
}
