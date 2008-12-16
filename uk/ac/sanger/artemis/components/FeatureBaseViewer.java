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

import uk.ac.sanger.artemis.*;
import java.awt.event.*;

/**
 *  A component for viewing the bases of a feature.  Once created this
 *  component listens for FeatureChange events to keep the the sequence up to
 *  date.
 *
 *  @author Kim Rutherford
 *  @version $Id: FeatureBaseViewer.java,v 1.2 2008-12-16 11:46:15 tjc Exp $
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
  
  /**
   *  Create a new FeatureBaseViewer component to display the bases of the
   *  given Feature.
   *  @param feature The feature to view.
   *  @param include_numbers If true then the sequence will be numbered
   *    (every second line of the display will be numbers rather than
   *    sequence).
   **/
  public FeatureBaseViewer (final Feature feature,
                            final boolean include_numbers) 
  {
    this.feature = feature;
    this.entry = feature.getEntry ();

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
    final StringBuffer header_buffer = new StringBuffer();

    header_buffer.append(getFeature().getSystematicName());
    header_buffer.append(" ");
    header_buffer.append(getFeature().getIDString());
    header_buffer.append(" ");

    final String product = getFeature().getProductString();

    if(product == null) 
      header_buffer.append("undefined product");
    else 
      header_buffer.append(product);
    
    header_buffer.append(" ").append(getFeature().getWriteRange());
    
    final String comment = ">" + header_buffer.toString();
      
    final String bases = getFeature ().getBases ().toUpperCase ();
    sequence_viewer.setSequence (comment, bases);
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
