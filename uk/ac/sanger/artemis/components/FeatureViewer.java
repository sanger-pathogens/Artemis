/* FeatureViewer.java
 *
 * created: Thu Nov 19 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/FeatureViewer.java,v 1.2 2005-01-11 16:07:35 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;

import java.io.*;
import java.awt.event.*;

/**
 *  A viewer for Feature objects.
 *
 *  @author Kim Rutherford
 *  @version $Id: FeatureViewer.java,v 1.2 2005-01-11 16:07:35 tjc Exp $
 *
 **/

public class FeatureViewer
    implements EntryChangeListener, FeatureChangeListener 
{

  /**
   *  The Feature this object is displaying.
   **/
  private Feature view_feature;

  /**
   *  The Entry that contains the Feature this object is displaying.
   **/
//private Entry entry;

  /**
   *  The FileViewer object that is displaying the feature.
   **/
  private FileViewer file_viewer;

  /**
   *  Create a new FeatureViewer object from the given Feature.
   **/
  public FeatureViewer(Feature view_feature) 
  {
    this.view_feature = view_feature;
//  this.entry = view_feature.getEntry();
    
    file_viewer = new FileViewer("Artemis Feature View: " +
                                 view_feature.getIDString());
    readFeature(view_feature);

    view_feature.getEntry().addEntryChangeListener(this);
    view_feature.addFeatureChangeListener(this);

    file_viewer.addWindowListener(new WindowAdapter() 
    {
      public void windowClosed(WindowEvent event)
      {
        stopListening();
      }
    });
  }

  /**
   *  Remove this object as a feature and entry change listener.
   **/
  public void stopListening()
  {
    view_feature.getEntry().removeEntryChangeListener(this);
    getFeature().removeFeatureChangeListener(this);
  }
  
  /**
   *  Implementation of the EntryChangeListener interface.  We listen to
   *  EntryChange events so we can delete this component if the feature gets
   *  deleted.
   **/
  public void entryChanged(EntryChangeEvent event) 
  {
    switch(event.getType()) 
    {
      case EntryChangeEvent.FEATURE_DELETED:
        if(event.getFeature() == view_feature) 
        {
          stopListening();
          file_viewer.dispose();
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
  public void featureChanged(FeatureChangeEvent event) 
  {
    // re-read the information from the feature
    readFeature(view_feature);
  }

  /**
   *  Read the given Feature into this FeatureViewer object.
   **/
  public void readFeature(Feature feature) 
  {
    try 
    {
      file_viewer.clear();
      file_viewer.appendFile(view_feature.toReader());
    }
    catch(uk.ac.sanger.artemis.io.ReadFormatException e) 
    {
      throw new Error("internal error - unexpected exception: " +
                      e.getMessage() +
                      (e.getLineNumber() > 1 ?
                       " at line " + e.getLineNumber() :
                       ""));
    } 
    catch(IOException e)
    {
      throw new Error("internal error - unexpected exception: " +
                      e.getMessage());
    }
  }

  /**
   *  Return the Feature we are viewing.
   **/
  public Feature getFeature() 
  {
    return view_feature;
  }
  
  /**
   *  Return the Entry that contains the Feature this object is displaying.
   **/
//private Entry getEntry()
//{
//  return entry;
//}

}

