/* FeaturePlotGroup.java
 *
 * created: Wed Dec 16 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/FeaturePlotGroup.java,v 1.2 2004-12-03 18:11:28 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.sequence.*;
import uk.ac.sanger.artemis.plot.*;

import java.awt.*;
import java.awt.event.*;
import java.util.Vector;

import javax.swing.*;

/**
 *  This is a super-component containing several FeaturePlot components, each
 *  of which can toggled off and on.  The component contains a row of toggle
 *  buttons and then the FeaturePlot components below.
 *
 *  @author Kim Rutherford
 *  @version $Id: FeaturePlotGroup.java,v 1.2 2004-12-03 18:11:28 tjc Exp $
 **/

public class FeaturePlotGroup extends JFrame
    implements EntryChangeListener, FeatureChangeListener {
  /**
   *  Create a new FeaturePlotGroup component for the given feature.
   **/
  public FeaturePlotGroup (Feature feature) {
    super ("Graphs for: " + feature.getIDString ());
    this.feature = feature;
    this.entry = feature.getEntry ();

    // don't repeat an algorithm in this array.
    final FeatureAlgorithm [] plot_value_producers = {
      new HydrophobicityAlgorithm (getFeature ()),
      new HydrophilicityAlgorithm (getFeature ()),
      new CoilFeatureAlgorithm (getFeature ())
    };

    final Font font = Options.getOptions ().getFont ();

    setFont (font);

    GridBagLayout gridbag = new GridBagLayout();
    getContentPane ().setLayout (gridbag);

    GridBagConstraints c = new GridBagConstraints();

    c.anchor = GridBagConstraints.NORTH;
    c.gridwidth = GridBagConstraints.REMAINDER;
    c.gridheight = 1;
    c.weightx = 1;
    c.fill = GridBagConstraints.BOTH;
    c.weighty = 1;
//    c.insets = new Insets (0,0,5,0);

    for (int i = 0 ; i < plot_value_producers.length ; ++i) {
      final FeatureAlgorithm this_algorithm = plot_value_producers[i];

      final FeaturePlot new_feature_plot = new FeaturePlot (this_algorithm);

      gridbag.setConstraints (new_feature_plot, c);
      getContentPane ().add (new_feature_plot);
      new_feature_plot.setVisible (true);
    }

    getFeature ().getEntry ().addEntryChangeListener (this);
    getFeature ().addFeatureChangeListener (this);

    addWindowListener (new WindowAdapter () {
      public void windowClosing (WindowEvent event) {
        stopListening ();
        FeaturePlotGroup.this.dispose ();
      }
    });

    addComponentListener (new ComponentAdapter () {
      public void componentResized (ComponentEvent event) {
        fixScrollbar ();
        fireAdjustmentEvent (scrollbar.getValue ());
      }
    });

    int new_x_size = feature.getTranslation ().length () + SCROLL_BAR_SIZE;

    if (new_x_size > 1000) {
      new_x_size = 1000;
    }

    if (new_x_size < 300) {
      new_x_size = 300;
    }

    scrollbar = new JScrollBar (Scrollbar.HORIZONTAL);
    c.fill = GridBagConstraints.HORIZONTAL;
    c.weighty = 0;
    gridbag.setConstraints (scrollbar, c);
    scrollbar.addAdjustmentListener (new AdjustmentListener () {
      public void adjustmentValueChanged(AdjustmentEvent e) {
        fireAdjustmentEvent (e.getValue ());
      }
    });
    getContentPane ().add (scrollbar);

    final Component [] children = getContentPane ().getComponents ();

    for (int i = 0 ; i < children.length ; ++i) {
      if (children[i] instanceof FeaturePlot) {
        addDisplayAdjustmentListener ((FeaturePlot)children[i]);
      }
    }


    bottom_button_panel = new JPanel ();
    c.fill = GridBagConstraints.HORIZONTAL;
    c.weighty = 0;
    gridbag.setConstraints (bottom_button_panel, c);
    getContentPane ().add (bottom_button_panel);

    close_button = new JButton ("Close");
    close_button.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent e) {
        stopListening ();
        FeaturePlotGroup.this.dispose ();
      }
    });

    bottom_button_panel.add (close_button);

    pack ();

    // give each FeaturePlot component a height of 200
    setSize (new_x_size, getSize ().height);

    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
    setLocation (new Point ((screen.width - getSize ().width) / 2,
                            (screen.height - getSize ().height) / 2));

    setVisible(true);

    fixScrollbar ();

    fireAdjustmentEvent (1);
  }

  /**
   *  Fudge factor to take the scroll bar width into account when creating
   *  windows and scrollbars.
   **/
  private final int SCROLL_BAR_SIZE = 100;

  /**
   *  Remove this object as a entry change listener and then call
   *  stopListening () on each FeaturePlot component.
   **/
  private void stopListening () {
    getEntry ().removeEntryChangeListener (this);

    final Component [] children = getComponents ();

    for (int i = 0 ; i < children.length ; ++i) {
      if (children[i] instanceof FeaturePlot) {
        ((FeaturePlot)children[i]).stopListening ();
      }
    }

    getFeature ().getEntry ().removeEntryChangeListener (this);
    getFeature ().removeFeatureChangeListener (this);
  }

  /**
   *  Set the scrollbar maximum, minimum and value.
   **/
  private void fixScrollbar () {
    final int feature_length = getFeature ().getTranslation ().length ();

    int scroll_max = feature_length;

    if (scroll_max < 1) {
      scroll_max = 1;
    }

    scrollbar.setValues (1, getSize ().width,
                         1, scroll_max + SCROLL_BAR_SIZE);

    scrollbar.setBlockIncrement (getSize ().width);
  }

  /**
   *  Adds the specified event adjustemnt listener to receive adjustment
   *  change events from this object.
   *  @param l the event change listener.
   **/
  private  void addDisplayAdjustmentListener (DisplayAdjustmentListener l) {
    adjustment_listener_list.addElement (l);
  }

  /**
   *  Removes the specified event listener so that it no longer receives
   *  adjustment change events from this object.
   *  @param l the event change listener.
   **/
  private  void removeDisplayAdjustmentListener (DisplayAdjustmentListener l) {
    adjustment_listener_list.removeElement (l);
  }

  /**
   *  Implementation of the FeatureChangeListener interface.
   *  @param event The change event.
   **/
  public void featureChanged (FeatureChangeEvent event) {
    fixScrollbar ();
    fireAdjustmentEvent (scrollbar.getValue ());
  }

  /**
   *  Implementation of the EntryChangeListener interface.  We listen to
   *  EntryChange events so we can delete this component if the feature gets
   *  deleted.
   **/
  public void entryChanged (EntryChangeEvent event) {
    switch (event.getType ()) {
    case EntryChangeEvent.FEATURE_DELETED:
      if (event.getFeature () == getFeature ()) {
        stopListening ();
        dispose ();
      }
      break;
    default:
      // do nothing;
      break;
    }
  }

  /**
   *  Send a DisplayAdjustmentEvent to the objects that are listening for it.
   **/
  private void fireAdjustmentEvent (final int scroll_value) {
    final int feature_length = getFeature ().getTranslation ().length ();

    int end_value = scroll_value + getSize ().width;

    if (end_value > feature_length) {
      end_value = feature_length;
    }

    final DisplayAdjustmentEvent event =
      new DisplayAdjustmentEvent (this,
                                  scroll_value,
                                  end_value,
                                  getSize ().width,
                                  1,
                                  0,// this arg will be ignored
                                  false,
                                  DisplayAdjustmentEvent.SCALE_ADJUST_EVENT);

    final Vector targets;

    // copied from a book - synchronising the whole method might cause a
    // deadlock
    synchronized (this) {
      targets = (Vector) adjustment_listener_list.clone ();
    }

    for ( int i = 0 ; i < targets.size () ; ++i ) {
      DisplayAdjustmentListener target =
        (DisplayAdjustmentListener) targets.elementAt (i);

      target.displayAdjustmentValueChanged (event);
    }
  }

  /**
   *  Return the feature that this component is plotting.
   **/
  private Feature getFeature () {
    return feature;
  }

  /**
   *  Return the Entry that contains the Feature this object is displaying.
   **/
  private Entry getEntry () {
    return entry;
  }

  /**
   *  The feature we are plotting.
   **/
  private Feature feature;

  /**
   *  The Entry that contains the Feature this object is displaying.
   **/
  private Entry entry;

  /**
   *  Pressing this button will distroy the JFrame.
   **/
  private JButton close_button;

  /**
   *  A JPanel to hold the buttons.
   **/
  private JPanel bottom_button_panel;

  private JScrollBar scrollbar;

  /**
   *  A vector of those objects listening for adjustment events.
   **/
  final private Vector adjustment_listener_list = new Vector ();
}
