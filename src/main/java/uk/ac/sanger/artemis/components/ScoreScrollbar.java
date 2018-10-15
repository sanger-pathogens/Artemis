/* ScoreScrollbar.java
 *
 * created: Thu Oct 21 1999
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 1999  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/ScoreScrollbar.java,v 1.1 2004-06-09 09:47:31 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import java.awt.*;
import java.awt.event.*;

import javax.swing.*;

/**
 *  This component is a Scrollbar that generates ScoreChange events when the
 *  value changes.
 *
 *  @author Kim Rutherford
 *  @version $Id: ScoreScrollbar.java,v 1.1 2004-06-09 09:47:31 tjc Exp $
 **/

public class ScoreScrollbar extends JScrollBar
  implements AdjustmentListener {
  /**
   *  Constructs a new horizontal score scroll bar with an initial value of 0.
   *  @param minimum_value The minimum allowable value for the scroll bar.
   *  @param maximum_value The maximum allowable value for the scroll bar.
   **/
  public ScoreScrollbar (final int minimum_value, final int maximum_value) {
    this (Scrollbar.HORIZONTAL, minimum_value, minimum_value, maximum_value);
  }

  /**
   *  Constructs a new score scroll bar with the specified orientation.
   *
   *  The orientation argument must take one of the two values
   *  java.awt.Scrollbar.HORIZONTAL, or java.awt.Scrollbar.VERTICAL,
   *  indicating a horizontal or vertical scroll bar, respectively.
   *  @param orientation indicates the orientation of the scroll bar.
   *  @param value The initial value of the scrollbar.
   *  @param minimum_value The minimum allowable value for the scroll bar.
   *  @param maximum_value The maximum allowable value for the scroll bar.
   *  @exception IllegalArgumentException when an illegal value for the
   *    orientation argument is supplied or if the value parameter is less
   *    than minimum_value or greater than maximum_value.
   **/
  public ScoreScrollbar (final int orientation, final int value,
                         final int minimum_value, final int maximum_value)
      throws IllegalArgumentException {
    super (orientation,
           (value < minimum_value || value > maximum_value ?
            minimum_value :
            value),
           1, minimum_value, maximum_value + 1);
  }

  /**
   *  Add the given ScoreChangeListener as a listener for ScoreChange events
   *  from this components.
   **/
  public void addScoreChangeListener (final ScoreChangeListener l) {
    score_change_listeners.addElement (l);

    if (score_change_listeners.size () == 1) {
      addAdjustmentListener (this);
    }
  }

  /**
   *  Remove the given ScoreChangeListener as a listener for ScoreChange
   *  events from this components.
   **/
  public void removeScoreChangeListener (final ScoreChangeListener l) {
    score_change_listeners.addElement (l);

    if (score_change_listeners.size () == 0) {
      removeAdjustmentListener (this);
    }
  }

  /**
   *  Implementation of the AdjustmentListener interface.
   **/
  public void adjustmentValueChanged (AdjustmentEvent e) {
    for (int i = 0 ; i < score_change_listeners.size () ; ++i) {
      final ScoreChangeEvent event =
        new ScoreChangeEvent (this, getValue ());
      final ScoreChangeListener this_listener =
        (ScoreChangeListener)(score_change_listeners.elementAt (i));
      this_listener.scoreChanged (event);
    }
  }

  /**
   *  The ScoreChangeListener objects that have been added with
   *  addScoreChangeListener ().
   **/
  private final java.util.Vector score_change_listeners =
    new java.util.Vector ();
}

