/* ScoreChanger.java
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/ScoreChanger.java,v 1.1 2004-06-09 09:47:29 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import java.awt.*;
import java.awt.event.*;

import javax.swing.*;

/**
 *  This is a JFrame component that contains ScoreScrollbar components.
 *
 *  @author Kim Rutherford
 *  @version $Id: ScoreChanger.java,v 1.1 2004-06-09 09:47:29 tjc Exp $
 **/

public class ScoreChanger extends JFrame {
  /**
   *  Create a new ScoreChanger.
   *  @param minimum_value The minimum allowable value for the scroll bars.
   *  @param maximum_value The maximum allowable value for the scroll bars.
   **/
  public ScoreChanger (final String name,
                       final ScoreChangeListener minimum_listener,
                       final ScoreChangeListener maximum_listener,
                       final int minimum_value, final int maximum_value)
      throws IllegalArgumentException {

    super (name);

    this.minimum_listener = minimum_listener;
    this.maximum_listener = maximum_listener;
    this.minimum_value = minimum_value;
    this.maximum_value = maximum_value;

    getContentPane ().setLayout (new GridLayout (5, 1));

    minimum_label = new JLabel ();

    getContentPane ().add (minimum_label);

    minimum_score_scrollbar =
      new ScoreScrollbar (minimum_value, maximum_value);

    final ScoreChangeListener this_minimum_listener =
      new ScoreChangeListener () {
        public void scoreChanged (ScoreChangeEvent event) {
          minimum_label.setText ("Minimum Cutoff: " + event.getValue ());
        }
      };

    minimum_score_scrollbar.addScoreChangeListener (this_minimum_listener);

    getContentPane ().add (minimum_score_scrollbar);

    maximum_label = new JLabel ();

    getContentPane ().add (maximum_label);

    maximum_score_scrollbar =
      new ScoreScrollbar (minimum_value, maximum_value);

    final ScoreChangeListener this_maximum_listener =
      new ScoreChangeListener () {
        public void scoreChanged (ScoreChangeEvent event) {
          maximum_label.setText ("Maximum Cutoff: " + event.getValue ());
        }
      };

    maximum_score_scrollbar.addScoreChangeListener (this_maximum_listener);

    getContentPane ().add (maximum_score_scrollbar);

    final FlowLayout flow_layout =
      new FlowLayout (FlowLayout.CENTER, 18, 0);

    final JPanel button_panel = new JPanel (flow_layout);

    final JButton reset_button = new JButton ("Reset");

    reset_button.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent e) {
        reset ();
      }
    });

    button_panel.add (reset_button);


    final JButton close_button = new JButton ("Close");

    close_button.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent e) {
        dispose ();
      }
    });


    button_panel.add (close_button);

    getContentPane ().add (button_panel);

    addWindowListener (new WindowAdapter () {
      public void windowClosing (WindowEvent event) {
        dispose ();
      }
    });

    pack ();

    setSize (270, getSize ().height);

    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
    setLocation (new Point ((screen.width - getSize ().width) / 2,
                            (screen.height - getSize ().height) / 2));

    reset ();
    
    minimum_score_scrollbar.addScoreChangeListener (minimum_listener);
    maximum_score_scrollbar.addScoreChangeListener (maximum_listener);
  }

  /**
   *  Call reset (), clean up the listeners then call super.dispose ().
   **/
  public void dispose () {
    reset ();

    minimum_score_scrollbar.removeScoreChangeListener (minimum_listener);
    maximum_score_scrollbar.removeScoreChangeListener (maximum_listener);

    super.dispose ();
  }

  /**
   *  Reset the minimum_score_scrollbar to 0 and the maximum_score_scrollbar
   *  to 100.
   **/
  private void reset () {
    minimum_score_scrollbar.setValue (minimum_value);
    maximum_score_scrollbar.setValue (maximum_value);

    minimum_label.setText ("Minimum Cutoff: " + minimum_value);
    maximum_label.setText ("Maximum Cutoff: " + maximum_value);

    final ScoreChangeEvent minimum_event =
      new ScoreChangeEvent (this, minimum_value);

    minimum_listener.scoreChanged (minimum_event);

    final ScoreChangeEvent maximum_event =
      new ScoreChangeEvent (this, maximum_value);

    maximum_listener.scoreChanged (maximum_event);
  }

  /**
   *  The ScoreChangeListener for the minimum score that was passed to the
   *  constructor.
   **/
  final ScoreChangeListener minimum_listener;

  /**
   *  The minimum score scrollbar the created by the constructor.
   **/
  final ScoreScrollbar minimum_score_scrollbar;

  /**
   *  The ScoreChangeListener for the maximum score that was passed to the
   *  constructor.
   **/
  final ScoreChangeListener maximum_listener;

  /**
   *  The maximum score scrollbar the created by the constructor.
   **/
  final ScoreScrollbar maximum_score_scrollbar;

  /**
   *  The minimum allowable value for the scroll bars.
   **/
  final int minimum_value;

  /**
   *  The maximum allowable value for the scroll bars.
   **/
  final int maximum_value;

  /**
   *  A Label that shows something like this: "Minimum Cutoff: 0"
   **/
  final JLabel minimum_label;

  /**
   *  A Label that shows something like this: "Maximum Cutoff: 100"
   **/
  final JLabel maximum_label;
}

