/* GotoMenu.java
 *
 * created: Thu Jan  7 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/GotoMenu.java,v 1.3 2006-01-17 16:05:05 tjc Exp $
 **/

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.sequence.Marker;
import uk.ac.sanger.artemis.sequence.MarkerRange;
import uk.ac.sanger.artemis.util.OutOfRangeException;

import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import javax.swing.JFrame;
import javax.swing.JMenuItem;
import javax.swing.KeyStroke;

/**
 *  A JMenu with commands for moving around the entries.
 *
 *  @author Kim Rutherford
 *  @version $Id: GotoMenu.java,v 1.3 2006-01-17 16:05:05 tjc Exp $
 **/

public class GotoMenu extends SelectionMenu {
  
  private static final long serialVersionUID = 1L;

  /** GotoEventSource object that was passed to the constructor */
  private GotoEventSource goto_event_source;

  /**
   *  Set by the constructor to the reference of the first feature in
   *  selection_features or null if selection_features is empty.
   **/
  private Feature selection_feature;

  /**
   *  Set by the constructor to the reference of the first segment in
   *  selection_segments or null if selection_segments is empty.
   **/
  private FeatureSegment selection_segment;
  
  /** shortcut for the Navigator */
  private final static KeyStroke NAVIGATOR_KEY =
    KeyStroke.getKeyStroke (KeyEvent.VK_G, 
                            Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()); 

  /** shortcut to go to the start of the selection */
  private final static KeyStroke START_OF_SELECTION_KEY =
    KeyStroke.getKeyStroke (KeyEvent.VK_LEFT, 
                            Toolkit.getDefaultToolkit().getMenuShortcutKeyMask());

  /** shortcut to go to the end of the selection */
  private final static KeyStroke END_OF_SELECTION_KEY =
    KeyStroke.getKeyStroke (KeyEvent.VK_RIGHT, 
                            Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()); 

  /** shortcut to go to the start of the sequence */
  private final static KeyStroke START_OF_SEQUENCE_KEY =
    KeyStroke.getKeyStroke (KeyEvent.VK_UP, 
                            Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()); 

  /** shortcut to go to the end of the sequence */
  private final static KeyStroke END_OF_SEQUENCE_KEY =
    KeyStroke.getKeyStroke (KeyEvent.VK_DOWN, 
                            Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()); 

  /** shortcut to go to the start of the feature */
  private final static KeyStroke START_OF_FEATURE_KEY =
    KeyStroke.getKeyStroke (KeyEvent.VK_COMMA, KeyEvent.SHIFT_MASK); 
  
  /** shortcut to go to the start of the feature */
  private final static KeyStroke END_OF_FEATURE_KEY =
    KeyStroke.getKeyStroke (KeyEvent.VK_PERIOD, KeyEvent.SHIFT_MASK); 
  
  /**
   *  Create a new GotoMenu object.
   *  @param frame The JFrame that owns this JMenu.
   *  @param selection The Selection that the commands in the menu will
   *    operate on.
   *  @param goto_event_source The object the we will call makeBaseVisible ()
   *    on.
   *  @param entry_group The EntryGroup object used to create the Navigator
   *    component.
   *  @param menu_name The name of the new menu.
   **/
  public GotoMenu (final JFrame frame,
                   final Selection selection,
                   final GotoEventSource goto_event_source,
                   final EntryGroup entry_group,
                   final String menu_name) {
    super (frame, menu_name, selection);
    this.goto_event_source = goto_event_source;

    final JMenuItem navigator_item = new JMenuItem ("Navigator ...");
    navigator_item.setAccelerator (NAVIGATOR_KEY);
    navigator_item.addActionListener(new ActionListener () {
      private Navigator navigator = null;
      public void actionPerformed (ActionEvent event) {
        if(navigator == null || navigator.isVisible())
        {
          if(navigator != null)
            navigator.setDisposeOnClose(true);
          navigator = new Navigator (getSelection (),
                       GotoMenu.this.goto_event_source,
                       entry_group);
        }
        else
          navigator.setVisible(true);
      }
    });

    final JMenuItem goto_selection_start_item = new JMenuItem ("Start of Selection");
    goto_selection_start_item.setAccelerator (START_OF_SELECTION_KEY);
    goto_selection_start_item.addActionListener(new ActionListener () {
      public void actionPerformed (ActionEvent event) {
        makeSelectionStartVisible ();
      }
    });

    final JMenuItem goto_selection_end_item = new JMenuItem ("End of Selection");
    goto_selection_end_item.setAccelerator (END_OF_SELECTION_KEY);
    goto_selection_end_item.addActionListener(new ActionListener () {
      public void actionPerformed (ActionEvent event) {
        makeSelectionEndVisible ();
      }
    });

    final JMenuItem goto_feature_start_item = new JMenuItem ("Feature Start");
    goto_feature_start_item.setAccelerator(START_OF_FEATURE_KEY);
    goto_feature_start_item.addActionListener(new ActionListener () {
      public void actionPerformed (ActionEvent event) {
        gotoFeatureStart ();
      }
    });

    final JMenuItem goto_feature_end_item = new JMenuItem ("Feature End");
    goto_feature_end_item.setAccelerator(END_OF_FEATURE_KEY);
    goto_feature_end_item.addActionListener(new ActionListener () {
      public void actionPerformed (ActionEvent event) {
        gotoFeatureEnd ();
      }
    });

    final JMenuItem goto_feature_base = new JMenuItem ("Feature Base Position ...");
    goto_feature_base.addActionListener(new ActionListener () {
      public void actionPerformed (ActionEvent event) {
        gotoFeaturePosition (false);
      }
    });

    final JMenuItem goto_feature_aa_position = new JMenuItem ("Feature Amino Acid ...");
    goto_feature_aa_position.addActionListener(new ActionListener () {
      public void actionPerformed (ActionEvent event) {
        gotoFeaturePosition (true);
      }
    });


    final JMenuItem goto_first_item = new JMenuItem ("Start of Sequence");
    goto_first_item.setAccelerator (START_OF_SEQUENCE_KEY);
    goto_first_item.addActionListener(new ActionListener () {
      public void actionPerformed (ActionEvent event) {
        goto_event_source.gotoFirstBase ();
      }
    });

    final JMenuItem goto_last_item = new JMenuItem ("End of Sequence");
    goto_last_item.setAccelerator (END_OF_SEQUENCE_KEY);
    goto_last_item.addActionListener(new ActionListener () {
      public void actionPerformed (ActionEvent event) {
        goto_event_source.gotoLastBase ();
      }
    });

    add (navigator_item);
    addSeparator ();
    add (goto_selection_start_item);
    add (goto_selection_end_item);
    add (goto_feature_start_item);
    add (goto_feature_end_item);
    add (goto_first_item);
    add (goto_last_item);
    add (goto_feature_base);
    add (goto_feature_aa_position);
  }
  
  /**
   *  Create a new GotoMenu object.
   *  @param frame The JFrame that owns this JMenu.
   *  @param selection The Selection that the commands in the menu will
   *    operate on.
   *  @param goto_event_source The object the we will call makeBaseVisible ()
   *    on.
   *  @param entry_group The EntryGroup object used to create the Navigator
   *    component.
   **/
  public GotoMenu (final JFrame frame,
                   final Selection selection,
                   final GotoEventSource goto_event_source,
                   final EntryGroup entry_group) {
    this (frame, selection, goto_event_source, entry_group, "Goto");
  }

  /**
   *  This method sends an GotoEvent to all the GotoEvent listeners that will
   *  make the first base of the selection visible.
   **/
  private void makeSelectionStartVisible () {
    final GotoEvent new_event =
      new GotoEvent (this, getSelection ().getStartBaseOfSelection ());

    goto_event_source.sendGotoEvent (new_event);
  }

  /**
   *  This method sends an GotoEvent to all the GotoEvent listeners that will
   *  make the last base of the selection visible.
   **/
  private void makeSelectionEndVisible () {
    final GotoEvent new_event =
      new GotoEvent (this, getSelection ().getEndBaseOfSelection ());

    goto_event_source.sendGotoEvent (new_event);
  }

  /**
   *  Goto the start of the (first) selected feature.  (FeatureDisplay only.)
   **/
  private void gotoFeatureStart () {
    setInternalVariables ();
    if (selection_feature == null) {
      if (selection_segment == null) {
        // do nothing
      } else {
        // go to the start of the parent feature of the first selected segment
        final Feature segment_feature = selection_segment.getFeature ();
        makeFeatureStartVisible (segment_feature);
      }
    } else {
      makeFeatureStartVisible (selection_feature);
    }
  }

  /**
   *  Goto the end of the (first) selected feature.  (FeatureDisplay only.)
   **/
  private void gotoFeatureEnd () {
    setInternalVariables ();
    if (selection_feature == null) {
      if (selection_segment == null) {
        // do nothing
      } else {
        // go to the end of the parent feature of the first selected segment
        final Feature segment_feature = selection_segment.getFeature ();
        makeFeatureEndVisible (segment_feature);
      }
    } else {
      makeFeatureEndVisible (selection_feature);
    }
  }

  /**
   *  This method will ask the user for a base or amino acid position within
   *  the first selected feature (using a TextRequester component) and will
   *  then goto that position.
   *  @param goto_aa_position If true goto an amino acid position, otherwise
   *    goto a base position
   **/
  private void gotoFeaturePosition (final boolean goto_aa_position) {
    if (!checkForSelectionFeatures ()) {
      return;
    }

    final FeatureVector selected_features = getSelection ().getAllFeatures ();

    if (selected_features.size () > 1) {
      new MessageDialog (getParentFrame (), "select only one feature");
      return;
    }
    
    final Feature first_feature = selected_features.elementAt (0);

    final TextRequester text_requester =
      goto_aa_position ?
      new TextRequester ("amino acid position within selected feature:",
                         18, "") :
      new TextRequester ("base position within selected feature:",
                         18, "");

    text_requester.addTextRequesterListener (new TextRequesterListener () {
      public void actionPerformed (final TextRequesterEvent event) {
        final String position_string = event.getRequesterText ().trim ();

        if (position_string.length () == 0) {
          return;
        }

        try {
          final int feature_position =
            Integer.valueOf (position_string).intValue ();

          final Marker sequence_marker;

          if (goto_aa_position) {
            // base and aa positions are numbered from 1
            final int aa_position = (feature_position - 1) * 3 + 1;

            sequence_marker =
              first_feature.getPositionInSequence (aa_position);
          } else {
            sequence_marker =
              first_feature.getPositionInSequence (feature_position);
          }

          final MarkerRange range = new MarkerRange (sequence_marker);

          goto_event_source.gotoBase (sequence_marker);

          /*final*/ MarkerRange selection_range;
          
          if (goto_aa_position) {
            // we want to select the whole codon if we are going to a amino
            // acid position
            try { 
              final MarkerRange codon_end_marker_range =
                new MarkerRange (sequence_marker.moveBy (2));
              selection_range =
                range.combineRanges (codon_end_marker_range, false);
            } catch (OutOfRangeException e) {
              // just select one base
              selection_range = range;
            }
          } else {
            selection_range = range;
          }

          // select that base (or range) 
          getSelection ().setMarkerRange (selection_range);

        } catch (NumberFormatException e) {
          new MessageDialog (getParentFrame (),
                             "this is not a number: " + position_string);
        } catch (OutOfRangeException e) {
          new MessageDialog (getParentFrame (), "the base position is not " +
                             "within the selection feature: " +
                             position_string);
        }
      }
    });

    text_requester.setVisible(true);
  }

  /**
   *  This method sends an GotoEvent to all the GotoEvent listeners that will
   *  make the given base visible.
   **/
  private void makeBaseVisible (final Marker base_marker) {
    goto_event_source.gotoBase (base_marker);
  }

  /**
   *  Scroll the display so that the start of the given feature is visible
   *  (it's start will be in the middle of the screen after the call).
   *  @param feature The feature to make visible.
   **/
  private void makeFeatureStartVisible (Feature feature) {
    if (feature == null)
      return;
    makeBaseVisible (feature.getFirstBaseMarker ());
  }

  /**
   *  Scroll the display so that the start of the given feature is visible
   *  (it's start will be in the middle of the screen after the call).
   *  @param feature The feature to make visible.
   **/
  private void makeFeatureEndVisible (Feature feature) {
    if (feature == null)
      return;
    makeBaseVisible (feature.getLastBaseMarker ());
  }

  /**
   *  Initialise selection_feature and selection_segment.
   **/
  private void setInternalVariables () {
    final FeatureVector selection_features =
      getSelection ().getSelectedFeatures ();
    final FeatureSegmentVector selection_segments =
      getSelection ().getSelectedSegments ();

    if (selection_features.size () > 0) {
      selection_feature = selection_features.elementAt (0);
    } else {
      selection_feature = null;
    }

    if (selection_segments.size () > 0) {
      selection_segment = selection_segments.elementAt (0);
    } else {
      selection_segment = null;
    }
  }
}
