/* FeatureList.java
 *
 * created: Fri Oct  9 1998
 *
 * This file is part of Artemis
 *
 * Copyright (C) 1998,1999,2000,2001,2002  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/FeatureList.java,v 1.1 2004-06-09 09:46:40 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.sequence.*;
import uk.ac.sanger.artemis.plot.*;

import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.InvalidRelationException;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.QualifierInfo;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.io.StreamQualifier;
import uk.ac.sanger.artemis.util.StringVector;

import java.awt.event.*;
import java.awt.*;
import java.text.*;
import java.util.*;

import javax.swing.*;

/**
 *  This component gives the user a list containing the details the current
 *  Features.
 *
 *  @author Kim Rutherford
 *  @version $Id: FeatureList.java,v 1.1 2004-06-09 09:46:40 tjc Exp $
 *
 **/

public class FeatureList extends EntryGroupPanel
  implements EntryGroupChangeListener,
             EntryChangeListener, FeatureChangeListener,
             SelectionChangeListener, DisplayComponent
{
  /**
   *  Create a new FeatureList with the default number of rows.
   *  @param entry_group The EntryGroup that this component will display.
   *  @param selection The Selection object for this component.  Selected
   *    objects will be highlighted.
   *  @param goto_event_source The object to use when we need to call
   *    gotoBase ().
   **/
  public FeatureList (final EntryGroup entry_group,
                      final Selection selection,
                      final GotoEventSource goto_event_source,
                      final BasePlotGroup base_plot_group) {
    super (entry_group, selection, goto_event_source, base_plot_group);

    getCanvas ().addMouseListener (new MouseAdapter () {
      /**
       *  Listen for mouse press events so that we can do popup menus and
       *  selection.
       **/
      public void mousePressed (MouseEvent event) {
        if (isMenuTrigger (event)) {
          final FeaturePopup popup =
            new FeaturePopup (FeatureList.this,
                              getEntryGroup (),
                              getSelection (),
                              getGotoEventSource (),
                              getBasePlotGroup ());
          final JComponent parent = (JComponent) event.getSource ();

          popup.show (parent, event.getX (), event.getY ());
        } else {
          handleCanvasMousePress (event);
        }
      }
    });

    createScrollbars ();

    addComponentListener (new ComponentAdapter () {
      public void componentShown (ComponentEvent e) {
        repaintCanvas ();
      }
      public void componentResized (ComponentEvent e) {
        repaintCanvas ();
      }
    });

    getSelection ().addSelectionChangeListener (this);

    // changes to the EntryGroup will be noticed by listening for EntryChange
    // and FeatureChange events.

    getEntryGroup ().addEntryGroupChangeListener (this);
    getEntryGroup ().addEntryChangeListener (this);
    getEntryGroup ().addFeatureChangeListener (this);


    // find the maximum posible width for the high and low positions
    final int sequence_length = getEntryGroup ().getSequenceLength ();
    max_base_pos_width = (int)(Math.log (sequence_length)/Math.log (10)) + 1;

    if (max_base_pos_width < 4) {
      max_base_pos_width = 4;
    }

    repaintCanvas ();
  }

  /**
   *  Remove this component from all the listener lists it is on.
   **/
  void stopListening () {
    getSelection ().removeSelectionChangeListener (this);

    getEntryGroup ().removeEntryGroupChangeListener (this);
    getEntryGroup ().removeEntryChangeListener (this);
    getEntryGroup ().removeFeatureChangeListener (this);
  }

  /**
   *  Returns the value of a flag that indicates whether this component can be
   *  traversed using Tab or Shift-Tab keyboard focus traversal - returns true
   *  for FeatureDisplay components
   **/
// tjc - deprecated replaced by isFocusable()
//public boolean isFocusTraversable () 
//{
//  return true;
//}

  /**
   *  Set value of the show correlation scores flag.
   *  @param show_correlation_scores Show correlation scores in the list if
   *    and only if this argument is true.
   **/
  public void setCorrelationScores (final boolean show_correlation_scores) {
    if (this.show_correlation_scores != show_correlation_scores) {
      this.show_correlation_scores = show_correlation_scores;
      repaintCanvas ();
    } else {
      // do nothing
    }
  }

  /**
   *  Get the value of the "show correlation scores" flag.
   **/
  public boolean getCorrelationScores () {
    return show_correlation_scores;
  }

  /**
   *  Set value of the show /gene flag.
   *  @param show_gene_names If true this component will show the /gene (really
   *    Feature.getIDString ()) instead of the key.
   **/
  public void setShowGenes (final boolean show_gene_names) {
    if (this.show_gene_names != show_gene_names) {
      this.show_gene_names = show_gene_names;
      repaintCanvas ();
    } else {
      // do nothing
    }
  }

  /**
   *  Get the value of the "show genes" flag.
   **/
  public boolean getShowGenes () {
    return show_gene_names;
  }

  /**
   *  Set value of the show qualifiers flag.
   *  @param show_quailfiers If true this component will show all the
   *    qualifiers after the note.
   **/
  public void setShowQualifiers (final boolean show_qualifiers) {
    if (this.show_qualifiers != show_qualifiers) {
      this.show_qualifiers = show_qualifiers;
      repaintCanvas ();
    } else {
      // do nothing
    }
  }

  /**
   *  Get the value of the "show qualifiers" flag.
   **/
  public boolean getShowQualifiers () {
    return show_qualifiers;
  }

  /**
   *  Set value of the show /product flag.
   *  @param show_products If true this component will show the /product
   *    qualifier instead of the /note.
   **/
  public void setShowProducts (final boolean show_products) {
    if (this.show_products != show_products) {
      this.show_products = show_products;
      repaintCanvas ();
    } else {
      // do nothing
    }
  }

  /**
   *  Get the value of the "show products" flag.
   **/
  public boolean getShowProducts () {
    return show_products;
  }

  /**
   *  Implementation of the EntryGroupChangeListener interface.  We listen to
   *  EntryGroupChange events so that we can update the display if entries
   *  are added or deleted.
   **/
  public void entryGroupChanged (EntryGroupChangeEvent event) {
    switch (event.getType ()) {
    case EntryGroupChangeEvent.ENTRY_ADDED:
    case EntryGroupChangeEvent.ENTRY_ACTIVE:
    case EntryGroupChangeEvent.ENTRY_DELETED:
    case EntryGroupChangeEvent.ENTRY_INACTIVE:
      repaintCanvas ();
      break;
    }
  }

  /**
   *  Implementation of the FeatureChangeListener interface.
   **/
  public void featureChanged (FeatureChangeEvent event) {
    if (!isVisible ()) {
      return;
    }

    repaintCanvas ();
  }

  /**
   *  Implementation of the EntryChangeListener interface.  We listen to
   *  EntryChange events so that we can update the list if features are added
   *  or deleted.
   **/
  public void entryChanged (EntryChangeEvent event) {
    if (!isVisible ()) {
      return;
    }

    repaintCanvas ();
  }

  /**
   *  Implementation of the SelectionChangeListener interface.  We listen to
   *  SelectionChange events so that we can update the list to reflect the
   *  current selection.
   **/
  public void selectionChanged (SelectionChangeEvent event) {
    if (!isVisible ()) {
      return;
    }

    if (event.getSource () == this) {
      // don't bother with events we sent ourself
      return;
    }

    if (getSelection ().getMarkerRange () != null &&
        event.getType () == SelectionChangeEvent.OBJECT_CHANGED) {
      // if the selected range changes we don't care
      return;
    }

    selection_changed_flag = true;

    repaintCanvas ();
  }

  /**
   *  Return a vector containing the text that is shown in the list - one
   *  String per line.
   **/
  public StringVector getListStrings () {
    final StringVector return_vector = new StringVector ();

    final FeatureEnumeration test_enumerator = getEntryGroup ().features ();

    while (test_enumerator.hasMoreFeatures ()) {
      final Feature this_feature = test_enumerator.nextFeature ();

      return_vector.add (makeFeatureString (this_feature, true));
    }

    return return_vector;
  }

  /**
   *  Create the scroll bar.
   **/
  private void createScrollbars () {
    scrollbar = new JScrollBar (Scrollbar.VERTICAL);
    scrollbar.setValues (0, 1, 0,
                         getEntryGroup ().getAllFeaturesCount ());
    scrollbar.setUnitIncrement (1);
    scrollbar.setBlockIncrement (1);
    scrollbar.addAdjustmentListener (new AdjustmentListener () {
      public void adjustmentValueChanged(AdjustmentEvent e) {
        setFirstIndex (e.getValue ());
      }
    });

    horiz_scrollbar = new JScrollBar (Scrollbar.HORIZONTAL);
    horiz_scrollbar.setValues (0, 1, 0,
                               getEntryGroup ().getAllFeaturesCount () - 1);
    horiz_scrollbar.setUnitIncrement (getFontWidth ());
    horiz_scrollbar.setBlockIncrement (100);
    horiz_scrollbar.addAdjustmentListener (new AdjustmentListener () {
      public void adjustmentValueChanged (AdjustmentEvent e) {
        repaintCanvas ();
      }
    });

    getMidPanel ().add (horiz_scrollbar, "South");

    add (scrollbar, "East");
  }

  /**
   *  Set the extent, max and value of the horizontal scrollbar 
   **/
  private void fixHorizScrollBar (final int max_width) {
    int old_value = horiz_scrollbar.getValue ();

    horiz_scrollbar.setValues (horiz_scrollbar.getValue (),
                               getCanvas ().getSize ().width,
                               0, max_width * getFontWidth ());
    horiz_scrollbar.setBlockIncrement (getCanvas ().getSize ().width);
  }

  /**
   *  Set the first visible index.
   **/
  public void setFirstIndex (final int first_index) {
    this.first_index = first_index;
    need_to_fix_horiz_scrollbar = true;
    repaintCanvas ();
  }

  /**
   *  Handle a mouse press event on the drawing canvas - select on click,
   *  select and broadcast it on double click.
   **/
  private void handleCanvasMousePress (final MouseEvent event) {
    if (event.getID() != MouseEvent.MOUSE_PRESSED) {
      return;
    }

    getCanvas ().requestFocus ();

    if (!event.isShiftDown ()) {
      getSelection ().clear ();
    }

    final int clicked_feature_index =
      scrollbar.getValue () + event.getY () / getLineHeight ();

    if (clicked_feature_index < getEntryGroup ().getAllFeaturesCount ()) {
      final FeatureVector selected_features =
        getSelection ().getAllFeatures ();

      final Feature clicked_feature =
        getEntryGroup ().featureAt (clicked_feature_index);

      if (selected_features.contains (clicked_feature)) {
        getSelection ().remove (clicked_feature);
        getSelection ().removeSegmentsOf (clicked_feature);
      } else {
        getSelection ().add (clicked_feature);
      }

      if (event.getClickCount () == 2) {
        makeSelectionVisible ();

        if ((event.getModifiers () & InputEvent.BUTTON2_MASK) != 0 ||
            event.isAltDown ()) {
          if (Options.readWritePossible ()) {
            new FeatureEdit (clicked_feature, getEntryGroup (),
                             getSelection (),
                             getGotoEventSource ()).show ();
          }
        }
      }
    }
  }

  /**
   *  The main paint function for the canvas.  An off screen image used for
   *  double buffering when drawing the canvas.
   *  @param g The Graphics object of the canvas.
   **/
  protected void paintCanvas (Graphics g) {
    refresh ();
    
    if (!isVisible ()) {
      return;
    }

    if (selection_changed_flag) {
      selection_changed_flag = false;

      final FeatureVector selected_features =
        getSelection ().getAllFeatures ();

      if (selected_features.size () > 0) {
        // set to true if any of the selected features is visible
        boolean a_selected_feature_is_visible = false;

        int first_line_in_view = scrollbar.getValue ();

        if (first_line_in_view == -1) {
          first_line_in_view = 0;
        }

        final int feature_count = getEntryGroup ().getAllFeaturesCount ();

        for (int i = first_line_in_view ;
             i < feature_count && i < first_line_in_view + linesInView () ;
             ++i) {
          final Feature this_feature = getEntryGroup ().featureAt (i);
          if (selected_features.contains (this_feature)) {
            a_selected_feature_is_visible = true;
            break;
          }
        }

        if (!a_selected_feature_is_visible) {
          // make the first selected feature visible
          final Feature first_selected_feature =
            selected_features.elementAt (0);

          final int index_of_first_selected_feature =
            getEntryGroup ().indexOf (first_selected_feature);

          if (index_of_first_selected_feature < scrollbar.getValue () ||
              index_of_first_selected_feature >=
              scrollbar.getValue () + linesInView ()) {

            scrollbar.setValue (index_of_first_selected_feature);
          }
        }
      }
    }

    g.setColor (background_colour);

    g.fillRect (0, 0, getCanvasWidth (), getCanvasHeight ());

    g.setColor (Color.black);

    final int all_feature_count = getEntryGroup ().getAllFeaturesCount ();
    
    if (all_feature_count == 0) {
      fixHorizScrollBar (0);
    } else {
      final int lines_in_view = linesInView ();
      int first_index_in_view = scrollbar.getValue ();

      if (first_index_in_view == -1) {
        first_index_in_view = 0;
      }

      final int feature_count = getEntryGroup ().getAllFeaturesCount ();

      /* jikes 1.15 bug  final */  int last_index_in_view;

      if (lines_in_view < feature_count - first_index_in_view) {
        last_index_in_view = first_index_in_view + lines_in_view;
      } else {
        last_index_in_view = feature_count - 1;
      }

      final FeatureVector features_in_view =
        getEntryGroup ().getFeaturesInIndexRange (first_index_in_view,
                                                  last_index_in_view);

      /**
       *  The maximum width of the strings we have seen - used to set the
       *  horiz_scrollbar maximum.
       **/
      int max_width = -1;

      for (int i = 0 ; i <= last_index_in_view - first_index_in_view ; ++i) {
        final Feature this_feature = features_in_view.elementAt (i);
        final String feature_string = makeFeatureString (this_feature, false);
        drawFeatureLine (g, this_feature, feature_string,i);

        if (feature_string.length () > max_width) {
          max_width = feature_string.length ();
        }
      }

      if (need_to_fix_horiz_scrollbar) {
        fixHorizScrollBar (max_width);
      }
    }
  }

  /**
   *  Return the number of visible text lines on canvas.
   **/
  private int linesInView () {
    return getCanvas ().getSize ().height / getLineHeight ();
  }

  /**
   *  Update the scrollbar.
   **/
  private void refresh () {
    scrollbar.setMaximum (getEntryGroup ().getAllFeaturesCount ());
    final int lines_in_view = linesInView ();
    scrollbar.setBlockIncrement (lines_in_view > 0 ? lines_in_view : 1);
    scrollbar.setUnitIncrement (1);
    scrollbar.setVisibleAmount (linesInView ());
  }

  /**
   *  Draw the given Feature at the given line of the list, taking the
   *  selection into account.
   **/
  private void drawFeatureLine (final Graphics g,
                                final Feature feature,
                                final String feature_string,
                                final int line) {
    final int y_pos = line * getLineHeight ();

    // the width of the coloured blob at the left of the text
    final int BOX_WIDTH = getLineHeight ();

    final Color feature_colour = feature.getColour ();

    if (feature_colour == null) {
      // default colour is white
      g.setColor (Color.white);
    } else {
      g.setColor (feature_colour);
    }
    g.fillRect (1 - horiz_scrollbar.getValue (), y_pos + 1,
                BOX_WIDTH, getLineHeight () - 1);

    if (getSelection ().contains (feature)) {
      // draw in reverse
      g.setColor (Color.black);
      g.fillRect (BOX_WIDTH + 4 - horiz_scrollbar.getValue (), y_pos,
                  getCanvas ().getSize ().width + horiz_scrollbar.getValue (),
                  getLineHeight ());
      g.setColor (background_colour);
    } else {
      g.setColor (Color.black);
    }

    g.setFont (getFont ());

    g.drawString (feature_string,
                  BOX_WIDTH + 5 -
                  horiz_scrollbar.getValue (),
                  y_pos + getFontAscent ());

    g.setPaintMode ();
  }

  /**
   *  Return the list index of a feature.
   **/
  private int indexOf (Feature feature) {
    return getEntryGroup ().indexOf (feature);
  }

  /**
   *  Return a String object suitable for displaying in the list of features.
   *  @param dont_truncate if true the gene name / key field won't be
   *    truncated if it is longer than the field width
   **/
  private String makeFeatureString (final Feature feature,
                                    final boolean dont_truncate) {
    String key_string;

    final int KEY_FIELD_WIDTH = 15;

    if (show_gene_names) {
      key_string = feature.getIDString ();

      if (key_string.length () > KEY_FIELD_WIDTH && !dont_truncate) {
        key_string = key_string.substring (0, KEY_FIELD_WIDTH);
      }
    } else {
      key_string = feature.getKey ().toString ();
    }

    final Marker low_marker = feature.getFirstBaseMarker ();
    final Marker high_marker = feature.getLastBaseMarker ();

    final StringBuffer description_string_buffer = new StringBuffer ();

    if (show_products) {
      final String product_string = feature.getProductString ();

      if (product_string == null) {
        if (feature.isCDS ()) {
          description_string_buffer.append ("[no /product]");
        } else {
          // description is blank
        }
      } else {
        description_string_buffer.append (product_string);
      }
    } else {
      final String note = feature.getNote ();

      if (note != null && note.length () != 0) {
        final int QUALIFIER_COLUMN = 10;

        final String note_string =
          padRightWithSpaces (feature.getNote (), QUALIFIER_COLUMN);

        description_string_buffer.append (note_string);
        description_string_buffer.append ("   ");
      }

      if (show_qualifiers) {
        description_string_buffer.append (getQualifierString (feature));
      }
    }

    final String low_pos;
    final String high_pos;

    if (low_marker == null || high_marker == null) {
      low_pos = "unknown";
      high_pos = "unknown";
    } else {
      if (low_marker.getRawPosition () < high_marker.getRawPosition ()) {
        low_pos = String.valueOf (low_marker.getRawPosition ());
        high_pos = String.valueOf (high_marker.getRawPosition ());
      } else {
        low_pos = String.valueOf (high_marker.getRawPosition ());
        high_pos = String.valueOf (low_marker.getRawPosition ());
      }
    }

    StringBuffer new_list_line = new StringBuffer ();

    new_list_line.append (padRightWithSpaces (key_string, KEY_FIELD_WIDTH));
    new_list_line.append (" ");

    new_list_line.append (padLeftWithSpaces (low_pos, max_base_pos_width));
    new_list_line.append (" ");
    new_list_line.append (padLeftWithSpaces (high_pos, max_base_pos_width));
    new_list_line.append (" ");

    if (feature.isForwardFeature ()) {
      new_list_line.append ("   ");
    } else {
      new_list_line.append ("c  ");
    }

    if (show_correlation_scores) {
      if (feature.isCDS ()) {
        new_list_line.append (getScoresString (feature));
        new_list_line.append ("  ");
      } else {
        new_list_line.append ("                         ");
        if (getBasePlotGroup ().getCodonUsageAlgorithm () != null) {
          new_list_line.append ("      ");
        }
      }
    }

    new_list_line.append (description_string_buffer);

    return new_list_line.toString ();
  }

  /**
   *  Return the characters width of the Canvas.
   **/
  private int getWidthInChars () {
    return getCanvas ().getSize ().width / getFontWidth ();
  }

  /**
   *  Return a String containing the given Qualifier and it's values (in EMBL
   *  format).
   *  @param start_index ignore the values before this index
   **/
  private String formatQualifier (final String qualifier_name,
                                  final Feature feature,
                                  final int start_index) {
    final StringBuffer buffer = new StringBuffer ();

    try {
      final Qualifier qualifier = feature.getQualifierByName (qualifier_name);

      if (qualifier != null) {
        final EntryInformation entry_information =
          feature.getEntry ().getEntryInformation ();

        final QualifierInfo qualifier_info =
          entry_information.getQualifierInfo (qualifier_name);

        final StringVector qualifier_strings =
          StreamQualifier.toStringVector (qualifier_info,
                                          qualifier);

        for (int i = start_index ; i < qualifier_strings.size () ; ++i) {

          final String qualifier_string = qualifier_strings.elementAt (i);

          buffer.append (qualifier_string + " ");
        }
      }
    } catch (InvalidRelationException e) {
      // ignore
    }

    return buffer.toString ();
  }

  /**
   *  Return a String containing all the qualifiers of the given Feature
   *  (except /note) in EMBL format.  Any /similarity qualifier will come
   *  first.
   **/
  private String getQualifierString (final Feature feature) {
    final StringBuffer buffer = new StringBuffer ();

    final QualifierVector qualifiers = feature.getQualifiers ();

    // if there is a /note and it has more than one value put it next (without
    // the first value)
    final Qualifier note_qualifier =
      qualifiers.getQualifierByName ("note");

    if (note_qualifier != null && note_qualifier.getValues ().size () > 1) {
      buffer.append (formatQualifier ("note", feature, 1));
      buffer.append (" ");
    }

    // put /similarity before all but the /note qualifier
    final Qualifier similarity_qualifier =
      qualifiers.getQualifierByName ("similarity");

    if (similarity_qualifier != null) {
      buffer.append (formatQualifier ("similarity", feature, 0));
      buffer.append (" ");
    }

    for (int i = 0 ; i < qualifiers.size () ; ++i) {
      final Qualifier this_qualifier = qualifiers.elementAt (i);

      final String this_qualifier_name = this_qualifier.getName ();

      if (!this_qualifier_name.equals ("note") &&
          !this_qualifier_name.equals ("similarity")) {
        buffer.append (formatQualifier (this_qualifier_name, feature, 0));
        buffer.append (" ");
      }
    }

    return buffer.toString ();
  }

  /**
   *  Return a String containing the correlation scores.
   **/
  public String getScoresString (final Feature feature) {
    final int base_total = feature.getTranslationBases ().length ();

    final int c_total = feature.getBaseCount (Bases.getIndexOfBase ('c'));
    final int g_total = feature.getBaseCount (Bases.getIndexOfBase ('g'));

    final int g1_count =
      feature.getPositionalBaseCount (0, Bases.getIndexOfBase ('g'));

    final int c3_count =
      feature.getPositionalBaseCount (2, Bases.getIndexOfBase ('c'));
    final int g3_count =
      feature.getPositionalBaseCount (2, Bases.getIndexOfBase ('g'));

    final double c3_score = 100.0 * (3 * c3_count - c_total) / c_total;
    final double g1_score = 100.0 * (3 * g1_count - g_total) / g_total;
    final double g3_score = 100.0 * (3 * g3_count - g_total) / g_total;

    final double cor1_2_score = feature.get12CorrelationScore ();

    final NumberFormat number_format = NumberFormat.getNumberInstance ();

    number_format.setMaximumFractionDigits (1);
    number_format.setMinimumFractionDigits (1);

    final String cor1_2_score_string = number_format.format (cor1_2_score);
    final String c3_score_string;
    final String g1_score_string;
    final String g3_score_string;


    if (c_total == 0) {
      c3_score_string = "ALL";
    } else {
      c3_score_string = number_format.format (c3_score);
    }

    if (g_total == 0) {
      g1_score_string = "ALL";
    } else {
      g1_score_string = number_format.format (g1_score);
    }

    if (g_total == 0) {
      g3_score_string = "ALL";
    } else {
      g3_score_string = number_format.format (g3_score);
    }

    String codon_usage_score_string = "";

    final CodonUsageAlgorithm codon_usage_alg =
      getBasePlotGroup ().getCodonUsageAlgorithm ();

    if (codon_usage_alg != null) {
      number_format.setMaximumFractionDigits (3);
      number_format.setMinimumFractionDigits (3);

      codon_usage_score_string =
        number_format.format (codon_usage_alg.getFeatureScore (feature)) + " ";
    }

    return
      codon_usage_score_string +
      padRightWithSpaces (cor1_2_score_string, 5) + " " +
      padRightWithSpaces (c3_score_string, 5) + " " +
      padRightWithSpaces (g1_score_string, 5) + " " +
      padRightWithSpaces (g3_score_string, 5);
  }

  /**
   *  Return the given string padded with spaces to the given width.  The
   *  spaces are added on the right of the string.
   **/
  private String padRightWithSpaces (final String string, final int width) {
    if (string.length () == width) {
      return string;
    }

    final StringBuffer buffer = new StringBuffer (string);

    for (int i = 0 ; i < width - string.length () ; ++i) {
      buffer.append (' ');
    }

    return buffer.toString ();
  }

  /**
   *  Return the given string padded with spaces to the given width.  The
   *  spaces are added on the left of the string.
   **/
  private String padLeftWithSpaces (final String string, final int width) {
    if (string.length () == width) {
      return string;
    }

    final StringBuffer buffer = new StringBuffer ();

    for (int i = 0 ; i < width - string.length () ; ++i) {
      buffer.append (' ');
    }

    buffer.append (string);

    return buffer.toString ();
  }

  /**
   *  Return the height each line of the display should be.  Each feature will
   *  be drawn into one line.
   **/
  private int getLineHeight () {
    return getFontAscent () + 2;
  }

  /**
   *  This variable is true if correlation scores should be shown in the list.
   **/
  private boolean show_correlation_scores = false;

  /**
   *  The JScrollBar for this FeatureList object.
   **/
  private JScrollBar scrollbar = null;

  /**
   *  The JScrollBar for horizontal scrolling in this FeatureList object.
   **/
  private JScrollBar horiz_scrollbar = null;

  /**
   *  The index of the first visible feature in the list.
   **/
  private int first_index;

  /**
   *  This is set to true by selectionChanged () and used by paintCanvas ().
   **/
  private boolean selection_changed_flag = false;

  /**
   *  The colour used to draw the background.
   **/
  private Color background_colour = Color.white;

  /**
   *  If true this component will show Feature.getIDString () (ie /gene or
   *  /label) instead of the key.
   **/
  private boolean show_gene_names = false;

  /**
   *  If true this component will show the /product qualifier instead of the
   *  /note field.
   **/
  private boolean show_products = false;

  /**
   * If true this component will show all the qualifiers after the note.
   **/
  private boolean show_qualifiers = false;

  /**
   *  The is the maximum width of the strings containing the feature start and
   *  stop positions.  Set in the constructor.
   **/
  private int max_base_pos_width;

  /**
   *  Set to true when paintCanvas() needs to call fixHorizScrollBar().
   **/
  private boolean need_to_fix_horiz_scrollbar = true;
}
