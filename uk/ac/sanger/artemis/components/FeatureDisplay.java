/* FeatureDisplay.java
 *
 * created: Fri Oct  9 1998
 *
 * This file is part of Artemis
 *
 * Copyright(C) 1998,1999,2000,2001,2002  Genome Research Limited
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or(at your option) any later version.
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/FeatureDisplay.java,v 1.8 2004-11-09 14:24:41 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.util.StringVector;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.sequence.*;

import java.awt.event.*;
import java.awt.*;
import java.lang.Math;
import java.util.Vector;
import java.util.Comparator;

import javax.swing.JScrollBar;
import javax.swing.JComponent;

/**
 *  This component is used for displaying an Entry.
 *
 *  @author Kim Rutherford
 *  @version $Id: FeatureDisplay.java,v 1.8 2004-11-09 14:24:41 tjc Exp $
 **/

public class FeatureDisplay extends EntryGroupPanel
  implements EntryGroupChangeListener,
             EntryChangeListener, FeatureChangeListener,
             SelectionChangeListener, GotoListener, SequenceChangeListener,
             DisplayComponent, OptionChangeListener, DisplayAdjustmentListener
{

  /** Key code for calling zoomToSelection(). */
  final static public int ZOOM_TO_SELECTION_KEY = KeyEvent.VK_Z;

  final static public int SCROLLBAR_AT_TOP = 1;
  final static public int SCROLLBAR_AT_BOTTOM = 2;
  final static public int NO_SCROLLBAR = 3;

  public final static int FORWARD = Bases.FORWARD;
  public final static int REVERSE = Bases.REVERSE;

  public final static int NO_FRAME        = FeatureSegment.NO_FRAME;
  public final static int FORWARD_STRAND  = FeatureSegment.FORWARD_STRAND;
  public final static int REVERSE_STRAND  = FeatureSegment.REVERSE_STRAND;
  public final static int FORWARD_FRAME_1 = FeatureSegment.FORWARD_FRAME_1;
  public final static int FORWARD_FRAME_2 = FeatureSegment.FORWARD_FRAME_2;
  public final static int FORWARD_FRAME_3 = FeatureSegment.FORWARD_FRAME_3;
  public final static int REVERSE_FRAME_3 = FeatureSegment.REVERSE_FRAME_3;
  public final static int REVERSE_FRAME_2 = FeatureSegment.REVERSE_FRAME_2;
  public final static int REVERSE_FRAME_1 = FeatureSegment.REVERSE_FRAME_1;
  public final static int SCALE_LINE      = FeatureSegment.SCALE_LINE;

  /**
   *  The JScrollBar for this FeatureDisplay object.  We create the scrollbar
   *  as part of this object rather than in the EntryEdit component because we
   *  may need to change the parameters of the scrollbar later.
   **/
  private JScrollBar scrollbar = null;

  /** A scroll bar for changing the viewing scale. */
  private JScrollBar scale_changer = null;

  /** Used to colour the frames. */
  private Color light_grey = new Color(240, 240, 240);

  /** Used to colour sequence line. */
  private Color not_so_light_grey = new Color(200, 200, 200);

  /**
   *  The colour used for the active entry line when 
   *  one_line_per_entry is set.
   **/
  private Color active_entry_colour = new Color(255, 255, 140);

  /**
   *  This Vector containing the references of those features that are
   *  currently visible.
   **/
  private FeatureVector visible_features = new FeatureVector();

  /**
   *  If true updateVisibleFeatureVector() will be called by paint().
   *  updateVisibleFeatureVector() sets this to false,
   *  needVisibleFeatureVectorUpdate() sets this to true.
   **/
  private boolean update_visible_features = true;

  /** Contains those objects listening for adjustment events. */
  final private Vector adjustment_listener_list = new Vector();

  /**
   *  The index of the first base that we are displaying.  
   *  Can be negative if hard_left_edge is true.
   **/
  private int left_edge_base = 1;

  /** See getScaleFactor(). */
  private int scale_factor = 3;

  /** See getScaleFactor(). */
  private float scale_value = 1;

  /** true if labels should be shown. */
  private boolean show_labels = true;

  /**
   *  This variable is true if the forward frame lines(or forward entry lines
   *  - see one_line_per_entry) should be drawn.
   **/
  private boolean show_forward_lines = true;

  /**
   *  This variable is true if the reverse frame lines(or reverse entry lines
   *  - see one_line_per_entry) should be drawn.
   **/
  private boolean show_reverse_lines = true;

  /**
   *  If true draw all features, sequence and scale lines reverse complemented.
   **/
  private boolean rev_comp_display = false;

  /**
   *  This variable is true if(for each strand) each entry should be on a
   *  separate line.
   **/
  private boolean one_line_per_entry = false;

  /**
   *  If true the there will never be a gap between the left edge of the
   *  screen and the first visible base.
   **/
  private boolean hard_left_edge = true;

  /** true if source features should be shown. */
  private boolean show_source_features = false;

  /** true if stop codons should be shown.     */
  private boolean show_stop_codons = true;

  /** true if start codons should be shown.    */
  private boolean show_start_codons = false;

  /** true if directional arrows should be shown on features. */
  private boolean show_feature_arrows;

  /** true a black border will be drawn around each feature.  */
  private boolean show_feature_borders;

  /**
   *  This variable is true if each base should be drawn in a different colour
   *  at scale feature 1.
   **/
  private boolean show_base_colours = false;

  /**
   *  All features(not just CDS features) will be drawn on the frame lines if
   *  and only if this variable is true.  See setFrameFeaturesFlag() and
   *  getFrameFeaturesFlag().
   **/
  private boolean frame_features_flag = false;

  /**
   *  The position(s) of the last mouse click on the dna line.  The
   *  MarkerRange contains a reference to the appropriate Strand and contains
   *  the base positions.  See getMarkerRangeFromPosition() to understand why
   *  one click can give multiple bases.
   **/
  private MarkerRange click_range = null;

  /**
   *  The last(FeatureSegment) Marker that the user clicked on.  This is used
   *  for dragging the ends of segments.
   **/
  private Marker click_segment_marker = null;

  /**
   *  This is true if click_segment_marker is the Marker at the start of
   *  segment false otherwise.  The value is only useful if
   *  click_segment_marker is set.
   **/
  private boolean click_segment_marker_is_start_marker = false;

  /**
   *  When a FeatureSegment Marker drag starts, this is set to the Marker at
   *  the other end of the segment.  This is used to check that the drag has
   *  not move the Marker too far(past the end of the segment).
   **/
  private Marker other_end_of_segment_marker = null;

  /**
   *  Features with a /score qualifier less than this value will not be shown.
   **/
  private int current_min_score = 0;

  /**
   *  Features with a /score qualifier greater than this value will not be
   *  shown.
   **/
  private int current_max_score = 100;

  private MouseEvent last_mouse_press_event;

  /**
   *  If set no DisplayAdjustment events will be sent.  This is set by
   *  displayAdjustmentValueChanged() to prevent an event we send from
   *  returning to us(a FeatureDisplay can listen for DisplayAdjustment
   *  events from another FeatureDisplay).
   **/
  private boolean disable_display_events = false;

  /**
   *  Set to true by selectionChanged() to tell updateVisibleFeatureVector()
   *  to raise the contents of the select before updating.
   **/
  private boolean raise_selection_flag = false;
 
  /** the minimum distance in pixels between the labels. */
  private final static int MINIMUM_LABEL_SPACING = 80;
 
  /** colour used for A. */
  private static Color dark_green = new Color(0, 150, 0);

  /**
   *  Used by drawOneLetter() - declared here so that we we don't need to
   *  allocated a whole array each time drawOneLetter() is called.
   **/
  private final char [] draw_one_char_temp_array = new char [1];

  private int scrollbar_style;

  /**
   *  Create a new FeatureDisplay object with the horizontal scrollbar at the
   *  bottom of the component.
   *  @param entry_group The EntryGroup that this component will display.
   *  @param selection The Selection object for this component.  Selected
   *    objects will be highlighted.
   *  @param goto_event_source The object to use when we need to call
   *    gotoBase().
   *  @param base_plot_group The BasePlotGroup associated with this JMenu -
   *    needed to call getCodonUsageAlgorithm()
   **/
  public FeatureDisplay(final EntryGroup entry_group,
                        final Selection selection,
                        final GotoEventSource goto_event_source,
                        final BasePlotGroup base_plot_group) 
  {
    this(entry_group, selection, goto_event_source,
          base_plot_group, SCROLLBAR_AT_BOTTOM);
  }

  /**
   *  Create a new FeatureDisplay object.
   *  @param entry_group The EntryGroup that this component will display.
   *  @param owning_component The EntryEdit object that contains the selection
   *    that this component uses.
   *  @param scrollbar_at_top If true the horizontal scrollbar will be at the
   *    top of component.
   *  @param base_plot_group The BasePlotGroup associated with this JMenu -
   *    needed to call getCodonUsageAlgorithm()
   *  @param scrollbar_style Controls the type of horizontal scrollbar.  Must
   *    be one of SCROLLBAR_AT_TOP, SCROLLBAR_AT_BOTTOM or NO_SCROLLBAR.
   **/
  public FeatureDisplay(final EntryGroup entry_group,
                        final Selection selection,
                        final GotoEventSource goto_event_source,
                        final BasePlotGroup base_plot_group,
                        final int scrollbar_style) 
  {
    super(entry_group, selection, goto_event_source, base_plot_group);

    this.scrollbar_style = scrollbar_style;

    show_feature_arrows =
      Options.getOptions().getPropertyTruthValue("draw_feature_arrows");

    show_feature_borders =
      Options.getOptions().getPropertyTruthValue("draw_feature_borders");

    frame_features_flag =
      Options.getOptions().getPropertyTruthValue("features_on_frame_lines");

    one_line_per_entry =
      Options.getOptions().getPropertyTruthValue("one_line_per_entry");

    show_labels =
      Options.getOptions().getPropertyTruthValue("feature_labels");

    addComponentListener(new ComponentAdapter()
    {
      public void componentResized(ComponentEvent e) 
      {
        // update the scroll bar as soon as we know the size of the canvas
        fixScrollbar();
        needVisibleFeatureVectorUpdate();
        fireAdjustmentEvent(DisplayAdjustmentEvent.ALL_CHANGE_ADJUST_EVENT);
      }

      public void componentShown(ComponentEvent e) 
      {
        // update the scroll bar as soon as we know the size of the canvas
        fixScrollbar();
        needVisibleFeatureVectorUpdate();
        fireAdjustmentEvent(DisplayAdjustmentEvent.ALL_CHANGE_ADJUST_EVENT);
      }
    });

    setScaleValue();

    if(scrollbar_style == SCROLLBAR_AT_TOP) 
      createScrollbar(true);
    else if(scrollbar_style == SCROLLBAR_AT_BOTTOM)
      createScrollbar(false);

    createScaleScrollbar();
//  fixCanvasSize();
//  fixScrollbar();
    addListeners();

    needVisibleFeatureVectorUpdate();

//  fireAdjustmentEvent(DisplayAdjustmentEvent.ALL_CHANGE_ADJUST_EVENT);

    getSelection().addSelectionChangeListener(this);
    getGotoEventSource().addGotoListener(this);

    getEntryGroup().addEntryGroupChangeListener(this);
    getEntryGroup().addEntryChangeListener(this);
    getEntryGroup().addFeatureChangeListener(this);

    getBases().addSequenceChangeListener(this, Bases.MIN_PRIORITY);

    Options.getOptions().addOptionChangeListener(this);
    setBackground(Color.white);
  }

  /**
   *  Returns the value of a flag that indicates whether this component can be
   *  traversed using Tab or Shift-Tab keyboard focus traversal - returns true
   *  for FeatureDisplay components
   **/
// tjc - deprecated replaced by isFocusable().
//public boolean isFocusTraversable()
//{
//  return true;
//}

  /**
   *  Overriden to call fixCanvasSize()
   **/
  public void setVisible(final boolean visible) 
  {
    super.setVisible(visible);
    fixCanvasSize();
  }

  /**
   *  Set value of the show label flag.
   *  @param show_label Show labels if and only if this argument is true.
   **/
  public void setShowLabels(boolean show_labels) 
  {
    if(this.show_labels != show_labels) 
    {
      this.show_labels = show_labels;
      fixCanvasSize();
    } 
  }

  /**
   *  Get the value of the "show label" flag.
   **/
  public boolean getShowLabels() 
  {
    return show_labels;
  }

  /**
   *  Set value of the "show forward frame lines" flag.
   *  @param show_forward_lines Show forward frame lines if and only if
   *    this argument is true.
   **/
  public void setShowForwardFrameLines(boolean show_forward_lines) 
  {
    if(this.show_forward_lines != show_forward_lines) 
    {
      this.show_forward_lines = show_forward_lines;
      fixCanvasSize();
    } 
  }

  /**
   *  Get the value of the "show forward frame lines" flag.
   **/
  public boolean getShowForwardFrameLines() 
  {
    return show_forward_lines;
  }

  /**
   *  Set value of the "show reverse frame lines" flag.
   *  @param show_reverse_lines Show frame lines if and only if this
   *    argument is true.
   **/
  public void setShowReverseFrameLines(boolean show_reverse_lines) 
  {
    if(this.show_reverse_lines != show_reverse_lines) 
    {
      this.show_reverse_lines = show_reverse_lines;
      fixCanvasSize();
    }
  }

  /**
   *  Get the value of the "show source features" flag.
   **/
  public boolean getShowSourceFeatures() 
  {
    return show_source_features;
  }

  /**
   *  Set value of the "show source features" flag.
   *  @param show_source_features Show features with a "source" key if and
   *    only if this argument is true.
   **/
  public void setShowSourceFeatures(boolean show_source_features) 
  {
    if(this.show_source_features != show_source_features) 
    {
      this.show_source_features = show_source_features;
      needVisibleFeatureVectorUpdate();
      repaint();
    } 
  }

  /**
   *  Get the value of the "show frame lines" flag.
   **/
  public boolean getShowReverseFrameLines()
  {
    return show_reverse_lines;
  }

  /**
   *  Set value of the show base colours flag.
   *  @param show_base_colours At scale_factor less than two show each base in
   *    a different colour if and only if this argument is true.
   **/
  public void setShowBaseColours(boolean show_base_colours) 
  {
    if(this.show_base_colours != show_base_colours) 
    {
      this.show_base_colours = show_base_colours;
      if(getScaleFactor() > 1) 
        setScaleFactor(1);
      repaint();
    } 
  }

  /**
   *  Get the value of the "show base colours" flag.
   **/
  public boolean getShowBaseColours() 
  {
    return show_base_colours;
  }

  /**
   *  Set value of the "one line per entry" flag.
   *  @param one_line_per_entry If true then each entry will be shown on a
   *    different line, instead of showing frame lines.
   **/
  public void setOneLinePerEntry(final boolean one_line_per_entry)
  {
    if(this.one_line_per_entry != one_line_per_entry) 
    {
      this.one_line_per_entry = one_line_per_entry;
      fixCanvasSize();
    }
  }

  /**
   *  Get the value of the "one line per entry" flag.
   **/
  public boolean getOneLinePerEntryFlag() 
  {
    return one_line_per_entry;
  }

  /**
   *  Set value of the "hard left edge" flag.
   *  @param hard_left_edge If true the there will never be a gap between the
   *    left edge of the screen and the first visible base.  If false base one
   *    can be moved to the centre of the display.
   **/
  public void setHardLeftEdge(final boolean hard_left_edge) 
  {
    if(this.hard_left_edge != hard_left_edge) 
    {
      this.hard_left_edge = hard_left_edge;
      if(hard_left_edge && getForwardBaseAtLeftEdge() < 1) 
        setFirstVisibleForwardBase(1);
      
      fixScrollbar();
    } 
  }

  /**
   *  Get the value of the "hard left edge" flag.
   **/
  public boolean getHardLeftEdgeFlag() 
  {
    return hard_left_edge;
  }

  /**
   *  Set value of the show stop codons flag.
   *  @param show_stop_codons Show stop codons if and only if this argument is
   *    true.
   **/
  public void setShowStopCodons(boolean show_stop_codons) 
  {
    if(this.show_stop_codons != show_stop_codons) 
    {
      this.show_stop_codons = show_stop_codons;
      repaint();
    } 
  }

  /**
   *  Return the value of the "show stop codons" flag.
   **/
  public boolean getShowStopCodons() 
  {
    return show_stop_codons;
  }

  /**
   *  Set value of the show start codons flag.
   *  @param show_start_codons Show start codons if and only if this argument
   *    is true.
   **/
  public void setShowStartCodons(boolean show_start_codons) 
  {
    if(this.show_start_codons != show_start_codons) 
    {
      this.show_start_codons = show_start_codons;
      repaint();
    } 
  }

  /**
   *  Return the value of the "show start codons" flag.
   **/
  public boolean getShowStartCodons() 
  {
    return show_start_codons;
  }

  /**
   *  Set value of the reverse complement display flag.
   *  @param show_start_codons Draw all features and sequence reverse
   *    complemented if and only if this argument is true.
   **/
  public void setRevCompDisplay(boolean rev_comp_display) 
  {
    if(this.rev_comp_display != rev_comp_display) 
    {
      this.rev_comp_display = rev_comp_display;
      int remember_position = getCentreForwardBase();

      // we want to keep the selection visible after the flip, so
      // that will override the centre position
      final Marker first_base_marker =
        getSelection().getStartBaseOfSelection();

      if(first_base_marker != null && baseVisible(first_base_marker)) 
        remember_position = first_base_marker.getRawPosition();

      final Marker last_base_marker =
        getSelection().getStartBaseOfSelection();

      if(last_base_marker != null && baseVisible(last_base_marker)) 
        remember_position = last_base_marker.getRawPosition();

      fireAdjustmentEvent(DisplayAdjustmentEvent.REV_COMP_EVENT);

      makeBaseVisibleInternal(remember_position, isRevCompDisplay(), true);

      needVisibleFeatureVectorUpdate();
      fixScrollbar();
      repaint();
    }
  }

  /**
   *  Return the value of the "reverse complement display" flag.
   **/
  public boolean isRevCompDisplay() 
  {
    return rev_comp_display;
  }

  /**
   *  Set value of the show feature arrows flag.
   *  @param show_feature_arrows Show directional arrows if and only if this
   *    argument is true.
   **/
  public void setShowFeatureArrows(boolean show_feature_arrows)
  {
    if(this.show_feature_arrows != show_feature_arrows)
    {
      this.show_feature_arrows = show_feature_arrows;
      repaint();
    } 
  }

  /**
   *  Return the value of the "show feature arrows" flag.
   **/
  public boolean getShowFeatureArrows() 
  {
    return show_feature_arrows;
  }

  /**
   *  Set value of the show feature borders flag.
   *  @param show_feature_borders Draw a border around each feature if and
   *    only if this argument is true.
   **/
  public void setShowFeatureBorders(boolean show_feature_borders) 
  {
    if(this.show_feature_borders != show_feature_borders) 
    {
      this.show_feature_borders = show_feature_borders;
      repaint();
    } 
  }

  /**
   *  Return the value of the "show feature borders" flag.
   **/
  public boolean getShowFeatureBorders() 
  {
    return show_feature_borders;
  }

  /**
   *  Set value of the show frame features flag.
   *  @param frame_features_flag All features(not just CDS features) will be
   *    drawn on the frame lines if and only if this argument is true.
   **/
  public void setFrameFeaturesFlag(boolean frame_features_flag) 
  {
    if(this.frame_features_flag != frame_features_flag) 
    {
      this.frame_features_flag = frame_features_flag;
      repaint();
    } 
  }

  /**
   *  Return the value of the "show frame features" flag.
   **/
  public boolean getFrameFeaturesFlag() 
  {
    return frame_features_flag;
  }

  /**
   *  Set the value of the minimum score for this FeatureDisplay - features
   *  that have a /score lower than this value are never shown.
   **/
  public void setMinimumScore(final int minimum_score) 
  {
    current_min_score = minimum_score;
    needVisibleFeatureVectorUpdate();
    repaint();
  }

  /**
   *  Return the value of the minimum score for this FeatureDisplay - see
   *  setMinimumScore().
   **/
  public int getMinimumScore() 
  {
    return current_min_score;
  }

  /**
   *  Set the value of the maximum score for this FeatureDisplay - features
   *  that have a /score higher than this value are never shown.
   **/
  public void setMaximumScore(final int maximum_score) 
  {
    current_max_score = maximum_score;
    needVisibleFeatureVectorUpdate();
    repaint();
  }

  /**
   *  Return the value of the maximum score for this FeatureDisplay - see
   *  setMaximumScore().
   **/
  public int getMaximumScore() 
  {
    return current_max_score;
  }

  /**
   *  Redraw this component.  This method is public so that other classes can
   *  force an update.  for example this is called when the options files is
   *  re-read.
   **/
  public void redisplay() 
  {
    repaint();
  }

  /**
   *  Adds the specified event adjustment listener to receive adjustment
   *  change events from this object.
   *  @param l the event change listener.
   **/
  public void addDisplayAdjustmentListener(DisplayAdjustmentListener l) 
  {
    adjustment_listener_list.addElement(l);
  }

  /**
   *  Removes the specified event listener so that it no longer receives
   *  adjustment change events from this object.
   *  @param l the event change listener.
   **/
  public void removeDisplayAdjustmentListener(DisplayAdjustmentListener l) 
  {
    adjustment_listener_list.removeElement(l);
  }

  /**
   *  Handle key press events.  This is static because making it non-static
   *  triggered a java.lang.VerifyError
   **/
  private static void handleKeyPress(final FeatureDisplay feature_display,
                                      final KeyEvent event) 
  {
    // this is done so that menu shortcuts don't cause each action to be
    // performed twice
    if(event.getModifiers() != 0) 
      return;

    switch(event.getKeyCode()) 
    {
      case ZOOM_TO_SELECTION_KEY:
        FeaturePopup.zoomToSelection(feature_display);
        break;
      default:
        break;
    }
  }

  /**
   *  Set the scale factor and update the display if the scale factor has
   *  changed.  A factor of zero means the full translation will be visible.
   *  At higher scale factors only stop codons are visible, and a bigger
   *  number will mean more bases are visible.
   **/
  public void setScaleFactor(int scale_factor) 
  {
    if(this.scale_factor != scale_factor) 
    {
      // we will try to keep the base in the centre of the view to stay where
      // it is, so we save it's position in remember_position.
      int remember_position = getCentreForwardBase();

      // if the first base is visible then keep it visible
      if(hard_left_edge && getFirstVisibleForwardBase() == 1) 
        remember_position = 1;

      if(!getSelection().isEmpty()) 
      {
        // but, we want to keep the selection visible after a scale change, so
        // that will override the centre position
        final Marker first_base_marker =
          getSelection().getStartBaseOfSelection();

        final int first_base_marker_raw_position =
                  first_base_marker.getRawPosition();
        final int first_base_marker_position;

        if(isRevCompDisplay()) 
          first_base_marker_position =
              getBases().getComplementPosition(first_base_marker_raw_position);
        else 
          first_base_marker_position = first_base_marker_raw_position;

        final Marker last_base_marker =
                              getSelection().getEndBaseOfSelection();

        final int last_base_marker_raw_position =
                                   last_base_marker.getRawPosition();

        final int last_base_marker_position;

        if(isRevCompDisplay()) 
          last_base_marker_position =
            getBases().getComplementPosition(last_base_marker_raw_position);
        else 
          last_base_marker_position = last_base_marker_raw_position;

        final int lowest_visible_base = getFirstVisibleForwardBase();
        final int highest_visible_base = getLastVisibleForwardBase();

        // first selected base or first visible base, whichever is greater
        int restricted_first_selected_base = lowest_visible_base;

        // last selected base or last visible base, whichever is smaller
        int restricted_last_selected_base = highest_visible_base;

        if(first_base_marker != null) 
        {
          if(first_base_marker_position > lowest_visible_base &&
              first_base_marker_position < highest_visible_base) 
            restricted_first_selected_base = first_base_marker_position;
        }

        if(last_base_marker != null) 
        {
          if(last_base_marker_position < highest_visible_base &&
              last_base_marker_position > lowest_visible_base) 
            restricted_last_selected_base = last_base_marker_position;
        }

        if(getSelection().getMarkerRange() == null) 
          remember_position = restricted_first_selected_base;
        else 
        {
          // keep the centre of the selection in the middle of the display if
          // a range of bases is selected
          remember_position = restricted_first_selected_base +
                             (restricted_last_selected_base -
                               restricted_first_selected_base) / 2;
        }
      }

      this.scale_factor = scale_factor;

      setScaleValue();
      scale_changer.setValue(scale_factor);
      setCentreVisibleForwardBase(remember_position);
      fixScrollbar();
      fireAdjustmentEvent(DisplayAdjustmentEvent.SCALE_ADJUST_EVENT);
      needVisibleFeatureVectorUpdate();

      repaint();
    }
  }

  /**
   *  Implementation of the FeatureChangeListener interface.  We listen to
   *  FeatureChange events so that we can update the display if qualifiers
   *  change.
   **/
  public void featureChanged(final FeatureChangeEvent event) 
  {
    final Feature event_feature = event.getFeature();

    // the feature isn't in an active entry
    if(!getEntryGroup().contains(event_feature)) 
      return;

    // if the feature is visible now or is in the list of visible features
    //(ie. it was visible previously) then redisplay.
    if(featureVisible(event_feature) ||
       getVisibleFeatures().contains(event_feature)) 
    {
      // update the visible_features vector
      if(getVisibleFeatures().contains(event_feature) &&
         !featureVisible(event_feature)) 
        getVisibleFeatures().remove(event_feature);
      else 
      {
        // the visibility of the feature has changed
        if(!getVisibleFeatures().contains(event_feature) &&
            featureVisible(event_feature)) 
          getVisibleFeatures().add(event_feature);
      }

      repaint();
    }
  }

  /**
   *  Implementation of the EntryGroupChangeListener interface.  We listen to
   *  EntryGroupChange events so that we can update the display if entries
   *  are added or deleted.
   **/
  public void entryGroupChanged(final EntryGroupChangeEvent event) 
  {
    switch(event.getType()) 
    {
      case EntryGroupChangeEvent.ENTRY_ADDED:
      case EntryGroupChangeEvent.ENTRY_ACTIVE:
      case EntryGroupChangeEvent.ENTRY_DELETED:
      case EntryGroupChangeEvent.ENTRY_INACTIVE:
      if(getOneLinePerEntryFlag()) 
        fixCanvasSize();
      
      needVisibleFeatureVectorUpdate();
      break;
    }

    repaint();
  }

  /**
   *  Implementation of the EntryChangeListener interface.  We listen to
   *  EntryChange events so that we can update the display if features are
   *  added or deleted.
   **/
  public void entryChanged(final EntryChangeEvent event) 
  {
    switch(event.getType()) 
    {
      case EntryChangeEvent.FEATURE_DELETED:
        remove(event.getFeature());
        break;
      case EntryChangeEvent.FEATURE_ADDED:
        add(event.getFeature());
        break;
    }
  }

  /**
   *  Implementation of the SelectionChangeListener interface.  We listen to
   *  SelectionChange events so that we can update the list to reflect the
   *  current selection.
   **/
  public void selectionChanged(final SelectionChangeEvent event) 
  {
    // don't bother with events we sent ourself
    if(event.getSource() == this) 
      return;

    needVisibleFeatureVectorUpdate();

    if(event.getType() == SelectionChangeEvent.SELECTION_CHANGED) 
      raise_selection_flag = true;

    repaint();
  }

  /**
   *  Implementation of the SequenceChangeListener interface.
   **/
  public void sequenceChanged(final SequenceChangeEvent event) 
  {
    visible_features = new FeatureVector();

    if(event.getType() == SequenceChangeEvent.REVERSE_COMPLEMENT) 
    {
      final int old_centre_position = getCentreForwardBase();

      final int new_centre_position =
                   getBases().getComplementPosition(old_centre_position);

      makeBaseVisibleInternal(new_centre_position, true, false);
    } 
    else 
      makeBaseVisibleInternal(event.getPosition(), true, false);

    fixScrollbar();
    needVisibleFeatureVectorUpdate();
    repaint();

    if(event.getType() == SequenceChangeEvent.REVERSE_COMPLEMENT) 
      fireAdjustmentEvent(DisplayAdjustmentEvent.REV_COMP_EVENT);
    else 
      fireAdjustmentEvent(DisplayAdjustmentEvent.SCROLL_ADJUST_EVENT);
  }

  /**
   *  Invoked when an Option is changed.
   **/
  public void optionChanged(OptionChangeEvent event) 
  {
    repaint();
  }

  /**
   *  Implementation of the GotoListener interface.  Invoked when the listener
   *  should goto to the given base.
   **/
  public void performGoto(final GotoEvent event) 
  {
    makeBaseVisible(event.getMarker());
  }

  /**
   *  Implementation of the DisplayAdjustmentListener interface.  Invoked when
   *  a component(FeatureDisplay) scrolls or changes the scale.  Set the
   *  position and scale of this FeatureDisplay to match the event
   **/
  public void displayAdjustmentValueChanged(DisplayAdjustmentEvent event) 
  {
    disable_display_events = true;
    try 
    {
      setScaleFactor(event.getScaleFactor());
      setFirstVisibleForwardBase(event.getStart());
    }
    finally 
    {
      disable_display_events = false;
    }
  }

  /**
   *  Return a MarkerRange that exactly covers the visible bases
   **/
  public MarkerRange getVisibleMarkerRange() 
  {
    final int first_base = getFirstVisibleForwardBase();
    final int last_base  = getLastVisibleForwardBase();

    final Strand strand = getBases().getForwardStrand();

    try 
    {
      return strand.makeMarkerRangeFromPositions(first_base, last_base);
    }
    catch(OutOfRangeException e) 
    {
      throw new Error("internal error - unexpected exception: " + e);
    }
  }

  /**
   *  Arrange for the given feature to be drawn last - ie. so that it will
   *  appear at the front of the other features.
   **/
  void raiseFeature(Feature feature) 
  {
    if(getVisibleFeatures().remove(feature)) 
    {
      getVisibleFeatures().addElementAtEnd(feature);
      repaint();
    }
  }

  /**
   *  Arrange for the given feature to be drawn first - ie. so that it will
   *  appear at the back of the other features.
   **/
  void lowerFeature(Feature feature) 
  {
    if(getVisibleFeatures().remove(feature))
    {
      getVisibleFeatures().insertElementAt(feature, 0);
      repaint();
    }
  }

  /**
   *  Redraw the display with the smallest features on top.  Selected features
   *  will always be on top.
   **/
  void smallestToFront() 
  {
    visible_features = new FeatureVector();
    needVisibleFeatureVectorUpdate();
    repaint();
  }

  /**
   *  Scroll the display so the given base on the Strand currently display on
   *  top will be in the centre of the display.  No events are sent.
   *  @param base_position The new centre base position.
   **/
  private void setCentreVisibleForwardBase(final int base_position) 
  {
    final int max_visible_bases = getMaxVisibleBases();
    final int possible_base_position = base_position - max_visible_bases / 2;
    int real_base_position;

    if(possible_base_position < 1 && hard_left_edge) 
      real_base_position = 1;
    else 
      real_base_position = possible_base_position;

    if(real_base_position > getSequenceLength()) 
      real_base_position = getSequenceLength();

    setFirstVisibleForwardBase(real_base_position);
  }

  /**
   *  Scroll the display so that the given base is in the middle of the
   *  screen.  This method obeys the rev_comp_display flag.
   *  @param base_position The base position to make visible.
   *  @param forward true means FORWARD - the base_position refers to a
   *    position on the forward strand, false means REVERSE - the
   *    base_position refers to the reverse strand.
   *  @param send_event Send a DisplayAdjustmentEvent if and only if this is
   *    true
   **/
  private void makeBaseVisibleInternal(final int base_position,
                                        final boolean forward,
                                        final boolean send_event) 
  {
    int forward_base_position = base_position;

    if(!forward ^ isRevCompDisplay()) 
      forward_base_position =
        getBases().getComplementPosition(forward_base_position);

    setCentreVisibleForwardBase(forward_base_position);

    if(send_event) 
      fireAdjustmentEvent(DisplayAdjustmentEvent.SCROLL_ADJUST_EVENT);
  }

  /**
   *  Scroll the display so that the given base is in the middle of the screen
   *  @param base_marker The Marker of base to make visible.
   **/
  private void makeBaseVisible(final Marker base_marker) 
  {
    makeBaseVisibleInternal(base_marker.getPosition(),
                             base_marker.getStrand().isForwardStrand(),
                             true);
  }

  /**
   *  Scroll the display so that the given base is in the middle of the screen.
   *  @param base The base to scroll to.
   **/
  public void makeBaseVisible(final int base) 
  {
    makeBaseVisibleInternal(base, true, true);
  }

  /**
   *  This method is called to add a single feature to the display by adding it
   *  to the vector of visible features(if it is visible that is).
   *  @param feature The object to add
   **/
  private void add(Feature feature) 
  {
    if(getEntryGroup().isActive(feature.getEntry()) &&
        featureVisible(feature)) 
    {
      if(visible_features.contains(feature)) 
      {
//        throw new Error("internal error - feature added a second time");
      } 
      else 
        getVisibleFeatures().addElementAtEnd(feature);
    }

    repaint();
  }


  /**
   *  The method is called when a Feature should be removed from the display.
   *  The feature is removed from the vector of visible features.
   *  @param feature The object to remove
   **/
  private void remove(Feature feature) 
  {
    if(visible_features != null) 
      visible_features.remove(feature);

    repaint();
  }

  /**
   *  Returns the vector containing those features that are currently visible.
   **/
  private FeatureVector getVisibleFeatures() 
  {
    return visible_features;
  }

  /**
   *  Returns a copy of the vector containing those features that are
   *  currently visible.
   **/
  public FeatureVector getCurrentVisibleFeatures() 
  {
    return(FeatureVector)visible_features.clone();
  }

  /**
   *  Returns a Range that starts at the first visible base and ends at the
   *  last visible base.
   **/
  private Range getVisibleRange()
  {
    final int first_visible_base = getFirstVisibleForwardBase();
    final int last_visible_base = getLastVisibleForwardBase();

    if(first_visible_base <= last_visible_base) 
      return newRange(first_visible_base, last_visible_base);
    else 
      return null;
  }

  /**
   *  This is called after scrolling to add all features that have become
   *  visible and remove those that have become invisible.  The method changes
   *  visible_features as little as possible so that features that were at the
   *  end(and hence drawn last/on top) stay there.  Selected features and
   *  segments will always be on top of unselected ones.
   **/
  private void updateVisibleFeatureVector() 
  {
    final Range visible_range;

    if(getSize().width == 0)
    {
      // don't bother doing any thinking
      visible_features = new FeatureVector();
      return;
    }

    if(raise_selection_flag)
    {
      final FeatureVector all_features = getSelection().getAllFeatures();

      for(int i = 0 ; i < all_features.size() ; ++i) 
        raiseFeature(all_features.elementAt(i));

      raise_selection_flag = false;
    }

    if(isRevCompDisplay())
    {
      final int first_visible_base = getFirstVisibleReverseBase();
      final int last_visible_base = getLastVisibleReverseBase();
      visible_range = newRange(first_visible_base, last_visible_base);
    } 
    else 
      visible_range = getVisibleRange();

    if(visible_range == null) 
    {
      visible_features = new FeatureVector();
      return;
    }

    final FeatureVector real_visible_features =
      getSortedFeaturesInRange(visible_range);

    final FeatureVector new_visible_features = new FeatureVector();

    // add features that are in visible_features and
    // real_visible_features - ie features that are still visible
    for(int i = 0 ; i < visible_features.size() ; ++i) 
    {
      final Feature new_feature = visible_features.elementAt(i);
      if(real_visible_features.contains(new_feature)) 
        new_visible_features.addElementAtEnd(new_feature);
    }

    // add features that are in real_visible_features and not currently
    // in visible_features and are not selected(selected features will be
    // added last so that they stay on top).
    for(int i = 0 ; i < real_visible_features.size() ; ++i) 
    {
      final Feature new_feature = real_visible_features.elementAt(i);

      if(!visible_features.contains(new_feature) &&
          !getSelection().contains(new_feature)) {
        new_visible_features.addElementAtEnd(new_feature);
      }
    }

    final FeatureVector selection_features = getSelection().getAllFeatures();

    // now add features that are in real_visible_features, are not in
    // visible_features and are selected(selected features are added last so
    // that they stay on top).
    for(int i = 0 ; i < real_visible_features.size() ; ++i) 
    {
      final Feature new_feature = real_visible_features.elementAt(i);
      if(!visible_features.contains(new_feature) &&
          selection_features.contains(new_feature)) 
        new_visible_features.addElementAtEnd(new_feature);
    }

    visible_features = new_visible_features;
    update_visible_features = false;
  }

  /**
   *  This is used by getSortedFeaturesInRange().
   **/
  final private static Comparator feature_comparator = new Comparator() 
  {
    /**
     *  Compare two Objects with respect to ordering.
     *  @return a negative number if feature1_object is less than
     *    feature2_object ; a positive number if feature1_object is greater
     *    than feature2_object; else 0
     **/
    public int compare(final Object feature1_object,
                       final Object feature2_object) 
    {
      final Feature feature1 =(Feature) feature1_object;
      final Feature feature2 =(Feature) feature2_object;

      final int feature1_size = feature1.getBaseCount();
      final int feature2_size = feature2.getBaseCount();

      if(feature1_size > feature2_size) 
        return -1;
      else 
      {
        if(feature1_size < feature2_size) 
          return 1;
        else
        {
          // use hash value as a last resort
          if(feature1.hashCode() < feature2.hashCode()) 
            return -1;
          else 
          {
            if(feature1.hashCode() == feature2.hashCode())
              return 0;
            else
              return 1;
          }
        }
      }
    }
  };

  /**
   *  Return a vector containing the references of the Feature objects in the
   *  EntryGroup that are within the given range are which should be
   *  displayed.  source features are not returned unless show_source_features
   *  is true.  features with a score less than getMinimumScore() or higher
   *  than getMaximumScore() aren't returned.
   *  @param range Return features that overlap this range - ie the start of
   *    the feature is less than or equal to the end of the range and the end
   *    of the feature is greater than or equal to the start of the range.
   *  @return The features the are within the given range.  The returned
   *    object is a copy - changes will not effect the EntryGroup object
   *    itself.  The Features are sorted so that the features with the least
   *    bases come last, which is the reverse of Entry.getFeaturesInRange().
   **/
  private FeatureVector getSortedFeaturesInRange(final Range range) 
  {
    try
    {
      FeatureVector features_from_entry =
        getEntryGroup().getFeaturesInRange(range);

      final int min_score = getMinimumScore();
      final int max_score = getMaximumScore();

      final FeatureVector filtered_features = new FeatureVector();

      // filter out low and high scoring features and(possibly) source
      // features
      for(int i = features_from_entry.size() - 1 ; i >= 0 ; --i)
      {
        final Feature this_feature = features_from_entry.elementAt(i);

        if(this_feature.getKey().equals("source") &&
            !getShowSourceFeatures() &&
            !getSelection().contains(this_feature)) 
          continue;

        if(min_score > 0 || max_score < 100) 
        {
          final int this_feature_score = this_feature.getScore();

          // features with no /score are always shown
          if(this_feature_score != -1 &&
             (this_feature_score < getMinimumScore() ||
               this_feature_score > getMaximumScore()) &&
              !getSelection().contains(this_feature)) 
            continue;
        }

        filtered_features.add(this_feature);
      }

      features_from_entry = filtered_features;

      final FeatureVector sorted_features =
        features_from_entry.sort(feature_comparator);

      return sorted_features;
    }
    catch(OutOfRangeException e) 
    {
      throw new Error("internal error - unexpected exception: " + e);
    }
  }

  /**
   *  The main paint function for the canvas.  An off screen image 
   *  used for double buffering when drawing the canvas.
   *  @param g The Graphics object of the canvas.
   **/
  protected void paintComponent(Graphics g) 
  {
    super.paintComponent(g);
    if(!isVisible()) 
      return;

    int scrollbar_hgt = 0;
    if(scrollbar_style == SCROLLBAR_AT_TOP)
    {
      scrollbar_hgt = scrollbar.getPreferredSize().height;
      ((Graphics2D)g).translate(0,scrollbar_hgt);
    }

    if(update_visible_features) 
      updateVisibleFeatureVector();

    fillBackground(g);

    final Selection selection = getSelection();
    final FeatureVector selected_features = selection.getAllFeatures();

    final FeatureSegmentVector selected_segments =
      selection.getSelectedSegments();

    // we draw the feature backgrounds first then the visible indication
    // that there is a MarkerRange selected, then the feature outlines.
    for(int i = 0 ; i < getVisibleFeatures().size() ; ++i) 
      drawFeature(g, getVisibleFeatures().elementAt(i),
                  true, selected_features, selected_segments);

    drawBaseSelection(g);
    drawScale(g);
    drawCodons(g);
    drawBases(g);

    // now draw the feature outlines
    for(int i = 0 ; i < getVisibleFeatures().size() ; ++i) 
      drawFeature(g, getVisibleFeatures().elementAt(i),
                  false, selected_features, selected_segments);

    if(scrollbar_style == SCROLLBAR_AT_TOP)
      ((Graphics2D)g).translate(0,-scrollbar_hgt);
    
    Thread.yield();
  }

  /**
   *  Draw the background colour of the frames.
   **/
  private void fillBackground(Graphics g) 
  {
    final int first_visible_base_coord =
      getLowXPositionOfBase(getFirstVisibleForwardBase());

    final int last_visible_base_coord =
      getHighXPositionOfBase(getLastVisibleForwardBase());

    if(getOneLinePerEntryFlag()) 
    {
      final int group_size = getEntryGroup().size();
      for(int i = 0 ; i < group_size; ++i) 
      {
        final int forward_entry_line = getDisplayLineOfEntryIndex(i, true);
        final int reverse_entry_line = getDisplayLineOfEntryIndex(i, false);

        final Entry current_entry = getEntryGroup().elementAt(i);

        if(getEntryGroup().getDefaultEntry() == current_entry &&
            Options.getOptions().highlightActiveEntryFlag())
        {
          g.setColor(active_entry_colour);
          fillLane(g, forward_entry_line,
                   first_visible_base_coord, last_visible_base_coord);
          fillLane(g, reverse_entry_line,
                   first_visible_base_coord, last_visible_base_coord);
        } 
        else
        {
          g.setColor(light_grey);
          fillLane(g, forward_entry_line,
                   first_visible_base_coord, last_visible_base_coord);
          fillLane(g, reverse_entry_line,
                   first_visible_base_coord, last_visible_base_coord);
        }
      }
    } 
    else
    {
      g.setColor(light_grey);
      if(show_forward_lines) 
      {
        fillLane(g, getFrameDisplayLine(FORWARD_FRAME_1),
                 first_visible_base_coord, last_visible_base_coord);
        fillLane(g, getFrameDisplayLine(FORWARD_FRAME_2),
                 first_visible_base_coord, last_visible_base_coord);
        fillLane(g, getFrameDisplayLine(FORWARD_FRAME_3),
                 first_visible_base_coord, last_visible_base_coord);
      }

      if(show_reverse_lines) 
      {
        fillLane(g, getFrameDisplayLine(REVERSE_FRAME_1),
                 first_visible_base_coord, last_visible_base_coord);
        fillLane(g, getFrameDisplayLine(REVERSE_FRAME_2),
                 first_visible_base_coord, last_visible_base_coord);
        fillLane(g, getFrameDisplayLine(REVERSE_FRAME_3), 
                 first_visible_base_coord, last_visible_base_coord);
      }
    }

    g.setColor(not_so_light_grey);
    fillLane(g, getFrameDisplayLine(FORWARD_STRAND), 
             first_visible_base_coord, last_visible_base_coord);
    fillLane(g, getFrameDisplayLine(REVERSE_STRAND), 
             first_visible_base_coord, last_visible_base_coord);
  }

  /**
   *  Fill one lane/frame line with the given colour.
   **/
  private void fillLane(Graphics g, int fill_line_number, 
                        int first_visible_base_coord, int last_visible_base_coord) 
  {
    int fill_line_top = fill_line_number*getLineHeight() + 1;

//  if(scrollbar_style == SCROLLBAR_AT_TOP)
//    fill_line_top += scrollbar.getPreferredSize().height;

//  g.setColor(colour);

//  final int first_visible_base_coord =
//    getLowXPositionOfBase(getFirstVisibleForwardBase());

//  final int last_visible_base_coord =
//    getHighXPositionOfBase(getLastVisibleForwardBase());

    g.fillRect(first_visible_base_coord, fill_line_top,
               last_visible_base_coord - first_visible_base_coord + 1,
               getFeatureHeight());
  }


  /**
   *  Draw the line showing the base numbers in the middle of the canvas. 
   *  The smallest numbers will be to the left.
   **/
  private void drawScale(Graphics g) 
  {
    g.setColor(Color.black);

    final int scale_line = getFrameDisplayLine(SCALE_LINE);
    final int scale_number_y_pos = scale_line * getLineHeight();

    final float bases_per_pixel =
     (float)getMaxVisibleBases() / getWidth();

    final int base_label_spacing;

    if(getScaleFactor() == 0) 
    {
      // set the spacing so that the labels are at multiples of ten at the
      // lowest scale factor
      base_label_spacing =
          (int) Math.ceil(MINIMUM_LABEL_SPACING * bases_per_pixel / 10) * 10;
    } 
    else 
    {
      // set the spacing so that the labels are at multiples of 100 otherwise
      base_label_spacing =
          (int) Math.ceil(MINIMUM_LABEL_SPACING * bases_per_pixel / 100) * 100;
    }

    final int label_spacing =(int)(base_label_spacing / bases_per_pixel);

    final int possible_index_of_first_label;
    final int seq_length = getSequenceLength();

    if(isRevCompDisplay()) 
      possible_index_of_first_label = (seq_length -
                         getLastVisibleForwardBase() + 1) / base_label_spacing;
    else 
      possible_index_of_first_label =
                         getFirstVisibleForwardBase() / base_label_spacing;

    final int index_of_first_label;

    if(possible_index_of_first_label <= 0) 
      index_of_first_label = 1;
    else 
      index_of_first_label = possible_index_of_first_label;

    final int index_of_last_label;

    if(isRevCompDisplay()) 
      index_of_last_label = (seq_length -
                   getFirstVisibleForwardBase() + 1) / base_label_spacing;
    else 
      index_of_last_label =
                   getLastVisibleForwardBase() / base_label_spacing;

    g.setFont(getFont());
    final int font_ascent   = getFontAscent();
    final int font_width    = getFontWidth();
    final int font_line_hgt = getLineHeight();

    for(int i = index_of_first_label; i <= index_of_last_label; ++i)
    {
      final String label_string=
          String.valueOf((int)(i * base_label_spacing));
      final int scale_number_x_pos;

      if(isRevCompDisplay()) 
        scale_number_x_pos =
          getLowXPositionOfBase(seq_length -
                                 i * base_label_spacing + 1);
      else 
        scale_number_x_pos = 
          getLowXPositionOfBase(i * base_label_spacing);

      g.drawString(label_string,
                    scale_number_x_pos + 2,
                    scale_number_y_pos + font_ascent + 1);

      g.drawLine(scale_number_x_pos, scale_number_y_pos,
                  scale_number_x_pos, scale_number_y_pos + font_line_hgt);

      if(isRevCompDisplay())
        g.drawLine(scale_number_x_pos, scale_number_y_pos,
                    scale_number_x_pos + font_width, scale_number_y_pos);
      else
        g.drawLine(scale_number_x_pos,
                    scale_number_y_pos + font_line_hgt,
                    scale_number_x_pos + font_width,
                    scale_number_y_pos + font_line_hgt);
    }
  }

  /**
   *  Draw the bases of the forward and reverse strands into a Graphics
   *  object.
   **/
  private void drawBases(Graphics g) 
  {
    if(getScaleFactor() > 1 ||
        getScaleFactor() == 1 && !show_base_colours) 
      return;

    final Strand forward_strand;
    final Strand reverse_strand;

    if(isRevCompDisplay()) 
    {
      reverse_strand = getBases().getForwardStrand();
      forward_strand = getBases().getReverseStrand();
    } 
    else 
    {
      forward_strand = getBases().getForwardStrand();
      reverse_strand = getBases().getReverseStrand();
    }

    final Range forward_range =
      newRange(getFirstVisibleForwardBase(),
                getLastVisibleForwardBase());

    String forward_visible_bases =
      forward_strand.getSubSequence(forward_range).toUpperCase();

    final int forward_frame_line = getFrameDisplayLine(FORWARD_STRAND);

    final int forward_sequence_length = forward_visible_bases.length();

    final int offset;

    if(getForwardBaseAtLeftEdge() < 1) 
      offset = 1 - getForwardBaseAtLeftEdge();
    else 
      offset = 0;
    

    for(int base_index = 0; base_index < forward_sequence_length;
        ++base_index)
    {
      final char char_to_draw =
        forward_visible_bases.charAt(base_index);

      if(getScaleFactor() == 0) 
        drawOneBaseLetter(g,
                           char_to_draw,
                          (offset + base_index) * getFontWidth(),
                           forward_frame_line * getLineHeight());
      else 
        drawOnePixelBase(g,
                          char_to_draw,
                          offset + base_index,
                          forward_frame_line * getLineHeight());
    }

    final Range reverse_range = newRange(getFirstVisibleReverseBase(),
                                          getLastVisibleReverseBase());

    String reverse_visible_bases =
      reverse_strand.getSubSequence(reverse_range).toUpperCase();

    final int reverse_frame_line = getFrameDisplayLine(REVERSE_STRAND);

    final int reverse_sequence_length = reverse_visible_bases.length();

    for(int base_index = 0; base_index < reverse_sequence_length;
        ++base_index) 
    {
      final char char_to_draw =
        reverse_visible_bases.charAt(reverse_sequence_length -
                                      base_index - 1);

      if(getScaleFactor() == 0) 
        drawOneBaseLetter(g,
                           char_to_draw,
                          (offset + base_index) * getFontWidth(),
                           reverse_frame_line * getLineHeight());
      else 
        drawOnePixelBase(g,
                          char_to_draw,
                          offset + base_index,
                          reverse_frame_line * getLineHeight());
    }
  }


  /**
   *  Draw the codons of the sequence into a graphics object.  If the scale is
   *  0 then the codon letters will be drawn, otherwise just the stop codons
   *  will marked.
   *  @param g The object to draw into.
   **/
  private void drawCodons(Graphics g) 
  {
    g.setColor(Color.black);

    if(getScaleFactor() == 0)  
    {
      if(!getOneLinePerEntryFlag()) 
      {
        if(show_forward_lines) 
          drawForwardCodonLetters(g);
        
        if(show_reverse_lines) 
          drawReverseCodonLetters(g);
      }
    }
    else 
    {
      final int MAX_STOP_CODON_SCALE_FACTOR = 7;
      if((show_stop_codons || show_start_codons) &&
          getScaleFactor() <= MAX_STOP_CODON_SCALE_FACTOR) 
      {
        if(!getOneLinePerEntryFlag()) 
        {
          if(show_forward_lines) 
            drawForwardCodons(g);

          if(show_reverse_lines) 
            drawReverseCodons(g);
        }
      }
    }
  }

  /**
   *  Mark the start and stop codons into the three forward frame lines.
   *  @param g The object to draw into.
   **/
  private void drawForwardCodons(Graphics g) 
  {
    final Strand strand;

    if(isRevCompDisplay()) 
      strand = getBases().getReverseStrand();
    else 
      strand = getBases().getForwardStrand();

    final int first_visible_base = getForwardBaseAtLeftEdge();

    final int frame_shift =(first_visible_base - 1) % 3;

    // base to start translation at - we start slightly off the left of the
    // screen
    final int start_base = first_visible_base - frame_shift;

    // base to end translation at
    // we + 3 to the upper bound because partial codons do not get translated
    // by getTranslation()
    final int end_base = getLastVisibleForwardBase() + 3;

    // not used if show_stop_codons is false
    int [][] forward_stop_codons = null;

    if(show_stop_codons) 
    {
      forward_stop_codons = new int [][] {
        strand.getStopCodons(newRange(start_base, end_base)),
        strand.getStopCodons(newRange(start_base + 1, end_base)),
        strand.getStopCodons(newRange(start_base + 2, end_base))
      };
    }

    // not used if show_start_codons is false
    int [][] forward_start_codons = null;

    if(show_start_codons) 
    {
      final StringVector start_codons;

      if(Options.getOptions().isEukaryoticMode()) 
        start_codons = Options.getOptions().getEukaryoticStartCodons();
      else 
        start_codons = Options.getOptions().getProkaryoticStartCodons();

      forward_start_codons = new int [][] {
        strand.getMatchingCodons(newRange(start_base, end_base),
                                  start_codons),
        strand.getMatchingCodons(newRange(start_base + 1, end_base),
                                  start_codons),
        strand.getMatchingCodons(newRange(start_base + 2, end_base),
                                  start_codons)
      };
    }

    for(int i = 0 ; i < 3 ; ++i) 
    {
      final int frame_line = getFrameDisplayLine(FORWARD_FRAME_1 + i);

      final int frame_x_start = i;

      if(show_start_codons) 
      {
        final int [] this_frame_start_codons = forward_start_codons [i];

        drawCodonMarkLine(g,
                           frame_line,
                           this_frame_start_codons,
                           FORWARD,
                           80);
      }

      if(show_stop_codons) 
      {
        final int [] this_frame_stop_codons = forward_stop_codons [i];

        drawCodonMarkLine(g,
                           frame_line,
                           this_frame_stop_codons,
                           FORWARD,
                           100);
      }
    }
  }

  /**
   *  Mark the stop codons into the three reverse frame lines.
   *  @param g The object to draw into.
   **/
  private void drawReverseCodons(Graphics g) 
  {
    final Strand strand;

    if(isRevCompDisplay()) 
      strand = getBases().getForwardStrand();
    else 
      strand = getBases().getReverseStrand();

    final int first_visible_base = getFirstVisibleReverseBase();

    final int frame_shift = (first_visible_base - 1) % 3;

    // base to start translation at
    final int start_base = first_visible_base - frame_shift;

    // base to end translation at
    // we + 3 to the upper bound because partial codons do not get translated
    // by getTranslation()
    final int end_base = getLastVisibleReverseBase() + 3;

    // not used if show_stop_codons is false
    int [][] reverse_stop_codons = null;

    if(show_stop_codons) 
    {
      reverse_stop_codons = new int [][] {
        strand.getStopCodons(newRange(start_base, end_base)),
        strand.getStopCodons(newRange(start_base + 1, end_base)),
        strand.getStopCodons(newRange(start_base + 2, end_base))
      };
    }

    // not used if show_start_codons is false
    int [][] reverse_start_codons = null;

    if(show_start_codons) 
    {
      final StringVector start_codons;

      if(Options.getOptions().isEukaryoticMode()) 
        start_codons = Options.getOptions().getEukaryoticStartCodons();
      else 
        start_codons = Options.getOptions().getProkaryoticStartCodons();

      reverse_start_codons = new int [][] {
        strand.getMatchingCodons(newRange(start_base, end_base),
                                  start_codons),
        strand.getMatchingCodons(newRange(start_base + 1, end_base),
                                  start_codons),
        strand.getMatchingCodons(newRange(start_base + 2, end_base),
                                  start_codons)
      };
    }

    for(int i = 0 ; i < 3 ; ++i) 
    {
      final int frame_line = getFrameDisplayLine(REVERSE_FRAME_1 - i);
      final int frame_x_start = i;

      if(show_start_codons) 
      {
        final int [] this_frame_start_codons = reverse_start_codons [i];
        drawCodonMarkLine(g,
                          frame_line,
                          this_frame_start_codons,
                          REVERSE,
                          80);
      }

      if(show_stop_codons) 
      {
        final int [] this_frame_stop_codons = reverse_stop_codons [i];
        drawCodonMarkLine(g,
                          frame_line,
                          this_frame_stop_codons,
                          REVERSE,
                          100);
      }
    }
  }

  /**
   *  Draw the codon letters into the three forward frame lines.
   *  @param g The object to draw into.
   **/
  private void drawForwardCodonLetters(Graphics g) 
  {
    final Strand strand;

    if(isRevCompDisplay()) 
      strand = getBases().getReverseStrand();
    else 
      strand = getBases().getForwardStrand();

    final int first_visible_base = getFirstVisibleForwardBase();

    for(int i = 0 ; i < 3 ; ++i) 
    {
      final int frame_shift = 1 -(first_visible_base + 3 - i) % 3;

      // base to start translation at - we start slightly off the left of the
      // screen so that the first partial codon is translated as '.'
      final int start_base = first_visible_base + frame_shift;

      // base to end translation at - we end slightly off the right of the
      // screen
      final int end_base = getLastVisibleForwardBase() + 1;

      final int frame_line = getFrameDisplayLine(FORWARD_FRAME_1 + i);

      int frame_x_start = frame_shift;

      final AminoAcidSequence this_frame_translation =
        strand.getTranslation(newRange(start_base, end_base),
                               false);

      final String this_frame_translation_string =
        this_frame_translation.toString();

      drawCodonLine(g,
                     frame_x_start,
                     frame_line,
                     this_frame_translation_string,
                     FORWARD);
    }
  }

  /**
   *  Draw the codon letters into the three reverse frame lines.
   *  @param g The object to draw into.
   **/
  private void drawReverseCodonLetters(Graphics g) 
  {
    final Strand strand;

    if(isRevCompDisplay()) 
      strand = getBases().getForwardStrand();
    else 
      strand = getBases().getReverseStrand();

    final int first_visible_base = getFirstVisibleReverseBase();

    for(int i = 0 ; i < 3 ; ++i) 
    {
      final int frame_shift = 1 -(first_visible_base + 3 - i) % 3;

      // base to start translation at - we start slightly off the right of the
      // screen so that the first partial codon is translated as '.'
      final int start_base = first_visible_base + frame_shift;

      // base to end translation at - we end slightly off the left of the
      // screen
      final int end_base = getLastVisibleReverseBase() + 1;

      final int frame_line = getFrameDisplayLine(REVERSE_FRAME_1 - i);

      int frame_x_start = frame_shift;

      final AminoAcidSequence this_frame_translation =
        strand.getTranslation(newRange(start_base, end_base),
                               false);

      final String this_frame_translation_string =
        this_frame_translation.toString();

      drawCodonLine(g,
                     frame_x_start,
                     frame_line,
                     this_frame_translation_string,
                     REVERSE);
    }
  }


  /**
   *  Draw one line of codons.
   *  @param g The object to draw into.
   *  @param frame_start When drawing in the FORWARD direction this is the
   *    offset from the left of the screen at which to start drawing the line.
   *    For REVERSE it is the offset from the last visible base.
   *  @param line_number The frame line to draw into. (see
   *    getFrameDisplayLine() for more details about these line numbers).
   *  @param codons A String containing the codon letters to draw.
   *  @param direction The direction to draw the letters in(see the
   *    frame_start parameter).
   **/
  private void drawCodonLine(Graphics g,
                              int frame_start,
                              int line_number,
                              String codons,
                              int direction) 
  {
    final int offset;

    if(getForwardBaseAtLeftEdge() < 1 && direction != REVERSE) 
      offset = 1 - getForwardBaseAtLeftEdge();
    else 
      offset = 0;

    final String upper_case_codons = codons.toUpperCase();
    final int draw_y_position = line_number * getLineHeight();

    // draw each AA letter
    for(int i = 0 ; i < upper_case_codons.length() ; ++i) 
    {
      final char aa_char = upper_case_codons.charAt(i);

      int draw_x_position =
       (int)((offset + frame_start + 3 * i + 1) * getScaleValue());

      if(direction == REVERSE) 
      {
        // draw in the opposite direction
        draw_x_position =
          getLowXPositionOfBase(getLastVisibleForwardBase()) -
          draw_x_position;
      }

      drawOneLetter(g, aa_char, draw_x_position, draw_y_position);
    }
  }

  /**
   *  Draw one line of codons marks(stop codons or start codons).
   *  @param g The object to draw into.
   *  @param line_number The frame line to draw into. (see
   *    getFrameDisplayLine() for more details about these line numbers).
   *  @param codon_positions A vector containing the base positions of the
   *    codons to mark.
   *  @param direction The direction to draw the letters in(see the
   *    frame_start parameter).
   *  @param height_percent The height of the mark as a percentage of the
   *    frame line height.
   **/
  private void drawCodonMarkLine(final Graphics g,
                                 final int line_number,
                                 final int [] codon_positions,
                                 final int direction,
                                 final int height_percent) 
  {

    final int first_visible_forward_base = getForwardBaseAtLeftEdge();
    final int first_visible_reverse_base = getFirstVisibleReverseBase();
    final int last_visible_forward_base = getLastVisibleForwardBase();
    final double scale_value = getScaleValue();

    final int line_height = getLineHeight();

    final Integer colour_number =
      Options.getOptions().getIntegerProperty("colour_of_start_codon");

    final Color start_codon_colour;

    if(colour_number != null) 
    {
      final int int_colour_number = colour_number.intValue();
      start_codon_colour =
        Options.getOptions().getColorFromColourNumber(int_colour_number);
    } 
    else 
      start_codon_colour = Color.blue;

    final int draw_y_position = line_number * line_height +
        (int)(1.0 * line_height *(100 - height_percent) / 100 / 2 + 0.5);

    final int mark_height = height_percent * line_height / 100;

    // used to remember the position of the previously drawn codon
    int last_x_position = -1;

    for(int i = 0 ; i < codon_positions.length ; ++i) 
    {
      final int codon_base_start = codon_positions[i];

      // zero is the end of data marker
      if(codon_base_start == 0) 
        break;

      int draw_x_position;

      if(direction == FORWARD) 
        draw_x_position = getLowXPositionOfBase(codon_base_start);
      else
      { 
        final int raw_base_position =
          getBases().getComplementPosition(codon_base_start);
        draw_x_position = getLowXPositionOfBase(raw_base_position);
      }

      // draw only if we haven't drawn on the position already
      if(last_x_position == -1 || draw_x_position != last_x_position) 
      {
        if(height_percent < 100) 
          drawOneCodonMark(g, draw_x_position, draw_y_position,
                            direction, mark_height, start_codon_colour);
        else 
          drawOneCodonMark(g, draw_x_position, draw_y_position,
                            direction, mark_height, Color.black);

        last_x_position = draw_x_position;
      }
    }
  }

  /**
   *  Return the colour in which the given base should be drawn: blue for C,
   *  red for T, green for A and black for G.
   **/
  private static Color getColourOfBase(final char base_char) 
  {
    switch(base_char) 
    {
      case 'c': case 'C':
        return Color.blue;
      case 't': case 'T':
        return Color.red;
      case 'a': case 'A':
        return dark_green;
      case 'g': case 'G':
        return Color.black;
      default:
        return Color.white;
    }
  }


  /**
   *  Draw a letter at a specified position in a Graphics object.
   *  @param g The object to draw into.
   *  @param letter The character to draw.
   *  @param x_pos The x position at which to draw the letter.
   *  @param y_pos The y position at which to draw the letter.
   **/
  private void drawOneLetter(final Graphics g, final char letter,
                              final int x_pos, final int y_pos) 
  {
    draw_one_char_temp_array[0] = letter;

    g.setFont(getFont());
    g.drawChars(draw_one_char_temp_array, 0, 1, x_pos,
                y_pos + getFontAscent() + 1);
  }

  /**
   *  Draw a letter at a specified position in a Graphics object.
   *  @param g The object to draw into.
   *  @param base_char The character of the base to draw(A, T, G or C).
   *  @param x_pos The x position at which to draw the letter.
   *  @param y_pos The y position at which to draw the letter.
   **/
  private void drawOneBaseLetter(final Graphics g, final char base_char,
                                  final int x_pos, final int y_pos) 
  {
    if(show_base_colours) 
      g.setColor(getColourOfBase(base_char));

    drawOneLetter(g, base_char, x_pos, y_pos);
  }


  /**
   *  Draw a codon mark at the given position.  The codon is represented
   *  by a one pixel wide line for scale factor greater than 1 and 3 pixels
   *  wide if the scale factor is one.
   *  @param g The object to draw into.
   *  @param x_pos The x position at which to draw the codon.
   *  @param y_pos The y position at which to draw the codon.  The line
   *    is drawn down from the given x,y position.
   *  @param direction If this is FORWARD the line is draw from the y position
   *    to the right, otherwise it is drawn to the left.
   *  @param height The height in pixels of the codon mark.
   *  @param colour The colour to use when drawing the mark.
   **/
  private void drawOneCodonMark(final Graphics g,
                                 final int x_pos, final int y_pos,
                                 final int direction,
                                 final int height,
                                 final Color colour) 
  {
    g.setColor(colour);

    if(getScaleFactor() == 1) 
    {
      // draw a line three pixels/bases wide
      if(direction == FORWARD) 
      {
        g.drawLine(x_pos + 1, y_pos,
                    x_pos + 1, y_pos + height - 1);
        g.drawLine(x_pos + 2, y_pos,
                    x_pos + 2, y_pos + height - 1);
      } 
      else
      {
        g.drawLine(x_pos - 1, y_pos,
                    x_pos - 1, y_pos + height - 1);
        g.drawLine(x_pos - 2, y_pos,
                    x_pos - 2, y_pos + height - 1);
      }
    }
    g.drawLine(x_pos, y_pos, x_pos, y_pos + height - 1);
  }

  /**
   *  Draw a one pixel wide line representing a single base on the forward or
   *  reverse strand.  The base is coloured according to which base is passed
   *  to the method: blue for C, red for T, green for A and black for G.
   *  @param g The object to draw into.
   *  @param aa_char The character to draw.
   *  @param x_pos The x position at which to draw the line representing the
   *    base.
   *  @param y_pos The y position at which to draw the top of the line
   *    representing the base.
   **/
  private void drawOnePixelBase(Graphics g, char base_char,
                                 int x_pos, int y_pos) 
  {
    g.setColor(getColourOfBase(base_char));
    g.drawLine(x_pos, y_pos, x_pos, y_pos + getLineHeight() - 1);
  }

  /**
   *  Return the line on the canvas where this feature segment should be drawn.
   *  @param segment The segment in question.
   *  @return The line to draw into.
   **/
  private int getSegmentDisplayLine(FeatureSegment segment)
  {
    if(getOneLinePerEntryFlag()) 
    {
      final Feature feature = segment.getFeature();
      if((feature.isProteinFeature() || frame_features_flag) &&
         (show_forward_lines &&(feature.isForwardFeature() ^
                                  isRevCompDisplay()) ||
          show_reverse_lines &&(!feature.isForwardFeature() ^
                                  isRevCompDisplay()))) 
      {
        return getFeatureDisplayLine(feature);
      } 
      else
      {
        if(feature.isForwardFeature() ^ isRevCompDisplay()) 
          return getFrameDisplayLine(FORWARD_STRAND);
        else 
          return getFrameDisplayLine(REVERSE_STRAND);
      }
    } 
    else
    {
      final int frame_id =
        maybeFlipFrameDirection(getSegmentFrameID(segment));

      return getFrameDisplayLine(frame_id);
    }
  }

  /**
   *  Return the frame ID of to use when drawing the given segment.
   **/
  private int getSegmentFrameID(final FeatureSegment segment) 
  {
    final int frame_id = segment.getFrameID();

    final Feature feature = segment.getFeature();

    if((feature.isProteinFeature() || frame_features_flag) &&
       (show_forward_lines &&(feature.isForwardFeature() ^
                                isRevCompDisplay())||
         show_reverse_lines &&(! feature.isForwardFeature() ^
                                isRevCompDisplay()))) 
    {
      return frame_id;
    }
    else 
    {
      switch(frame_id) 
      {
        case FORWARD_FRAME_1:
          // fall through
        case FORWARD_FRAME_2:
          // fall through
        case FORWARD_FRAME_3:
          // fall through
        case FORWARD_STRAND:
          return FORWARD_STRAND;
        case REVERSE_FRAME_1:
          // fall through
        case REVERSE_FRAME_2:
          // fall through
        case REVERSE_FRAME_3:
          // fall through
        case REVERSE_STRAND:
          return REVERSE_STRAND;
        default:
          return frame_id;
      }
    }
  }

  /**
   *  Return the line on the canvas where the given feature should be drawn
   *  when one_line_per_entry is true.  Call getLineCount() to find out the
   *  total number of lines.  Each line is getLineHeight() pixels high.
   *  @return The line to draw into.
   **/
  private int getFeatureDisplayLine(final Feature feature) 
  {
    final int entry_index = getEntryGroup().indexOf(feature.getEntry());
    return getDisplayLineOfEntryIndex(entry_index,
                                      feature.isForwardFeature() ^
                                      isRevCompDisplay());
  }

  /**
   *  Return the line in the display where features from the entry specified
   *  by the given entry_index should be drawn.
   *  @param is_forward_feature indicates whether the feature of interest is
   *    on the forward or reverse strand.
   **/
  private int getDisplayLineOfEntryIndex(final int entry_index,
                                         final boolean is_forward_feature) 
  {
    if(is_forward_feature) 
    {
      if(getShowForwardFrameLines()) 
      {
        if(getShowLabels()) 
          return entry_index * 2;
        else 
          return entry_index;
      }
      else
        return 0;
    } 
    else
    {
      int return_value = 0;
      if(getShowLabels()) 
      {
        return_value = getEntryGroup().size() * 2 + 3 - entry_index * 2;

        if(getShowForwardFrameLines()) 
          return_value += getEntryGroup().size() * 2;
      }
      else 
      {
        return_value = getEntryGroup().size() + 2 - entry_index;

        if(getShowForwardFrameLines()) 
          return_value += getEntryGroup().size();
      }
      return return_value;
    }
  }

  /**
   *  If rev_comp_display is true return the corresponding frame_id from the
   *  opposite strand(eg. FORWARD_FRAME_1 gives REVERSE_FRAME_3) otherwise
   *  return the frame_id unchanged.
   **/
  private int maybeFlipFrameDirection(final int frame_id) 
  {
    if(isRevCompDisplay()) 
    {
      // flip the frame so that forward becomes reverse and reverse becomes
      // forward
      switch(frame_id) 
      {
        case FORWARD_FRAME_1:
          return REVERSE_FRAME_1;
        case FORWARD_FRAME_2:
          return REVERSE_FRAME_2;
        case FORWARD_FRAME_3:
          return REVERSE_FRAME_3;
        case FORWARD_STRAND:
          return REVERSE_STRAND;
        case REVERSE_FRAME_1:
          return FORWARD_FRAME_1;
        case REVERSE_FRAME_2:
          return FORWARD_FRAME_2;
        case REVERSE_FRAME_3:
          return FORWARD_FRAME_3;
        case REVERSE_STRAND:
          return FORWARD_STRAND;
        default:
          return frame_id;
      }
    }
    return frame_id;
  }

  /**
   *  Return the line on the canvas where this frame ID should be drawn.  Call
   *  getLineCount() to find out the total number of lines.  Each line is
   *  getLineHeight() pixels high.
   *  @param frame_id The frame ID.
   *  @return The line to draw into.
   **/
  private int getFrameDisplayLine(int frame_id) 
  {
    if(getOneLinePerEntryFlag()) 
    {
      int return_value;
      switch(frame_id) 
      {
        case FORWARD_STRAND:
          return_value = 0;
          break;
        case REVERSE_STRAND:
          if(show_labels) 
            return_value = 3;
          else 
            return_value = 2;
          break;
        case SCALE_LINE:
          if(show_labels) 
            return_value = 2;
          else 
            return_value = 1;
          break;
        default:
          throw new Error("internal error - unexpected value: " + frame_id);
      }

      if(show_forward_lines) 
      {
        if(show_labels)
          return_value += getEntryGroup().size();

        return return_value + getEntryGroup().size();
      } 
      else 
        return return_value;
    }

    final int line_number;

    switch(frame_id) 
    {
      case FORWARD_FRAME_1:
        line_number = 0;
        break;
      case FORWARD_FRAME_2:
        if(show_labels) 
          line_number = 2;
        else 
          line_number = 1;
        break;
      case FORWARD_FRAME_3:
        if(show_labels) 
          line_number = 4;
        else 
          line_number = 2;
        break;
      case FORWARD_STRAND:
        if(show_forward_lines)
        {
          if(show_labels)
            line_number = 6;
          else
            line_number = 3;
        }
        else
          line_number = 0;
        break;
      case REVERSE_STRAND:
        if(show_forward_lines) 
        {
          if(show_labels) 
            line_number = 9;
          else
            line_number = 5;
        }
        else
        {
          if(show_labels)
            line_number = 3;
          else 
            line_number = 2;
        }
        break;
      case REVERSE_FRAME_3:
        if(show_forward_lines) 
        {
          if(show_labels) 
            line_number = 11;
          else 
            line_number = 6;
        } 
        else
        {
          if(show_labels) 
            line_number = 5;
          else 
            line_number = 3;
        }
        break;
      case REVERSE_FRAME_2:
        if(show_forward_lines) 
        {
          if(show_labels) 
            line_number = 13;
          else 
            line_number = 7;
        }
        else
        {
          if(show_labels) 
            line_number = 7;
          else 
            line_number = 4;
        }
        break;
      case REVERSE_FRAME_1:
        if(show_forward_lines) 
        {
          if(show_labels) 
            line_number = 15;
          else 
            line_number = 8;
        }
        else
        {
          if(show_labels) 
            line_number = 9;
          else 
            line_number = 5;
        }
        break;
      default:
        if(show_forward_lines) 
        {
          if(show_labels)
            line_number = 8;
          else
            line_number = 4;
        }
        else
        {
          if(show_labels) 
            line_number = 2;
          else 
            line_number = 1;
        }
        break;
    }

    return line_number;
  }

  /**
   *  Return the number of lines of text we need to fit on the canvas.  This
   *  is used to set the height of the canvas.
   **/
  private int getLineCount() 
  {
    int line_count;

    if(show_labels) 
      line_count = 5;
    else 
      line_count = 3;

    int extra_line_count;

    if(getOneLinePerEntryFlag())
    {
      // some number of entry lines
      extra_line_count = getEntryGroup().size();
    } 
    else
    {
      // three frame line
      extra_line_count = 3;
    }

    if(show_labels) 
      extra_line_count *= 2;

    if(show_forward_lines) 
      line_count += extra_line_count;

    if(show_reverse_lines) 
      line_count += extra_line_count;

    return line_count;
  }

  /**
   *  Return the vertical offset from the top of the canvas for this feature.
   *  The returned value will be the y-coordinate of the top of the lane that
   *  this feature should be draw into.
   **/
  private int getSegmentVerticalOffset(FeatureSegment segment) 
  {
    // one "line" is font_height pixels high,(the number of lines times the
    // font_height is the height of the height of canvas)
    final int line = getSegmentDisplayLine(segment);

    return getLineOffset(line);
  }

  /**
   *  Given a line return the vertical offset(in pixels) from the top of the
   *  canvas.
   **/
  private int getLineOffset(int line) 
  {
    return line * getLineHeight();
  }

  /**
   *  Return the lowest on screen(with respect to the canvas) x coordinate of
   *  a base.  If the scale_factor is zero then one base will be font_width
   *  wide and this method will return a different value than
   *  getHighXPositionOfBase().  If scale_factor is greater than one then
   *  the two methods will return the same thing.
   *  @param base_number The(forward) base to calculate the position of.
   **/
  private int getLowXPositionOfBase(int base_number) 
  {
    return(int)((base_number - getForwardBaseAtLeftEdge()) *
                  getScaleValue());
  }

  /**
   *  Return the highest on screen(ie with respect to the canvas) x
   *  coordinate of a base.  See comment on getLowXPositionOfBase().
   *  @param base_number The(forward) base to calculate the position of.
   **/
  private int getHighXPositionOfBase(int base_number) 
  {
    if(getScaleFactor() == 0) 
      return getLowXPositionOfBase(base_number) + getFontWidth() - 1;
    else
      return getLowXPositionOfBase(base_number);
  }

  /**
   *  Return the low on screen x coordinate of the base at the given Marker
   *  position.  If the scale_factor is zero then one base will be font_width
   *  wide so that this method will return a different position than
   *  getHighXPositionOfMarker().
   *  @param marker The Marker position to return the screen position of.
   **/
  private int getLowXPositionOfMarker(Marker marker) 
  {
    int position = marker.getRawPosition();

    if(isRevCompDisplay()) 
      position = getSequenceLength() - position + 1;

    if(marker.getStrand().isForwardStrand() ^ isRevCompDisplay()) 
      return getLowXPositionOfBase(position);
    else 
      return getHighXPositionOfBase(position);
  }

  /**
   *  Return the high on screen x coordinate of the base at the given Marker
   *  position.  If the scale_factor is zero then one base will be font_width
   *  wide so that this method will return a different position than
   *  getLowXPositionOfMarker().
   *  @param marker The Marker position to return the screen position of.
   **/
  private int getHighXPositionOfMarker(Marker marker) 
  {
    int position = marker.getRawPosition();

    if(isRevCompDisplay()) 
      position = getSequenceLength() - position + 1;

    if(marker.getStrand().isForwardStrand() ^ isRevCompDisplay()) 
      return getHighXPositionOfBase(position);
    else 
      return getLowXPositionOfBase(position);
  }


  /**
   *  This is the approximate opposite of getLowXPositionOfBase() - it returns
   *  MarkerRange corresponding to the given screen position.  It returns a
   *  MarkerRange rather than a Marker because if the user clicks on a FRAME
   *  line then we want to select a whole codon not just one base.
   *  @param position A position on the canvas.
   *  @return A MarkerRange covering the clicked on bases.
   **/
  private MarkerRange getMarkerRangeFromPosition(Point position)
  {
    return getMarkerRangeFromPosition(position, true);
  }

  /**
   *  Given a Point and a direction(either FORWARD_STRAND or REVERSE_STRAND),
   *  return a the corresponding base position on that Strand.
   **/
  private int getBasePositionOfPoint(final Point position,
                                      final int direction) 
  {
    if(direction == FORWARD_STRAND) 
    {
      return(int)(1.0 * position.x / getScaleValue() +
                    getForwardBaseAtLeftEdge());
    } 
    else
    {
      if(direction == REVERSE_STRAND) 
      {
        final int raw_base_position =
         (int)(1.0 * position.x / getScaleValue() +
                 getForwardBaseAtLeftEdge());
        return getBases().getComplementPosition(raw_base_position);
      } 
      else 
        throw new Error("internal error - unexpected value: " + direction);
    }
  }

  /**
   *  This is the approximate opposite of getLowXPositionOfBase() - it returns
   *  MarkerRange corresponding to the given screen position.  It returns a
   *  MarkerRange rather than a Marker because if the user clicks on a FRAME
   *  line then we want to select a whole codon not just one base.
   *  @param position A position on the canvas.
   *  @param whole_codon If true then return a range that covers a whole codon
   *    if the click was on a frame line.
   *  @return A MarkerRange covering the clicked on bases.
   **/
  private MarkerRange getMarkerRangeFromPosition(final Point position,
                                                  final boolean whole_codon) 
  {

    final int frame_id;

    int base_position;
    final Strand strand;

    if(scrollbar_style == SCROLLBAR_AT_TOP)
      position.x += scrollbar.getPreferredSize().height;

    if(getOneLinePerEntryFlag()) 
    {
      final int scale_line = getFrameDisplayLine(SCALE_LINE);
      final int position_line = getLineFromPoint(position);

      if(position_line < scale_line) 
      {
        if(isRevCompDisplay()) 
          strand = getBases().getReverseStrand();
        else 
          strand = getBases().getForwardStrand();
        
        frame_id = FORWARD_STRAND;
      } 
      else 
      {
        if(position_line > scale_line)
        {
          if(isRevCompDisplay()) 
            strand = getBases().getForwardStrand();
          else 
            strand = getBases().getReverseStrand();
          
          frame_id = REVERSE_STRAND;
        } 
        else 
          return null;
      }

      base_position = getBasePositionOfPoint(position, frame_id);
    } 
    else
    {
      frame_id = getFrameFromPoint(position);

      // calculate the frame/strand line and base position
      switch(frame_id) 
      {
        case FORWARD_FRAME_1:
        case FORWARD_FRAME_2:
        case FORWARD_FRAME_3:
        case FORWARD_STRAND:
          if(isRevCompDisplay()) 
            strand = getBases().getReverseStrand();
          else 
            strand = getBases().getForwardStrand();
          
          base_position = getBasePositionOfPoint(position, FORWARD_STRAND);
          break;
        case REVERSE_FRAME_1:
        case REVERSE_FRAME_2:
        case REVERSE_FRAME_3:
        case REVERSE_STRAND:
          if(isRevCompDisplay()) 
            strand = getBases().getForwardStrand();
          else 
            strand = getBases().getReverseStrand();
          
          base_position = getBasePositionOfPoint(position, REVERSE_STRAND);
          break;
        default:
          base_position = 0;
          strand = null;
      }
    }

    final int start_base_position;
    final int end_base_position;

    if(whole_codon) 
    {
      start_base_position =
        adjustBasePositionForFrame(frame_id, base_position, true);
      end_base_position =
        adjustBasePositionForFrame(frame_id, base_position, false);
    }
    else
    {
      start_base_position = base_position;
      end_base_position = base_position;
    }

    if(strand == null) 
      return null;
    else
    {
      try 
      {
        return strand.makeMarkerRangeFromPositions(start_base_position,
                                                    end_base_position);
      }
      catch(OutOfRangeException e) 
      {
        // XXX
        return null;
      }
    }
  }


  /**
   *  Helper method for getMarkerRangeFromPosition().  Returns a base
   *  position offset to the start/end of the codon if the given frame id on a
   *  FRAME.
   *  @param get_start If true then the base position returned is base
   *    position of the start of the codon, otherwise the end base is
   *    returned.
   **/
  private int adjustBasePositionForFrame(int frame_id,
                                         final int base_position,
                                         final boolean get_start) 
  {
    final int return_base_position;
    final int end_offset;

    if(get_start) 
      end_offset = 0;
    else 
    {
      // the end is two bases ahead of the start
      end_offset = 2;
    }

    switch(frame_id) 
    {
      case FORWARD_FRAME_1:
      case REVERSE_FRAME_1:
        {
          final int base_pos_mod3 =(base_position - 1) % 3;
          return_base_position = base_position + end_offset - base_pos_mod3;
        }
        break;
      case FORWARD_FRAME_2:
      case REVERSE_FRAME_2:
        {
          final int base_pos_mod3 =(base_position - 2) % 3;
          return_base_position = base_position + end_offset - base_pos_mod3;
        }
        break;
      case FORWARD_FRAME_3:
      case REVERSE_FRAME_3:
        {
          final int base_pos_mod3 =(base_position - 3) % 3;
          return_base_position = base_position + end_offset - base_pos_mod3;
        }
        break;
      default:
        return_base_position = base_position;
        // do nothing
    }

    return return_base_position;
  }

  /**
   *  Draw one feature at the correct place in a Graphics object.
   *  @param g The Graphics object on which to draw.
   *  @param draw_feature_fill If true then draw just the solid block of colour
   *    inside the segments.  If false then draw just the outline.
   *  @param feature The feature to draw.
   **/
  private void drawFeature(final Graphics g,
                           final Feature feature,
                           final boolean draw_feature_fill,
                           final FeatureVector selected_features,
                           final FeatureSegmentVector selected_segments) 
  {
    final FeatureSegmentVector this_feature_segments = feature.getSegments();

    // don't try to draw a feature with no segments
    if(this_feature_segments.size() == 0) 
      return;

    // set to true if and only if the whole of this feature should be
    // highlighted
    final boolean highlight_feature_flag;

    if(selected_features.contains(feature))
    {
      // ignore the possibility that a feature and a segment from the same
      // feature could be in the selection vector at the same time
      highlight_feature_flag = true;
    }
    else
    {
      // if the feature border flag is off and the feature is not selected
      // then don't draw the border
      if(!show_feature_borders && !draw_feature_fill) 
        return;

      highlight_feature_flag = false;
    }

    // draw each segment/exon
    for(int i = 0 ; i < this_feature_segments.size() ; ++i) 
    {
      final FeatureSegment current_segment =
                           this_feature_segments.elementAt(i);
      final boolean highlight_segment_flag;

      if(selected_segments.indexOf(current_segment) == -1) 
        highlight_segment_flag = false;
      else 
        highlight_segment_flag = true;

      final boolean draw_direction_arrow_flag;

      // draw an arrow only on the last segment
      if(i == this_feature_segments.size() - 1 && show_feature_arrows) 
        draw_direction_arrow_flag = true;
      else 
        draw_direction_arrow_flag = false;

      drawSegment(g, current_segment,
                   highlight_feature_flag, highlight_segment_flag,
                   draw_feature_fill,
                   draw_direction_arrow_flag);

      if(i + 1 < this_feature_segments.size()) 
      {
          // draw a line between the segments

        final FeatureSegment next_segment =
          this_feature_segments.elementAt(i + 1);

        drawSegmentConnection(g, current_segment, next_segment);
      }
    }

    // draw the label last if the is no label line because in this case the
    // label is draw on top of the feature segments
    if(draw_feature_fill) 
      drawFeatureLabel(g, feature);

//  Thread.yield();
  }

  /**
   *  Draw a bent line between two segments which represents the connection
   *  between exons in feature.
   *  @param g The Graphics object on which to draw.
   *  @param lower_segment The reference of the segment that is closest to the
   *    start of the Strand.  The connection line will start at the beginning
   *    of this segment.
   *  @param upper_segment The reference of the segment that is furthest from
   *    the start of the Strand.  The connection line will finish at the end
   *    of this segment.
   **/
  private void drawSegmentConnection(Graphics g,
                                     FeatureSegment lower_segment,
                                     FeatureSegment upper_segment) 
  {
    final Marker upper_segment_start_marker = upper_segment.getStart();
    final Marker lower_segment_end_marker = lower_segment.getEnd();

    int next_segment_start_coord =
      getLowXPositionOfMarker(upper_segment_start_marker);

    // make sure we don't wrap around when drawing
    if(next_segment_start_coord > 16000) 
      next_segment_start_coord = 16000;

    // make sure we don't wrap around when drawing
    if(next_segment_start_coord < -16000) 
      next_segment_start_coord = -16000;

    int this_segment_end_coord =
      getHighXPositionOfMarker(lower_segment_end_marker);

    // make sure we don't wrap around when drawing
    if(this_segment_end_coord > 16000) 
      this_segment_end_coord = 16000;

    // make sure we don't wrap around when drawing
    if(this_segment_end_coord < -16000) 
      this_segment_end_coord = -16000;

    final int this_segment_vertical_offset =
      getSegmentVerticalOffset(lower_segment);
    final int next_segment_vertical_offset =
      getSegmentVerticalOffset(upper_segment);

    // we draw the line with a bend in the middle - this is the vertical
    // position of the bend
    final int line_y_position_for_centre;

    if(this_segment_vertical_offset < next_segment_vertical_offset) 
      line_y_position_for_centre = this_segment_vertical_offset;
    else 
      line_y_position_for_centre = next_segment_vertical_offset;

    final int line_y_position_for_this =
      this_segment_vertical_offset + getFeatureHeight() / 2;
    final int line_y_position_for_next =
      next_segment_vertical_offset + getFeatureHeight() / 2;

    final int horizontal_position_of_centre =
     (this_segment_end_coord + next_segment_start_coord) / 2;

    final Feature segment_feature = lower_segment.getFeature();
    final Color feature_colour = segment_feature.getColour();

    // draw in black if no colour is specified or if the feature is selected
    if(feature_colour == null ||
       getSelection().contains(lower_segment.getFeature())) 
      g.setColor(Color.black);
    else 
      g.setColor(feature_colour);

    g.drawLine(this_segment_end_coord,
                line_y_position_for_this,
                horizontal_position_of_centre,
                line_y_position_for_centre);

    g.drawLine(horizontal_position_of_centre,
                line_y_position_for_centre,
                next_segment_start_coord,
                line_y_position_for_next);
  }

  /**
   *  Draw a label a feature.  This is a helper function for drawFeature().
   *  If show_labels is true the labels will be drawn below the features,
   *  otherwise they will be drawn within the features.
   *  @param g The Graphics object on which to draw.
   *  @param feature The feature to draw the label for.
   **/
  private void drawFeatureLabel(Graphics g, Feature feature) 
  {

    // the three frame translation is visible when the scale factor is 0,
    // don't draw labels over it
    if(!show_labels && getScaleFactor() == 0) 
      return;

    String label_or_gene = feature.getIDString();
    final String label = feature.getLabel();

    // special case - don't display a label if the label qualifier is "*"
    if(label != null && label.equals("*")) 
      label_or_gene = "";

    // don't waste time drawing nothing
    if(label_or_gene.length() == 0) 
      return;

    final FontMetrics fm = g.getFontMetrics();
    final int string_width = fm.stringWidth(label_or_gene);

    final FeatureSegment first_segment = feature.getSegments().elementAt(0);
    final int label_x_coord;

    if(feature.isForwardFeature() ^ isRevCompDisplay())
    {
      int segment_start_pos =
        first_segment.getStart().getRawPosition();

      if(isRevCompDisplay()) 
        segment_start_pos = getSequenceLength() - segment_start_pos + 1;

      label_x_coord = getLowXPositionOfBase(segment_start_pos);
    } 
    else
    {
      int segment_end_pos =
        first_segment.getEnd().getRawPosition();

      if(isRevCompDisplay()) 
        segment_end_pos = getSequenceLength() - segment_end_pos + 1;

      // reverse the start and end on the complementary strand
      label_x_coord = getLowXPositionOfBase(segment_end_pos);
    }

    if(label_x_coord >= getSize().width) 
      return;

//  if(label_x_coord + string_width <= 0) {
      // don't draw the label if it is not visible on screen
//  }


    int vertical_offset =
      getSegmentVerticalOffset(first_segment);

    if(show_labels) 
      vertical_offset += getLineHeight(); // move to the label line

    // save this so we can restore it later
    final Shape saved_clip = g.getClip();

    // if we have a labels line draw a white background behind the
    // label
    if(show_labels) 
    {
      g.setColor(Color.white);

      g.fillRect(label_x_coord - getFontWidth(),
                 vertical_offset,
                 string_width + getFontWidth() * 2,
                 getLineHeight());
    } 
    else
    {
      // if there is no label line clip to the size of the first segment
      // and draw in there
      final int segment_start_coord =
        getSegmentStartCoord(first_segment);
      final int segment_end_coord =
        getSegmentEndCoord(first_segment);

      if(Math.abs(segment_end_coord - segment_start_coord) > 5) 
      {
        if(segment_end_coord > segment_start_coord)
          g.setClip(segment_start_coord, vertical_offset,
                     segment_end_coord - segment_start_coord,
                     getFeatureHeight());
        else 
          g.setClip(segment_end_coord, vertical_offset,
                     segment_start_coord - segment_end_coord,
                     getFeatureHeight());
      }
      else 
        return;  // don't draw small labels if there is no room
    }

    g.setColor(Color.black);
    g.setFont(getFont());
    g.drawString(label_or_gene, label_x_coord + 1,
                  vertical_offset + getFontAscent() + 1);

    if(!show_labels)
      g.setClip(saved_clip);
  }

  /**
   *  Return the position on the canvas where this segment starts.  The
   *  returned value will be -1 if the position is off the left of the screen
   *  and will be(width of the canvas) if the position is off the right of
   *  the screen.
   **/
  private int getSegmentStartCoord(FeatureSegment segment) 
  {
    final Marker segment_start_marker = segment.getStart();
    int segment_start_coord =
      getLowXPositionOfMarker(segment_start_marker);

    // make sure we don't wrap around when drawing
    if(segment_start_coord > getSize().width) 
      segment_start_coord = getSize().width;

    // make sure we don't wrap around when drawing
    if(segment_start_coord < 0) 
      segment_start_coord = -1;

    return segment_start_coord;
  }

  /**
   *  Return the position on the canvas where this segment ends.  The returned
   *  value will be -1 if the position is off the left of the screen and will
   *  be(width of the canvas) if the position is off the right of the screen.
   **/
  private int getSegmentEndCoord(FeatureSegment segment) 
  {
    final Marker segment_end_marker = segment.getEnd();
    int segment_end_coord = getHighXPositionOfMarker(segment_end_marker);

    // make sure we don't wrap around when drawing
    if(segment_end_coord > getSize().width) 
      segment_end_coord = getSize().width;

    // make sure we don't wrap around when drawing
    if(segment_end_coord < 0) 
      segment_end_coord = -1;

    return segment_end_coord;
  }

  /**
   *  Return true if the base given by the Marker is visible.
   **/
  private boolean baseVisible(Marker marker) 
  {
    final int first_visible_base = getForwardBaseAtLeftEdge();
    final int last_visible_base = getLastVisibleForwardBase();

    int marker_position = marker.getRawPosition();

    if(isRevCompDisplay()) 
      marker_position = getBases().getComplementPosition(marker_position);

    if(marker_position < first_visible_base ||
        marker_position > last_visible_base) 
      return false;
    else 
      return true;
  }

  /**
   *  Return if and only if the segment is(partly) within the range of bases
   *  we are currently displaying.
   *  @param segment The FeatureSegment to check.
   **/
  private boolean segmentVisible(FeatureSegment segment) 
  {
    int segment_start_coord = getSegmentStartCoord(segment);
    int segment_end_coord = getSegmentEndCoord(segment);

    if(segment_end_coord < 0 && segment_start_coord < 0 ||
        segment_start_coord >= getSize().width &&
        segment_end_coord >= getSize().width) 
      return false;
    else 
      return true;
  }

  /**
   *  Return if and only if some part of a feature is(partly) within the
   *  range of bases we are currently displaying.
   *  @param feature The Feature to check.
   **/
  private boolean featureVisible(Feature feature) 
  {
    if(getMinimumScore() > 0) 
    {
      final int feature_score = feature.getScore();

      if(feature_score != -1 && feature_score < getMinimumScore()) 
        return false;
    }

    final FeatureSegmentVector segments = feature.getSegments();

    for(int i = 0 ; i < segments.size() ; ++i) 
    {
      if(segmentVisible(segments.elementAt(i))) 
        return true;
    }

    return false;
  }

  /**
   *  Draw one FeatureSegment into a Graphics object.
   *  @param g The Graphics object on which to draw.
   *  @param segment The feature segment to draw
   *  @param highlight_feature If true draw an extra thick line
   *  @param highlight_segment If true draw the segment with a doubly thick
   *    line.
   *  @param draw_feature_fill If true then just the solid block of colour
   *    inside the segments will be drawn.  If false then only the feature
   *    outline is drawn.
   *  @param draw_arrow If true draw a direction arrow at the end of the
   *    segment.
   **/
  private void drawSegment(Graphics g, FeatureSegment segment,
                            boolean highlight_feature,
                            boolean highlight_segment,
                            boolean draw_feature_fill,
                            boolean draw_arrow) 
  {
    // not on screen
    if(!segmentVisible(segment)) 
      return;

    final Feature segment_feature = segment.getFeature();

    final int vertical_offset = getSegmentVerticalOffset(segment) + 1;

    int segment_start_coord = getSegmentStartCoord(segment);
    int segment_end_coord = getSegmentEndCoord(segment);

    final int segment_height = getFeatureHeight() - 1;

    // this is 1 if the feature is on the forward strand or on a forward frame
    // and -1 otherwise.  this used to draw the feature arrow in the right
    // direction.
    final int feature_direction;
    if(segment.getFeature().isForwardFeature() ^ isRevCompDisplay()) 
      feature_direction = 1;
    else 
      feature_direction = -1;

    final int [] x_points = {
      segment_start_coord, segment_end_coord,
      segment_end_coord, segment_start_coord
    };

    final int [] y_points = {
      vertical_offset, vertical_offset,
      vertical_offset + segment_height, vertical_offset + segment_height
    };

    if(draw_feature_fill) 
    {
      final Color feature_colour = segment_feature.getColour();

      // no colour means draw in white
      if(feature_colour == null) 
        g.setColor(Color.white);
      else 
        g.setColor(feature_colour);

      if(segment_feature.isForwardFeature() ^ isRevCompDisplay()) 
      {
        final int segment_width = segment_end_coord - segment_start_coord + 1;

        g.fillRect(segment_start_coord, vertical_offset,
                    segment_width, segment_height + 1);
      }
      else
      {
        final int segment_width = segment_start_coord - segment_end_coord + 1;

        g.fillRect(segment_end_coord, vertical_offset,
                    segment_width, segment_height + 1);
      }
    }
    else 
    {
      g.setColor(Color.black);
      g.drawPolygon(x_points, y_points, 4);

      if(highlight_feature)  // highlight selected features
      {
        // selected - highlight by drawing a thicker line

        x_points[0] -= feature_direction; x_points[1] += feature_direction;
        x_points[2] += feature_direction; x_points[3] -= feature_direction;
        y_points[0] -= 1; y_points[1] -= 1;
        y_points[2] += 1; y_points[3] += 1;
        g.drawPolygon(x_points, y_points, 4);

        x_points[0] -= feature_direction; x_points[1] += feature_direction;
        x_points[2] += feature_direction; x_points[3] -= feature_direction;
        y_points[0] -= 1; y_points[1] -= 1;
        y_points[2] += 1; y_points[3] += 1;
        g.drawPolygon(x_points, y_points, 4);

        if(highlight_segment) {
          x_points[0] -= feature_direction; x_points[1] += feature_direction;
          x_points[2] += feature_direction; x_points[3] -= feature_direction;
          y_points[0] -= 1; y_points[1] -= 1;
          y_points[2] += 1; y_points[3] += 1;
          g.drawPolygon(x_points, y_points, 4);

          x_points[0] -= feature_direction; x_points[1] += feature_direction;
          x_points[2] += feature_direction; x_points[3] -= feature_direction;
          y_points[0] -= 1; y_points[1] -= 1;
          y_points[2] += 1; y_points[3] += 1;
          g.drawPolygon(x_points, y_points, 4);
        }
      }

      if(draw_arrow)
      {
        // now draw the arrow point
        final int arrow_tip_x =
          x_points[1] + feature_direction * getFontWidth() * 8 / 10;
        final int arrow_tip_y =(y_points[1] + y_points[2]) / 2;

        g.drawLine(x_points[1], y_points[1], arrow_tip_x, arrow_tip_y);
        g.drawLine(arrow_tip_x, arrow_tip_y, x_points[2], y_points[2]);
      }
    }
  }

  /**
   *  Draw/highlight the selected range of bases on the FORWARD_STRAND or
   *  REVERSE_STRAND lines.
   **/
  private void drawBaseSelection(Graphics g) 
  {
    final MarkerRange selection_range = getSelection().getMarkerRange();

    if(selection_range == null) 
      return;

    final Marker raw_start_base;
    final Marker raw_end_base;

    if(selection_range.getStrand().isForwardStrand()) 
    {
      raw_start_base = selection_range.getRawStart();
      raw_end_base = selection_range.getRawEnd();
    }
    else
    {
      raw_start_base = selection_range.getRawEnd();
      raw_end_base = selection_range.getRawStart();
    }

    final int start_coord = getLowXPositionOfMarker(raw_start_base);
    final int end_coord = getHighXPositionOfMarker(raw_end_base);

    // not on screen
    if(start_coord > getSize().width &&
        end_coord > getSize().width) 
      return;

    // not on screen
    if(start_coord < 0 && end_coord < 0) 
      return;

    // the ID of the strand on which to draw the highlighting - this will be
    // FORWARD_STRAND or REVERSE_STRAND.
    final int strand_id;

    if(selection_range.getStrand().isForwardStrand() ^
        isRevCompDisplay()) 
      strand_id = FORWARD_STRAND;
    else
      strand_id = REVERSE_STRAND;

    if(!getOneLinePerEntryFlag()) 
    {
      // the ID of the frame on which to draw the highlighting - this will be
      // FORWARD_FRAME_1, 2 or 3 or REVERSE_FRAME_1, 2 or 3.
      final int frame_id;

      if(selection_range.getStrand().isForwardStrand() ^
          isRevCompDisplay()) 
      {
        switch((selection_range.getStart().getPosition() - 1) % 3) 
        {
          case 0:
            frame_id = FORWARD_FRAME_1;
            break;
          case 1:
            frame_id = FORWARD_FRAME_2;
            break;
          case 2:
            frame_id = FORWARD_FRAME_3;
            break;
          default:
            frame_id = NO_FRAME;
        }
      } 
      else 
      {
        switch((selection_range.getStart().getPosition() - 1) % 3)
        {
          case 0:
            frame_id = REVERSE_FRAME_1;
            break;
          case 1:
            frame_id = REVERSE_FRAME_2;
            break;
          case 2:
            frame_id = REVERSE_FRAME_3;
            break;
          default:
            frame_id = NO_FRAME;
        }
      }

      if(show_forward_lines && strand_id == FORWARD_STRAND ||
          show_reverse_lines && strand_id == REVERSE_STRAND) 
        drawBaseRange(g, start_coord, end_coord, frame_id, Color.pink);
    }

    drawBaseRange(g, start_coord, end_coord, strand_id, Color.yellow);
  }

  /**
   *  Draw a rectangle representing the currently selected bases.
   *  @param g The Graphics object on which to draw.
   *  @param start_coord The start x coordinate.
   *  @param start_coord The end x coordinate.
   *  @param frame_id The ID of the frame line where the box should be drawn.
   *  @param fill_colour The colour to use to draw the box(the border will be
   *    black).
   **/
  private void drawBaseRange(Graphics g,
                             int start_coord, int end_coord,
                             int frame_id,
                             Color fill_colour) 
  {
    if(start_coord < -1) 
      start_coord = -1;

    if(end_coord < -1) 
      end_coord = -1;

    if(end_coord > getWidth()) 
      end_coord = getWidth();

    if(start_coord > getWidth()) 
      start_coord = getWidth();

    final int frame_line = getFrameDisplayLine(frame_id);

    final int frame_line_y_coord = getLineOffset(frame_line);

    final int [] x_points = {
      start_coord, end_coord, end_coord, start_coord
    };

    final int [] y_points = {
      frame_line_y_coord, frame_line_y_coord,
      frame_line_y_coord + getFeatureHeight() + 1,
      frame_line_y_coord + getFeatureHeight() + 1
    };

    g.setColor(fill_colour);

    if(start_coord > end_coord) 
      g.fillRect(end_coord, frame_line_y_coord + 1,
                 start_coord - end_coord + 1, getFeatureHeight());
    else 
      g.fillRect(start_coord, frame_line_y_coord + 1,
                 end_coord - start_coord + 1, getFeatureHeight());

    g.setColor(Color.black);
    g.drawPolygon(x_points, y_points, 4);
  }

  /**
   *  Return the number of button the user has down.
   **/
  private int buttonDownCount(final MouseEvent event) 
  {
    int count = 0;
    if((event.getModifiers() & InputEvent.BUTTON1_MASK) > 0) 
      ++count;
    
    if((event.getModifiers() & InputEvent.BUTTON2_MASK) > 0) 
      ++count;
    
    if((event.getModifiers() & InputEvent.BUTTON3_MASK) > 0) 
      ++count;
    
    return count;
  }

  /**
   *  Add mouse and key listeners to the canvas.
   **/
  private void addListeners() 
  {

    addMouseListener(new MouseAdapter() 
    {
      /**
       *  Listen for mouse press events so that we can do popup menus and
       *  selection.
       **/
      public void mousePressed(MouseEvent event)
      {
        // adjust for scrollbar 
        if(scrollbar_style == SCROLLBAR_AT_TOP)
          event.translatePoint(0,-scrollbar.getPreferredSize().height);

        // finish dragging if the user presses any other button because
        // we may not get a mouseReleased() call on some platforms
        if(click_segment_marker != null) 
        {
          getEntryGroup().getActionController().endAction();
          click_segment_marker = null;
        }

        //ignore
        if(buttonDownCount(event) > 1) 
          return;

        last_mouse_press_event = event;
        handleCanvasMousePress(event);

        if(isMenuTrigger(event)) 
        {
          final FeaturePopup popup =
            new FeaturePopup(FeatureDisplay.this,
                             getEntryGroup(),
                             getSelection(),
                             getGotoEventSource(),
                             getBasePlotGroup());
          final JComponent parent = (JComponent)event.getSource();

          parent.add(popup);
          popup.show(parent, event.getX(), event.getY());
        }
      }

      /**
       *  Listen for mouse release events so that we can call endAction().
       **/
      public void mouseReleased(MouseEvent event) 
      {
        // adjust for scrollbar
        if(scrollbar_style == SCROLLBAR_AT_TOP)
          event.translatePoint(0,-scrollbar.getPreferredSize().height);

        last_mouse_press_event = null;

        if(click_segment_marker != null) 
        {
          getEntryGroup().getActionController().endAction();
          click_segment_marker = null;
        }
      }
    });

    // Listen for mouse motion events so that we can select ranges of bases.
    addMouseMotionListener(new MouseMotionAdapter() 
    {
      public void mouseDragged(MouseEvent event) 
      {
        // adjust for scrollbar
        if(scrollbar_style == SCROLLBAR_AT_TOP)
          event.translatePoint(0,-scrollbar.getPreferredSize().height);

        if(last_mouse_press_event != null &&
           !isMenuTrigger(last_mouse_press_event)) 
          handleCanvasMouseDrag(event);
      }
    });

    addKeyListener(new KeyAdapter() 
    {
      public void keyPressed(final KeyEvent event) 
      {
        handleKeyPress(FeatureDisplay.this, event);
      }
    });
  }

  /**
   *  Handle a mouse press event on the drawing canvas - select on click,
   *  select and broadcast it on double click.
   **/
  private void handleCanvasMousePress(MouseEvent event) 
  {
    if(event.getID() != MouseEvent.MOUSE_PRESSED) 
      return;

    requestFocus();

    if(event.getClickCount() == 2) 
      handleCanvasDoubleClick(event);
    else
      handleCanvasSingleClick(event);

    repaint();
  }

  /**
   *  Handle a double click on the canvas.
   **/
  private void handleCanvasDoubleClick(MouseEvent event) 
  {
    if(isMenuTrigger(event)) 
      return;

    if((event.getModifiers() & InputEvent.BUTTON2_MASK) != 0 ||
        event.isAltDown()) 
    {

      final Selectable clicked_thing = getThingAtPoint(event.getPoint());

      if(clicked_thing == null) 
      {
        // select the ORF around the click position
        final MarkerRange new_click_range =
          getMarkerRangeFromPosition(event.getPoint());

        if(new_click_range == null) 
          return;

        final MarkerRange orf_range =
          Strand.getORFAroundMarker(new_click_range.getStart(), true);

        // user clicked on a stop codon
        if(orf_range == null) 
          return;
        else 
          getSelection().setMarkerRange(orf_range);
      }
      else 
      {
        // edit the selected Feature

        final Feature clicked_feature = getFeatureOf(clicked_thing);

        getSelection().set(clicked_feature);

        if(Options.readWritePossible()) 
          new FeatureEdit(clicked_feature, getEntryGroup(),
                           getSelection(),
                           getGotoEventSource()).show();
      }
    }

    makeSelectionVisible();
  }

  /**
   *  Handle a single click on the canvas.
   **/
  private void handleCanvasSingleClick(MouseEvent event) 
  {
    click_segment_marker = null;

    final Selectable clicked_thing = getThingAtPoint(event.getPoint());

    // treate alt modifier like pressing button 2
    if(clicked_thing == null ||
       (event.getModifiers() & InputEvent.BUTTON2_MASK) != 0) 
    {

      // if the user presses the mouse button 2 on feature or segment we treat
      // it like the feature/segment isn't there

      // if the user hasn't clicked on anything then don't change the
      // selection
      if(isMenuTrigger(event)) 
        return;

      // something is selected but the user clicked on nothing with the
      // shift key down - do nothing
      if(event.isShiftDown() &&
          getSelection().getAllFeatures().size() > 0) 
        return;

      final MarkerRange old_selected_range = getSelection().getMarkerRange();

      final MarkerRange new_click_range =
        getMarkerRangeFromPosition(event.getPoint());

      if(new_click_range == null) 
      {
        click_range = null;
        getSelection().clear();
        return;
      }

      final MarkerRange new_selection_range;

      if(!event.isShiftDown() || old_selected_range == null) 
        new_selection_range = new_click_range;
      else 
      {
        if(old_selected_range.getStrand() ==
            new_click_range.getStrand()) 
        {
          // extend the selection
          new_selection_range =
            old_selected_range.combineRanges(new_click_range, true);
        }
        else 
          new_selection_range = old_selected_range;
      }

      getSelection().setMarkerRange(new_selection_range);
      click_range = new_selection_range;

    }  
    else
    {
      // can't select a feature and a MarkerRange
      if(getSelection().getMarkerRange() != null &&
         event.isShiftDown()) 
        return;

      // clear the MarkerRange because this isn't the start of a range drag
      click_range = null;

      getSelection().setMarkerRange(null);
      final Feature clicked_feature = getFeatureOf(clicked_thing);
      raiseFeature(clicked_feature);

      if(clicked_thing instanceof Feature) 
      {
        // toggle the feature in or out of the selection unless this event
        // is a popup trigger in which case we should just make sure the
        // feature is in the selection
        if(getSelection().contains(clicked_feature)) 
        {
          if(! isMenuTrigger(event)) 
          {
            // change the selection only if the user isn't popping up a menu
            if(event.isShiftDown()) 
              getSelection().remove(clicked_feature);
            else 
              getSelection().set(clicked_feature);
          }
        } 
        else
        {
          if(event.isShiftDown()) 
            getSelection().add(clicked_feature);
          else 
            getSelection().set(clicked_feature);
        }
      } 
      else
      {
        // must be a FeatureSegment

        // toggle the feature segment in or out of the selection unless
        // this event is a popup trigger in which case we should just make
        // sure the segment is in the selection
        final FeatureSegment clicked_segment =
         (FeatureSegment) clicked_thing;

        final FeatureSegmentVector selected_feature_segments =
          getSelection().getSelectedSegments();

        if(selected_feature_segments.contains(clicked_segment)) 
        {
          if(! isMenuTrigger(event)) 
          {
            if(event.isShiftDown()) 
              getSelection().remove(clicked_segment);
            else 
              getSelection().set(clicked_segment);
          }
        } 
        else
        {
          if(event.isShiftDown()) 
          {
            final FeatureVector selected_features =
              getSelection().getSelectedFeatures();

            if(selected_features.contains(clicked_feature)) 
              getSelection().remove(clicked_feature);
            else 
              getSelection().add(clicked_segment);
          }
          else 
            getSelection().set(clicked_segment);
        }
      }

      if(Options.getOptions().canDirectEdit() &&
          !clicked_feature.isReadOnly()) 
      {
        // check to see if the click was on a marker, the second argument is
        // false because we want the range to cover just one base

        final MarkerRange new_click_range =
          getMarkerRangeFromPosition(event.getPoint(), false);

        // search the feature to find if the start Marker of any of the
        // segments matches the start Marker of new_click_range or if an end
        // Marker matches the end Marker of new_click_range

        final FeatureSegmentVector segments = clicked_feature.getSegments();

        for(int i = 0 ; i < segments.size() ; ++i) 
        {
          final FeatureSegment this_segment = segments.elementAt(i);

          if(new_click_range.getStart().equals(this_segment.getStart()) &&
              this_segment.canDirectEdit()) 
          {
            click_segment_marker = this_segment.getStart();
            click_segment_marker_is_start_marker = true;
            other_end_of_segment_marker = this_segment.getEnd();
            getEntryGroup().getActionController().startAction();
            break;
          }

          if(new_click_range.getEnd().equals(this_segment.getEnd()) &&
              this_segment.canDirectEdit()) 
          {
            click_segment_marker = this_segment.getEnd();
            click_segment_marker_is_start_marker = false;
            other_end_of_segment_marker = this_segment.getStart();
            getEntryGroup().getActionController().startAction();
            break;
          }
        }
      }
    }
  }

  /**
   *  This is a helper method for handleCanvasSingleClick() and
   *  handleCanvasDoubleClick().
   *  @param object The should be an instance of Feature or FeatureSegment.
   *  @return If the the argument is a Feature it is returned, if it is a
   *    FeatureSegment the owning Feature is returned.
   **/
  private Feature getFeatureOf(final Object object) 
  {
    if(object instanceof Feature) 
      return(Feature) object;
    else 
      // object must be a FeatureSegment
      return((FeatureSegment) object).getFeature();
  }

  /**
   *  Handle a mouse drag that started on a Marker from a FeatureSegment.
   **/
  private void handleMarkerDrag(final MouseEvent event) 
  {
    MarkerRange drag_range =
      getMarkerRangeFromPosition(event.getPoint(), false);

    if(drag_range == null) 
      return;

    final int new_position;

    if(click_segment_marker_is_start_marker) 
    {
      // the Marker is at the start of the segment
      new_position = drag_range.getStart().getPosition();

      // don't go past the other end of the segment
      if(new_position > other_end_of_segment_marker.getPosition()) 
        return;
    }
    else
    {
      new_position = drag_range.getEnd().getPosition();

      // don't go past the other end of the segment
      if(new_position < other_end_of_segment_marker.getPosition()) 
        return;
    }

    try
    {
      click_segment_marker.setPosition(new_position);

      if(!baseVisible(click_segment_marker)) 
        makeBaseVisible(click_segment_marker);
    }
    catch(OutOfRangeException e) 
    {
      throw new Error("internal error - unexpected OutOfRangeException");
    }
  }

  /**
   *  Handle a mouse drag event on the drawing canvas.
   **/
  private void handleCanvasMouseDrag(final MouseEvent event)
  {
    if(event.isShiftDown() &&
       getSelection().getAllFeatures().size() > 0) 
      return;

    if(click_segment_marker != null) 
    {
      handleMarkerDrag(event);
      return;
    }

    MarkerRange drag_range = getMarkerRangeFromPosition(event.getPoint());

    if(drag_range == null) 
      return;

    final MarkerRange selected_range = getSelection().getMarkerRange();

    // if the start and current positions of the drag are not on the
    // same Strand then ignore this event
    if(selected_range != null &&
       drag_range.getStrand() != selected_range.getStrand()) 
      return;

    try 
    {
      // user has dragged off the screen so use a marker at position 1
      if(drag_range.getStart().getPosition() < 1) 
        drag_range = new MarkerRange(drag_range.getStrand(), 1, 1);

      // user has dragged off the screen so use a marker at position the
      // end of the sequence
      if(drag_range.getEnd().getPosition() > getSequenceLength()) 
        drag_range = new MarkerRange(drag_range.getStrand(),
                                      getSequenceLength(),
                                      getSequenceLength());
    }
    catch(OutOfRangeException e) 
    {
      throw new Error("internal error - unexpected OutOfRangeException");
    }

    final boolean start_base_is_visible = baseVisible(drag_range.getStart());
    final boolean end_base_is_visible   = baseVisible(drag_range.getEnd());

    if(!start_base_is_visible) 
      makeBaseVisible(drag_range.getStart());

    if(!end_base_is_visible) 
      makeBaseVisible(drag_range.getEnd());

    // scrolling was necessary so update the visible features vector
    if(!start_base_is_visible || !end_base_is_visible) 
      needVisibleFeatureVectorUpdate();

    final MarkerRange new_marker_range;

    if(click_range == null) 
      new_marker_range = drag_range;
    else 
      new_marker_range = selected_range.combineRanges(drag_range, true);

    getSelection().setMarkerRange(new_marker_range);

    repaint();
  }

  /**
   *  Return the reference of the Object at p on the canvas or null if there
   *  is none.  The returned Object will be a Feature or a FeatureSegment.
   **/
  private Selectable getThingAtPoint(Point p) 
  {
    final int line_of_point = getLineFromPoint(p);

    // point is not on a forward or reverse frame or strand
    if(line_of_point == -1) 
      return null;

    // go through the features in reverse order because the feature that is
    // drawn last will be on top
    for(int i = getVisibleFeatures().size() - 1 ; i >= 0 ; --i) 
    {
      final Feature current_feature = getVisibleFeatures().elementAt(i);
      final FeatureSegmentVector segments = current_feature.getSegments();

      for(int segment_index = 0; segment_index < segments.size();
          ++segment_index) 
      {
        final FeatureSegment current_segment =
                                         segments.elementAt(segment_index);

        final int line_of_segment = getSegmentDisplayLine(current_segment);

        if(line_of_segment == line_of_point ||
            show_labels && line_of_segment + 1 == line_of_point) 
        {
          // this segment is in the right frame so check that it is between
          // the start and end positions of the segment

          if(p.x >= getSegmentStartCoord(current_segment) &&
              p.x <= getSegmentEndCoord(current_segment) ||
              p.x <= getSegmentStartCoord(current_segment) &&
              p.x >= getSegmentEndCoord(current_segment)) 
          {
            final Feature segment_feature = current_segment.getFeature();

            // special case - if there is only one segment then select the
            // whole feature
            if(segment_feature.getSegments().size() == 1) 
              return segment_feature;
            else
              return current_segment;
          }
        }
      }
    }

    return null;
  }

  /**
   *  Create and add a new scroll bar to this FeatureDisplay.
   **/
  private void createScrollbar(final boolean scrollbar_at_top) 
  {
    scrollbar = new JScrollBar(Scrollbar.HORIZONTAL);
     
    scrollbar.addAdjustmentListener(new AdjustmentListener() 
    {
      public void adjustmentValueChanged(AdjustmentEvent e) 
      {
        if(left_edge_base != e.getValue()) 
        {
          if(e.getValue() < getSequenceLength()) 
            left_edge_base = e.getValue();
          else
          {
            if(left_edge_base == getSequenceLength()) 
              return;
            else 
              left_edge_base = getSequenceLength();
          }

          fireAdjustmentEvent(DisplayAdjustmentEvent.SCROLL_ADJUST_EVENT);
          needVisibleFeatureVectorUpdate();
          repaint();
        }
      }
    });


    if(scrollbar_at_top) 
      add(scrollbar, "North");
    else 
      add(scrollbar, "South");
    
//  fixScrollbar();
//  fireAdjustmentEvent(DisplayAdjustmentEvent.ALL_CHANGE_ADJUST_EVENT);
  }

  /**
   *  Send a DisplayAdjustmentEvent with the current start and end base to
   *  each of the listeners.  This has package scope because EntryEdit
   *  components need to send this event when a new BasePlot component is
   *  added.  The BasePlot needs to know the first visible base, last visible
   *  base and the maximum number of visible bases before the plot can be
   *  drawn.
   *  @param type the type of event: DisplayAdjustmentEvent.SCALE_ADJUST_EVENT,
   *    DisplayAdjustmentEvent.SCROLL_ADJUST_EVENT or
   *    DisplayAdjustmentEvent.ALL_CHANGE_ADJUST_EVENT.
   **/
  void fireAdjustmentEvent(final int type) 
  {
    if(disable_display_events) 
      return;

    final DisplayAdjustmentEvent event =
      new DisplayAdjustmentEvent(this,
                                 getForwardBaseAtLeftEdge(),
                                 getLastVisibleForwardBase(),
                                 getMaxVisibleBases(),
                                 getScaleValue(), getScaleFactor(),
                                 isRevCompDisplay(), type);

    fireAction(adjustment_listener_list, event);
  }

  /**
   *  Send an event to those object listening for it.
   *  @param listeners A Vector of the objects that the event should be sent
   *    to.
   *  @param event The event to send
   **/
  private void fireAction(Vector listeners, ChangeEvent event) 
  {
    final Vector targets;
    // copied from a book - synchronizing the whole method might cause a
    // deadlock
    synchronized(this) 
    {
      targets = (Vector)listeners.clone();
    }

    for(int i = 0; i < targets.size(); ++i) 
    {
      ChangeListener change_listener = (ChangeListener)targets.elementAt(i);

      final DisplayAdjustmentListener target =
                               (DisplayAdjustmentListener)change_listener;
      target.displayAdjustmentValueChanged((DisplayAdjustmentEvent)event);
    }
  }

  /**
   *  Create the scroll bar used for changing the scale and add it to the
   *  FeatureDisplay.
   **/
  private void createScaleScrollbar() 
  {
    scale_changer = new JScrollBar(Scrollbar.VERTICAL);
    // try to arrange for the scrollbar to have a maximum value big enough
    // that the whole sequence can be visible at once
    final int MAX_FACTOR =
      (int)Math.round(Math.log(getSequenceLength()/20) /  Math.log(3));
    scale_changer.setValues(getScaleFactor(), 1, 0, MAX_FACTOR);
    scale_changer.setBlockIncrement(1);
    scale_changer.setUnitIncrement(1);
    scale_changer.addAdjustmentListener(new AdjustmentListener() 
    {
      public void adjustmentValueChanged(AdjustmentEvent e) 
      {
        setScaleFactor(e.getValue());
      }
    });

    add(scale_changer, "East");

    if(scale_factor >= MAX_FACTOR) 
      setScaleFactor(MAX_FACTOR - 1);
  }

  /**
   *  Update the parameters of the scrollbar taking changes to the entry_group
   *  into account.
   **/
  private void fixScrollbar() 
  {
    if(scrollbar == null)
      return;

    final int sequence_length = getSequenceLength();

    final int max_visible_bases;

    // make sure max_visible_bases is at least 1 to stop setBlockIncrement()
    // from complaining
    if(getMaxVisibleBases() > 0) 
      max_visible_bases = getMaxVisibleBases();
    else 
      max_visible_bases = 1;

    // initial_value, visible, minimum, maximum
    scrollbar.setValues(getForwardBaseAtLeftEdge(),
                         max_visible_bases,
                         hard_left_edge ? 1 : -max_visible_bases / 2,
                         sequence_length + max_visible_bases / 2);

    scrollbar.setBlockIncrement(max_visible_bases);

    if(max_visible_bases >= 10 && getScaleFactor() > 0)
      scrollbar.setUnitIncrement(max_visible_bases / 10);
    else
      scrollbar.setUnitIncrement(1);
  }

  /**
   *  Set the size of the canvas, taking the value of font_height and
   *  show_label into account.
   **/
  private void fixCanvasSize() 
  {
    int new_width  = getSize().width;
    int new_height = getLineHeight() * getLineCount();

    if(scale_changer != null)
      new_width  += scale_changer.getPreferredSize().height;
    if(scrollbar != null)
      new_height += scrollbar.getPreferredSize().height;

    if(new_height != getSize().width ||
        new_width != getSize().height) 
    {
      final Dimension preferred_size =
        new Dimension(getSize().width, new_height);
      setPreferredSize(preferred_size);

      final Dimension minimum_size = new Dimension(1, new_height);
      setMinimumSize(minimum_size);

      revalidate();
      repaint();
    }
  }

  /**
   *  Returns the base length of the sequence we are displaying in this
   *  component.
   **/
  public int getSequenceLength() 
  {
    return getEntryGroup().getSequenceLength();
  }

  /**
   *  Set the first visible base in the forward direction (set the base on the
   *  forward strand that is on the left of the canvas).  This method will
   *  scroll the display so that this base is is at the very left hand side of
   *  the canvas.
   **/
  private void setFirstVisibleForwardBase(int new_position) 
  {
    if(left_edge_base != new_position) 
    {
      left_edge_base = new_position;
      if(scrollbar != null) 
        scrollbar.setValue(new_position);
      
      needVisibleFeatureVectorUpdate();
      repaint();
    }
  }

  /**
   *  Return the height each line of the display should be.  Each frame will
   *  be drawn into one line.
   **/
  private int getLineHeight() 
  {
    return getFontHeight();
  }

  /**
   *  Return the current scale factor.  The scale factor is a number greater
   *  than or equal to zero than controls the number of bases that can appear
   *  on screen.  See getScaleValue().
   **/
  public int getScaleFactor() 
  {
    return scale_factor;
  }

  /**
   *  Return the ID of the first(ie. top) line of the display.
   **/
  private int getFirstLineID() 
  {
    if(show_forward_lines) 
      return FORWARD_FRAME_1;
    else 
      return FORWARD_STRAND;
  }

  /**
   *  Return the ID of the last(ie. bottom) line of the display.
   **/
  private int getLastLineID() 
  {  
    if(show_reverse_lines) 
      return REVERSE_FRAME_1;
    else 
      return REVERSE_STRAND;
  }


  private int getDisplayHeight()
  {
    return getHeight() - scrollbar.getHeight();
  }


  protected int getDisplayWidth()
  {
    return getWidth() - scale_changer.getWidth();
  }


  /**
   *  Return the display line that contains the given point, or -1 if
   *  the point is off screen.
   **/
  private int getLineFromPoint(final Point point) 
  {
    if(point.y >= getDisplayHeight()) 
      return -1;

    final int return_value = point.y / getLineHeight();

    if(return_value < 0)
      return -1;
    else 
      return return_value;
  }

  /**
   *  Return the ID of the frame at the given Point in the canvas or NO_FRAME
   *  if the Point is not within a frame. If labels are on then points on the
   *  label lines get returned as the frame id of the STRAND or FRAME line
   *  above.
   **/
  private int getFrameFromPoint(final Point point) 
  {
    if(getOneLinePerEntryFlag()) 
    {
      // there are no frame lines so just look at the strand lines
      final int line_point = getLineFromPoint(point);
      if(line_point == getFrameDisplayLine(FORWARD_STRAND) ^
         isRevCompDisplay()) 
        return FORWARD_STRAND;
      
      if(line_point == getFrameDisplayLine(REVERSE_STRAND) ^
          isRevCompDisplay()) 
        return REVERSE_STRAND;
      
      // fall through
    } 
    else 
    {
      final int start_frame = getFirstLineID();
      final int end_frame = getLastLineID();

      for(int i = start_frame ; i <= end_frame ; ++i) {
        final int top_of_line = getLineOffset(getFrameDisplayLine(i));
        final int line_height;

        if(show_labels) 
          line_height = getLineHeight() * 2;
        else 
          line_height = getLineHeight();

        if(point.y >= top_of_line &&
           point.y < top_of_line + line_height) 
          return i;
      }
    }
    // the point isn't in any or the frames
    return NO_FRAME;
  }

  /**
   *  Return the base number of the first visible base displayed in the
   *  forward direction, ie the base on the forward(or top strand) that is on
   *  the left each of the canvas.  Note that the strand that is currently
   *  being display in the forward direction will be the reverse strand if
   *  rev_comp_display is true.  The number returned will always be > 1 and <
   *  the sequence length.  If hard_left_edge is true this method will always
   *  return the same as getForwardBaseAtLeftEdge*(.
   **/
  public int getFirstVisibleForwardBase() 
  {
    if(left_edge_base < 1) 
      return 1;
    else 
      return left_edge_base;
  }

  /**
   *  Return the base number of the first visible base displayed in the
   *  forward direction, ie the base on the forward(or top strand) that is on
   *  the left each of the canvas.  Note that the strand that is currently
   *  being display in the forward direction will be the reverse strand if
   *  rev_comp_display is true.
   **/
  public int getForwardBaseAtLeftEdge() 
  {
    return left_edge_base;
  }

  /**
   *  Return a base marker of the first visible base in the forward
   *  direction, ie the base on the forward strand that is on the left
   *  each of the canvas.  See comments on getFirstVisibleForwardBase().
   **/
  private Marker getFirstVisibleForwardBaseMarker() 
  {
    try 
    {
      final Strand strand;
      if(isRevCompDisplay()) 
        strand = getBases().getReverseStrand();
      else 
        strand = getBases().getForwardStrand();
      
      return strand.makeMarker(getFirstVisibleForwardBase());
    } 
    catch(OutOfRangeException e) 
    {
      throw new Error("internal error - unexpected OutOfRangeException");
    }
  }

  /**
   *  Return the base number of the last visible base in the forward
   *  direction, ie the base on the forward strand that is on the right
   *  edge of the canvas.  The number returned will always be > 1 and < the
   *  sequence length.
   **/
  public int getLastVisibleForwardBase() 
  {
    final int possible_last_base =
     (int)(getForwardBaseAtLeftEdge() + getMaxVisibleBases());

    if(possible_last_base > getSequenceLength()) 
      return getSequenceLength();
    else 
    {
      if(possible_last_base > 0) 
        return possible_last_base;
      else 
        return 1;
    }
  }

  /**
   *  Return the base number of the first visible base in the reverse
   *  direction , ie the base on the reverse strand that is on the left
   *  each of the canvas.
   **/
  private int getFirstVisibleReverseBase() 
  {
    final int last_forward_base = getLastVisibleForwardBase();
    return getBases().getComplementPosition(last_forward_base);
  }

  /**
   *  Return the base number of the last visible base in the reverse
   *  direction, ie the base on the reverse strand that is on the right
   *  edge of the canvas.
   **/
  private int getLastVisibleReverseBase() 
  {
    // XXX
    final int first_forward_base = getFirstVisibleForwardBase();
    return getBases().getComplementPosition(first_forward_base);
  }

  /**
   *  Return the base of the forward Strand of the sequence that is closest to
   *  the centre of the view.
   **/
  private int getCentreForwardBase() 
  {
    int possible_return_position =
      getForwardBaseAtLeftEdge() + getMaxVisibleBases() / 2;

    int return_position;

    if(possible_return_position < getSequenceLength()) 
      return_position = possible_return_position;
    else 
      return_position = getSequenceLength();

    return return_position;
  }

  /**
   *  Return the amount display is scaled.  The value is 3 to the power of
   *  -(scale_factor - 1), except for a scale_factor of zero which gives a
   *  scale value of font_width.
   **/
  float getScaleValue() 
  {
    return scale_value;
  }

  /**
   *  Return the height in pixels we should use for drawing features.
   **/
  private int getFeatureHeight() 
  {
    // don't use getLineHeight() because we want a nice space between
    // frames/lines
    return getFontAscent() + 2;
  }

  /**
   *  Return the number of bases we can fit on screen at once, ie the number
   *  that will fit side by side on the canvas.
   **/
  public int getMaxVisibleBases() 
  {
    return(int)(getWidth() / getScaleValue());
  }

  /**
   *  Return the Bases object of the EntryGroup that this FeatureDisplay is
   *  displaying.
   **/
  public Bases getBases() 
  {
    return getEntryGroup().getBases();
  }

  /**
   *  Update the value of scale_value to reflect the current value of
   *  scale_factor.
   **/
  private void setScaleValue() 
  {
    final int scale_factor = getScaleFactor();

    if(scale_factor > 0) 
      scale_value =(float)(1 / Math.pow(3, scale_factor - 1));
    else 
      scale_value = getFontWidth();
  }

  /**
   *  Scroll and scale the display so that the given first base is at the left
   *  edge of the screen and the given last base is at the right edge.
   **/
  public void setFirstAndLastBase(final int first, final int last) 
  {
    left_edge_base = first;
    setScaleValue(1.0F * getWidth() /(last - first + 1));
  }

  /**
   *  Scroll the display so that the given first base is at the left edge of
   *  the screen.
   **/
  public void setFirstBase(int base_position) 
  {
    if(base_position > getSequenceLength()) 
      base_position = getSequenceLength();

    if(base_position < 1 && hard_left_edge) 
      base_position = 1;

    setFirstVisibleForwardBase(base_position);
    fireAdjustmentEvent(DisplayAdjustmentEvent.SCROLL_ADJUST_EVENT);
  }

  /**
   *  Set the scale value to use for this FeatureDisplay.  Note that it is
   *  better to call setScaleFactor() which will also set the scale value.
   **/
  private void setScaleValue(final float scale_value) 
  {
    if(scale_value <= 0) 
      throw new Error("internal error in FeatureDisplay.setScaleValue() - " +
                       "scale value must be positive");

    this.scale_value = scale_value;

    if(scale_value > 1) 
      scale_factor = 1;
    else 
      scale_factor =
       (int) Math.round(Math.log(1/scale_value) /  Math.log(3)) + 1;

    scale_changer.setValue(scale_factor);

    fixScrollbar();
    needVisibleFeatureVectorUpdate();
    repaint();
    fireAdjustmentEvent(DisplayAdjustmentEvent.SCALE_ADJUST_EVENT);
  }

  /**
   *  Make a new Range and throw a Error is an OutOfRangeException occurs.
   **/
  private Range newRange(final int start, final int end) 
  {
    try 
    {
      return new Range(start, end);
    } 
    catch(OutOfRangeException e) 
    {
      throw new Error("internal error - unexpected exception: " + e);
    }
  }

  /**
   *  Arrange for updateVisibleFeatureVector() to be called in the next call
   *  to paint()
   **/
  private void needVisibleFeatureVectorUpdate() 
  {
    update_visible_features = true;
  }

}
