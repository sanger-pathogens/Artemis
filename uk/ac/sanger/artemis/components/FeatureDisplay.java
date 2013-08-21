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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/FeatureDisplay.java,v 1.66 2009-06-12 13:50:35 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.components.filetree.*;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.RemoteFileDocument;
import uk.ac.sanger.artemis.util.ReadOnlyException;
import uk.ac.sanger.artemis.util.FileDocument;
import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.util.StringVector;
import uk.ac.sanger.artemis.io.ChadoCanonicalGene;
import uk.ac.sanger.artemis.io.InvalidRelationException;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.SimpleEntryInformation;
import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.RawStreamSequence;
import uk.ac.sanger.artemis.io.FastaStreamSequence;
import uk.ac.sanger.artemis.io.Sequence;
import uk.ac.sanger.artemis.io.Location;
import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.sequence.*;
import uk.ac.sanger.artemis.components.genebuilder.*;

import java.io.IOException;
import java.io.File;
import java.awt.event.*;
import java.awt.*;
import java.lang.Math;
import java.util.Vector;
import java.util.Comparator;
import javax.swing.border.Border;
import javax.swing.border.BevelBorder;

import java.awt.datatransfer.*;
import java.awt.dnd.*;
import javax.swing.Box;
import javax.swing.JOptionPane;
import javax.swing.JScrollBar;
import javax.swing.JComponent;
import javax.swing.UIManager;
import javax.swing.ImageIcon;
import javax.swing.JFrame;

import org.apache.batik.svggen.SVGGraphics2D;

/**
 *  This component is used for displaying an Entry.
 *
 *  @author Kim Rutherford
 *  @version $Id: FeatureDisplay.java,v 1.66 2009-06-12 13:50:35 tjc Exp $
 **/

public class FeatureDisplay extends EntryGroupPanel
  implements EntryGroupChangeListener,
             EntryChangeListener, FeatureChangeListener,
             SelectionChangeListener, GotoListener, SequenceChangeListener,
             DisplayComponent, OptionChangeListener, DisplayAdjustmentListener,
             DragGestureListener, DropTargetListener,
             DragSourceListener
{

  /** */
  private static final long serialVersionUID = 1L;

  private int highlight_drop_base = -1;

  /** Key code for calling zoomToSelection(). */
  private final static int ZOOM_TO_SELECTION_KEY = KeyEvent.VK_Z;
  private final static int ARROW_LEFT = KeyEvent.VK_LEFT;
  private final static int ARROW_RIGHT = KeyEvent.VK_RIGHT;

  protected final static int SCROLLBAR_AT_TOP = 1;
  protected final static int SCROLLBAR_AT_BOTTOM = 2;

  private final static int FORWARD = Bases.FORWARD;
  private final static int REVERSE = Bases.REVERSE;

  private final static int NO_FRAME        = FeatureSegment.NO_FRAME;
  private final static int FORWARD_STRAND  = FeatureSegment.FORWARD_STRAND;
  private final static int REVERSE_STRAND  = FeatureSegment.REVERSE_STRAND;
  private final static int FORWARD_FRAME_1 = FeatureSegment.FORWARD_FRAME_1;
  private final static int FORWARD_FRAME_2 = FeatureSegment.FORWARD_FRAME_2;
  private final static int FORWARD_FRAME_3 = FeatureSegment.FORWARD_FRAME_3;
  private final static int REVERSE_FRAME_3 = FeatureSegment.REVERSE_FRAME_3;
  private final static int REVERSE_FRAME_2 = FeatureSegment.REVERSE_FRAME_2;
  private final static int REVERSE_FRAME_1 = FeatureSegment.REVERSE_FRAME_1;
  private final static int SCALE_LINE      = FeatureSegment.SCALE_LINE;


  /**
   *  The JScrollBar for this FeatureDisplay object.  We create the scrollbar
   *  as part of this object rather than in the EntryEdit component because we
   *  may need to change the parameters of the scrollbar later.
   **/
  private JScrollBar scrollbar = null;

  /** A scroll bar for changing the viewing scale. */
  private ZoomScrollBar scale_changer = null;

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
  final private Vector<DisplayAdjustmentListener> adjustment_listener_list = 
      new Vector<DisplayAdjustmentListener>();

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
  
  private boolean feature_stack_view = false;
  
  private short MAX_LINES_FEATURE_STACK = 10;
  
  /** expand / collapse stack view, 2 or 1 respectively */
  private int STACK_EXPAND_FACTOR = 1;
  
  /** visible features sorted by size */
  private FeatureVector visibleFeaturesSortBySize;
  /** visible features sorted by position */
  private FeatureVector visibleFeaturesSortByPosition;

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

  private final int scrollbar_style;
  
  private static Vector<String> contigKeys;
  private static Vector<String> allPossibleContigKeys;

  private Object[] protein_keys = { "CDS", 
                                    "BLASTCDS",
                                    "polypeptide",
                                    "pseudogenic_exon"};

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
    
    feature_stack_view =
      Options.getOptions().getPropertyTruthValue("feature_stack_view");

    show_labels =
      Options.getOptions().getPropertyTruthValue("feature_labels");

    show_reverse_lines =
      Options.getOptions().getPropertyTruthValue("show_reverse_lines");
    
    show_forward_lines =
      Options.getOptions().getPropertyTruthValue("show_forward_lines");
    
    final StringVector frame_line_features = 
      Options.getOptions().getOptionValues("frame_line_features");
    if(frame_line_features != null)
      protein_keys = frame_line_features.toArray();
    
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
    addListeners();

    needVisibleFeatureVectorUpdate();

    getSelection().addSelectionChangeListener(this);
    getGotoEventSource().addGotoListener(this);

    getEntryGroup().addEntryGroupChangeListener(this);
    getEntryGroup().addEntryChangeListener(this);
    getEntryGroup().addFeatureChangeListener(this);

    getBases().addSequenceChangeListener(this, Bases.MIN_PRIORITY);

    Options.getOptions().addOptionChangeListener(this);
    setBackground(Color.white);

    DragSource dragSource = DragSource.getDefaultDragSource();

    dragSource.createDefaultDragGestureRecognizer(
       this,                             // component where drag originates
       DnDConstants.ACTION_COPY_OR_MOVE, // actions
       this);                            // drag gesture recognizer

    setDropTarget(new DropTarget(this,this));
  }


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
  protected void setShowLabels(boolean show_labels) 
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
  protected boolean getShowLabels() 
  {
    return show_labels;
  }

  /**
   *  Set value of the "show forward frame lines" flag.
   *  @param show_forward_lines Show forward frame lines if and only if
   *    this argument is true.
   **/
  protected void setShowForwardFrameLines(boolean show_forward_lines) 
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
  protected boolean getShowForwardFrameLines() 
  {
    return show_forward_lines;
  }

  /**
   *  Set value of the "show reverse frame lines" flag.
   *  @param show_reverse_lines Show frame lines if and only if this
   *    argument is true.
   **/
  protected void setShowReverseFrameLines(boolean show_reverse_lines) 
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
  protected boolean getShowSourceFeatures() 
  {
    return show_source_features;
  }

  /**
   *  Set value of the "show source features" flag.
   *  @param show_source_features Show features with a "source" key if and
   *    only if this argument is true.
   **/
  protected void setShowSourceFeatures(boolean show_source_features) 
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
  protected boolean getShowReverseFrameLines()
  {
    return show_reverse_lines;
  }

  /**
   *  Set value of the show base colours flag.
   *  @param show_base_colours At scale_factor less than two show each base in
   *    a different colour if and only if this argument is true.
   **/
  protected void setShowBaseColours(boolean show_base_colours) 
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
  protected boolean getShowBaseColours() 
  {
    return show_base_colours;
  }

  /**
   *  Set value of the "one line per entry" flag.
   *  @param one_line_per_entry If true then each entry will be shown on a
   *    different line, instead of showing frame lines.
   **/
  protected void setOneLinePerEntry(final boolean one_line_per_entry)
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
  protected boolean getOneLinePerEntryFlag() 
  {
    return one_line_per_entry;
  }
  
  
  /**
   *  Set value of the feature stack view flag.
   *  @param feature_stack_view If true then each CDS will be shown on a
   *    different line, instead of showing frame lines.
   **/
  protected void setFeatureStackViewFlag(final boolean feature_stack_view)
  {
    if(this.feature_stack_view != feature_stack_view) 
    {
      this.feature_stack_view = feature_stack_view;
      fixCanvasSize();
    }
  }

  /**
   *  Get the value of the "one line per entry" flag.
   **/
  protected boolean getFeatureStackViewFlag() 
  {
    return feature_stack_view;
  }

  /**
   *  Set value of the "hard left edge" flag.
   *  @param hard_left_edge If true the there will never be a gap between the
   *    left edge of the screen and the first visible base.  If false base one
   *    can be moved to the centre of the display.
   **/
  protected void setHardLeftEdge(final boolean hard_left_edge) 
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
   *  Set value of the show stop codons flag.
   *  @param show_stop_codons Show stop codons if and only if this argument is
   *    true.
   **/
  protected void setShowStopCodons(boolean show_stop_codons) 
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
  protected boolean getShowStopCodons() 
  { 
    return show_stop_codons;
  }

  /**
   *  Set value of the show start codons flag.
   *  @param show_start_codons Show start codons if and only if this argument
   *    is true.
   **/
  protected void setShowStartCodons(boolean show_start_codons) 
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
  protected boolean getShowStartCodons() 
  {
    return show_start_codons;
  }

  /**
   *  Set value of the reverse complement display flag.
   *  @param show_start_codons Draw all features and sequence reverse
   *    complemented if and only if this argument is true.
   **/
  protected void setRevCompDisplay(boolean rev_comp_display) 
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
  protected boolean isRevCompDisplay() 
  {
    return rev_comp_display;
  }

  /**
   *  Set value of the show feature arrows flag.
   *  @param show_feature_arrows Show directional arrows if and only if this
   *    argument is true.
   **/
  protected void setShowFeatureArrows(boolean show_feature_arrows)
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
  protected boolean getShowFeatureArrows() 
  {
    return show_feature_arrows;
  }

  /**
   *  Set value of the show feature borders flag.
   *  @param show_feature_borders Draw a border around each feature if and
   *    only if this argument is true.
   **/
  protected void setShowFeatureBorders(boolean show_feature_borders) 
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
  protected boolean getShowFeatureBorders() 
  {
    return show_feature_borders;
  }

  /**
   *  Set value of the show frame features flag.
   *  @param frame_features_flag All features(not just CDS features) will be
   *    drawn on the frame lines if and only if this argument is true.
   **/
  protected void setFrameFeaturesFlag(boolean frame_features_flag) 
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
  protected boolean getFrameFeaturesFlag() 
  {
    return frame_features_flag;
  }

  /**
   *  Set the value of the minimum score for this FeatureDisplay - features
   *  that have a /score lower than this value are never shown.
   **/
  protected void setMinimumScore(final int minimum_score) 
  {
    current_min_score = minimum_score;
    needVisibleFeatureVectorUpdate();
    repaint();
  }

  /**
   *  Return the value of the minimum score for this FeatureDisplay - see
   *  setMinimumScore().
   **/
  private int getMinimumScore() 
  {
    return current_min_score;
  }

  /**
   *  Set the value of the maximum score for this FeatureDisplay - features
   *  that have a /score higher than this value are never shown.
   **/
  protected void setMaximumScore(final int maximum_score) 
  {
    current_max_score = maximum_score;
    needVisibleFeatureVectorUpdate();
    repaint();
  }

  /**
   *  Return the value of the maximum score for this FeatureDisplay - see
   *  setMaximumScore().
   **/
  private int getMaximumScore() 
  {
    return current_max_score;
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
      case ARROW_LEFT:
        feature_display.setFirstBase(
            feature_display.getFirstVisibleForwardBase()-
            feature_display.scrollbar.getUnitIncrement());
        break;
      case ARROW_RIGHT:
        feature_display.setFirstBase(
            feature_display.getFirstVisibleForwardBase()+
            feature_display.scrollbar.getUnitIncrement());
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
  protected void setScaleFactor(int scale_factor) 
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
      
      if(scale_changer != null)
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
    else if(event.getType() != SequenceChangeEvent.CONTIG_REVERSE_COMPLEMENT)
      makeBaseVisibleInternal(event.getPosition(), true, false);

    fixScrollbar();
    needVisibleFeatureVectorUpdate();
    repaint();

    if(event.getType() == SequenceChangeEvent.REVERSE_COMPLEMENT ) 
      fireAdjustmentEvent(DisplayAdjustmentEvent.REV_COMP_EVENT);
    else if(event.getType() == SequenceChangeEvent.CONTIG_REVERSE_COMPLEMENT ) 
    {
      final Range range = event.getRange();
      
      fireAction(adjustment_listener_list, 
                  new DisplayAdjustmentEvent(FeatureDisplay.this,
                           range.getStart(),
                           range.getEnd(),
                           getMaxVisibleBases(),
                           getScaleValue(), getScaleFactor(),
                           isRevCompDisplay(), 
                           DisplayAdjustmentEvent.CONTIG_REV_COMP_EVENT));
    }
    else if(event.getType() == SequenceChangeEvent.CONTIG_REORDER)
    {
      final Range range  = event.getRange();
      final int drop_position = event.getPosition();

      fireAction(adjustment_listener_list,
                  new DisplayAdjustmentEvent(
                           FeatureDisplay.this,
                           range.getStart(), range.getEnd(),
                           drop_position,
                           DisplayAdjustmentEvent.CONTIG_REORDER));
    }   
    else 
      fireAdjustmentEvent(DisplayAdjustmentEvent.SCROLL_ADJUST_EVENT);
  }

  /**
   *  Invoked when an Option is changed.
   **/
  public void optionChanged(OptionChangeEvent event) 
  {
    AminoAcidSequence.setGeneCode();
    getBases().clearCodonCache();
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
    if(event.getType() == DisplayAdjustmentEvent.CONTIG_REORDER)
      return;

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
  protected MarkerRange getVisibleMarkerRange() 
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
  protected void raiseFeature(Feature feature) 
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
  protected void lowerFeature(Feature feature) 
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
  protected void smallestToFront() 
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
    //int max_visible_bases = getVisibleRange().getCount()-1; 
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
   *  to the vector of visible features (if it is visible that is).
   *  @param feature The object to add
   **/
  private void add(Feature feature) 
  {
    if(getEntryGroup().isActive(feature.getEntry()) &&
       featureVisible(feature)) 
    {
      if(!visible_features.contains(feature)) 
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
  protected FeatureVector getCurrentVisibleFeatures() 
  {
    return (FeatureVector)visible_features.clone();
  }

  /**
   *  Returns a Range that starts at the first visible base and ends at the
   *  last visible base.
   **/
  private Range getVisibleRange()
  {
    final int first_visible_base = getFirstVisibleForwardBase();
    final int last_visible_base  = getLastVisibleForwardBase();

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

      final int all_features_size = all_features.size();
      for(int i = 0 ; i < all_features_size; ++i) 
        raiseFeature(all_features.elementAt(i));

      raise_selection_flag = false;
    }

    if(isRevCompDisplay())
    {
      final int first_visible_base = getFirstVisibleReverseBase();
      final int last_visible_base  = getLastVisibleReverseBase();
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
    final int visible_features_size = visible_features.size();
    for(int i = 0 ; i < visible_features_size; ++i) 
    {
      final Feature new_feature = visible_features.elementAt(i);
      if(real_visible_features.contains(new_feature)) 
        new_visible_features.addElementAtEnd(new_feature);
    }

    // add features that are in real_visible_features and not currently
    // in visible_features and are not selected(selected features will be
    // added last so that they stay on top).
    final int real_visible_features_size = real_visible_features.size();
    for(int i = 0 ; i < real_visible_features_size; ++i) 
    {
      final Feature new_feature = real_visible_features.elementAt(i);

      if(!visible_features.contains(new_feature) &&
         !getSelection().contains(new_feature)) 
        new_visible_features.addElementAtEnd(new_feature);
    }

    final FeatureVector selection_features = getSelection().getAllFeatures();

    // now add features that are in real_visible_features, are not in
    // visible_features and are selected (selected features are added last so
    // that they stay on top).
    for(int i = 0 ; i < real_visible_features_size; ++i) 
    {
      final Feature new_feature = real_visible_features.elementAt(i);
      if(!visible_features.contains(new_feature) &&
         selection_features.contains(new_feature))
      {
        try
        {
          new_visible_features.addElementAtEnd(new_feature);
        }catch(Error e){}
      }
    }

    visible_features = new_visible_features;
    update_visible_features = false;
  }

  /**
   *  This is used by getSortedFeaturesInRange().
   **/
  final private static Comparator<Feature> feature_comparator = new Comparator<Feature>() 
  {
    /**
     *  Compare two Objects with respect to ordering.
     *  @return a negative number if feature1_object is less than
     *    feature2_object ; a positive number if feature1_object is greater
     *    than feature2_object; else 0
     **/
    public int compare(final Feature feature1,
                       final Feature feature2) 
    {
      final int feature1_size = feature1.getBaseCount();
      final int feature2_size = feature2.getBaseCount();

      if(feature1_size > feature2_size) 
        return -1;
      else if(feature1_size < feature2_size)
        return 1;
      else if(feature1.hashCode() < feature2.hashCode())  // use hash value as a last resort
        return -1;
      else if(feature1.hashCode() == feature2.hashCode())
        return 0;
      else
        return 1;
    }
  };
  
  final private static Comparator<Feature> feature_position_comparator = new Comparator<Feature>() 
  {
    /**
     *  Compare two Objects with respect to ordering.
     *  @return a negative number if feature1_object is less than
     *    feature2_object ; a positive number if feature1_object is greater
     *    than feature2_object; else 0
     **/
    public int compare(final Feature f1,
                       final Feature f2) 
    {
      final int pos1 = f1.getFirstBase();
      final int pos2 = f2.getFirstBase();

      if(pos1 > pos2) 
        return -1;
      else if(pos1 < pos2)
        return 1;
      else if(f1.hashCode() < f2.hashCode())  // use hash value as a last resort
        return -1;
      else if(f1.hashCode() == f2.hashCode())
        return 0;
      else
        return 1;
    }
  };

  /**
   *  Return a vector containing the references of the Feature objects in the
   *  EntryGroup that are within the given range are which should be
   *  displayed.  source features are not returned unless show_source_features
   *  is true.  features with a score less than getMinimumScore() or higher
   *  than getMaximumScore() aren't returned.
   *  
   *  GFFStreamFeature features check whether they are set to be visible.
   *  
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

      // filter out low and high scoring features and (possibly) source
      // features
      for(int i = features_from_entry.size() - 1; i >= 0; --i)
      {
        final Feature this_feature = features_from_entry.elementAt(i);

        if(this_feature.getKey().equals("source") &&
            !getShowSourceFeatures() &&
            !getSelection().contains(this_feature)) 
          continue;

        if(this_feature.getEmblFeature() instanceof GFFStreamFeature &&
           !((GFFStreamFeature)this_feature.getEmblFeature()).isVisible())
          continue;
          
        if(min_score > 0 || max_score < 100) 
        {
          final int this_feature_score = this_feature.getScore();

          // features with no /score are always shown
          if(this_feature_score != -1 &&
             (this_feature_score < min_score ||
              this_feature_score > max_score) &&
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

//  System.out.println("1 "+ System.currentTimeMillis());
    int scrollbar_hgt = 0;
    if(scrollbar_style == SCROLLBAR_AT_TOP)
    {
      scrollbar_hgt = scrollbar.getPreferredSize().height;
      ((Graphics2D)g).translate(0,scrollbar_hgt);
    }

    if(update_visible_features) 
    {
      updateVisibleFeatureVector();
      visibleFeaturesSortBySize = null;
      visibleFeaturesSortByPosition = null;
    }

    fillBackground(g);

    g.setFont(getFont());
//  System.out.println("2 "+ System.currentTimeMillis());
    final Selection selection = getSelection();
    final FeatureVector selected_features = selection.getAllFeatures();

    final FeatureSegmentVector selected_segments =
      selection.getSelectedSegments();

    int num_visible_features = getVisibleFeatures().size();
    final int segment_height = getFeatureHeight() - 1;
    final int seq_length = getSequenceLength();

//  System.out.println("3 "+ System.currentTimeMillis());
    FontMetrics fm = g.getFontMetrics();

    // we draw the feature backgrounds first then the visible indication
    // that there is a MarkerRange selected, then the feature outlines.
    final Vector<SegmentBorder> segment_borders = 
        new Vector<SegmentBorder>(num_visible_features);

    for(int i = 0; i < num_visible_features; ++i)
      drawFeature(g, segment_borders, getVisibleFeatures().elementAt(i),
                  true, selected_features, selected_segments,
                  segment_height, seq_length, fm);

//  System.out.println("4 "+ System.currentTimeMillis());
    drawBaseSelection(g);
    drawScale(g);
    drawCodons(g);
    drawBases(g);

//  System.out.println("5 "+ System.currentTimeMillis());

    // draw the segment borders
    g.setColor(Color.black);
    final int arrowWidth = getFontWidth() * 8 / 10;
    for(SegmentBorder sb: segment_borders)
      sb.drawSegmentBorder(g, segment_height, arrowWidth);

//  System.out.println("6 "+ System.currentTimeMillis());

    if(scrollbar_style == SCROLLBAR_AT_TOP)
      ((Graphics2D)g).translate(0,-scrollbar_hgt);
 
// draw drag and drop line
    if(highlight_drop_base > 0)
    { 
      g.setColor(Color.red);
      final int draw_x_position = getLowXPositionOfBase(highlight_drop_base);
      int nlines = 16;
       
      if(!show_forward_lines)
        nlines -= 6; 
      if(!show_reverse_lines)
        nlines -= 6;

      g.drawLine(draw_x_position, 0,
                 draw_x_position, (nlines*getFontHeight()));
    }
//  Thread.yield();

    if(getFeatureStackViewFlag())
    {
      int w = fm.charWidth('+');
      g.setColor(Color.RED);
      if(STACK_EXPAND_FACTOR == 1)
        g.drawString("+", getDisplayWidth()-w-2, 10);
      else
        g.drawString("-", getDisplayWidth()-w-2, 10);
    }
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

    if(getFeatureStackViewFlag())
    {
      
    }
    else if(getOneLinePerEntryFlag()) 
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
    int fill_line_top = fill_line_number*getFontHeight() + 1;

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
    final int scale_number_y_pos = scale_line * getFontHeight();

    final float bases_per_pixel =
     (float)getMaxVisibleBases()/getWidth();

    final int base_label_spacing;

    if(getScaleFactor() == 0) 
    {
      // set the spacing so that the labels are at multiples of ten at the
      // lowest scale factor
      base_label_spacing =
          (int)Math.ceil(MINIMUM_LABEL_SPACING * bases_per_pixel / 10) * 10;
    } 
    else 
    {
      // set the spacing so that the labels are at multiples of 100 otherwise
      base_label_spacing =
          (int)Math.ceil(MINIMUM_LABEL_SPACING * bases_per_pixel / 100) * 100;
    }

    //final int label_spacing = (int)(base_label_spacing / bases_per_pixel);

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

    final int font_ascent   = getFontAscent();
    final int font_width    = getFontWidth();
    final int font_line_hgt = getFontHeight();

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
    
    if(getFeatureStackViewFlag())
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
    
    g.setFont(getFont());
    int yposition = forward_frame_line * getFontHeight();

    // draw fwd bases
    if(getScaleFactor() == 0)
    {
      if(!(g instanceof SVGGraphics2D))
        g.drawString(forward_visible_bases, offset * getFontWidth(),
            yposition + getFontAscent() + 1);
      else
      {
        // for svg graphics
        for(int i=0;i<forward_visible_bases.length();i++)
          g.drawString(String.valueOf(forward_visible_bases.charAt(i)), 
              (offset+i)*getFontWidth(), 
              yposition + getFontAscent() + 1);
      }
    }
    else
    {
      for(int base_index = 0; base_index < forward_sequence_length;
          ++base_index)
        drawOnePixelBase(g, forward_visible_bases.charAt(base_index),
                         offset + base_index, yposition);
    }

    final Range reverse_range = newRange(getFirstVisibleReverseBase(),
                                         getLastVisibleReverseBase());

    String reverse_visible_bases =
      reverse_strand.getSubSequence(reverse_range).toUpperCase();

    reverse_visible_bases = reverse(reverse_visible_bases);

    final int reverse_frame_line = getFrameDisplayLine(REVERSE_STRAND);
    final int reverse_sequence_length = reverse_visible_bases.length();
    yposition = reverse_frame_line * getFontHeight();

    // draw bwd bases
    if(getScaleFactor() == 0)
    {  
      if(!(g instanceof SVGGraphics2D))
        g.drawString(reverse_visible_bases, offset * getFontWidth(),
            yposition + getFontAscent() + 1);
      else
      {
        // for svg graphics
        for(int i=0;i<reverse_visible_bases.length();i++)
          g.drawString(String.valueOf(reverse_visible_bases.charAt(i)), 
              (offset+i)*getFontWidth(), 
              yposition + getFontAscent() + 1);
      }
    }
    else
    {
      for(int base_index = 0; base_index < reverse_sequence_length;
          ++base_index) 
        drawOnePixelBase(g, reverse_visible_bases.charAt(base_index),
                         offset + base_index, yposition);
    }
  }


  /** 
  *
  * Reverse the String.
  *
  */
  private String reverse(String string) 
  {
    StringBuffer sb = new StringBuffer(string);
    return sb.reverse().toString();
  }


  /**
   *  Draw the codons of the sequence into a graphics object.  If the scale is
   *  0 then the codon letters will be drawn, otherwise just the stop codons
   *  will marked.
   *  @param g The object to draw into.
   **/
  private void drawCodons(Graphics g) 
  {
    if(getOneLinePerEntryFlag() || getFeatureStackViewFlag()) 
      return;
    
    g.setColor(Color.black);

    if(getScaleFactor() == 0)  
    {
      if(show_forward_lines) 
        drawForwardCodonLetters(g);
        
      if(show_reverse_lines) 
        drawReverseCodonLetters(g);
    }
    else 
    {
      final int MAX_STOP_CODON_SCALE_FACTOR = 7;
      if((show_stop_codons || show_start_codons) &&
          getScaleFactor() <= MAX_STOP_CODON_SCALE_FACTOR) 
      {
        if(show_forward_lines) 
          drawCodons(g,true);

        if(show_reverse_lines) 
          drawCodons(g,false);
      }
    }
  }


  /**
   *  Mark the start and stop codons on the three forward frame lines.
   *  @param g The object to draw into.
   **/
  private void drawCodons(Graphics g, boolean fwd) 
  {
    final Strand strand;
    final int first_visible_base;
    final int end_base;
    final int direction;

    if(fwd)
    {
      direction = FORWARD;
      if(isRevCompDisplay()) 
        strand = getBases().getReverseStrand();
      else 
        strand = getBases().getForwardStrand();
      first_visible_base = getForwardBaseAtLeftEdge();
      // base to end translation at
      // we + 3 to the upper bound because partial codons do not get translated
      // by getTranslation()
      end_base   = getLastVisibleForwardBase() + 3;
    }
    else
    {
      direction = REVERSE;
      if(isRevCompDisplay()) 
        strand = getBases().getForwardStrand();
      else 
        strand = getBases().getReverseStrand();
      first_visible_base = getFirstVisibleReverseBase();
      end_base   = getLastVisibleReverseBase() + 3;
    }
     
    final int frame_shift = (first_visible_base - 1) % 3;

    // base to start translation at - we start slightly off the 
    // left of the screen
    int start_base = first_visible_base - frame_shift;

    // not used if show_stop_codons is false
    int [][] stop_codons = null;

    if(show_stop_codons) 
    {
      if(fwd && start_base < 1)
        start_base = 1;
 
      stop_codons = strand.getStopOrStartCodons(newRange(start_base, end_base), 
                                                null);
//    stop_codons = new int [][] 
//    {
//      strand.getStopCodons(newRange(start_base, end_base)),
//      strand.getStopCodons(newRange(start_base + 1, end_base)),
//      strand.getStopCodons(newRange(start_base + 2, end_base))
//    };
    }

    // not used if show_start_codons is false
    int [][] start_codons = null;

    if(show_start_codons) 
    {
      final StringVector starts = Options.getOptions().getStartCodons();
      
      start_codons = strand.getStopOrStartCodons(newRange(start_base, end_base), 
                                                 starts);
//      start_codons = new int [][] 
//      {
//        strand.getMatchingCodons(newRange(start_base, end_base),
//                                 starts),
//        strand.getMatchingCodons(newRange(start_base + 1, end_base),
//                                 starts),
//        strand.getMatchingCodons(newRange(start_base + 2, end_base),
//                                 starts)
//      };
    }

    for(int i = 0 ; i < 3 ; ++i) 
    {
      final int this_frame_line;
      if(fwd)
        this_frame_line = getFrameDisplayLine(FORWARD_FRAME_1 + i);
      else
        this_frame_line = getFrameDisplayLine(REVERSE_FRAME_1 - i);

      if(show_start_codons) 
      {
        final int[] this_frame_start_codons = start_codons[i];

        drawCodonMarkLine(g, this_frame_line,
                          this_frame_start_codons,
                          direction, 80);
      }

      if(show_stop_codons) 
      {
        final int[] this_frame_stop_codons = stop_codons[i];

        drawCodonMarkLine(g, this_frame_line,
                          this_frame_stop_codons,
                          direction, 100);
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

    // base to end translation at - we end slightly off the 
    // right of the screen
    final int end_base   = getLastVisibleForwardBase() + 1;

    for(int i = 0 ; i < 3 ; ++i) 
    {
      final int frame_shift = 1 - (first_visible_base + 3 - i) % 3;

      // base to start translation at - we start slightly off the left of the
      // screen so that the first partial codon is translated as '.'
      final int start_base = first_visible_base + frame_shift;
      final int frame_line = getFrameDisplayLine(FORWARD_FRAME_1 + i);

      final AminoAcidSequence this_frame_translation =
                   strand.getSpacedTranslation(newRange(start_base, end_base), false);

      final String this_frame_translation_string =
                   this_frame_translation.toString();

      drawCodonLine(g, frame_shift, frame_line,
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

    // base to end translation at - we end slightly off the 
    // left of the screen
    final int first_visible_base = getFirstVisibleReverseBase();
    final int end_base = getLastVisibleReverseBase() + 1;

    for(int i = 0 ; i < 3 ; ++i) 
    {
      final int frame_shift = 1 - (first_visible_base + 3 - i) % 3;

      // base to start translation at - we start slightly off the right of the
      // screen so that the first partial codon is translated as '.'
      final int start_base = first_visible_base + frame_shift;
      final int frame_line = getFrameDisplayLine(REVERSE_FRAME_1 - i);

      final AminoAcidSequence this_frame_translation =
                    strand.getSpacedTranslation(newRange(start_base, end_base), false);

      final String this_frame_translation_string =
                    this_frame_translation.toString();

      drawCodonLine(g, frame_shift, frame_line,
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
  private void drawCodonLine(Graphics g, int frame_start,
                             int line_number, String codons,
                             int direction) 
  {
    final int offset;

    if(getForwardBaseAtLeftEdge() < 1 && direction != REVERSE) 
      offset = 1 - getForwardBaseAtLeftEdge();
    else 
      offset = 0;

//  codons = codons.toUpperCase();
    final int draw_y_position = line_number * getFontHeight();
    final int draw_x_position;
    
    if(direction == REVERSE)
    {
      draw_x_position = getLowXPositionOfBase(getLastVisibleForwardBase()) -
           (int)((offset+frame_start+codons.length()) * getScaleValue());

      codons = reverse(codons);
    }
    else
       draw_x_position = (int)((offset + frame_start + 1) * getScaleValue());

    if(!(g instanceof SVGGraphics2D))
      g.drawString(codons, draw_x_position,
                 draw_y_position + getFontAscent() + 1);
    else
    {
      for(int i=0;i<codons.length();i++)
        g.drawString(String.valueOf(codons.charAt(i)), draw_x_position+(i*getFontWidth() ), 
          draw_y_position + getFontAscent() + 1);
    }
  }

  /**
   *  Draw one line of codons marks (stop codons or start codons).
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

    //final int first_visible_forward_base = getForwardBaseAtLeftEdge();
    //final int first_visible_reverse_base = getFirstVisibleReverseBase();
    //final int last_visible_forward_base = getLastVisibleForwardBase();
    //final double scale_value = getScaleValue();

    final int line_height = getFontHeight();

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

    if(height_percent < 100)
      g.setColor(start_codon_colour);
    else
      g.setColor(Color.black);


    int length = 0;
    if(direction != FORWARD)
      length = getBases().getLength();

    final int codon_positions_length = codon_positions.length;
    for(int i = 0; i < codon_positions_length; ++i) 
    {
      final int codon_base_start = codon_positions[i];

      // zero is the end of data marker
      if(codon_base_start == 0) 
        break;

      int draw_x_position;

      if(direction == FORWARD) 
        draw_x_position = getLowXPositionOfBase(codon_base_start+2);
      else
      { 
        final int raw_base_position = length - codon_base_start + 1;
//      final int raw_base_position =
//        getBases().getComplementPosition(codon_base_start);
        draw_x_position = getLowXPositionOfBase(raw_base_position);
      }

      // draw only if we haven't drawn on the position already
      if(draw_x_position != last_x_position || last_x_position == -1) 
      {
        drawOneCodonMark(g, draw_x_position, draw_y_position,
                         direction, mark_height);

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
   **/
  private void drawOneCodonMark(final Graphics g,
                                 final int x_pos, final int y_pos,
                                 final int direction,
                                 final int height)
  {
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
    g.drawLine(x_pos, y_pos, x_pos, y_pos + getFontHeight() - 1);
  }
  
  /**
   * Get the 'feature stack view' line number for a feature.
   * @param thisCDSFeature
   * @param parentId
   * @param predicate
   * @param sortByPosition
   * @return
   */
  private int getFeatureStackLineNumber(final Feature thisCDSFeature, 
                                        final String parentId, 
                                        final FeaturePredicate predicate, 
                                        final boolean sortByPosition)
  {
    final String sysName = thisCDSFeature.getSystematicName();
    final Location loc = thisCDSFeature.getLocation();
    final Range range = loc.getTotalRange();
    final FeatureVector features;
    
    if(sortByPosition)
    {
      if(visibleFeaturesSortByPosition == null)
        visibleFeaturesSortByPosition = getVisibleFeatures().sort(feature_position_comparator);
      features = visibleFeaturesSortByPosition;
    }
    else
    {
      if(visibleFeaturesSortBySize == null)
        visibleFeaturesSortBySize = getVisibleFeatures().sort(feature_comparator);
      features = visibleFeaturesSortBySize;
    }
    
    int cnt = 0;
    for(int i=0; i<features.size(); i++)
    {
      final Feature f = features.elementAt(i);

      if( range.getStart() == f.getLocation().getTotalRange().getStart() )
      {
        if(f.getSystematicName().equals(sysName))
        {
          if(f.getKey().equals(thisCDSFeature.getKey()))
          {
            if(parentId != null && parentId.equals(getParentQualifier(f)))
              break;
            else if(thisCDSFeature == f)
              break; 
          }
        }
      }

      if(predicate.testPredicate(f) && 
         range.fuzzyOverlaps(f.getLocation().getTotalRange(), 10))
        cnt++;

      if(cnt > MAX_LINES_FEATURE_STACK)
        break;
    }

    if(((cnt*STACK_EXPAND_FACTOR)+8)  > MAX_LINES_FEATURE_STACK)
      cnt = (MAX_LINES_FEATURE_STACK-8)/STACK_EXPAND_FACTOR;
    return cnt;
  }
  
  /**
   * Return locus_tag, Parent (GFF) or transcript_id (GTF) qualifier or
   * null.
   * @param f
   * @return
   */
  public static String getParentQualifier(final Feature f)
  {
    try
    {
      if(f.getQualifierByName("locus_tag") != null)
        return (String) f.getQualifierByName("locus_tag").getValues().get(0);
      else if(f.getQualifierByName("Parent") != null)
        return (String) f.getQualifierByName("Parent").getValues().get(0);
      else if(f.getQualifierByName("transcript_id") != null)
        return (String) f.getQualifierByName("transcript_id").getValues().get(0);
    }
    catch (InvalidRelationException e){}
    return null;
  }
  
  /**
   * Return true if the feature is a CDS, pseudogenic_exon or an exon
   * that has been added as a frame line feature. This defines what
   * feature types get stacked in the feature stack view.
   * @param f       feature
   * @param keyStr  feature key
   * @return
   */
  private boolean isStackingFeature(final Feature f)
  {
    if( isProteinFeature(f) )
      return true;
    return false;
  }
  
  /**
   *  Return the line on the canvas where this feature segment should be drawn.
   *  @param segment The segment in question.
   *  @return The line to draw into.
   **/
  private int getSegmentDisplayLine(FeatureSegment segment)
  {
    if(getFeatureStackViewFlag())
    {
      final FeaturePredicate predicate = new FeaturePredicate(){
        public boolean testPredicate(Feature f)
        {
          if(isStackingFeature(f))
            return true;
          return false;
        }
      };
      
      String keyStr = segment.getFeature().getKey().toString();
      if(keyStr.equals("gene") || keyStr.equals("pseudogene"))
        return getLineCount()-2;
      else if( isStackingFeature(segment.getFeature()) )
      {
        final String parentId = getParentQualifier(segment.getFeature());
        return getLineCount()-4-( getFeatureStackLineNumber(
            segment.getFeature(), parentId, predicate, false) * STACK_EXPAND_FACTOR);
      }
      else if(keyStr.indexOf("utr") > -1 || keyStr.indexOf("UTR") > -1 )
      {
        try
        {
          // try to identify with associated CDS and return same line
          final String parentId = getParentQualifier(segment.getFeature());
          final FeatureVector visibleFeatures = getVisibleFeatures().sort(feature_comparator);
          for(int i=0; i<visibleFeatures.size(); i++)
          {
            final Feature f = visibleFeatures.elementAt(i);
            if(isStackingFeature(f))
            {
              if(parentId.equals( getParentQualifier(f) ))
              {
                int ln = getLineCount()-4-( getFeatureStackLineNumber(
                    f, parentId, predicate, false) * STACK_EXPAND_FACTOR );

                if(ln < 1)
                  ln = STACK_EXPAND_FACTOR;
                if(STACK_EXPAND_FACTOR == 2)
                  ln--;
                return ln;
              }
            }
          }
        }
        catch (Exception e){}
        return MAX_LINES_FEATURE_STACK-5;
      }
      else if(getAllPossibleContigKeys().contains(keyStr))
      {
        if(keyStr.equals("gap"))
          return 1;
        return 2;
      }
 
      ///
      ///
      final FeaturePredicate predicate2 = new FeaturePredicate(){
        public boolean testPredicate(Feature f)
        {
          final String thisKeyStr = f.getKey().getKeyString();
          if( !thisKeyStr.equals("gene") && !thisKeyStr.equals("pseudogene") &&
              !isStackingFeature(f) &&
              (thisKeyStr.indexOf("utr") == -1 || thisKeyStr.indexOf("UTR") == -1 ) &&
              !getAllPossibleContigKeys().contains(thisKeyStr) )
            return true;
          return false;
        }
      };
      int ln = getFeatureStackLineNumber(segment.getFeature(), null, predicate2, true);
      if(ln % 2 == 0)
        return 3;

      return 4;
    }
    else if(getOneLinePerEntryFlag()) 
    {
      final Feature feature = segment.getFeature();
      if((isProteinFeature(feature) || frame_features_flag) &&
          (show_forward_lines && (segment.isForwardSegment() ^
                                  isRevCompDisplay()) ||
           show_reverse_lines && (!segment.isForwardSegment() ^
                                  isRevCompDisplay())))
      {
        return getFeatureDisplayLine(feature,segment);
      } 
      else
      {
        if(segment.isForwardSegment() ^ isRevCompDisplay())
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

  
  protected Object[] getProteinKeys()
  {
    return protein_keys;
  }

  protected void setProteinKeys(Object[] protein_keys)
  {
    this.protein_keys = protein_keys; 
    repaint();
  }

  private boolean isProteinFeature(Feature feature)
  {
    final String key = feature.getKey().toString();
    
    if(isExonOfNonCodingTranscript(feature, key))
      return false;
     
    for(int i=0; i<protein_keys.length; i++)
    {
      if(key.equals((String)protein_keys[i]))
      {
        return true;
      }
    }
    return false;
  }
  
  /**
   * Check if this feature is an exon and is a child of a non-coding transcript
   * and is a GFF3 feature.
   * @param feature
   * @param key
   * @return
   */
  private boolean isExonOfNonCodingTranscript(final Feature feature, final String key)
  {
    if(key.equals(DatabaseDocument.EXONMODEL) && 
       feature.getEmblFeature() instanceof GFFStreamFeature)
    {
      final String nonCodingTranscripts[] = GeneUtils.getNonCodingTranscripts();
      try
      {
        Qualifier qualifier = feature.getQualifierByName("Parent");
        if(qualifier != null)
        {
          final ChadoCanonicalGene chadoGene = 
            ((GFFStreamFeature)feature.getEmblFeature()).getChadoGene();
          final String transcriptName = qualifier.getValues().get(0);
          final GFFStreamFeature transcript = 
            (GFFStreamFeature)chadoGene.getFeatureFromId(transcriptName);
          final String transcriptKey = transcript.getKey().getKeyString();

          for(int i=0; i<nonCodingTranscripts.length; i++)
            if(nonCodingTranscripts[i].equals(transcriptKey))
              return true;
        }
      }
      catch(Exception e){}
    }
    return false;
  }

  /**
   *  Return the frame ID of to use when drawing the given segment.
   **/
  private int getSegmentFrameID(final FeatureSegment segment) 
  {
    final int frame_id = segment.getFrameID();
    final Feature feature = segment.getFeature();

    if((isProteinFeature(feature) || frame_features_flag) &&
       (show_forward_lines && (segment.isForwardSegment() ^
                                isRevCompDisplay())||
         show_reverse_lines && (!segment.isForwardSegment() ^
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
   *  total number of lines.  Each line is getFontHeight() pixels high.
   *  @return The line to draw into.
   **/
  private int getFeatureDisplayLine(final Feature feature,
                                    final FeatureSegment segment) 
  {
    final int entry_index = getEntryGroup().indexOf(feature.getEntry());
    return getDisplayLineOfEntryIndex(entry_index,
                                      segment.isForwardSegment() ^
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
        if(show_labels) 
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
      if(show_labels) 
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
   *  getFontHeight() pixels high.
   *  @param frame_id The frame ID.
   *  @return The line to draw into.
   **/
  private int getFrameDisplayLine(int frame_id) 
  {
    if(getFeatureStackViewFlag()) 
    {
      return 0;
    }
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

    if(getFeatureStackViewFlag())
      return MAX_LINES_FEATURE_STACK;
    else if(getOneLinePerEntryFlag())  // some number of entry lines
      extra_line_count = getEntryGroup().size();
    else    // three frame line
      extra_line_count = 3;

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
  private int getLineOffset(final int line) 
  {
    return line * getFontHeight();
  }

  /**
   *  Return the lowest on screen(with respect to the canvas) x coordinate of
   *  a base.  If the scale_factor is zero then one base will be font_width
   *  wide and this method will return a different value than
   *  getHighXPositionOfBase().  If scale_factor is greater than one then
   *  the two methods will return the same thing.
   *  @param base_number The(forward) base to calculate the position of.
   **/
  private int getLowXPositionOfBase(final int base_number) 
  {
    return (int)((base_number - getForwardBaseAtLeftEdge()) *
                  getScaleValue());
  }

  /**
   *  Return the highest on screen(ie with respect to the canvas) x
   *  coordinate of a base.  See comment on getLowXPositionOfBase().
   *  @param base_number The(forward) base to calculate the position of.
   **/
  private int getHighXPositionOfBase(final int base_number) 
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

//  if(scrollbar_style == SCROLLBAR_AT_TOP)
//    position.x += scrollbar.getPreferredSize().height;

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
                           final Vector<SegmentBorder> segment_borders,
                           final Feature feature,
                           final boolean draw_feature_fill,
                           final FeatureVector selected_features,
                           final FeatureSegmentVector selected_segments, 
                           final int segment_height,
                           final int seq_length,
                           final FontMetrics fm) 
  {
    final FeatureSegmentVector this_feature_segments = feature.getSegments();
    final int num_segs = this_feature_segments.size();

    // don't try to draw a feature with no segments
    if(num_segs == 0) 
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
      highlight_feature_flag = false;
    }

    boolean highlight_segment_flag;
    boolean draw_direction_arrow_flag;
    FeatureSegment current_segment;
     
    if(show_labels)
      drawFeatureLabel(g, feature, seq_length, fm);

    // draw each segment/exon
    boolean trans_spliced = false;
    for(int i = 0; i < num_segs; ++i) 
    {
      current_segment = this_feature_segments.elementAt(i);

      if(selected_segments.indexOf(current_segment) == -1) 
        highlight_segment_flag = false;
      else 
        highlight_segment_flag = true;

      // draw an arrow only on the last segment
      if(i == num_segs - 1 && show_feature_arrows) 
        draw_direction_arrow_flag = true;
      else 
        draw_direction_arrow_flag = false;

      SegmentBorder fb = drawSegment(g, current_segment,
                                     highlight_feature_flag, highlight_segment_flag,
                                     draw_direction_arrow_flag, segment_height);
      if(fb != null)
        segment_borders.add(fb);

      // draw a line between the segments
      if(i + 1 < num_segs) 
      {
        final FeatureSegment next_segment =
          this_feature_segments.elementAt(i + 1);

        trans_spliced = drawSegmentConnection(g, current_segment, 
                                    next_segment, trans_spliced);
      }
    }

    // draw the label last if the is no label line because in this case the
    // label is draw on top of the feature segments
    if(!show_labels) 
      drawFeatureLabel(g, feature, seq_length, fm);
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
  private boolean drawSegmentConnection(Graphics g,
                                     FeatureSegment lower_segment,
                                     FeatureSegment upper_segment,
                                     boolean trans_spliced) 
  {
    Marker upper_segment_start_marker = upper_segment.getStart();
    Marker lower_segment_end_marker   = lower_segment.getEnd();

    // trans-spliced
    if((upper_segment.isForwardSegment() ^ lower_segment.isForwardSegment()) ||
       trans_spliced)
    {
      trans_spliced = true;
      int start = upper_segment.getStart().getRawPosition();
      int end   = upper_segment.getEnd().getRawPosition();
      if(end < start)
        upper_segment_start_marker = upper_segment.getEnd();

      start = lower_segment.getStart().getRawPosition();
      end   = lower_segment.getEnd().getRawPosition();

      if(end < start)
        lower_segment_end_marker = lower_segment.getStart(); 
    }

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

    return trans_spliced;
  }

  /**
   *  Draw a label a feature.  This is a helper function for drawFeature().
   *  If show_labels is true the labels will be drawn below the features,
   *  otherwise they will be drawn within the features.
   *  @param g The Graphics object on which to draw.
   *  @param feature The feature to draw the label for.
   **/
  private void drawFeatureLabel(Graphics g, final Feature feature, 
                                final int seq_length,
                                final FontMetrics fm) 
  {
    // the three frame translation is visible when the scale factor is 0,
    // don't draw labels over it
    if(!show_labels && getScaleFactor() == 0) 
      return;

    final String label = feature.getLabel();

    // special case - don't display a label if the label qualifier is "*"
    if(label != null && label.equals("*"))
      return;
  
    final String label_or_gene = feature.getIDString(); 

    // don't waste time drawing nothing
    if(label_or_gene.length() == 0)
      return;

    final int string_width = fm.stringWidth(label_or_gene); 
    final FeatureSegment first_segment = feature.getSegments().elementAt(0);
    final int label_x_coord;

    if(first_segment.isForwardSegment() ^ isRevCompDisplay())
    {
      int segment_start_pos =
        first_segment.getStart().getRawPosition();

      if(isRevCompDisplay()) 
        segment_start_pos = seq_length - segment_start_pos + 1;

      label_x_coord = getLowXPositionOfBase(segment_start_pos);
    } 
    else
    {
      int segment_end_pos =
        first_segment.getEnd().getRawPosition();

      if(isRevCompDisplay()) 
        segment_end_pos = seq_length - segment_end_pos + 1;

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
      vertical_offset += getFontHeight(); // move to the label line

    // save this so we can restore it later
    Shape saved_clip = null;
    if(!show_labels)
      saved_clip = g.getClip();
 
    // if we have a labels line draw a white background behind the
    // label
    if(show_labels) 
    {
      g.setColor(Color.white);

      g.fillRect(label_x_coord - getFontWidth(),
                 vertical_offset+2,
                 string_width + getFontWidth() * 2,
                 getFontHeight()-1);
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
    g.drawString(label_or_gene, label_x_coord + 1,
                 vertical_offset + getFontAscent() + 1);

    if(!show_labels)
      g.setClip(saved_clip);

  }

  /**
   *  Return the position on the canvas where this segment starts.  The
   *  returned value will be -1 if the position is off the left of the screen
   *  and will be (width of the canvas) if the position is off the right of
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
    else if(segment_start_coord < 0) 
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
   *  Return if and only if the segment is (partly) within the range of bases
   *  we are currently displaying.
   *  @param segment The FeatureSegment to check.
   **/
  private boolean segmentVisible(FeatureSegment segment) 
  {
    final int segment_start_coord = getSegmentStartCoord(segment);
    final int segment_end_coord = getSegmentEndCoord(segment);
    final int width = getSize().width;
    
    if(segment_end_coord < 0 && segment_start_coord < 0 ||
        segment_start_coord >= width &&
        segment_end_coord >= width) 
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
  private SegmentBorder drawSegment(Graphics g, FeatureSegment segment,
                                    final boolean highlight_feature,
                                    final boolean highlight_segment,
                                    final boolean draw_arrow,
                                    final int segment_height) 
  {
    // not on screen
    if(!segmentVisible(segment)) 
      return null;

    final Feature segment_feature = segment.getFeature();
    final int vertical_offset = getSegmentVerticalOffset(segment) + 1;

    int segment_start_coord = getSegmentStartCoord(segment);
    int segment_end_coord   = getSegmentEndCoord(segment);

    // this is 1 if the feature is on the forward strand or on a forward frame
    // and -1 otherwise.  this used to draw the feature arrow in the right
    // direction.
    final int feature_direction;
    if(segment.isForwardSegment() ^ isRevCompDisplay())
    {
      feature_direction = 1;
      if(segment_end_coord < segment_start_coord)
      {
        int tmp = segment_start_coord;
        segment_start_coord = segment_end_coord;
        segment_end_coord   = tmp;
      }
    }
    else 
    {
      feature_direction = -1;
      if(segment_end_coord > segment_start_coord)
      {
        int tmp = segment_start_coord;
        segment_start_coord = segment_end_coord;
        segment_end_coord   = tmp;
      }
    }


    final Color feature_colour = segment_feature.getColour();

    // no colour means draw in white
    if(feature_colour == null) 
      g.setColor(Color.white);
    else 
      g.setColor(feature_colour);

    final int segment_width;
    if(feature_direction == 1)
    {
      segment_width = segment_end_coord - segment_start_coord + 1;

      g.fillRect(segment_start_coord, vertical_offset,
                 segment_width, segment_height + 1);
    }
    else
    {
      segment_width = segment_start_coord - segment_end_coord + 1;

      g.fillRect(segment_end_coord, vertical_offset,
                 segment_width, segment_height + 1);
    }

    if(!show_feature_borders && !highlight_feature)
      return null;

    if(feature_direction == 1)
      return new SegmentBorder(highlight_feature, highlight_segment, draw_arrow,
                               segment_start_coord, vertical_offset,
                               segment_width, 
                               feature_direction);
    else
      return new SegmentBorder(highlight_feature, highlight_segment, draw_arrow,
                               segment_end_coord, vertical_offset,
                               segment_width, 
                               feature_direction);
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

    if((event.getModifiers() & InputEvent.BUTTON2_DOWN_MASK) > 0) 
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
      private FeaturePopup popup = null;

      public void mouseEntered(MouseEvent event)
      {
        // grab focus to enable scrolling with ARROW_LEFT/RIGHT
        requestFocusInWindow(); 
      }
      
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
          if(popup == null)
            popup = new FeaturePopup(FeatureDisplay.this,
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
        {        
          if(clicked_feature.getEmblFeature() instanceof GFFStreamFeature &&
             ((GFFStreamFeature)clicked_feature.getEmblFeature()).getChadoGene() != null)
            new GeneBuilderFrame(clicked_feature, getEntryGroup(),
                            getSelection(), getGotoEventSource());
          else
          {
            final JFrame frame = new JFrame("Artemis Feature Edit: " + 
                clicked_feature.getIDString() +
                (clicked_feature.isReadOnly() ?
                    "  -  (read only)" :
                    ""));
            
            final FeatureEdit fe = new FeatureEdit(clicked_feature, getEntryGroup(),
                                       getSelection(), getGotoEventSource(), frame);
            
            frame.addWindowListener(new WindowAdapter() 
            {
              public void windowClosing(WindowEvent event) 
              {
                fe.stopListening();
                frame.dispose();
              }
            });
            
            frame.getContentPane().add(fe);
            frame.pack();

            final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
            frame.setLocation(new Point((screen.width - getSize().width)/2,
                                        (screen.height - getSize().height)/2));
            frame.setVisible(true);
          }
        }
      }
    }

    makeSelectionVisible();
  }

  /**
   *  Handle a single click on the canvas.
   **/
  private void handleCanvasSingleClick(MouseEvent event) 
  {
    if(getFeatureStackViewFlag())
    {
      if(event.getPoint().y < 10 &&
         event.getPoint().x > getDisplayWidth()-10)
      {
        if(STACK_EXPAND_FACTOR == 1) // expand/collapse
          STACK_EXPAND_FACTOR = 2;
        else
          STACK_EXPAND_FACTOR = 1;
        
        updateOneLinePerFeatureFlag();
        return;
      }
    }
    
    click_segment_marker = null;

    final Selectable clicked_thing = getThingAtPoint(event.getPoint());

    // note: ALT_MASK and BUTTON2_MASK are the same value
    // so check this is not BUTTON1 to allow BUTTON1+ALT to drag
    if( clicked_thing == null ||
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
        if(getFeatureStackViewFlag())
          return;
        
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

          if(this_segment.getStart().getStrand() ==  new_click_range.getStrand() &&
             this_segment.canDirectEdit())
          {
            if(event.isShiftDown() && getSelection().getSelectedSegments().size() == 0)
            {
              int baseGrab = 15;
              if(getScaleFactor() > 2)
                baseGrab = 30;

              if( (new_click_range.getStart().getPosition() >= this_segment.getStart().getPosition() &&
                   new_click_range.getStart().getPosition() <= this_segment.getStart().getPosition()+baseGrab) ||
                  (new_click_range.getEnd().getPosition() <= this_segment.getEnd().getPosition() &&
                   new_click_range.getEnd().getPosition() >= this_segment.getEnd().getPosition()-baseGrab) )
              {
                int distFromBeg = new_click_range.getStart().getPosition()-this_segment.getStart().getPosition();
                int distFromEnd = this_segment.getEnd().getPosition()-new_click_range.getEnd().getPosition();
                if(distFromBeg < distFromEnd)
                {
                  click_segment_marker = this_segment.getStart();
                  click_segment_marker_is_start_marker = true;
                  other_end_of_segment_marker = this_segment.getEnd();
                }
                else
                {
                  click_segment_marker = this_segment.getEnd();
                  click_segment_marker_is_start_marker = false;
                  other_end_of_segment_marker = this_segment.getStart();
                }

                getSelection().add(this_segment);
                getEntryGroup().getActionController().startAction();
                break;
              }
            }
            else
            {
              if (new_click_range.getStart().equals(this_segment.getStart()))
              {
                click_segment_marker = this_segment.getStart();
                click_segment_marker_is_start_marker = true;
                other_end_of_segment_marker = this_segment.getEnd();
                getEntryGroup().getActionController().startAction();
                break;
              }

              if (new_click_range.getEnd().equals(this_segment.getEnd()))
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
      return (Feature)object;
    else 
      // object must be a FeatureSegment
      return ((FeatureSegment)object).getFeature();
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
    if( event.isShiftDown() &&
       (getSelection().getAllFeatures().size() > 1 ||
        getSelection().getSelectedSegments().size() > 1)) 
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
            {
              //
              // exon-model returns the segment - this is a fix for merging chado 
              if(segment_feature.getKey().getKeyString().equals(DatabaseDocument.EXONMODEL))
                return current_segment;
              return segment_feature;
            }
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

          needVisibleFeatureVectorUpdate();
          updateOneLinePerFeatureFlag();
          
          fireAdjustmentEvent(DisplayAdjustmentEvent.SCROLL_ADJUST_EVENT);
          repaint();
        }
      }
    });

    Box box = Box.createHorizontalBox();
    box.add(scrollbar);

    if(UIManager.getLookAndFeel().getName().equals("Mac OS X Aqua"))
      box.add(Box.createHorizontalStrut(scrollbar.getPreferredSize().height));

    if(scrollbar_at_top) 
      add(box, "North");
    else 
      add(box, "South");
  }
  

  protected void updateOneLinePerFeatureFlag()
  {
    if(getFeatureStackViewFlag())
    {
      final Range visRange;
      if(isRevCompDisplay())
      {
        final int first_visible_base = getFirstVisibleReverseBase();
        final int last_visible_base  = getLastVisibleReverseBase();
        visRange = newRange(first_visible_base, last_visible_base);
      } 
      else 
        visRange = getVisibleRange();
      final FeatureVector visFeatures = getSortedFeaturesInRange(visRange);
      
      short maxOverlaps = 2;
      for(int i=0; i<visFeatures.size(); i++)
      {
        final Feature f1 = visFeatures.elementAt(i);
        if(!isProteinFeature(f1))
          continue;
        
        final Range r = f1.getMaxRawRange();
        short cnt = 0;
        for(int j=0; j<visFeatures.size(); j++)
        {
          final Feature f2 = visFeatures.elementAt(j);
          if( f1 == f2 )
            break;

          if(isProteinFeature(f2) && r.overlaps(f2.getMaxRawRange()))
            cnt++;
        }
        if(cnt > maxOverlaps)
          maxOverlaps = cnt;
      }

      if(((maxOverlaps*STACK_EXPAND_FACTOR)+10) != MAX_LINES_FEATURE_STACK)
      {
        MAX_LINES_FEATURE_STACK = (short) ((maxOverlaps*STACK_EXPAND_FACTOR)+10);
        if(MAX_LINES_FEATURE_STACK > 42)
          MAX_LINES_FEATURE_STACK = 42;
        fixCanvasSize();
      }
    }
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
   *  @param listeners A Vector of the objects that the event should be sent to.
   *  @param event The event to send
   **/
  private void fireAction(final Vector<DisplayAdjustmentListener> listeners, final ChangeEvent event) 
  {
    final Vector<DisplayAdjustmentListener> targets;
    // copied from a book - synchronizing the whole method might cause a
    // deadlock
    synchronized(this) 
    {
      //targets = (Vector)listeners.clone();
      targets = new Vector<DisplayAdjustmentListener>(listeners);
    }

    final int ntargets = targets.size();
    for(int i = 0; i < ntargets; ++i) 
    {
      DisplayAdjustmentListener target = targets.elementAt(i);
      target.displayAdjustmentValueChanged((DisplayAdjustmentEvent)event);
    }
  }

  /**
   *  Create the scroll bar used for changing the scale and add it to the
   *  FeatureDisplay.
   **/
  private void createScaleScrollbar() 
  {
    scale_changer = new ZoomScrollBar(this);
    
    if(System.getProperty("autohide") == null ||
       System.getProperty("autohide").equals("false") )
      add(scale_changer, "East");
    else
      scale_changer.addMouseMotionListenerToFeatureDisplay();
  }

  /**
   *  Update the parameters of the scrollbar taking changes to the entry_group
   *  into account.
   **/
  protected void fixScrollbar() 
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
    int new_height = getFontHeight() * getLineCount();

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
//    repaint();
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
      updateOneLinePerFeatureFlag();
      repaint();
    }
  }

  /**
   *  Return the current scale factor.  The scale factor is a number greater
   *  than or equal to zero than controls the number of bases that can appear
   *  on screen.  See getScaleValue().
   **/
  protected int getScaleFactor() 
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
    if(point.y >= getHeight() - scrollbar.getHeight()) 
      return -1;

    final int return_value = point.y / getFontHeight();

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
          line_height = getFontHeight() * 2;
        else 
          line_height = getFontHeight();

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
  protected int getFirstVisibleForwardBase() 
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
  /*private Marker getFirstVisibleForwardBaseMarker() 
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
  }*/

  /**
   *  Return the base number of the last visible base in the forward
   *  direction, ie the base on the forward strand that is on the right
   *  edge of the canvas.  The number returned will always be > 1 and < the
   *  sequence length.
   **/
  protected int getLastVisibleForwardBase() 
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
  protected float getScaleValue() 
  {
    return scale_value;
  }

  /**
   *  Return the height in pixels we should use for drawing features.
   **/
  private int getFeatureHeight() 
  {
    // we want a nice space between
    // frames/lines
    return getFontAscent() + 2;
  }

  /**
   *  Return the number of bases we can fit on screen at once, ie the number
   *  that will fit side by side on the canvas.
   **/
  public int getMaxVisibleBases() 
  {
    return (int)(getWidth()/getScaleValue());
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
  protected void setFirstAndLastBase(final int first, final int last) 
  {
    left_edge_base = first;
    setScaleValue(1.0F * getWidth() /(last - first + 1));
  }

  /**
   *  Scroll the display so that the given first base is at the left edge of
   *  the screen.
   **/
  protected void setFirstBase(int base_position) 
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
  protected void needVisibleFeatureVectorUpdate() 
  {
    update_visible_features = true;
  }


////////////////////
// DRAG AND DROP
////////////////////

  /**
  *
  *  Read an entry from a remote file node (ssh)
  *
  **/
  private void readAnEntryFromRemoteFileNode(final RemoteFileNode node)
  {
    SwingWorker entryWorker = new SwingWorker()
    {
      public Object construct()
      {
        try
        {
          EntryInformation new_entry_information =
             new SimpleEntryInformation(Options.getArtemisEntryInformation());

          final Entry entry =  new Entry(getEntryGroup().getBases(),
                           EntryFileDialog.getEntryFromFile(null,
                           new RemoteFileDocument(node),
                           new_entry_information, true));
          if(entry != null)
            getEntryGroup().add(entry);
        }
        catch(final OutOfRangeException e)
        {
          new MessageDialog(null,
                         "read failed: one of the features " +
                         "in the entry has an out of " +
                         "range location: " +
                         e.getMessage());
        }
        return null;
      }
    };
    entryWorker.start();
  }


  /**
   *  Read an entry
   **/
  private void readAnEntryFromFile(final File file,
                                   final EntryGroup entry_group)
  {
    SwingWorker entryWorker = new SwingWorker()
    {
      public Object construct()
      {
        try
        {
          EntryInformation new_entry_information =
             new SimpleEntryInformation(Options.getArtemisEntryInformation());

          final Entry new_entry =  new Entry(entry_group.getBases(),
                         EntryFileDialog.getEntryFromFile(null,
                          new FileDocument(file),
                          new_entry_information, true));

          if(new_entry != null)
            entry_group.add(new_entry);
        }
        catch(final OutOfRangeException e)
        {
          new MessageDialog(null,
                         "read failed: one of the features " +
                         "in the entry has an out of " +
                         "range location: " +
                         e.getMessage());
        }
        return null;
      }
    };
    entryWorker.start();
  }

  protected static Vector<String> getContigKeys()
  {
    if(contigKeys == null)
    {
      contigKeys = new Vector<String>(3);
      contigKeys.add("fasta_record");
      contigKeys.add("contig");
      contigKeys.add("insertion_gap");
    }
    
    return contigKeys;
  }
  
  protected static Vector<String> getAllPossibleContigKeys()
  {
    if(allPossibleContigKeys == null)
    {
      allPossibleContigKeys = new Vector<String>();
      allPossibleContigKeys.add("fasta_record");
      allPossibleContigKeys.add("contig");
      allPossibleContigKeys.add("insertion_gap");
      allPossibleContigKeys.add("gap");
      allPossibleContigKeys.add("scaffold");
      allPossibleContigKeys.add("source");
    }
    return allPossibleContigKeys;
  }

  /**
  *
  */
  protected FeatureVector getContigs()
  {
    final FeatureVector tmpContigFeatures = new FeatureVector();
    // find all fasta_record features
    final Vector<String> contigKeys = getContigKeys();
    final FeaturePredicate key_predicate_contig[]
           =  new FeatureKeyQualifierPredicate[contigKeys.size()];

    for(int i=0; i<contigKeys.size(); i++)
      key_predicate_contig[i] = new FeatureKeyQualifierPredicate(
                                             new Key(contigKeys.get(i)),
                                             null, // match any qialifier
                                             false);

    final FeatureEnumeration test_enumerator = getEntryGroup().features();

    while(test_enumerator.hasMoreFeatures())
    {
      final Feature this_feature = test_enumerator.nextFeature();

      for(int i=0; i<contigKeys.size(); i++)
      {
        if(key_predicate_contig[i].testPredicate(this_feature))
          tmpContigFeatures.add(this_feature);
      }
    }
    
    
    final FeaturePredicate subsetPredicate = new FeaturePredicate()
    {
      public boolean testPredicate(final Feature feature)
      {
        final Location loc = feature.getLocation();
        final int start = loc.getFirstBase();
        final int end   = loc.getLastBase();

        for(int i=0; i<tmpContigFeatures.size(); i++)
        {
          Feature contig = tmpContigFeatures.elementAt(i);
          if(contig.getSystematicName().equals(feature.getSystematicName()))
            continue;
          
          int this_start = contig.getLocation().getFirstBase();
          int this_end   = contig.getLocation().getLastBase();

          if( (this_start >= start && this_start < end) &&
              (this_end > start && this_end <= end) )
            return true;
        }
        return false;
      }
    };
    
    FeatureVector contigFeatures = new FeatureVector();
    for(int i=0; i<tmpContigFeatures.size(); i++)
    {
      Feature contig = tmpContigFeatures.elementAt(i);
      if(!subsetPredicate.testPredicate(contig))
        contigFeatures.add(contig);
    }
    
    tmpContigFeatures.removeAllElements();
    
    return contigFeatures;
  }

  /**
  *
  * Contig reordering. The current selection is moved to highlight_drop_base
  * base position in the sequence.
  *
  */
  protected void reorder(int drop_base, Feature selected_feature)
  {
    // rearrange contigs
    try
    {
      Sequence sequence = getBases().getSequence();
      FeatureVector contig_features = null;
      int old_pos[] = null;

      if(sequence instanceof FastaStreamSequence)
      {
        contig_features = getContigs();

        // get fasta_record old positions
        old_pos = new int[contig_features.size()];
        for(int i=0; i<old_pos.length; i++)
          old_pos[i] = contig_features.elementAt(i).getMaxRawRange().getStart()-1;
      }

      //
      // Test that there are no features overlapping contig boundaries
      //
      final FeaturePredicate overlapPredicate = new FeaturePredicate()
      {
        private FeatureVector contig_features = getContigs();
        public boolean testPredicate(final Feature feature)
        {
          //if(getContigKeys().contains(feature.getKey().getKeyString()))
          //  return false;
          
          final Location loc = feature.getLocation();

          final int start = loc.getFirstBase();
          final int end   = loc.getLastBase();

          for(int i=0; i<contig_features.size(); i++)
          {
            Feature contig = contig_features.elementAt(i);
            int this_start = contig.getLocation().getFirstBase();
            int this_end   = contig.getLocation().getLastBase();
            
            if( (this_start > start && this_start < end) ||
                (this_end > start && this_end < end) )
              return true;
          }
          return false;
        }
      };

      
      final FilteredEntryGroup filtered_entry_group =
        new FilteredEntryGroup(getEntryGroup(), overlapPredicate, "Overlapping Features");

      if(filtered_entry_group.getAllFeaturesCount() > 0)
      {
        JOptionPane.showMessageDialog(null,
            "There is/are "+filtered_entry_group.getAllFeaturesCount()+
            " feature(s) overlapping the contig boundaries."+
            "\nA list of these features will appear."+
            "\nThis needs fixing before contigs can be reordered.",
            "Warning",
            JOptionPane.WARNING_MESSAGE);
        
        FeatureListFrame list = new FeatureListFrame("Overlapping Features (to be fixed before reordering)",
            getSelection(), getGotoEventSource(), filtered_entry_group, getBasePlotGroup());
        list.setVisible(true);
        
        return;
      }
      
      // rearrange contig order
      getEntryGroup().getBases().contigRearrange(selected_feature,
                                                 drop_base);

      // get fasta_record new positions
      if(sequence instanceof FastaStreamSequence &&
         old_pos != null)
      {
        int new_pos[] = new int[old_pos.length];
        for(int i=0; i<new_pos.length; i++)
          new_pos[i] = contig_features.elementAt(i).getMaxRawRange().getStart()-1;

        // update header record
        ((RawStreamSequence)sequence).setFastaHeaderPosition(old_pos,new_pos);
      }
    }
    catch(ReadOnlyException roe)
    {
      final String message =
        "one or more of the features is read-only or is in a " +
        "read-only entry - cannot continue";
      new MessageDialog(null, message);
      highlight_drop_base = -1;

      return;
    }

  }
// drop
  protected static Border dropBorder = new BevelBorder(BevelBorder.LOWERED);

  private void getNearestFeatureEnd(Point loc)
  {
    final int base_pos = getBasePositionOfPoint(loc, FORWARD_STRAND);
    FeatureVector features_from_entry = getVisibleFeatures();   
    int first;
    int last;
    Vector<String> contig_keys = getContigKeys();

    for(int i = 0; i < features_from_entry.size(); i++)
    {
      final Feature this_feature = features_from_entry.elementAt(i);

      if(contig_keys.contains(this_feature.getKey()))
      {
        first = this_feature.getRawFirstBase();
        last  = this_feature.getRawLastBase();
 
        if( Math.abs(first - base_pos) < Math.abs(base_pos - highlight_drop_base) )
          highlight_drop_base = first;
        if( Math.abs(last - base_pos) < Math.abs(base_pos - highlight_drop_base) )
          highlight_drop_base = last+1;
      }
    }
  }

  public void drop(DropTargetDropEvent e)
  {
    Transferable t = e.getTransferable();
    try
    {
      if(t.isDataFlavorSupported(FileNode.FILENODE))
      {
        FileNode fn = (FileNode)t.getTransferData(FileNode.FILENODE);
        readAnEntryFromFile(fn.getFile(), getEntryGroup());
      }
      else if(t.isDataFlavorSupported(RemoteFileNode.REMOTEFILENODE))
      {
        final RemoteFileNode node =
            (RemoteFileNode)t.getTransferData(RemoteFileNode.REMOTEFILENODE);
        readAnEntryFromRemoteFileNode(node);
      }
      else if(e.isDataFlavorSupported(DataFlavor.stringFlavor))
      {
        final FeatureVector selected_features = getSelection().getAllFeatures();
        if(selected_features.size() > 1)
          return;
 
        reorder(highlight_drop_base, selected_features.elementAt(0));       // rearrange contigs
      }
      else
        e.rejectDrop();
    }
    catch(UnsupportedFlavorException ufe)
    {
      ufe.printStackTrace();
    }
    catch(IOException ioe)
    {
      ioe.printStackTrace();
    }
    finally
    {
      highlight_drop_base = -1;
      setBorder(null);
    }
  }

  public void dragExit(DropTargetEvent e)
  {
    highlight_drop_base = -1;
    setBorder(null);
    repaint();
  }

  public void dropActionChanged(DropTargetDragEvent e) {}

  public void dragOver(DropTargetDragEvent e)
  {
    if(e.isDataFlavorSupported(FileNode.FILENODE))
    {
      setBorder(dropBorder);
      e.acceptDrag(DnDConstants.ACTION_COPY_OR_MOVE);
      return;
    }
    else if(e.isDataFlavorSupported(DataFlavor.stringFlavor))
    {
      Point ploc = e.getLocation();
      getNearestFeatureEnd(ploc);
      repaint();
    }
    else
      e.rejectDrag();
  }

  public void dragEnter(DropTargetDragEvent e)
  {
    if(e.isDataFlavorSupported(FileNode.FILENODE))
      e.acceptDrag(DnDConstants.ACTION_COPY_OR_MOVE);
  }


// drag source
  public void dragGestureRecognized(DragGestureEvent e)
  {
    // ignore if mouse popup trigger
    InputEvent ie = e.getTriggerEvent();
    if(ie instanceof MouseEvent)
      if(((MouseEvent)ie).isPopupTrigger())
        return;

    if(getFeatureStackViewFlag() || getEntryGroup().isReadOnly())
      return;
    
    final Vector<String> contig_keys = getContigKeys();
    final FeatureVector selected_features = getSelection().getAllFeatures();
    if(selected_features.size() == 1 &&
       contig_keys.contains(selected_features.elementAt(0).getKey()))
    {
      ClassLoader cl = this.getClass().getClassLoader();
      ImageIcon icon = new ImageIcon(cl.getResource("images/icon.gif"));
      final Image icon_image = icon.getImage();

      //TransferableContig tcontig = new TransferableContig(selected_features.elementAt(0)); 
      StringSelection name = new StringSelection(selected_features.elementAt(0).getGeneName());

      e.startDrag(DragSource.DefaultCopyDrop,     // cursor
                  icon_image, new Point(-1, -1),
                 (Transferable)name,              // transferable data
                                       this);     // drag source listener
    }
  }

  public void dragDropEnd(DragSourceDropEvent e) {}
  public void dragEnter(DragSourceDragEvent e) 
  {
  }

  public void dragExit(DragSourceEvent e) 
  {
    highlight_drop_base = -1;
  }
  public void dragOver(DragSourceDragEvent e) {}
  public void dropActionChanged(DragSourceDragEvent e) {}

  protected ZoomScrollBar getScaleChanger()
  {
    return scale_changer;
  }

}
