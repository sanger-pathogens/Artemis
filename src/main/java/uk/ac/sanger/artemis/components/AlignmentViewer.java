/* AlignmentViewer.java
 *
 * created: Mon Jul 12 1999
 *
 * This file is part of Artemis
 *
 * Copyright (C) 1999,2000,2001  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/AlignmentViewer.java,v 1.43 2008-11-28 17:51:09 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.sequence.Strand;
import uk.ac.sanger.artemis.sequence.Bases;
import uk.ac.sanger.artemis.sequence.SequenceChangeListener;
import uk.ac.sanger.artemis.sequence.SequenceChangeEvent;
import uk.ac.sanger.artemis.util.StringVector;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.io.RangeVector;

import java.awt.*;
import java.awt.event.*;
import java.util.Vector;
import java.util.Comparator;
import java.util.Arrays;
import javax.swing.*;

import org.apache.batik.svggen.SVGGraphics2D;

import java.io.FileWriter;
import java.io.IOException;

/**
 *  This component shows an alignment of two sequences using the data from a
 *  ComparisonData object.
 *
 *  @author Kim Rutherford
 *  @version $Id: AlignmentViewer.java,v 1.43 2008-11-28 17:51:09 tjc Exp $
 **/

public class AlignmentViewer extends CanvasPanel
    implements SequenceChangeListener 
{
  private static final long serialVersionUID = 1L;

  private Image offscreen;

  /** Comparison data that will be displayed in this component. */
  final private ComparisonData comparison_data;

  /** 
   *  All the AlignMatch objects from comparison_data (possibly in a
   *  different order.
   **/
  private AlignMatch[] all_matches = null;

  /**
   *  This is the last DisplayAdjustmentEvent reference that was passed to
   *  setSubjectSeqeuencePosition().
   **/
  private DisplayAdjustmentEvent last_subject_event;

  /**
   *  This is the last DisplayAdjustmentEvent reference that was passed to
   *  setQuerySeqeuencePosition().
   **/
  private DisplayAdjustmentEvent last_query_event;

  /** FeatureDisplay that is above this component.  (From the constructor). */
  private FeatureDisplay subject_feature_display;

  /** FeatureDisplay that is below this component.  (From the constructor). */
  private FeatureDisplay query_feature_display;

  /**
   *  Set by the constructor to be the original forward strand for the subject
   *  sequence.  This is use to determine whether to subject sequence has been
   *  reverse-complemented or not.
   **/
  private Strand orig_subject_forward_strand;

  /**
   *  Set by the constructor to be the original forward strand for the query
   *  sequence.  This is use to determine whether to query sequence has been
   *  reverse-complemented or not.
   **/
  private Strand orig_query_forward_strand;

  /** One of the two Entry objects that we are comparing. */
  final private EntryGroup subject_entry_group;

  /** One of the two Entry objects that we are comparing. */
  final private EntryGroup query_entry_group;

  /** Selected matches.  null means no matches are selected. */
  private AlignMatchVector selected_matches = null;

  /**
   *  The objects that are listening for AlignmentSelectionChangeEvents.
   **/
  private Vector<AlignmentSelectionChangeListener> selection_change_listeners = 
      new Vector<AlignmentSelectionChangeListener>();

  /**
   *  The number of shades of red and blue to use for percentage ID colouring.
   **/
  private static int NUMBER_OF_SHADES = 13;

  /** Reds used to display the percent identity of matches.  */
  private Color[] red_percent_id_colours;

  /** Blues used to display the percent identity of matches. */
  private Color[] blue_percent_id_colours;

  /** Scroll bar used to set the minimum length of the visible matches. */
  private JScrollBar scroll_bar = null;

  /** Matches with scores below this value will not be shown. */
  private int minimum_score = 0;

  /** Matches with scores above this value will not be shown. */
  private int maximum_score = 99999999;

  /**
   *  Matches with percent id values below this number will not be shown.
   **/
  private int minimum_percent_id = 0;

  /**
   *  Matches with percent id values above this number will not be shown.
   **/
  private int maximum_percent_id = 100;

  /**
   *  True if we should offer to flip the query sequence when the user
   *  double clicks on a flipped match.
   **/
  private boolean offer_to_flip_flag = false;

  /**
   *  If true ignore self matches (ie query start == subject start && query
   *  end == subject end)
   **/
  private boolean ignore_self_match_flag = false;

  /** Vector of those objects that are listening for AlignmentEvents */
  private Vector<AlignmentListener> alignment_event_listeners =  new Vector<AlignmentListener> ();

  /**
   *  If true then the FeatureDisplays above and below this AlignmentViewer
   *  should scroll together.
   **/
  private boolean displays_are_locked = true;

  /**
   *  Setting this to true will temporarily disable selectFromQueryRange() and
   *  selectFromSubjectRange() until enableSelection().  This is need to allow
   *  the selections of the top and bottom FeatureDisplays to be set without
   *  changing which AlignMatches are selected.
   **/
  private boolean disable_selection_from_ranges = false;

  /** user defined colours */
  private boolean reverseMatchColour = false;
  /** colour for reverse matches */
  private Color revMatchColour = Color.blue;
  /** colour for matches         */
  private Color matchColour    = Color.red;

  /**
   *  Create a new AlignmentViewer for the given entries.
   *  @param subject_feature_display The FeatureDisplay that is above this
   *    component.
   *  @param query_feature_display The FeatureDisplay that is below this
   *    component.
   *  @param comparison_data Provides the AlignMatch objects that will be
   *    displayed.
   **/
  public AlignmentViewer(final FeatureDisplay subject_feature_display,
                         final FeatureDisplay query_feature_display,
                         final ComparisonData comparison_data) 
  {
    this.subject_feature_display = subject_feature_display;
    this.query_feature_display   = query_feature_display;
    this.comparison_data         = comparison_data;
    this.all_matches             = getComparisonData().getMatches();

    subject_entry_group          = getSubjectDisplay().getEntryGroup();
    query_entry_group            = getQueryDisplay().getEntryGroup();

    final Bases subject_bases = getSubjectForwardStrand().getBases();
    final Bases query_bases   = getQueryForwardStrand().getBases();

    final Selection subject_selection = getSubjectDisplay().getSelection();
    final Selection query_selection = getQueryDisplay().getSelection();

    final SelectionChangeListener subject_listener =
      new SelectionChangeListener() 
    {
      JFrame frame = null;
      public void selectionChanged(SelectionChangeEvent event) 
      {
        if(frame == null)
          frame = subject_feature_display.getParentFrame();
        if(!frame.isVisible())
          return;

        final RangeVector ranges = subject_selection.getSelectionRanges();
        selectFromSubjectRanges(ranges);
      }
    };

    final SelectionChangeListener query_listener =
      new SelectionChangeListener() 
    {
      JFrame frame = null;
      public void selectionChanged (SelectionChangeEvent event) 
      {
        if(frame == null)
          frame = query_feature_display.getParentFrame();
        if(!frame.isVisible())
          return;

        final RangeVector ranges = query_selection.getSelectionRanges ();
        selectFromQueryRanges(ranges);
      }
    };

    makeColours();

    subject_selection.addSelectionChangeListener(subject_listener);
    query_selection.addSelectionChangeListener(query_listener);

    subject_bases.addSequenceChangeListener(this, 0);
    query_bases.addSequenceChangeListener(this, 0);

    orig_subject_forward_strand = getSubjectForwardStrand();
    orig_query_forward_strand   = getQueryForwardStrand();

    addMouseListener(new MouseAdapter() 
    {
      public void mousePressed(final MouseEvent event) 
      {
        // on windows we have to check isPopupTrigger in mouseReleased(),
        // but do it in mousePressed() on UNIX
        if(isMenuTrigger(event)) 
          popupMenu(event);
        else 
          handleCanvasMousePress(event);
      }
    });

    addMouseMotionListener(new MouseMotionAdapter() 
    {
      public void mouseDragged(final MouseEvent event) 
      {
        if(isMenuTrigger(event))
          return;

        if(!modifiersForLockToggle(event))
        {
          if(!event.isShiftDown()) 
          {
            selected_matches = null;
            toggleSelection (event.getPoint());
          }
          repaint();
        }
      }
    });

    scroll_bar = new JScrollBar(Scrollbar.VERTICAL);
    scroll_bar.setValues(1, 1, 1, 1000);
    scroll_bar.setBlockIncrement(10);

    scroll_bar.addAdjustmentListener(new AdjustmentListener() 
    {
      public void adjustmentValueChanged(AdjustmentEvent e) 
      {
        repaint();
      }
    });

    maximum_score = getComparisonData().getMaximumScore();

    add(scroll_bar, "East");
    setBackground(Color.white);
  }

  /**
   *  Returns true if and only if the given MouseEvent should toggle the lock
   *  displays toggle.
   **/
  private boolean modifiersForLockToggle(final MouseEvent event) 
  {
    return(event.getModifiers() & InputEvent.BUTTON2_MASK) != 0 ||
      event.isAltDown();
  }

  /**
   *  Select those matches that overlap the given range on the subject
   *  sequence.
   **/
  public void selectFromSubjectRanges(final RangeVector select_ranges) 
  {
    if(disable_selection_from_ranges) 
      return;

    selected_matches = null;
    final int all_matches_length = all_matches.length;
    final int select_ranges_size = select_ranges.size();

    final Strand current_subject_fwd_strand =
                              getSubjectForwardStrand();

    final int subject_length = current_subject_fwd_strand.getSequenceLength();

    for(int match_index = 0; match_index < all_matches_length; ++match_index)
    {
      final AlignMatch this_match = all_matches[match_index];

      if(!isVisible(this_match))
        continue;

      int subject_sequence_start = getRealSubjectSequenceStart(this_match,
          subject_length,
          (getOrigSubjectForwardStrand() != current_subject_fwd_strand));
      int subject_sequence_end = getRealSubjectSequenceEnd(this_match,
          subject_length,
          (getOrigSubjectForwardStrand() != current_subject_fwd_strand));

      if(subject_sequence_end < subject_sequence_start)
      {
        final int tmp = subject_sequence_start;
        subject_sequence_start = subject_sequence_end;
        subject_sequence_end = tmp;
      }

      for(int range_index = 0; range_index < select_ranges_size; ++range_index)
      {
        final Range select_range = (Range) select_ranges.elementAt(range_index);
        final int select_range_start = select_range.getStart();
        final int select_range_end = select_range.getEnd();

        if(select_range_start < subject_sequence_start
            && select_range_end < subject_sequence_start)
          continue;

        if(select_range_start > subject_sequence_end &&
           select_range_end   > subject_sequence_end) 
          continue;

        if(selected_matches == null) 
          selected_matches = new AlignMatchVector();

        //if(!selected_matches.contains(this_match)) 
        selected_matches.add(this_match);
        break;
      }
    }

    if(selected_matches != null)
      selectionChanged();
    else
      repaint();
  }

  /**
   *  Select those matches that overlap the given range on the query sequence.
   **/
  public void selectFromQueryRanges(final RangeVector select_ranges) 
  {
    if(disable_selection_from_ranges) 
      return;

    selected_matches = null;
    final int select_ranges_size = select_ranges.size();
    final int all_matches_length = all_matches.length;
    final Strand current_query_forward_strand = getQueryForwardStrand();
    final int query_length =
          current_query_forward_strand.getSequenceLength();

    for(int match_index = 0; match_index < all_matches_length; ++match_index)
    {
      final AlignMatch this_match = all_matches[match_index];

      if(!isVisible(this_match))
        continue;

      int query_sequence_start = getRealQuerySequenceStart(this_match,
          query_length,
          (getOrigQueryForwardStrand() != current_query_forward_strand));
      int query_sequence_end = getRealQuerySequenceEnd(this_match,
          query_length,
          (getOrigQueryForwardStrand() != current_query_forward_strand));

      if(query_sequence_end < query_sequence_start)
      {
        final int tmp = query_sequence_start;
        query_sequence_start = query_sequence_end;
        query_sequence_end = tmp;
      }

      for(int range_index = 0; range_index < select_ranges_size; ++range_index)
      {
        final Range select_range = (Range) select_ranges.elementAt(range_index);
        final int select_range_start = select_range.getStart();
        final int select_range_end = select_range.getEnd();
        if(select_range_start < query_sequence_start
            && select_range_end < query_sequence_start) 
          continue;

        if(select_range_start > query_sequence_end &&
           select_range_end > query_sequence_end) 
          continue;

        if(selected_matches == null) 
          selected_matches = new AlignMatchVector();

        //if(!selected_matches.contains(this_match)) 
        selected_matches.add(this_match);
        break;
      }
    }
    
    if(selected_matches != null)
      selectionChanged();
    else
      repaint();
  }

  /**
   *  Select the given match and move it to the top of the display.
   **/
  public void setSelection(final AlignMatch match) 
  {
    selected_matches = new AlignMatchVector();
    selected_matches.add(match);
    selectionChanged();
  }

  /**
   *  This method tells this AlignmentViewer component where the subject
   *  sequence is now.
   **/
  public void setSubjectSequencePosition(final DisplayAdjustmentEvent event) 
  {
    last_subject_event = event;
    repaint();
  }

  /**
   *  This method tells this AlignmentViewer component where the query
   *  sequence is now.
   **/
  public void setQuerySequencePosition(final DisplayAdjustmentEvent event) 
  {
    last_query_event = event;
    repaint();
  }

  /**
   *  Implementation of the SequenceChangeListener interface.  The display is
   *  redrawn if there is an event.
   **/
  public void sequenceChanged(final SequenceChangeEvent event) 
  {
    repaint();
  }

  /**
   *  Return true if and only if the given MouseEvent (a mouse press) should
   *  pop up a JPopupMenu.
   **/
  private boolean isMenuTrigger(final MouseEvent event) 
  {
    if( event.isPopupTrigger() ||
       (event.getModifiers() & InputEvent.BUTTON3_MASK) != 0) 
      return true;
    else 
      return false;
  }

  /**
   *  Popup a menu.
   **/
  private void popupMenu(final MouseEvent event) 
  {
    final JPopupMenu popup = new JPopupMenu();


    final JMenuItem save_matches = new JMenuItem("Save Comparison File...");
    save_matches.addActionListener(new ActionListener()
    {
      public void actionPerformed (ActionEvent actionEvent)
      {
        StickyFileChooser fc = new StickyFileChooser();

        int returnVal = fc.showSaveDialog(null);  
        if(returnVal != JFileChooser.APPROVE_OPTION) 
          return;
        else if(fc.getSelectedFile().exists())
        {
          Object[] possibleValues = { "YES", "NO" };
          int select = JOptionPane.showOptionDialog(null, 
                                 fc.getSelectedFile().getName()+"\n"+
                                 "exists. Overwrite?",
                                 "File Exists",
                                 JOptionPane.DEFAULT_OPTION,
                                 JOptionPane.QUESTION_MESSAGE,null,
                                 possibleValues, possibleValues[0]);
          if(select == 1)
            return;  
        }

        try
        {
/*
          if(!fc.getSelectedFile().canWrite())
          {
            JOptionPane.showMessageDialog(null,
                        "Cannot write to "+
                        fc.getSelectedFile().getCanonicalPath(),
                        "Warning",
                        JOptionPane.WARNING_MESSAGE);
            return;
          }
*/
          final FileWriter out_writer = new FileWriter(fc.getSelectedFile());
          final String query = getQueryEntryGroup().getDefaultEntry().getName();
          final String subject = getSubjectEntryGroup().getDefaultEntry().getName();

          for(int i = 0; i < all_matches.length; ++i)
            MSPcrunchComparisonData.writeMatchFromAlignMatch(all_matches[i],
                                               query, subject,
                                               out_writer);
          out_writer.close();
        }
        catch(IOException ioe)
        {
          JOptionPane.showMessageDialog(null,
                        "Error writing out comparison file.", 
                        "Warning",
                        JOptionPane.WARNING_MESSAGE);
          ioe.printStackTrace();
        }
      }
    });
    popup.add(save_matches);
    popup.add(new JSeparator());

    final JMenuItem alignmatch_list_item =
      new JMenuItem("View Selected Matches");

    popup.add(alignmatch_list_item);

    alignmatch_list_item.addActionListener(new ActionListener() 
    {
      public void actionPerformed (ActionEvent actionEvent) 
      {
        if(selected_matches == null) 
          new MessageFrame("No matches selected").setVisible (true);
        else 
        {
          final AlignMatchVector matches =
            (AlignMatchVector)selected_matches.clone();

          final AlignMatchViewer viewer =
            new AlignMatchViewer(AlignmentViewer.this, matches);

          viewer.setVisible(true);
        }
      }
    });

    final JMenuItem flip_subject_item =
      new JMenuItem("Flip Subject Sequence");

    popup.add (flip_subject_item);

    flip_subject_item.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent actionEvent) 
      {
        if(getSubjectDisplay().isRevCompDisplay()) 
          getSubjectDisplay().setRevCompDisplay(false);
        else 
          getSubjectDisplay().setRevCompDisplay(true);
      }
    });

    final JMenuItem flip_query_item =
      new JMenuItem ("Flip Query Sequence");

    popup.add (flip_query_item);

    flip_query_item.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent actionEvent) 
      {
        if(getQueryDisplay().isRevCompDisplay())
          getQueryDisplay().setRevCompDisplay(false);
        else
          getQueryDisplay().setRevCompDisplay(true);
      }
    });

    final JMenuItem cutoffs_item = new JMenuItem("Set Score Cutoffs ...");
    popup.add(cutoffs_item);

    cutoffs_item.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent actionEvent) 
      {
        final ScoreChangeListener minimum_listener =
          new ScoreChangeListener() 
        {
          public void scoreChanged(final ScoreChangeEvent event)
          {
            minimum_score = event.getValue();
            repaint();
          }
        };

        final ScoreChangeListener maximum_listener =
          new ScoreChangeListener()
        {
          public void scoreChanged(final ScoreChangeEvent event) 
          {
            maximum_score = event.getValue();
            repaint();
          }
        };

        final ScoreChanger score_changer =
          new ScoreChanger("Score Cutoffs",
                           minimum_listener, maximum_listener,
                           getComparisonData().getMinimumScore(),
                           getComparisonData().getMaximumScore());

        score_changer.setVisible (true);
      }
    });

    final JMenuItem percent_id_cutoffs_item =
      new JMenuItem("Set Percent ID Cutoffs ...");
    popup.add(percent_id_cutoffs_item);

    percent_id_cutoffs_item.addActionListener(new ActionListener () 
    {
      public void actionPerformed(ActionEvent actionEvent) 
      {
        final ScoreChangeListener minimum_listener =
          new ScoreChangeListener()
        {
          public void scoreChanged(final ScoreChangeEvent event) 
          {
            minimum_percent_id = event.getValue();
            repaint();
          }
        };

        final ScoreChangeListener maximum_listener =
          new ScoreChangeListener() 
        {
          public void scoreChanged(final ScoreChangeEvent event) 
          {
            maximum_percent_id = event.getValue();
            repaint();
          }
        };

        final ScoreChanger score_changer =
          new ScoreChanger("Percent Identity Cutoffs",
                           minimum_listener, maximum_listener,
                           0, 100);

        score_changer.setVisible(true);
      }
    });


    final JCheckBoxMenuItem lock_item = new JCheckBoxMenuItem("Lock Sequences");
    lock_item.setSelected(displaysAreLocked());
    popup.add(lock_item);

    lock_item.addItemListener(new ItemListener()
    {
      public void itemStateChanged(ItemEvent e)
      {
        if(lock_item.isSelected())
          lockDisplays();
        else
          unlockDisplays();
      }
    });

    popup.addSeparator();

    final JCheckBoxMenuItem sameColour = 
       new JCheckBoxMenuItem("Colour reverse & forward matches the same",reverseMatchColour);
    sameColour.addItemListener(new ItemListener()
    {
      public void itemStateChanged(ItemEvent e)
      {
        reverseMatchColour = sameColour.getState();
        repaint();
      }
    });
    popup.add(sameColour);

    JMenuItem colourMatches = new JMenuItem("Colour matches...");
    colourMatches.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)
      {
        ColorChooserShades shades = createColours("Colour Matches", 
                                                  red_percent_id_colours[NUMBER_OF_SHADES-1]);
        if(shades != null)
        {
          red_percent_id_colours = shades.getDefinedColour();
          repaint();
        }
      }
    });
    popup.add(colourMatches);

    if(!reverseMatchColour)
    {
      JMenuItem colourRevMatches = new JMenuItem("Colour reverse matches...");
      colourRevMatches.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent event)
        {
          ColorChooserShades shades = createColours("Colour Reverse Matches",
                                                   blue_percent_id_colours[NUMBER_OF_SHADES-1]);
          if(shades != null)
          {
            blue_percent_id_colours = shades.getDefinedColour();
            repaint();
          }
        }
      });
      popup.add(colourRevMatches);
    }

    popup.addSeparator();

    final JCheckBoxMenuItem offer_to_flip_item =
      new JCheckBoxMenuItem("Offer To RevComp", offer_to_flip_flag);

    offer_to_flip_item.addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent e) 
      {
        offer_to_flip_flag = !offer_to_flip_flag;
      }
    });

    popup.add(offer_to_flip_item);

    final JCheckBoxMenuItem ignore_self_match_item =
      new JCheckBoxMenuItem("Ignore Self Matches");

    ignore_self_match_item.addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent event)
      {
        ignore_self_match_flag = ignore_self_match_item.getState();
        repaint();
      }
    });

    ignore_self_match_item.setState(ignore_self_match_flag);

    popup.add(ignore_self_match_item);

    add(popup);
    popup.show(this, event.getX(), event.getY());
  }


  private ColorChooserShades createColours(String title, Color initialColour) 
  {
    //Make sure we have nice window decorations.
    JFrame.setDefaultLookAndFeelDecorated(true);

    //Create and set up the window.
    JFrame frame = new JFrame("ColorChooserDemo");
    frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

    //Create and set up the content pane.
    ColorChooserShades newContentPane = new ColorChooserShades(title,initialColour);

    Object[] possibleValues = { "OK", "CANCEL" };
    int select = JOptionPane.showOptionDialog(null, newContentPane,
                                 "Colour Selection",
                                 JOptionPane.DEFAULT_OPTION, 
                                 JOptionPane.QUESTION_MESSAGE,null,
                                 possibleValues, possibleValues[0]);
    if(select == 0)
      return newContentPane;

    return null;
  }


  /**
   *  Handle a mouse press event on the drawing canvas - select on click,
   *  select and broadcast it on double click.
   **/
  private void handleCanvasMousePress(final MouseEvent event) 
  {
    if(event.getID() != MouseEvent.MOUSE_PRESSED) 
      return;

    if(event.getClickCount() == 2) 
      handleCanvasDoubleClick(event);
    else 
      handleCanvasSingleClick(event);

    repaint();
  }

  /**
   *  Handle a double click on the canvas.
   **/
  private void handleCanvasDoubleClick(final MouseEvent event) 
  {
    // there should be only one match in the array
    if(selected_matches != null) 
      alignAt(selected_matches.elementAt(0));
  }

  /**
   *  Send an AlignmentEvent to all the AlignmentListeners.
   *  @param align_match The AlignMatch that we have just centred on.
   **/
  public void alignAt(final AlignMatch align_match) 
  {
    final Vector<AlignmentListener> targets;
    // copied from a book - synchronizing the whole method might cause a
    // deadlock
    synchronized(this) 
    {
      targets = new Vector<AlignmentListener>(alignment_event_listeners);
    }

    for(int i = 0; i < targets.size(); ++i) 
    {
      final AlignmentListener listener = targets.elementAt(i);
      listener.alignMatchChosen(new AlignmentEvent(align_match));
    }
  }

  /**
   *  Handle a single click on the canvas.
   **/
  private void handleCanvasSingleClick(final MouseEvent event) 
  {
    if(modifiersForLockToggle(event)) 
      toggleDisplayLock();
    else 
    {
      if(!event.isShiftDown()) 
        selected_matches = null;
      
      toggleSelection(event.getPoint());
    }
  }

  /**
   *  Add or remove the match at the given mouse position to the selection.
   **/
  private void toggleSelection(final Point point) 
  {
    final AlignMatch clicked_align_match =
      getAlignMatchFromPosition(point);

    if(clicked_align_match != null) 
    {
      if(selected_matches == null) 
      {
        selected_matches = new AlignMatchVector ();
        selected_matches.add (clicked_align_match);
      } 
      else 
      {
        if(selected_matches.contains(clicked_align_match)) 
        {
          selected_matches.remove(clicked_align_match);
          if(selected_matches.size() == 0)
            selected_matches = null;
        }
        else
          selected_matches.add(clicked_align_match);
      }
    }

    selectionChanged();
  }

  /**
   *  Return the AlignMatch at the given Point on screen or null if there is
   *  no match at that point.  The alignment_data_array is searched in reverse
   *  order.
   **/
  private AlignMatch getAlignMatchFromPosition(final Point click_point) 
  {
    final int canvas_height = getSize().height;
    final int canvas_width  = getSize().width;

    final int subject_length = getSubjectForwardStrand().getSequenceLength();
    final int query_length   = getQueryForwardStrand().getSequenceLength();

    final boolean subject_flipped = subjectIsRevComp();
    final boolean query_flipped   = queryIsRevComp();

    final int all_matches_length = all_matches.length;
    final float base_width       = last_subject_event.getBaseWidth();
    final float query_base_width = last_query_event.getBaseWidth();
    
    final int subject_start = last_subject_event.getStart();
    final int query_start   = last_query_event.getStart();
    final boolean subject_is_rev_comp = subjectIsRevComp();
    final boolean query_is_rev_comp   = queryIsRevComp();
    boolean is_rev_match;
    int[] match_x_positions;

    for(int i = all_matches_length - 1; i >= 0 ; --i) 
    {
      final AlignMatch this_match = all_matches [i];

      is_rev_match = this_match.isRevMatch();
      match_x_positions =
        getMatchCoords(canvas_width, this_match, subject_length, query_length,
                       subject_flipped, query_flipped, base_width, query_base_width, subject_start, 
                       query_start, subject_is_rev_comp, query_is_rev_comp, is_rev_match);

      if(match_x_positions == null) 
        continue;

      if(!isVisible(this_match)) 
        continue;

      final int subject_start_x = match_x_positions[0];
      final int subject_end_x   = match_x_positions[1];
      final int query_start_x   = match_x_positions[2];
      final int query_end_x     = match_x_positions[3];

      // this is the x coordinate of the point where the line y = click_point
      // hits the left edge of the match box
      final double match_left_x =
        subject_start_x +
        (1.0 * (query_start_x - subject_start_x)) *
        (1.0 * click_point.y / canvas_height);

      // this is the x coordinate of the point where the line y = click_point
      // hits the right edge of the match box
      final double match_right_x =
        subject_end_x +
        (1.0 * (query_end_x - subject_end_x)) *
        (1.0 * click_point.y / canvas_height);

      if(click_point.x >= match_left_x - 1 &&
         click_point.x <= match_right_x + 1 ||
         click_point.x <= match_left_x + 1 &&
         click_point.x >= match_right_x - 1) 
        return this_match;
    }

    return null;
  }

  /**
   *  This method is called by setSelection() and others whenever the list of
   *  selected/highlighted hits changes. Calls alignmentSelectionChanged()
   *  on all interested AlignmentSelectionChangeListener objects, moves the
   *  selected matches to the top of the display and then calls
   *  repaint
   **/
  private void selectionChanged()
  {
    for(int i = 0 ; i < selection_change_listeners.size() ; ++i) 
    {
      final AlignMatchVector matches =
        (AlignMatchVector) selected_matches.clone();

      final AlignmentSelectionChangeEvent ev =
        new AlignmentSelectionChangeEvent(this, matches);

      final AlignmentSelectionChangeListener listener = selection_change_listeners.elementAt(i);
      listener.alignmentSelectionChanged(ev);
    }

    if(selected_matches == null)
      return;

    final int selected_matches_size = selected_matches.size();
    if(selected_matches != null && selected_matches_size > 0) 
    {
      // a count of the number of selected matches seen so far
      int seen_and_selected_count = 0;
      final int all_matches_length = all_matches.length;

      for(int i = 0 ; i < all_matches_length ; ++i) 
      {
        final AlignMatch this_match = all_matches[i];

        if(selected_matches.contains(this_match)) 
          ++seen_and_selected_count;
        else
        { 
          if(seen_and_selected_count > 0) 
          {
            // move the matches down to fill the gap
            all_matches[i-seen_and_selected_count] = all_matches[i];
          }
        }
      }

      // put the selected_matches at the end of all_matches
      for(int i = 0; i < selected_matches_size; ++i) 
      {
        all_matches[all_matches_length - selected_matches_size + i] =
          selected_matches.elementAt(i);
      }
    }

    repaint();
  }

  /**
   *  Add the AlignmentSelectionChangeListener to the list of objects
   *  listening for AlignmentSelectionChangeEvents.
   **/
  public void
    addAlignmentSelectionChangeListener(final AlignmentSelectionChangeListener l) 
  {
    selection_change_listeners.addElement(l);
  }

  /**
   *  Remove the AlignmentSelectionChangeListener from the list of objects
   *  listening for AlignmentSelectionChangeEvents.
   **/
  public void
    removeAlignmentSelectionChangeListener(final AlignmentSelectionChangeListener l) 
  {
    selection_change_listeners.removeElement(l);
  }

  /**
   *  Adds the specified AlignmentEvent listener to receive events from this
   *  object.
   *  @param l the listener.
   **/
  public void addAlignmentListener(AlignmentListener l) 
  {
    alignment_event_listeners.addElement(l);
  }

  /**
   *  Removes the specified AlignmentEvent listener so that it no longer
   *  receives events from this object.
   *  @param l the listener.
   **/
  public void removeAlignmentListener(AlignmentListener l) 
  {
    alignment_event_listeners.removeElement(l);
  }

  /**
   *  Returns true if and only if we should offer to flip the query
   *  sequence when the user double clicks on a flipped match.
   **/
  public boolean offerToFlip() 
  {
    return offer_to_flip_flag;
  }


  /**
  *  Set the offscreen buffer to null as part of invalidation.
  **/
  public void invalidate()
  {
    super.invalidate();
    offscreen = null;
  }


  /**
   *  The main paint function for the canvas.  An off screen image used for
   *  double buffering when drawing the canvas.
   *  @param g The Graphics object of the canvas.
   **/
  protected void paintComponent(final Graphics g) 
  {
    super.paintComponent(g);
 
    if(last_subject_event != null && last_query_event != null) 
    {
      final int canvas_height = getSize().height;
      int canvas_width  = getSize().width;
      if(scroll_bar != null)
        canvas_width -= scroll_bar.getPreferredSize().width;

      // need to do this off-screen otherwise it
      // does not draw properly on windows
      if(offscreen == null)
        offscreen = createImage(canvas_width, canvas_height);

      Graphics og = offscreen.getGraphics();
      og.setClip(0,0,canvas_width,canvas_height);
      og.setColor(Color.white);
      og.fillRect(0, 0, canvas_width, canvas_height);
      drawAlignments(og);
      drawLabels(og);
      g.drawImage(offscreen, 0, 0, null);
      og.dispose();
    }
  }

 
  /**
   *  The paint function for printing.
   *  
   *  @param g The Graphics object of the canvas.
   *  @param drawLabel draw the labels
   **/
  protected void paintComponentForPrint(final Graphics g, final boolean drawLabel)
  {
    super.paintComponent(g);

    if(last_subject_event != null && last_query_event != null)
    {
      final int canvas_height = getSize().height;
      int canvas_width  = getSize().width;
      if(scroll_bar != null)
        canvas_width -= scroll_bar.getPreferredSize().width;

      if(offscreen == null)
        offscreen = createImage(canvas_width, canvas_height);

      final Graphics og;
      if(!(g instanceof SVGGraphics2D))
      {
        og = offscreen.getGraphics();
        og.setClip(0,0,canvas_width,canvas_height);
      }
      else
        og = g;
      
      og.setColor(Color.white);
      og.fillRect(0, 0, canvas_width, canvas_height);
      drawAlignments(og);
      drawLabels(og);
       
      if(!(g instanceof SVGGraphics2D))
      {
        g.drawImage(offscreen, 0, 0, null);
        og.dispose();
      }
    }
  }


  /**
   *  Draw the labels into the given Graphics object.  There is a label at the
   *  top left showing info about the current AlignMatch object,
   *
   *  XXX
   *  a label at
   *  the bottom left showing whether or not the query sequence is reverse
   *  complemented
   *  XXX
   *
   *  and a label beside the Scrollbar showing the current score
   *  cutoff.
   **/
  private void drawLabels(final Graphics g) 
  {
    final FontMetrics fm   = g.getFontMetrics();
    int canvas_width  = getSize().width;

    if(scroll_bar != null)
      canvas_width -= scroll_bar.getPreferredSize().width;

    final int canvas_height = getSize().height;

    final String cutoff_label =
      Integer.toString(scroll_bar.getValue());

    final int cutoff_label_width = fm.stringWidth(cutoff_label);

    int cutoff_label_position =
      (int)((scroll_bar.getValue() -
             scroll_bar.getMinimum()) / (1.0 *
            (scroll_bar.getMaximum() -
             scroll_bar.getMinimum())) * canvas_height);

    if(cutoff_label_position < getFontAscent()) 
      cutoff_label_position = getFontAscent();


    final int[] cutoff_x_points = 
    {
      canvas_width - cutoff_label_width,
      canvas_width - 2,
      canvas_width - 2,
      canvas_width - cutoff_label_width,
    };

    final int[] cutoff_y_points = 
    {
      cutoff_label_position + 1,
      cutoff_label_position + 1,
      cutoff_label_position - getFontAscent(),
      cutoff_label_position - getFontAscent(),
    };

    g.setColor(Color.white);
    g.fillPolygon(cutoff_x_points, cutoff_y_points, 4);

    g.setColor(Color.black);
    g.drawString(cutoff_label, canvas_width - cutoff_label_width,
                 cutoff_label_position);

    final int font_height = getFontAscent() + getFontDescent();

    if(selected_matches != null) 
    {
      final String match_string_1;

      if(selected_matches.size() > 1) 
        match_string_1 = selected_matches.size() + " matches selected";
      else 
      {
        final AlignMatch selected_align_match = selected_matches.elementAt(0);

        match_string_1 =
          selected_align_match.getQuerySequenceStart() + ".." +
          selected_align_match.getQuerySequenceEnd() + " -> " +
          selected_align_match.getSubjectSequenceStart() + ".." +
          selected_align_match.getSubjectSequenceEnd();
      }

      final int match_string_1_width = fm.stringWidth(match_string_1);

      final int[] match_1_x_points = 
      {
        0, 0, match_string_1_width, match_string_1_width
      };

      final int[] match_1_y_points = 
      {
        0, font_height, font_height, 0
      };

      g.setColor(Color.white);
      g.fillPolygon(match_1_x_points, match_1_y_points, 4);

      g.setColor(Color.black);
      g.drawString(match_string_1, 0, getFontAscent ());

      if(selected_matches.size() == 1) 
      {
        final AlignMatch selected_align_match = selected_matches.elementAt(0);

        final String match_string_2 = "score: " +
          selected_align_match.getScore() + "  percent id: " +
          selected_align_match.getPercentID() + "%";
        
        final int match_string_2_width = fm.stringWidth(match_string_2);
        
        final int[] match_2_x_points = 
        {
          0, 0, match_string_2_width, match_string_2_width
        };
        
        final int[] match_2_y_points = 
        {
          font_height, font_height * 2, font_height * 2, font_height
        };
        
        g.setColor(Color.white);
        g.fillPolygon(match_2_x_points, match_2_y_points, 4);
        
        g.setColor(Color.black);
        g.drawString(match_string_2, 0, getFontAscent() + font_height);
      }
    }

    final StringVector status_strings = new StringVector();

    if(displaysAreLocked()) 
      status_strings.add("LOCKED");

    if(getSubjectDisplay().isRevCompDisplay()) 
      status_strings.add("Subject: Flipped");

    if(getQueryDisplay().isRevCompDisplay()) 
      status_strings.add("Query: Flipped");

    if(getOrigSubjectForwardStrand() != getSubjectForwardStrand()) 
      status_strings.add("Subject: Reverse Complemented");

    if(getOrigQueryForwardStrand() != getQueryForwardStrand()) 
      status_strings.add("Query: Reverse Complemented");

    g.setColor(Color.white);

    for(int i = 0 ; i < status_strings.size() ; ++i) 
    {
      final String status_string = (String)status_strings.elementAt(i);

      final int status_string_width = fm.stringWidth(status_string);

      final int[] x_points = 
      {
        0, 0, status_string_width, status_string_width
      };

      final int string_offset = font_height * (status_strings.size () - i - 1);

      final int[] y_points = 
      {
        canvas_height - string_offset,
        canvas_height - font_height - string_offset,
        canvas_height - font_height - string_offset,
        canvas_height - string_offset
      };

      g.fillPolygon(x_points, y_points, 4);
    }

    g.setColor(Color.black);
    for(int i = 0; i < status_strings.size(); ++i) 
    {
      final String status_string = (String)status_strings.elementAt(i);
      final int string_offset = font_height * (status_strings.size() - i - 1);

      g.drawString(status_string, 0,
                   canvas_height - string_offset - getFontDescent());
    }
  }

  /**
   *  Draw the alignments into the given Graphics object.
   **/
  private void drawAlignments(final Graphics g) 
  {
    final int canvas_height = getSize().height;
    final int canvas_width  = getSize().width;
    final int OFFSCREEN     = 3000;

    final int subject_length = getSubjectForwardStrand().getSequenceLength();
    final int query_length   = getQueryForwardStrand().getSequenceLength();

    final boolean subject_flipped = subjectIsRevComp();
    final boolean query_flipped   = queryIsRevComp();

    final float base_width       = last_subject_event.getBaseWidth();
    final float query_base_width = last_query_event.getBaseWidth();
    
    final int subject_start = last_subject_event.getStart();
    final int query_start   = last_query_event.getStart();
    final boolean subject_is_rev_comp = subjectIsRevComp();
    final boolean query_is_rev_comp   = queryIsRevComp();
    boolean is_rev_match;
    int[] match_x_positions;
    AlignMatch this_match;
    final int all_matches_length = all_matches.length;    

    for(int i = 0 ; i < all_matches_length ; ++i) 
    {
      this_match = all_matches[i];

      is_rev_match = this_match.isRevMatch();
      match_x_positions =
        getMatchCoords(canvas_width, this_match, subject_length, query_length,
                       subject_flipped, query_flipped, base_width, query_base_width, subject_start, 
                       query_start, subject_is_rev_comp, query_is_rev_comp, is_rev_match);

      if(match_x_positions == null) 
        continue;

      if(!isVisible(this_match)) 
        continue;

      final int subject_start_x = match_x_positions[0];
      final int subject_end_x   = match_x_positions[1];
      final int query_start_x   = match_x_positions[2];
      final int query_end_x     = match_x_positions[3];

      final int[] x_coords = new int[4];
      final int[] y_coords = new int[4];

      x_coords[0] = subject_start_x;
      y_coords[0] = 0;
      x_coords[1] = query_start_x;
      y_coords[1] = canvas_height;
      x_coords[2] = query_end_x;
      y_coords[2] = canvas_height;
      x_coords[3] = subject_end_x;
      y_coords[3] = 0;

      final boolean highlight_this_match;

      if(selected_matches != null &&
         selected_matches.contains(this_match)) 
        highlight_this_match = true;
      else 
        highlight_this_match = false;

      final int percent_id = this_match.getPercentID();

      if(highlight_this_match)
        g.setColor (Color.yellow);
      else 
      {
        if(percent_id == -1) 
        {
          if(is_rev_match)
            g.setColor(revMatchColour);
          else 
            g.setColor(matchColour);
        } 
        else 
        {
          int colour_index = red_percent_id_colours.length - 1;

          if(maximum_percent_id > minimum_percent_id) 
          {
            colour_index =
              (int)(red_percent_id_colours.length * 0.999 *
                    (percent_id - minimum_percent_id) /
                    (maximum_percent_id - minimum_percent_id));
          }

          if(is_rev_match && !reverseMatchColour) 
            g.setColor(blue_percent_id_colours[colour_index]);
          else
            g.setColor(red_percent_id_colours[colour_index]);
        }
      }

      g.fillPolygon(x_coords, y_coords, x_coords.length);

      if(subject_end_x - subject_start_x < 5 &&
         subject_end_x - subject_start_x > -5 ||
         subject_start_x < -OFFSCREEN ||
         subject_end_x > OFFSCREEN ||
         query_start_x < -OFFSCREEN ||
         query_end_x > OFFSCREEN)
      {
        // match is (probably) narrow so draw the border to the same colour as
        // the polygon
      } 
      else 
      {
        // draw a black outline the make the match stand out
        g.setColor (Color.black);
      }

      g.drawLine(subject_start_x, 0, query_start_x, canvas_height);
      g.drawLine(subject_end_x, 0, query_end_x, canvas_height);
    }
  }

  /**
   *  Return true if and only if the given match is currently visible.
   **/
  private boolean isVisible(final AlignMatch match) 
  {
    if(ignore_self_match_flag && match.isSelfMatch()) 
      return false;

    final int score = match.getScore();
    if(score > -1) 
    {
      if(score < minimum_score || score > maximum_score) 
        return false;
    }

    final int percent_id = match.getPercentID();
    if(percent_id > -1) 
    {
      if(percent_id < minimum_percent_id || percent_id > maximum_percent_id) 
        return false;
    }

    final int match_length = match.getLength();
//    Math.abs(match.getSubjectSequenceStart() -
//             match.getSubjectSequenceEnd());

    if(match_length < scroll_bar.getValue()) 
      return false;

    return true;
  }

  /**
   *  Return the start position in the subject sequence of the given
   *  AlignMatch, taking into account the current orientation of the
   *  sequences. If the reverse_position argument is true reverse the
   *  complemented match coordinates will be returned.
   **/
  protected static int getRealSubjectSequenceStart(final AlignMatch match,
                                         final int sequence_length,
                                         final boolean reverse_position) 
  {
    if(reverse_position) 
      return sequence_length - match.getSubjectSequenceStart() + 1;
    else 
      return match.getSubjectSequenceStart();
  }

  /**
   *  Return the end position in the subject sequence of the given AlignMatch,
   *  taking into account the current orientation of the sequences.  If the
   *  reverse_position argument is true reverse the complemented match
   *  coordinates will be returned.
   **/
  protected static int getRealSubjectSequenceEnd(final AlignMatch match,
                                       final int sequence_length,
                                       final boolean reverse_position) 
  {
    if(reverse_position) 
      return sequence_length - match.getSubjectSequenceEnd() + 1;
    else 
      return match.getSubjectSequenceEnd();
  }

  /**
   *  Return the start position in the query sequence of the given AlignMatch,
   *  taking into account the current orientation of the sequences.    If the
   *  reverse_position argument is true reverse the complemented match
   *  coordinates will be returned.
   **/
  protected static int getRealQuerySequenceStart(final AlignMatch match,
                                       final int sequence_length,
                                       final boolean reverse_position) 
  {
    if(reverse_position) 
      return sequence_length - match.getQuerySequenceStart () + 1;
    else 
      return match.getQuerySequenceStart ();
  }

  /**
   *  Return the end position in the query sequence of the given AlignMatch,
   *  taking into account the current orientation of the sequences.  If the
   *  reverse_position argument is true reverse the complemented match
   *  coordinates will be returned.
   **/
  protected static int getRealQuerySequenceEnd(final AlignMatch match,
                                     final int sequence_length,
                                     final boolean reverse_position) 
  {
    if(reverse_position)
      return sequence_length - match.getQuerySequenceEnd() + 1;
    else 
      return match.getQuerySequenceEnd();
  }

  /**
   * Remove AlignMatch from the all_matches array
   * @param collection of indexes to be removed from the array
   */
  private void removeMatches(Vector<Integer> index)
  {
    AlignMatch tmp_matches[] = new AlignMatch[all_matches.length - index.size()];
    int start_old  = 0;
    int curr_index = 0;

    for(int i=0; i<index.size(); i++)
    {
      curr_index = ( index.get(i) ).intValue();

      if(curr_index !=0)
        System.arraycopy(all_matches, start_old, tmp_matches,
                         start_old-i, curr_index-start_old);
      start_old = curr_index+1;
    }

    if(start_old < all_matches.length)
      System.arraycopy(all_matches, start_old, tmp_matches,
                       start_old-index.size(), 
                       all_matches.length-start_old);

    this.all_matches = new AlignMatch[tmp_matches.length];
    this.all_matches = tmp_matches;
  }
  
  /**
   * Used in contig reordering.
   * Split alignment matches that go over boundaries of contig
   * that is being moved. If the match contains gaps (i.e. subject
   * and query length of match not same length) then we don't know
   * where to split so these are deleted.
   *
   * @param subject true if contig is on the subject sequence
   * @param start of the contig
   * @param end of the contig
   * @param drop_position position where the contig is moving to
   */
  private void splitMatches(final boolean subject,
                            final int start, final int end,
                            final int drop_position)
  {
    int curr_index;
    int match_start;
    int match_end;
    int split_at;
    int delete_overlaps = -1;

    Vector<Integer> matches_to_split = new Vector<Integer>();
    Vector<Integer> removals = new Vector<Integer>();

    for(int i=0; i<all_matches.length; i++)
    {
      if(subject)
      {
        match_start = all_matches[i].getSubjectSequenceStart();
        match_end   = all_matches[i].getSubjectSequenceEnd();
      }
      else
      {
        match_start = all_matches[i].getQuerySequenceStart();
        match_end   = all_matches[i].getQuerySequenceEnd();
      }

      // catch matches that span 2 contigs that need moving
      if( (match_start < start && match_end >= start) ||
          (match_start <= end   && match_end > end)   ||
          (match_start < drop_position && match_end >= drop_position) )
      {
        // check query and subject ranges same length
        //
        if( (all_matches[i].getQuerySequenceEnd()-
             all_matches[i].getQuerySequenceStart()) !=
            (all_matches[i].getSubjectSequenceStart()-
             all_matches[i].getSubjectSequenceEnd()) )
        {
          // this match extends past end of contig
          if(delete_overlaps == -1) 
          {
            Range q_range = all_matches[i].getQuerySequenceRange();
            Range s_range = all_matches[i].getSubjectSequenceRange();
            delete_overlaps = JOptionPane.showConfirmDialog(null,
                  "Found a match that extends past the boundary of a contig\n"+
                  "with query and subject ranges of different lengths:\n"+
                  q_range.toString()+
                  " ("+q_range.getCount()+")\n"+
                  s_range.toString()+
                  " ("+s_range.getCount()+")\n"+
                  "Delete all such matches?",
                  "Delete Overlapping Matches",
                  JOptionPane.YES_NO_OPTION);
         }

         if(delete_overlaps == JOptionPane.YES_OPTION)
           removals.add(new Integer(i));
        }
        else
          matches_to_split.add(new Integer(i));
      }
    }

    // now split the matches
    AlignMatch tmp_matches[] = new AlignMatch[all_matches.length+
                                              matches_to_split.size()];
    System.arraycopy(all_matches, 0, tmp_matches,
                     0, all_matches.length);
    int tmp_match_start;

    for(int i=0; i<matches_to_split.size(); i++)
    {
      curr_index = matches_to_split.get(i).intValue();

      //
      if(subject)
      {
        match_start = tmp_matches[curr_index].getSubjectSequenceStart();
        match_end   = tmp_matches[curr_index].getSubjectSequenceEnd();
      }
      else
      {
        match_start = tmp_matches[curr_index].getQuerySequenceStart();
        match_end   = tmp_matches[curr_index].getQuerySequenceEnd();
      }

      if(match_start <= start && match_end >= start)
        split_at = start-1;
      else if(match_start <= end   && match_end >= end)
        split_at = end;
      else
        split_at = drop_position-1;

      tmp_matches[curr_index].setRange(match_start, split_at, subject, false);
      tmp_matches[all_matches.length+i] = AlignMatch.copy(tmp_matches[curr_index]);
      tmp_matches[all_matches.length+i].setRange(split_at+1, match_end, subject, false);

      tmp_match_start = match_start;
      //
      if(!subject)
      {
        match_start = tmp_matches[curr_index].getSubjectSequenceStart();
        match_end   = tmp_matches[curr_index].getSubjectSequenceEnd();
      }
      else
      {
        match_start = tmp_matches[curr_index].getQuerySequenceStart();
        match_end   = tmp_matches[curr_index].getQuerySequenceEnd();
      }

      split_at = match_start+(tmp_match_start-split_at);
      tmp_matches[curr_index].setRange(match_start, split_at, !subject, false);

      if(tmp_matches[curr_index].isRevMatch())
        split_at--;
      else
        split_at++;

      tmp_matches[all_matches.length+i].setRange(split_at, match_end, !subject, false);
    }

    this.all_matches = new AlignMatch[tmp_matches.length];
    this.all_matches = tmp_matches;

    if(removals.size() > 0)
      removeMatches(removals);
  }

  /**
   * Reorder matches on reordering contigs
   * 
   * @param subject true if reordering the subject
   * @param start of the contig
   * @param end of the contig
   */
  protected void reorder(boolean subject, final int start, final int end,
                         final int drop_position)
  {
    // find matches that cross contig boundary
    splitMatches(subject,start,end,drop_position);

    int match_start;
    int match_end;
    //int delete_overlaps = -1;

    //Vector removals = new Vector();

    for(int i = 0; i < all_matches.length; ++i)
    {
      if(subject)
      {
        match_start = all_matches[i].getSubjectSequenceStart();
        match_end   = all_matches[i].getSubjectSequenceEnd();
      }
      else
      {
        match_start = all_matches[i].getQuerySequenceStart();
        match_end   = all_matches[i].getQuerySequenceEnd();
      }

      if(match_start >= start || match_start>=drop_position)
      {
        if(drop_position < start)
        {
          if(match_start <= start &&
             match_end < start)
          {
            match_start = match_start + (end-start+1);
            match_end   = match_end   + (end-start+1);
            all_matches[i].setRange(match_start, match_end, subject, false);
          } 
          else if(match_start >= start && match_start <= end &&  // within contig
                  match_end   >= start && match_end <= end)
          {
            match_start = match_start - (start - drop_position);
            match_end   = match_end   - (start - drop_position);
            all_matches[i].setRange(match_start, match_end, subject, false);
          }
        }
        else
        {
          if(match_start < end &&   // within contig
             match_end <= end)
          {
            match_start = match_start + (drop_position-end-1);
            match_end   = match_end   + (drop_position-end-1);
            all_matches[i].setRange(match_start, match_end, subject, false); 
          }
          else if(match_start > end && match_start < drop_position &&
                  match_end  > end && match_end < drop_position)
          {
            match_start = match_start - (end-start+1);
            match_end   = match_end   - (end-start+1);
            all_matches[i].setRange(match_start, match_end, subject, false);
          }
        }
      }
    }

  }

  /**
   * Flip the matches within a contig and deleted those that overlap
   * with other contigs
   * @param subject true if flipping the subject
   * @param start of the contig
   * @param end of the contig
   */
  protected void flippingContig(boolean subject, int start, int end)
  {
    int match_start;
    int match_end;
    int delete_overlaps = -1;
    Vector<Integer> removals = new Vector<Integer>();

    for(int i = 0; i < all_matches.length; ++i)
    {
      if(subject)
      {
        match_start = all_matches[i].getSubjectSequenceStart();
        match_end   = all_matches[i].getSubjectSequenceEnd();
			
       if(match_start >= start &&
          match_end <= end)
       {
         match_start = end - (match_start - start);
         match_end   = end - (match_end -start);
         all_matches[i].setRange(match_start, match_end, subject, true);
       }
       else if( (match_start >= start && match_start <= end) ||
                (match_end <= end && match_end >= start) )
       {
	  // this match extends past end of contig
	  if(delete_overlaps == -1)
	    delete_overlaps = JOptionPane.showConfirmDialog(null,
	                 "Found a match extending past the boundary of the contig:\n"+
	                 match_start+".."+match_end+
	                 "\nDelete all such matches?",
	                 "Delete Overlapping Matches",
	                 JOptionPane.YES_NO_OPTION);
		      
		      
          if(delete_overlaps == JOptionPane.YES_OPTION)
  	    removals.add(new Integer(i));
        }
      }
      else
      {
        match_start = all_matches[i].getQuerySequenceStart();
        match_end   = all_matches[i].getQuerySequenceEnd();

        if(match_start >= start &&
           match_end <= end)
        {
          match_start = end - (match_start - start);
          match_end   = end - (match_end -start);
          all_matches[i].setRange(match_start, match_end, subject, true);
        }
        else if( (match_start >= start && match_start <= end) ||
                 (match_end <= end && match_end >= start) )
        {
          // this match extends past end of contig
          if(delete_overlaps == -1)
            delete_overlaps = JOptionPane.showConfirmDialog(null,
	                 "Found a match extending past the boundary of the contig:\n"+
	                 match_start+".."+match_end+
	                 "\nDelete all such matches?",
	                 "Delete Overlapping Matches",
	                 JOptionPane.YES_NO_OPTION);
		      
		      
	  if(delete_overlaps == JOptionPane.YES_OPTION)
	    removals.add(new Integer(i));
        }
      }
    }

    if(removals.size() > 0)
      removeMatches(removals);
  }
  
  /**
  *
  * Find regions where there are no matches.
  * @return vector of coordinates of regions where there are 
  *         no matches.
  *
  */
  protected Vector<Integer[]> getDifferenceCoords(boolean subject)
  {
    final int length;
    final boolean flipped;

    if(subject)
    {
      flipped = subjectIsRevComp();
      length  = getSubjectForwardStrand().getSequenceLength();
    }
    else
    {
      flipped = queryIsRevComp();
      length  = getQueryForwardStrand().getSequenceLength();
    }

    AlignMatchComparator comparator = new AlignMatchComparator(subject, length,
                                                               flipped);
    int imatch = 0;
    final AlignMatch[] sorted_all_matches = new AlignMatch[all_matches.length];   
    for(int i = 0; i < sorted_all_matches.length; ++i)
    {
      if(isVisible(all_matches[i]))
        sorted_all_matches[imatch++] = AlignMatch.copy(all_matches[i]);
    }

    // sort alignment matches based on where they start
    Arrays.sort(sorted_all_matches, 0, imatch, comparator);

    int start = 1;
    Vector<Integer[]> differences = new Vector<Integer[]>();

    // find & record regions of no match
    for(int i = 0; i < imatch; ++i)
    {
      Integer coords[];
      final AlignMatch this_match = sorted_all_matches[i];

      if(this_match == null)
        continue;

      int this_start;
      int this_end;

      if(subject)
      {
        this_start = getRealSubjectSequenceStart(this_match,
                                                 length, flipped);
        this_end   = getRealSubjectSequenceEnd(this_match,
                                               length, flipped);
      }
      else
      {
        this_start = getRealQuerySequenceStart(this_match,
                                               length, flipped);
        this_end   = getRealQuerySequenceEnd(this_match,
                                             length, flipped);
      }

      if(this_start > this_end)
      {
        int tmp_this_start = this_start;
        this_start = this_end;
        this_end   = tmp_this_start;
      }

      if(i == 0 && this_start > 1)
      {
        coords = new Integer[2];
        coords[0] = new Integer(start);
        coords[1] = new Integer(this_start-1);
        differences.add(coords);
      }

      if(this_end > start)
        start = this_end;

      if(i < imatch-1)
      {
        int next_start = 0;

        // quicksort doesn't get it quite right so check
        // start position of several ahead
        for(int j=1; j<11; j++)
        {
          if(i+j > imatch-1)
            continue;

          final AlignMatch next_match = sorted_all_matches[i+j];
          int jstart;
          int jend;

          if(subject)
          {
            jstart = getRealSubjectSequenceStart(next_match, length, flipped);
            jend   = getRealSubjectSequenceEnd(next_match, length, flipped);
          }
          else
          {
            jstart = getRealQuerySequenceStart(next_match, length, flipped);
            jend   = getRealQuerySequenceEnd(next_match, length, flipped);
          }

          if(jend < jstart)
            jstart = jend;
  
          if(j == 1 || jstart < next_start)
            next_start = jstart;
        }

        if(next_start > start+1)
        {
          coords = new Integer[2];
          coords[0] = new Integer(start+1);
          coords[1] = new Integer(next_start-1);
          differences.add(coords);
        }
      }
      else if(i == imatch-1 &&
              this_end < length &&
              start+1 < length)
      {
        coords = new Integer[2];
        coords[0] = new Integer(start+1);
        coords[1] = new Integer(length);
        differences.add(coords);
      }
    }

    return differences;
  }

  /**
   *  Return the screen x positions of the corners of the given match.  The
   *  order is Top Left, Top Right, Bottom Left, Bottom Right, unless the
   *  match is an inversion, in which case it will be TL,TR,BR,BL.  Returns
   *  null if and only if the match is not currently visible.
   **/
  private int[] getMatchCoords(final int canvas_width, final AlignMatch this_match,
                               final int subject_length, final int query_length,
                               final boolean subject_flipped, final boolean query_flipped,
                               final float base_width, final float query_base_width,
                               final int subject_start, final int query_start,
                               final boolean subject_is_rev_comp, final boolean query_is_rev_comp,
                               final boolean is_rev_match)
  {
    int subject_sequence_start =
      getRealSubjectSequenceStart(this_match,
                                  subject_length, subject_flipped);
    int subject_sequence_end =
      getRealSubjectSequenceEnd(this_match,
                                subject_length, subject_flipped);
    int query_sequence_start =
      getRealQuerySequenceStart(this_match,
                                query_length, query_flipped);
    int query_sequence_end =
      getRealQuerySequenceEnd(this_match,
                              query_length, query_flipped);

    // add one base because we want to draw to the end of the base
    if(subject_is_rev_comp) 
      subject_sequence_start += 1;
    else 
      subject_sequence_end += 1;

    if(is_rev_match && !query_is_rev_comp ||
       !is_rev_match && query_is_rev_comp) 
      query_sequence_start += 1;
    else 
      query_sequence_end += 1;
 
    // this is the base that is at the left of the screen
    final int subject_start_x = getScreenPosition(base_width, subject_sequence_start,
                                                  subject_start);
    final int subject_end_x   = getScreenPosition(base_width, subject_sequence_end,
                                                  subject_start);

    final int query_start_x = getScreenPosition(query_base_width, query_sequence_start,
                                                query_start);
    final int query_end_x   = getScreenPosition(query_base_width, query_sequence_end,
                                                query_start);

    boolean subject_off_left  = false;
    boolean subject_off_right = false;
    boolean query_off_left    = false;
    boolean query_off_right   = false;

    if(subject_start_x < 0 && subject_end_x < 0) 
      subject_off_left = true;

    if(subject_start_x >= canvas_width && subject_end_x >= canvas_width) 
      subject_off_right = true;

    if(query_start_x < 0 && query_end_x < 0) 
      query_off_left = true;

    if(query_start_x >= canvas_width && query_end_x >= canvas_width) 
      query_off_right = true;

    if((subject_off_left ? 1 : 0) +
       (query_off_left ? 1 : 0) +
       (subject_off_right ? 1 : 0) +
       (query_off_right ? 1 : 0) == 2) 
      return null;
    else 
    {
      final int[] return_values = new int[4];

      return_values[0] = subject_start_x;
      return_values[1] = subject_end_x;
      return_values[2] = query_start_x;
      return_values[3] = query_end_x;

      return return_values;
    }
  }

  /**
   *  Return the current forward Strand of the subject EntryGroup.
   **/
  private Strand getSubjectForwardStrand() 
  {
    return getSubjectEntryGroup().getBases().getForwardStrand();
  }

  /**
   *  Return the current forward Strand of the query EntryGroup.
   **/
  private Strand getQueryForwardStrand()
  {
    return getQueryEntryGroup().getBases().getForwardStrand();
  }

  /**
   *  Return the subject EntryGroup that was passed to the constructor.
   **/
  protected EntryGroup getSubjectEntryGroup()
  {
    return subject_entry_group;
  }

  /**
   *  Return the subject EntryGroup that was passed to the constructor.
   **/
  protected EntryGroup getQueryEntryGroup() 
  {
    return query_entry_group;
  }

  /**
   *  Return the reference of the subject FeatureDisplay.
   **/
  private FeatureDisplay getSubjectDisplay() 
  {
    return subject_feature_display;
  }

  /**
   *  Return the reference of the query FeatureDisplay.
   **/
  private FeatureDisplay getQueryDisplay() 
  {
    return query_feature_display;
  }

  /**
   *  Returns true if and only if the subject sequence has been flipped since
   *  this object was created.
   **/
  protected boolean subjectIsRevComp() 
  {
    final Strand current_subject_forward_strand = getSubjectForwardStrand();

    if(getOrigSubjectForwardStrand() == current_subject_forward_strand ^
       getSubjectDisplay().isRevCompDisplay()) 
      return false;
    else 
      return true;
  }

  /**
   *  Returns true if and only if the query sequence has been flipped since
   *  this object was created.
   **/
  protected boolean queryIsRevComp() 
  {
    final Strand current_query_forward_strand = getQueryForwardStrand ();

    if(getOrigQueryForwardStrand() == current_query_forward_strand ^
       getQueryDisplay().isRevCompDisplay()) 
      return false;
    else 
      return true;
  }

  /**
   *  Return the forward Strand of the subject EntryGroup from when the
   *  Comparator was created.
   **/
  private Strand getOrigSubjectForwardStrand() 
  {
    return orig_subject_forward_strand;
  }

  /**
   *  Return the forward Strand of the query EntryGroup from when the
   *  Comparator was created.
   **/
  private Strand getOrigQueryForwardStrand() 
  {
    return orig_query_forward_strand;
  }

  /**
   *  Arrange for the two FeatureDisplay objects to scroll in parallel.
   **/
  protected void lockDisplays() 
  {
    displays_are_locked = true;
    repaint();
  }

  /**
   *  Arrange for the two FeatureDisplay objects to scroll independently.
   **/
  protected void unlockDisplays() 
  {
    displays_are_locked = false;
    repaint();
  }

  /**
   *  Return true if and only if the displays are currently locked.
   **/
  protected boolean displaysAreLocked() 
  {
    return displays_are_locked;
  }

  /**
   *  Toggle whether the two displays are locked.
   **/
  protected void toggleDisplayLock() 
  {
    displays_are_locked = !displays_are_locked;
    repaint();
  }

  /**
   *  Calling this method will temporarily disable selectFromQueryRange() and
   *  selectFromSubjectRange() until enableSelection().  This is need to allow
   *  the selections of the top and bottom FeatureDisplays to be set without
   *  changing which AlignMatches are selected.
   **/
  protected void disableSelection() 
  {
    disable_selection_from_ranges = true;
  }

  /**
   *  Enable selectFromQueryRange() and selectFromSubjectRange().
   **/
  protected void enableSelection() 
  {
    disable_selection_from_ranges = false;
  }

  /**
   *  Convert a base position into a screen x coordinate.
   **/
  private int getScreenPosition(final float base_width,
                                final int base_position,
                                final int screen_start_base) 
  {
    final float base_pos =  base_width *
                            (base_position - screen_start_base);
    if(base_pos > 30000) 
      return 30000;
    else if(base_pos < -30000) 
        return -30000;
    else 
      return (int)base_pos;
  }

  /**
   *  Return an array of colours that will be used for colouring the matches
   *  (depending on score).
   **/
  private void makeColours() 
  {
    red_percent_id_colours  = new Color[NUMBER_OF_SHADES];
    blue_percent_id_colours = new Color[NUMBER_OF_SHADES];

    for(int i = 0; i < blue_percent_id_colours.length; ++i) 
    {
      final int shade_value = 255 - (int) (256 * i / NUMBER_OF_SHADES);
      red_percent_id_colours[i] = new Color (255, shade_value, shade_value);
      blue_percent_id_colours[i] = new Color (shade_value, shade_value, 255);
    }
  }

  /**
   *  Return the ComparisonData object that was passed to the constructor.
   **/
  private ComparisonData getComparisonData () 
  {
    return comparison_data;
  }

  public class AlignMatchComparator implements Comparator<AlignMatch>
  {
    final private boolean subject;
    final private int length;
    final private boolean flipped;

    public AlignMatchComparator(boolean subject, int length,
                                boolean flipped)
    {
      this.subject = subject;
      this.length  = length;
      this.flipped = flipped;
    }

    public int compare(AlignMatch c1,AlignMatch c2) throws ClassCastException
    {        
      int start1;
      int start2;
    
      if(subject)
      {
        start1   = getRealSubjectSequenceStart(c1, length, flipped);
        int end1 = getRealSubjectSequenceEnd(c1, length, flipped);
        if(end1 < start1)
          start1 = end1;

        start2   = getRealSubjectSequenceStart(c2, length, flipped);
        int end2 = getRealSubjectSequenceEnd(c2, length, flipped);
        if(end2 < start2)
          start2 = end2;
      }
      else
      {
        start1   = getRealQuerySequenceStart(c1, length, flipped);
        int end1 = getRealQuerySequenceEnd(c1, length, flipped);
        if(end1 < start1)
          start1 = end1;

        start2   = getRealQuerySequenceStart(c2, length, flipped);
        int end2 = getRealQuerySequenceEnd(c2, length, flipped);
        if(end2 < start2)
          start2 = end2;
      }

      if(start1 < start2)
        return -1;
      else if(start1 > start2)
        return 1;
      else
        return 0;
    }
  }

  public class ColorChooserShades extends JPanel
                              implements javax.swing.event.ChangeListener
  {
    private static final long serialVersionUID = 1L;
    private JPanel bannerPanel = new JPanel(new BorderLayout());
    private JColorChooser tcc;
    private Color definedColour[] = new Color[NUMBER_OF_SHADES];
    private JLabel banner[]       = new JLabel[NUMBER_OF_SHADES];
    private JSlider scaleColour;
    private Box bacross = Box.createHorizontalBox();

    public ColorChooserShades(String title, Color initialColour)
    {
      super(new BorderLayout());

      Dimension d  = new Dimension(35, 35);
      double fract = (maximum_percent_id - minimum_percent_id)/
                    (NUMBER_OF_SHADES* 0.999);

      for(int i = 0; i < NUMBER_OF_SHADES; ++i)
      {
        int percent_id = (int)(i*fract)+minimum_percent_id;
        banner[i] = new JLabel(" "+percent_id+"  ",JLabel.CENTER);
        banner[i].setOpaque(true);
        banner[i].setPreferredSize(d);
        bacross.add(banner[i]);
      }

      //Set up color chooser for setting text color
      tcc = new JColorChooser(initialColour);
      tcc.getSelectionModel().addChangeListener(ColorChooserShades.this);
      tcc.setBorder(BorderFactory.createTitledBorder(title));

      //set scale
      scaleColour = new JSlider(0,20,3);
      scaleColour.addChangeListener(ColorChooserShades.this);

      makeColours();

      //Set up the banner at the top of the window
      colourBox();
      bannerPanel.add(bacross, BorderLayout.CENTER);
      bannerPanel.add(scaleColour, BorderLayout.SOUTH);
      bannerPanel.setBorder(BorderFactory.createTitledBorder("% ID Scale"));

      add(bannerPanel, BorderLayout.CENTER);
      add(tcc, BorderLayout.SOUTH);
    }

    public void stateChanged(javax.swing.event.ChangeEvent e) 
    {
      makeColours(); 
      colourBox();
      bannerPanel.repaint();
      repaint();
    } 

    private void colourBox()
    {
      for(int i = 0; i < NUMBER_OF_SHADES; ++i)
      {
        banner[i].setBackground(definedColour[i]);
        banner[i].repaint();
      }
    }

    /**
     *  Return an array of colours that will be used for colouring the matches
     *  (depending on score).
     **/
    private void makeColours()
    {
      Color newColour = tcc.getColor();

      for(int i = 0; i < NUMBER_OF_SHADES; ++i)
      {
        int R = newColour.getRed();
        int G = newColour.getGreen();
        int B = newColour.getBlue();

        int scale = (NUMBER_OF_SHADES-i)*scaleColour.getValue()*5;
        if((R+scale) < 255)
          R += scale;
        if((G+scale) < 255)
          G += scale;
        if((B+scale) < 255)
          B += scale;
         
        definedColour[i] = new Color(R,G,B);
      }
    }

    public Color[] getDefinedColour()
    {
      return definedColour;
    }
  }
}
