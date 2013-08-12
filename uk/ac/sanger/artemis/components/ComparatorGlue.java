/* ComparatorGlue.java
 *
 * created: Tue Sep 11 2001
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2001  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/ComparatorGlue.java,v 1.4 2005-12-02 14:58:57 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.sequence.*;

import uk.ac.sanger.artemis.util.OutOfRangeException;

import javax.swing.*;

/**
 *  This class contains the Event glue needed to combine two FeatureDisplays
 *  and an AlignmentViewer.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: ComparatorGlue.java,v 1.4 2005-12-02 14:58:57 tjc Exp $
 **/

public class ComparatorGlue {
  /**
   *  Create a new ComparatorGlue object that glues the given components
   *  togeather.
   **/
  public ComparatorGlue (final JFrame parent_frame,
                         final FeatureDisplay subject_feature_display,
                         final FeatureDisplay query_feature_display,
                         final AlignmentViewer alignment_viewer) {
    this.parent_frame = parent_frame;
    this.subject_feature_display = subject_feature_display;
    this.query_feature_display = query_feature_display;
    this.alignment_viewer = alignment_viewer;

    getSubjectEntryGroup ().addFeatureChangeListener (getSubjectSelection ());
    getQueryEntryGroup ().addEntryChangeListener (getQuerySelection ());

    addDisplayListeners (subject_feature_display, query_feature_display);

    orig_subject_forward_strand =
      getSubjectEntryGroup ().getBases ().getForwardStrand ();
    orig_query_forward_strand =
      getQueryEntryGroup ().getBases ().getForwardStrand ();

    makeAlignmentEventListener ();
  }

  /**
   *  Wire-up the two FeatureDisplay objects with DisplayAdjustmentListeners.
   **/
  private void addDisplayListeners (final FeatureDisplay subject_display,
                                    final FeatureDisplay query_display) 
  {

    subject_listener = new DisplayAdjustmentListener() 
    {
      public void displayAdjustmentValueChanged(DisplayAdjustmentEvent e) 
      {
        if(e.getType() == DisplayAdjustmentEvent.REV_COMP_EVENT) 
        {
          getAlignmentViewer().unlockDisplays();
          return;
        }
        else if(e.getType() == DisplayAdjustmentEvent.CONTIG_REV_COMP_EVENT) 
        {
          getAlignmentViewer().flippingContig(true, e.getStart(), e.getEnd());
          getAlignmentViewer().unlockDisplays();
          getSubjectDisplay().setFirstBase(e.getStart());
          return;
        }
        else if(e.getType() == DisplayAdjustmentEvent.CONTIG_REORDER)
        {
          getAlignmentViewer().reorder(true, e.getStart(), e.getEnd(), 
                                        e.getDropPosition());
          getAlignmentViewer().unlockDisplays();

          if(e.getStart() < e.getDropPosition())
            getSubjectDisplay().setFirstBase(e.getStart());
          else
            getSubjectDisplay().setFirstBase(e.getDropPosition());

          return;
        }

        subject_display.removeDisplayAdjustmentListener(subject_listener);
        query_display.removeDisplayAdjustmentListener(query_listener);

        query_display.setScaleFactor (e.getScaleFactor ());

        if (getAlignmentViewer ().displaysAreLocked () &&
            (e.getType () == DisplayAdjustmentEvent.SCROLL_ADJUST_EVENT ||
             e.getType () == DisplayAdjustmentEvent.ALL_CHANGE_ADJUST_EVENT)) {
          final int difference = e.getStart () - subject_first_base_position;

          query_first_base_position =
            query_display.getForwardBaseAtLeftEdge () + difference;

          query_display.setFirstBase (query_first_base_position);
        }

        subject_first_base_position = e.getStart ();

        subject_display.addDisplayAdjustmentListener (subject_listener);
        query_display.addDisplayAdjustmentListener (query_listener);
      }
    };

    query_listener = new DisplayAdjustmentListener () 
    {
      public void displayAdjustmentValueChanged (DisplayAdjustmentEvent e) 
      {
        if (e.getType () == DisplayAdjustmentEvent.REV_COMP_EVENT) 
        {
          getAlignmentViewer().unlockDisplays();
          return;
        }
        else if (e.getType () == DisplayAdjustmentEvent.CONTIG_REV_COMP_EVENT) 
        {
          getAlignmentViewer().flippingContig(false, e.getStart(), e.getEnd());
          getAlignmentViewer().unlockDisplays();
          getQueryDisplay().setFirstBase(e.getStart());
          return;
        }
        else if(e.getType() == DisplayAdjustmentEvent.CONTIG_REORDER)
        {
          getAlignmentViewer().reorder(false, e.getStart(), e.getEnd(),
                                       e.getDropPosition());
          getAlignmentViewer().unlockDisplays();

          if(e.getStart() < e.getDropPosition())
            getSubjectDisplay().setFirstBase(e.getStart());
          else
            getSubjectDisplay().setFirstBase(e.getDropPosition());

          return;
        }

        subject_display.removeDisplayAdjustmentListener (subject_listener);
        query_display.removeDisplayAdjustmentListener (query_listener);

        subject_display.setScaleFactor (e.getScaleFactor ());

        if (getAlignmentViewer ().displaysAreLocked () &&
            (e.getType () == DisplayAdjustmentEvent.SCROLL_ADJUST_EVENT ||
             e.getType () == DisplayAdjustmentEvent.ALL_CHANGE_ADJUST_EVENT)) 
        {
          final int difference = e.getStart () - query_first_base_position;

          subject_first_base_position =
            subject_display.getForwardBaseAtLeftEdge () + difference;

          subject_display.setFirstBase (subject_first_base_position);
        }

        query_first_base_position = e.getStart ();

        subject_display.addDisplayAdjustmentListener (subject_listener);
        query_display.addDisplayAdjustmentListener (query_listener);
      }
    };

    subject_display.addDisplayAdjustmentListener (subject_listener);
    query_display.addDisplayAdjustmentListener (query_listener);

    final DisplayAdjustmentListener subject_align_listener =
      new DisplayAdjustmentListener () 
      {
        public void displayAdjustmentValueChanged (DisplayAdjustmentEvent e) 
        {
          getAlignmentViewer ().setSubjectSequencePosition (e);
        }
      };

    subject_display.addDisplayAdjustmentListener (subject_align_listener);

    final DisplayAdjustmentListener query_align_listener =
      new DisplayAdjustmentListener () {
        public void displayAdjustmentValueChanged (DisplayAdjustmentEvent e) {
          getAlignmentViewer ().setQuerySequencePosition (e);
        }
      };

    query_display.addDisplayAdjustmentListener (query_align_listener);
  }

  /**
   *  Make an AlignmentListener which calls alignAt ().
   **/
  private void makeAlignmentEventListener () {
    final AlignmentListener listener =
      new AlignmentListener () {
        public void alignMatchChosen (AlignmentEvent e) {
          alignAt (e.getMatch ());
        }
      };

    getAlignmentViewer ().addAlignmentListener (listener);
  }

  /**
   *  Scroll both the subject and query so that they centre on the given
   *  AlignMatch.
   **/
  private void alignAt (final AlignMatch align_match) {
    getAlignmentViewer ().unlockDisplays ();
    getAlignmentViewer ().disableSelection ();

    maybeRevCompQuery (align_match);

    final int subject_length = getSubjectForwardStrand ().getSequenceLength ();
    final int query_length = getQueryForwardStrand ().getSequenceLength ();

    int subject_sequence_start =
      AlignmentViewer.getRealSubjectSequenceStart (align_match,
                                                   subject_length,
                                                   getAlignmentViewer ().subjectIsRevComp ());
    int subject_sequence_end =
      AlignmentViewer.getRealSubjectSequenceEnd (align_match,
                                                 subject_length,
                                                 getAlignmentViewer ().subjectIsRevComp ());
    int query_sequence_start =
      AlignmentViewer.getRealQuerySequenceStart (align_match,
                                                 query_length,
                                                 getAlignmentViewer ().queryIsRevComp ());
    int query_sequence_end =
      AlignmentViewer.getRealQuerySequenceEnd (align_match,
                                               query_length,
                                               getAlignmentViewer ().queryIsRevComp ());

    if (getSubjectDisplay ().isRevCompDisplay ()) {
      subject_sequence_start =
        getSubjectDisplay ().getSequenceLength () - subject_sequence_start + 1;
      subject_sequence_end =
        getSubjectDisplay ().getSequenceLength () - subject_sequence_end + 1;
    }

    if (getQueryDisplay ().isRevCompDisplay ()) {
      query_sequence_start =
        getQueryDisplay ().getSequenceLength () - query_sequence_start + 1;
      query_sequence_end =
        getQueryDisplay ().getSequenceLength () - query_sequence_end + 1;
    }

    final int new_subject_base =
      subject_sequence_start +
      (subject_sequence_end - subject_sequence_start) / 2;
    getSubjectDisplay ().makeBaseVisible (new_subject_base);

    try {
      final Strand subject_strand;

      if (align_match.isRevMatch () ^
          !getOrigSubjectForwardStrand ().isForwardStrand ()) {
        subject_strand = getSubjectReverseStrand ();
      } else {
        subject_strand = getSubjectForwardStrand ();
      }

      final MarkerRange new_subject_marker =
        subject_strand.makeMarkerRangeFromRawPositions (subject_sequence_start,
                                                        subject_sequence_end);
      getSubjectSelection ().setMarkerRange (new_subject_marker);
    } catch (OutOfRangeException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    }

    final int new_query_base =
      query_sequence_start +
      (query_sequence_end - query_sequence_start) / 2;
    getQueryDisplay ().makeBaseVisible (new_query_base);

    try {
      final Strand query_strand;

      if (getOrigQueryForwardStrand ().isForwardStrand ()) {
        query_strand = getQueryForwardStrand ();
      } else {
        query_strand = getQueryReverseStrand ();
      }

      final MarkerRange new_query_marker_range =
        query_strand.makeMarkerRangeFromRawPositions (query_sequence_start,
                                                      query_sequence_end);

      getQuerySelection ().setMarkerRange (new_query_marker_range);
    } catch (OutOfRangeException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    }

    getAlignmentViewer ().lockDisplays ();
    getAlignmentViewer ().enableSelection ();
  }

  /**
   *  Called by alignAt () to reverse and complement the query sequence if the
   *  given AlignMatch is currently a match from the subject sequence to the
   *  reverse complement of the query sequence.
   **/
  private void maybeRevCompQuery (final AlignMatch align_match) {
    if (!getAlignmentViewer ().offerToFlip ()) {
      return;
    }

    // true if and only if exactly one of the query and subject is rev-comped
    final boolean display_is_rev_comped;

    if (getOrigSubjectForwardStrand ().isForwardStrand () ^
        getOrigQueryForwardStrand ().isForwardStrand () ^
        getSubjectDisplay ().isRevCompDisplay () ^
        getQueryDisplay ().isRevCompDisplay ()) {
      display_is_rev_comped = true;
    } else {
      display_is_rev_comped = false;
    }

    if (align_match.isRevMatch () ^ display_is_rev_comped) {
      final YesNoDialog yes_no_dialog =
        new YesNoDialog (parent_frame,
                         "reverse and complement query sequence display?");
      if (yes_no_dialog.getResult ()) {
        if (getQueryDisplay ().isRevCompDisplay ()) {
          getQueryDisplay ().setRevCompDisplay (false);
        } else {
          getQueryDisplay ().setRevCompDisplay (true);
        }
      }
    }
  }

/*  *//**
   *  Send an event to those object listening for it.
   *  @param listeners A Vector of the objects that the event should be sent
   *    to.
   *  @param event The event to send
   **//*
  private void fireAction (final Vector listeners, final EventObject event) {
    final Vector targets;
    // copied from a book - synchronising the whole method might cause a
    // deadlock
    synchronized (this) {
      targets = (Vector) listeners.clone ();
    }

    for (int i = 0 ; i < targets.size () ; ++i) {
      GotoListener target = (GotoListener) targets.elementAt (i);

      if (event instanceof GotoEvent) {
        final GotoListener goto_listener = (GotoListener) target;
        goto_listener.performGoto ((GotoEvent) event);
      } else {
        throw new Error ("EntryEdit.fireAction () - unknown event");
      }
    }
  }*/

  /**
   *  Return the forward Strand of the subject EntryGroup from when the
   *  Comparator was created.
   **/
  public Strand getOrigSubjectForwardStrand () {
    return orig_subject_forward_strand;
  }

  /**
   *  Return the forward Strand of the query EntryGroup from when the
   *  Comparator was created.
   **/
  public Strand getOrigQueryForwardStrand () {
    return orig_query_forward_strand;
  }

  /**
   *  Return the GotoEventSource object used by the subject FeatureDisplay.
   **/
  public GotoEventSource getSubjectGotoEventSource () {
    return getSubjectDisplay ().getGotoEventSource ();
  }

  /**
   *  Return the GotoEventSource object used by the query FeatureDisplay.
   **/
  public GotoEventSource getQueryGotoEventSource () {
    return getQueryDisplay ().getGotoEventSource ();
  }

  /**
   *  Return the reference of the EntryGroup in view in the subject
   *  FeatureDisplay.
   **/
  public EntryGroup getSubjectEntryGroup () {
    return getSubjectDisplay ().getEntryGroup ();
  }

  /**
   *  Return the reference of the EntryGroup in view in the query
   *  FeatureDisplay.
   **/
  public EntryGroup getQueryEntryGroup () {
    return getQueryDisplay ().getEntryGroup ();
  }

  /**
   *  Return the Selection object used by the subject FeatureDisplay.
   **/
  public Selection getSubjectSelection () {
    return getSubjectDisplay ().getSelection ();
  }

  /**
   *  Return the Selection object used by the query FeatureDisplay.
   **/
  public Selection getQuerySelection () {
    return getQueryDisplay ().getSelection ();
  }

  /**
   *  Return the AlignmentViewer that was created in the constructor.
   **/
  private AlignmentViewer getAlignmentViewer () {
    return alignment_viewer;
  }

  /**
   *  Return the reference of the subject FeatureDisplay.
   **/
  public FeatureDisplay getSubjectDisplay () {
    return subject_feature_display;
  }

  /**
   *  Return the reference of the query FeatureDisplay.
   **/
  public FeatureDisplay getQueryDisplay () {
    return query_feature_display;
  }

  /**
   *  Return the current forward Strand of the subject EntryGroup.
   **/
  private Strand getSubjectForwardStrand () {
    return getSubjectEntryGroup ().getBases ().getForwardStrand ();
  }

  /**
   *  Return the current forward Strand of the query EntryGroup.
   **/
  private Strand getQueryForwardStrand () {
    return getQueryEntryGroup ().getBases ().getForwardStrand ();
  }

  /**
   *  Return the current reverse Strand of the subject EntryGroup.
   **/
  private Strand getSubjectReverseStrand () {
    return getSubjectEntryGroup ().getBases ().getReverseStrand ();
  }

  /**
   *  Return the current reverse Strand of the query EntryGroup.
   **/
  private Strand getQueryReverseStrand () {
    return getQueryEntryGroup ().getBases ().getReverseStrand ();
  }

  final private JFrame parent_frame;
  final private FeatureDisplay subject_feature_display;
  final private FeatureDisplay query_feature_display;
  final private AlignmentViewer alignment_viewer;

  private int subject_first_base_position = 1;
  private int query_first_base_position = 1;

  private DisplayAdjustmentListener subject_listener = null;
  private DisplayAdjustmentListener query_listener = null;

  /**
   *  The forward Strand of the subject EntryGroup when the Comparator was
   *  created.
   **/
  final private Strand orig_subject_forward_strand;

  /**
   *  The forward Strand of the query EntryGroup when the Comparator was
   *  created.
   **/
  final private Strand orig_query_forward_strand;

}
