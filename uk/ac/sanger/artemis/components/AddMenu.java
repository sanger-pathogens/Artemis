/* AddMenu.java
 *
 * created: Tue Dec 29 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/AddMenu.java,v 1.2 2004-12-03 18:11:28 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.sequence.*;
import uk.ac.sanger.artemis.plot.CodonUsageAlgorithm;

import uk.ac.sanger.artemis.util.*;
import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.io.RangeVector;
import uk.ac.sanger.artemis.io.Location;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierParseException;
import uk.ac.sanger.artemis.io.InvalidQualifierException;
import uk.ac.sanger.artemis.io.InvalidRelationException;
import uk.ac.sanger.artemis.io.InvalidKeyException;
import uk.ac.sanger.artemis.io.LocationParseException;
import uk.ac.sanger.artemis.io.EntryInformationException;

import java.awt.*;
import java.awt.event.*;

import javax.swing.*;

/**
 *  A Menu with commands that add new features/entries to an EntryGroup.  This
 *  should have been called CreateMenu.
 *
 *  @author Kim Rutherford
 *  @version $Id: AddMenu.java,v 1.2 2004-12-03 18:11:28 tjc Exp $
 **/

public class AddMenu extends SelectionMenu {
  /**
   *  The shortcut for "Create From Base Range".
   **/
  final static KeyStroke CREATE_FROM_BASE_RANGE_KEY =
    KeyStroke.getKeyStroke (KeyEvent.VK_C, InputEvent.CTRL_MASK);

  final static public int CREATE_FROM_BASE_RANGE_KEY_CODE = KeyEvent.VK_C;

  /**
   *  Create a new AddMenu object.
   *  @param frame The JFrame that owns this JMenu.
   *  @param selection The Selection that the commands in the menu will
   *    operate on.
   *  @param entry_group The EntryGroup object where new features/entries will
   *    be added.
   *  @param goto_event_source The object the we will call makeBaseVisible ()
   *    on.
   *  @param base_plot_group The BasePlotGroup associated with this JMenu -
   *    needed to call getCodonUsageAlgorithm()
   *  @param menu_name The name of the new menu.
   **/
  public AddMenu (final JFrame frame,
                  final Selection selection,
                  final EntryGroup entry_group,
                  final GotoEventSource goto_event_source,
                  final BasePlotGroup base_plot_group,
                  final String menu_name) {
    super (frame, menu_name, selection);

    this.entry_group = entry_group;
    this.base_plot_group = base_plot_group;

    new_feature_item = new JMenuItem ("New Feature");
    new_feature_item.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent event) {
        makeNewFeature ();
      }
    });

    add (new_feature_item);

    create_feature_from_range_item =
      new JMenuItem ("Create Feature From Base Range");
    create_feature_from_range_item.setAccelerator (CREATE_FROM_BASE_RANGE_KEY);
    create_feature_from_range_item.addActionListener(new ActionListener () {
      public void actionPerformed (ActionEvent event) {
        createFeatureFromBaseRange (getParentFrame (), getSelection (),
                                    entry_group, getGotoEventSource ());
      }
    });

    add (create_feature_from_range_item);

    create_intron_features_item =
      new JMenuItem ("Create Intron Features");
    create_intron_features_item.addActionListener(new ActionListener () {
      public void actionPerformed (ActionEvent event) {
        createIntronFeatures (getParentFrame (), getSelection (),
                              entry_group);
      }
    });

    add (create_intron_features_item);

    create_exon_features_item =
      new JMenuItem ("Create Exon Features");
    create_exon_features_item.addActionListener(new ActionListener () {
      public void actionPerformed (ActionEvent event) {
        createExonFeatures (getParentFrame (), getSelection (),
                            entry_group);
      }
    });

    add (create_exon_features_item);

    create_gene_features_item =
      new JMenuItem ("Create Gene Features");
    create_gene_features_item.addActionListener(new ActionListener () {
      public void actionPerformed (ActionEvent event) {
        createGeneFeatures (getParentFrame (), getSelection (),
                            entry_group);
      }
    });

    add (create_gene_features_item);

    addSeparator ();

    new_entry_item = new JMenuItem ("New Entry");
    new_entry_item.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent event) {
        makeNewEntry ();
      }
    });

    add (new_entry_item);

    addSeparator ();

    mark_orfs_with_size_item = new JMenuItem ("Mark Open Reading Frames ...");

    mark_orfs_with_size_item.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent event) {
        markORFSWithSize (false);
      }
    });

    add (mark_orfs_with_size_item);


    mark_empty_orfs_with_size_item = new JMenuItem ("Mark Empty ORFs ...");

    mark_empty_orfs_with_size_item.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent event) {
        markORFSWithSize (true);
      }
    });

    add (mark_empty_orfs_with_size_item);


    mark_orfs_range_item = new JMenuItem ("Mark ORFs In Range ...");
    mark_orfs_range_item.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent event) {
        markOpenReadingFramesInRange ();
      }
    });

    add (mark_orfs_range_item);


    mark_pattern_item = new JMenuItem ("Mark From Pattern ...");
    mark_pattern_item.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent event) {
        makeFeaturesFromPattern ();
      }
    });

    add (mark_pattern_item);

    mark_ambiguities_item = new JMenuItem ("Mark Ambiguities");
    mark_ambiguities_item.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent event) {
        markAmbiguities ();
      }
    });

    add (mark_ambiguities_item);
  }

  /**
   *  Create a new AddMenu object.
   *  @param frame The JFrame that owns this JMenu.
   *  @param selection The Selection that the commands in the menu will
   *    operate on.
   *  @param entry_group The EntryGroup object where new features/entries will
   *    be added.
   *  @param goto_event_source The object the we will call makeBaseVisible ()
   *    on.
   *  @param base_plot_group The BasePlotGroup associated with this JMenu -
   *    needed to call getCodonUsageAlgorithm()
   **/
  public AddMenu (final JFrame frame,
                  final Selection selection,
                  final EntryGroup entry_group,
                  final GotoEventSource goto_event_source,
                  final BasePlotGroup base_plot_group) {
    this (frame, selection, entry_group,
          goto_event_source, base_plot_group, "Create");
  }

  /**
   *  Create a new feature in the default Entry of the entry group.  See
   *  EntryGroup.createFeature () for details.
   **/
  private void makeNewFeature () {
    if (entry_group.size () > 0) {
      if (entry_group.getDefaultEntry () == null) {
        new MessageDialog (getParentFrame (), "There is no default entry");
      } else {

        try {
          entry_group.getActionController ().startAction ();

          final Feature new_feature = entry_group.createFeature ();

          final FeatureEdit feature_edit =
            new FeatureEdit (new_feature, entry_group, getSelection (),
                             getGotoEventSource ());

          final ActionListener cancel_listener =
            new ActionListener () {
              public void actionPerformed (ActionEvent e) {
                try {
                  new_feature.removeFromEntry ();
                } catch (ReadOnlyException exception) {
                  throw new Error ("internal error - unexpected exception: " +
                                   exception);
                }
              }
            };

          feature_edit.addCancelActionListener (cancel_listener);

          feature_edit.addApplyActionListener (new ActionListener () {
            public void actionPerformed (ActionEvent e) {
              // after apply is pressed cancel should not remove the new
              // feature
              feature_edit.removeCancelActionListener (cancel_listener);
            }
          });

          feature_edit.setVisible(true);
        } catch (ReadOnlyException e) {
          new MessageDialog (getParentFrame (), "feature not created: " +
                             "the default entry is read only");
        } finally {
          entry_group.getActionController ().endAction ();
        }
      }
    } else {
      new MessageDialog (getParentFrame (),
                         "Cannot make a feature without an existing entry");
    }
  }

  /**
   *  Create a new Entry in the first Entry of the entry group.
   **/
  private void makeNewEntry () {
    entry_group.createEntry ();
  }

  /**
   *  Create a new Feature in entry_group from the selected range of bases and
   *  then display a FeatureEdit component for the new Feature.
   *  @param frame The JFrame to use for MessageDialog components.
   *  @param selection The Selection containing the sequence to create the
   *    feature from.
   *  @param entry_group The EntryGroup to create the feature in.
   *  @param goto_event_source Needed to create a FeatureEdit component.
   **/
  static void createFeatureFromBaseRange (final JFrame frame,
                                          final Selection selection,
                                          final EntryGroup entry_group,
                                          final GotoEventSource
                                            goto_event_source) {
    try {
      entry_group.getActionController ().startAction ();

      if (!checkForSelectionRange (frame, selection)) {
        return;
      }

      final MarkerRange range = selection.getMarkerRange ();
      final Entry default_entry = entry_group.getDefaultEntry ();

      if (default_entry == null) {
        new MessageDialog (frame, "There is no default entry");
        return;
      }

      try {
        final Location new_location = range.createLocation ();

        /*final*/ Feature temp_feature;

        try {
          temp_feature = default_entry.createFeature (Key.CDS, new_location);
        } catch (EntryInformationException e) {
          // use the default key instead

          final Key default_key =
            default_entry.getEntryInformation ().getDefaultKey ();

          try {
            temp_feature =
              default_entry.createFeature (default_key, new_location);
          } catch (EntryInformationException ex) {
            throw new Error ("internal error - unexpected exception: " + ex);
          }
        }

        final Feature new_feature = temp_feature;

        selection.setMarkerRange (null);
        selection.set (new_feature);

        final FeatureEdit feature_edit =
          new FeatureEdit (new_feature, entry_group,
                           selection, goto_event_source);

        final ActionListener cancel_listener =
          new ActionListener () {
            public void actionPerformed (ActionEvent e) {
              try {
                new_feature.removeFromEntry ();
                selection.setMarkerRange (range);
              } catch (ReadOnlyException exception) {
                throw new Error ("internal error - unexpected exception: " +
                                 exception);
              }
            }
          };

        feature_edit.addCancelActionListener (cancel_listener);

        feature_edit.addApplyActionListener (new ActionListener () {
          public void actionPerformed (ActionEvent e) {
            // after apply is pressed cancel should not remove the new
            // feature
            feature_edit.removeCancelActionListener (cancel_listener);
          }
        });

        feature_edit.setVisible(true);
      } catch (ReadOnlyException e) {
        new MessageDialog (frame, "feature not created: " +
                           "the default entry is read only");
      } catch (OutOfRangeException e) {
        throw new Error ("internal error - unexpected exception: " + e);
      }
    } finally {
      entry_group.getActionController ().endAction ();
    }
  }

  /**
   *  Create a new intron between each pair of exons in the selected CDS
   *  features.  The introns are created in the Entry that contains the CDSs.
   *  @param frame The Frame to use for MessageDialog components.
   *  @param selection The Selection containing the CDS features to create the
   *    introns for.
   *  @param entry_group The EntryGroup to create the features in.
   **/
  static void createIntronFeatures (final JFrame frame,
                                    final Selection selection,
                                    final EntryGroup entry_group) {
    try {
      entry_group.getActionController ().startAction ();

      if (!checkForSelectionFeatures (frame, selection)) {
        return;
      }

      final FeatureVector selected_features = selection.getAllFeatures ();

      for (int feature_index = 0 ;
           feature_index < selected_features.size () ;
           ++feature_index) {

        final Feature selection_feature =
          selected_features.elementAt (feature_index);

        if (!selection_feature.isProteinFeature ()) {
          continue;
        }

        final Location cds_location = selection_feature.getLocation ();

        final RangeVector cds_ranges = cds_location.getRanges ();

        if (cds_ranges.size () < 2) {
          continue;
        }

        if (cds_location.isComplement ()) {
          cds_ranges.reverse ();
        }

        for (int range_index = 0 ;
             range_index < cds_ranges.size () - 1 ;
             ++range_index) {
          final int end_of_range_1 =
            cds_ranges.elementAt (range_index).getEnd ();
          final int start_of_range_2 =
            cds_ranges.elementAt (range_index + 1).getStart ();

          if (end_of_range_1 > start_of_range_2) {
            // ignore - the exons overlap so there is no room for an intron
            continue;
          }

          final Range new_range;

          try {
            new_range = new Range (end_of_range_1 + 1,
                                   start_of_range_2 - 1);
          } catch (OutOfRangeException e) {
            throw new Error ("internal error - unexpected exception: " + e);
          }

          final RangeVector intron_ranges = new RangeVector ();

          intron_ranges.add (new_range);

          final Key intron_key = new Key ("intron");
          final Location intron_location =
            new Location (intron_ranges, cds_location.isComplement ());
          final QualifierVector qualifiers = new QualifierVector ();

          try {
            selection_feature.getEntry ().createFeature (intron_key,
                                                         intron_location,
                                                         qualifiers);
          } catch (ReadOnlyException e) {
            throw new Error ("internal error - unexpected exception: " + e);
          } catch (EntryInformationException e) {
            throw new Error ("internal error - unexpected exception: " + e);
          } catch (OutOfRangeException e) {
            throw new Error ("internal error - unexpected exception: " + e);
          }
        }
      }
    } finally {
      entry_group.getActionController ().endAction ();
    }
  }

  /**
   *  Create a new exon for each FeatureSegment in the selected CDS features.
   *  The exons are created in the Entry that contains the CDSs.
   *  @param frame The Frame to use for MessageDialog components.
   *  @param selection The Selection containing the CDS features to create the
   *    exons for.
   *  @param entry_group The EntryGroup to create the features in.
   **/
  static void createExonFeatures (final JFrame frame,
                                  final Selection selection,
                                  final EntryGroup entry_group) {
    try {
      entry_group.getActionController ().startAction ();

      if (!checkForSelectionFeatures (frame, selection)) {
        return;
      }

      final FeatureVector selected_features = selection.getAllFeatures ();

      for (int feature_index = 0 ;
           feature_index < selected_features.size () ;
           ++feature_index) {

        final Feature selection_feature =
          selected_features.elementAt (feature_index);

        if (!selection_feature.isProteinFeature ()) {
          continue;
        }

        final Location cds_location = selection_feature.getLocation ();

        final RangeVector cds_ranges = cds_location.getRanges ();

        for (int range_index = 0 ;
             range_index < cds_ranges.size () ;
             ++range_index) {
          final Range this_range = cds_ranges.elementAt (range_index);

          final RangeVector exon_ranges = new RangeVector ();

          exon_ranges.add (this_range);

          final Key exon_key = new Key ("exon");
          final Location exon_location =
            new Location (exon_ranges, cds_location.isComplement ());
          final QualifierVector qualifiers = new QualifierVector ();

          try {
            selection_feature.getEntry ().createFeature (exon_key,
                                                         exon_location,
                                                         qualifiers);
          } catch (ReadOnlyException e) {
            throw new Error ("internal error - unexpected exception: " + e);
          } catch (EntryInformationException e) {
            throw new Error ("internal error - unexpected exception: " + e);
          } catch (OutOfRangeException e) {
            throw new Error ("internal error - unexpected exception: " + e);
          }
        }
      }
    } finally {
      entry_group.getActionController ().endAction ();
    }
  }

  /**
   *  Create a new gene for each of the selected CDS features.
   *  @param frame The Frame to use for MessageDialog components.
   *  @param selection The Selection containing the CDS features.
   *  @param entry_group The EntryGroup to create the features in.
   **/
  static void createGeneFeatures (final JFrame frame,
                                  final Selection selection,
                                  final EntryGroup entry_group) {
    /*
     *  XXX - FIXME - include 5'UTR at the start and 3'UTR at the end and if
     *  two (or more) CDSs have the same primary name create one gene feature
     *  that covers them.
     */    
    try {
      entry_group.getActionController ().startAction ();

      if (!checkForSelectionFeatures (frame, selection)) {
        return;
      }

      final FeatureVector selected_features = selection.getAllFeatures ();

      for (int feature_index = 0 ;
           feature_index < selected_features.size () ;
           ++feature_index) {

        final Feature selection_feature =
          selected_features.elementAt (feature_index);

        if (!selection_feature.isProteinFeature ()) {
          continue;
        }

        final Range max_range = selection_feature.getMaxRawRange ();
        final boolean complement_flag =
          selection_feature.getLocation ().isComplement ();

        final RangeVector ranges = new RangeVector ();
        ranges.add (max_range);

        final Key gene_key = new Key ("gene");
        final Location gene_location =
          new Location (ranges, complement_flag);
        final QualifierVector qualifiers = new QualifierVector ();

        try {
          selection_feature.getEntry ().createFeature (gene_key,
                                                       gene_location,
                                                       qualifiers);
        } catch (ReadOnlyException e) {
          throw new Error ("internal error - unexpected exception: " + e);
        } catch (EntryInformationException e) {
          throw new Error ("internal error - unexpected exception: " + e);
        } catch (OutOfRangeException e) {
          throw new Error ("internal error - unexpected exception: " + e);
        }
      }
    } finally {
      entry_group.getActionController ().endAction ();
    }
  }

  /**
   *  Open a TextRequester to ask the user for the minimum ORF size then call
   *  markOpenReadingFrames ().
   *  @param empty_only If true only those ORFS that don't already contain a
   *    segment will be marked.
   **/
  private void markORFSWithSize (final boolean empty_only) {
    final int default_minimum_orf_size =
      Options.getOptions ().getMinimumORFSize ();

    final TextRequester text_requester =
      new TextRequester ("minimum open reading frame size?",
                         18, String.valueOf (default_minimum_orf_size));

    text_requester.addTextRequesterListener (new TextRequesterListener () {
      public void actionPerformed (final TextRequesterEvent event) {
        if (event.getType () == TextRequesterEvent.CANCEL) {
          return;
        }

        final String requester_text = event.getRequesterText ().trim ();

        if (requester_text.length () == 0) {
          return;
        }

        try {
          final int minimum_orf_size =
            Integer.valueOf (requester_text).intValue ();

          markOpenReadingFrames (minimum_orf_size, empty_only);
        } catch (NumberFormatException e) {
          new MessageDialog (getParentFrame (),
                             "this is not a number: " + requester_text);
        }
      }
    });

    text_requester.setVisible (true);
  }

  /**
   *  Create a new Feature for each open reading frame.
   *  @param minimum_orf_size All the returned ORFs will be at least this many
   *    amino acids long.
   *  @param empty_only If true only those ORFS that don't already contain a
   *    segment will be marked.
   **/
  private void markOpenReadingFrames (final int minimum_orf_size,
                                      final boolean empty_only) {
    try {
      final Entry new_entry =
        entry_group.createEntry ("ORFS_" + minimum_orf_size + '+');

      final int sequence_length = entry_group.getSequenceLength ();

      final Strand forward_strand =
        entry_group.getBases ().getForwardStrand ();

      final MarkerRange forward_range =
        forward_strand.makeMarkerRangeFromPositions (1, sequence_length);

      markOpenReadingFrames (new_entry, forward_range, minimum_orf_size,
                             empty_only);

      final Strand backward_strand =
        entry_group.getBases ().getReverseStrand ();

      final MarkerRange backward_range =
        backward_strand.makeMarkerRangeFromPositions (1, sequence_length);

      markOpenReadingFrames (new_entry, backward_range, minimum_orf_size,
                             empty_only);
    } catch (OutOfRangeException e) {
      throw new Error ("internal error - unexpected OutOfRangeException");
    }
  }

  /**
   *  Create a new Feature for each open reading frame.  The minimum size of
   *  the ORFS is specified in the options file.
   **/
  private void markOpenReadingFramesInRange () {
    if (!checkForSelectionRange (getParentFrame (), getSelection ())) {
      return;
    }

    final int default_minimum_orf_size =
      Options.getOptions ().getMinimumORFSize ();

    final TextRequester text_requester =
      new TextRequester ("minimum open reading frame size?",
                         18, String.valueOf (default_minimum_orf_size));

    text_requester.addTextRequesterListener (new TextRequesterListener () {
      public void actionPerformed (final TextRequesterEvent event) {
        if (event.getType () == TextRequesterEvent.CANCEL) {
          return;
        }

        final String requester_text = event.getRequesterText ().trim ();

        if (requester_text.length () == 0) {
          return;
        }

        try {
          final int minimum_orf_size =
            Integer.valueOf (requester_text).intValue ();

          final Entry new_entry =
            entry_group.createEntry ("ORFS_" + minimum_orf_size + '+');

          final MarkerRange selection_range =
            getSelection ().getMarkerRange ();

          markOpenReadingFrames (new_entry, selection_range, minimum_orf_size,
                                 false);


        } catch (NumberFormatException e) {
          new MessageDialog (getParentFrame (),
                             "this is not a number: " + requester_text);
        }
      }
    });

    text_requester.setVisible (true);

  }

  /**
   *  Create a new Feature in the given Entry for each open reading frame that
   *  overlaps the given range.  The minimum size of the ORFS is specified in
   *  the options file.
   *  @param entry The new features are created in this entry.
   *  @param search_range The range of bases to search for ORFs.
   *  @param minimum_orf_size All the returned ORFs will be at least this many
   *    amino acids long.
   *  @param empty_only If true only those ORFS that don't already contain a
   *    segment will be marked.
   **/
  private void markOpenReadingFrames (final Entry entry,
                                      final MarkerRange search_range,
                                      final int minimum_orf_size,
                                      final boolean empty_only) {
    final MarkerRange [] forward_orf_ranges =
      Strand.getOpenReadingFrameRanges (search_range, minimum_orf_size);

    for (int i = 0 ; i < forward_orf_ranges.length ; ++i) {
      final MarkerRange this_range = forward_orf_ranges[i];

      final Feature new_feature;

      try {
        new_feature = makeFeatureFromMarkerRange (entry, this_range, Key.CDS);
      } catch (EntryInformationException e) {
        new MessageDialog (getParentFrame (), "cannot continue: " +
                           "the default entry does not support CDS features");
        return;
      } catch (ReadOnlyException e) {
        new MessageDialog (getParentFrame (), "cannot continue: " +
                           "the default entry is read only");
        return;
      }

      if (empty_only && overlapsAnActiveSegment (new_feature)) {
        try {
          new_feature.removeFromEntry ();
        } catch (ReadOnlyException exception) {
          throw new Error ("internal error - unexpected exception: " +
                           exception);
        }
      }
    }
  }

  /**
   *  Return true if and only if the given feature overlaps (and is in the
   *  same frame as) a segment in an active entry.
   **/
  private boolean overlapsAnActiveSegment (final Feature test_feature) {
    final Range test_feature_range = test_feature.getMaxRawRange ();

    FeatureVector overlapping_features;

    try {
      overlapping_features =
        entry_group.getFeaturesInRange (test_feature_range);
    } catch (OutOfRangeException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    }

    for (int feature_index = 0 ;
         feature_index < overlapping_features.size () ;
         ++feature_index ) {

      final Feature current_feature =
        overlapping_features.elementAt (feature_index);

      if (current_feature != test_feature && current_feature.isCDS ()) {
        final FeatureSegmentVector segments = current_feature.getSegments ();

        for (int segment_index = 0;
             segment_index < segments.size () ;
             ++segment_index) {
          final FeatureSegment this_segment =
            segments.elementAt (segment_index);

          if (test_feature_range.overlaps (this_segment.getRawRange ())) {
            final int test_feature_frame =
              test_feature.getSegments ().elementAt (0).getFrameID ();
            final int this_segment_frame = this_segment.getFrameID ();

            if (test_feature_frame == this_segment_frame) {
              return true;
            }
          }
        }
      }
    }

    return false;
  }

  /**
   *  Make a new Feature from the given MarkerRange in the given Entry.  The
   *  new feature will be given the key 'CDS' and it's location will match the
   *  MarkerRange.
   *  @param entry The new feature is created in this entry.
   *  @param range The location of the new feature.
   *  @param key The key give the new feature
   *  @exception EntryInformationException Thrown if this Entry does not
   *    support features with the given key.  Also thrown if any of these
   *    qualifiers aren't supported: note, label or gene.
   **/
  private Feature makeFeatureFromMarkerRange (final Entry entry,
                                              final MarkerRange range,
                                              final Key key)
      throws EntryInformationException, ReadOnlyException {
    try {
      final Location new_location = range.createLocation ();

      final QualifierVector qualifiers = new QualifierVector ();

      qualifiers.setQualifier (new Qualifier ("note", "none"));

      final Feature new_feature =
        entry.createFeature (key, new_location, qualifiers);

      final CodonUsageAlgorithm codon_usage_algorithm =
        base_plot_group.getCodonUsageAlgorithm ();

      if (codon_usage_algorithm != null) {
        int score =
          (int) (codon_usage_algorithm.getFeatureScore (new_feature) * 50);

        if (score < 0) {
          score = 0;
        }

        if (score > 100) {
          score = 100;
        }

        final String score_string = String.valueOf (score);
        new_feature.addQualifierValues (new Qualifier ("score",
                                                       score_string));

        final int var_colour = 255 - score * 5 / 2;
        final String colour_string = var_colour + " " + var_colour + " 255";
        new_feature.addQualifierValues (new Qualifier ("colour",
                                                       colour_string));
      }

      return new_feature;

    } catch (OutOfRangeException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    }
  }

  /**
   *  This method will ask the user for a BasePattern (using a TextRequester
   *  component) then search the sequence for the given pattern and make a new
   *  feature from each match.  The new features will created in an Entry
   *  called "matches: <pattern>".
   **/
  private void makeFeaturesFromPattern () {
    final TextRequester text_requester =
      new TextRequester ("create features from this pattern:", 18, "");

    text_requester.addTextRequesterListener (new TextRequesterListener () {
      public void actionPerformed (final TextRequesterEvent event) {
        final String pattern_string = event.getRequesterText ().trim ();

        try {
          if (pattern_string.length () == 0) {
            new MessageDialog (getParentFrame (), "the pattern is too short");
            return;
          }

          final BasePattern pattern = new BasePattern (pattern_string);

          makeFeaturesFromPattern (pattern);
        } catch (BasePatternFormatException e) {
          new MessageDialog (getParentFrame (),
                             "Illegal base pattern: " +
                             pattern_string);
        }
      }
    });

    text_requester.setVisible(true);
  }

  /**
   *  Search the sequence for the given pattern and make a new feature from
   *  each match.  The new features will created in an Entry called "matches:
   *  <pattern>".
   **/
  private void makeFeaturesFromPattern (final BasePattern pattern) {
    final MarkerRangeVector matches =
      pattern.findMatches (entry_group.getBases (),
                           null,        // search from start
                           entry_group.getSequenceLength ());

    if (matches.size () == 0) {
      new MessageDialog (getParentFrame (),
                         "no matches found for: " + pattern);
      return;
    }

    final int TOO_MANY_MATCHES = 100;

    if (matches.size () > TOO_MANY_MATCHES) {
      final YesNoDialog dialog =
        new YesNoDialog (getParentFrame (),
                         matches.size () + " matches, continue?");

      if (dialog.getResult ()) {
        // yes - continue
      } else {
        // no
        return;
      }
    }

    final Entry new_entry = entry_group.createEntry ("matches: " + pattern);

    final Key key = new_entry.getEntryInformation ().getDefaultKey ();

    for (int i = 0 ; i < matches.size () ; ++i) {
      try {
        final Feature new_feature =
          makeFeatureFromMarkerRange (new_entry, matches.elementAt (i), key);
        new_feature.setQualifier (new Qualifier ("note", pattern.toString ()));
      } catch (EntryInformationException e) {
        new MessageDialog (getParentFrame (), "cannot continue: " +
                           e.getMessage ());
        return;
      } catch (ReadOnlyException e) {
        new MessageDialog (getParentFrame (), "cannot continue: " +
                           "the default entry is read only");
        return;
      }
    }
  }

  /**
   *  Create a misc_feature for each block of ambiguous bases.  The new
   *  features will created in an Entry called "ambiguous bases".
   **/
  private void markAmbiguities () {
    Entry new_entry = null;

    final Bases bases = entry_group.getBases ();

    for (int i = 1 ; i <= bases.getLength () ; ++i) {
      try {
        if (! Bases.isLegalBase (bases.getBaseAt (i))) {
          final int start_index = i;

          while (i < bases.getLength () &&
                 ! Bases.isLegalBase (bases.getBaseAt (i + 1))) {
            ++i;
          }

          final int end_index = i;

          if (new_entry == null) {
            new_entry = entry_group.createEntry ("ambiguous bases");
          }

          final Range range = new Range (start_index, end_index);

          final String unsure_bases =
            bases.getSubSequence (range, Bases.FORWARD);

          final Location location = new Location (range);

          final QualifierVector qualifiers = new QualifierVector ();

          qualifiers.setQualifier (new Qualifier ("note", unsure_bases));

          final Feature feature =
            new_entry.createFeature (new Key ("unsure"), location, qualifiers);
        }
      } catch (ReadOnlyException e) {
        throw new Error ("internal error - unexpected exception: " + e);
      } catch (EntryInformationException e) {
        throw new Error ("internal error - unexpected exception: " + e);
      } catch (OutOfRangeException e) {
        throw new Error ("internal error - unexpected exception: " + e);
      }
    }

    if (new_entry == null) {
      new MessageDialog (getParentFrame (), "No ambiguities found");
    } else {
      if (new_entry.getFeatureCount () == 1) {
        new MessageDialog (getParentFrame (), "Created one feature");

      } else {
        new MessageDialog (getParentFrame (), "Created " +
                           new_entry.getFeatureCount () + " features");

      }
    }
  }

  /**
   *  Return the GotoEventSource object that was passed to the constructor.
   **/
  private GotoEventSource getGotoEventSource () {
    return goto_event_source;
  }

  /**
   *  The GotoEventSource object that was passed to the constructor.
   **/
  private GotoEventSource goto_event_source = null;

  /**
   *  The EntryGroup object that was passed to the constructor.
   **/
  private EntryGroup entry_group;

  private JMenuItem new_feature_item = null;
  private JMenuItem new_entry_item = null;
  private JMenuItem create_feature_from_range_item = null;
  private JMenuItem create_intron_features_item = null;
  private JMenuItem create_exon_features_item = null;
  private JMenuItem create_gene_features_item = null;
  private JMenuItem mark_orfs_with_size_item = null;
  private JMenuItem mark_empty_orfs_with_size_item = null;
  private JMenuItem mark_orfs_range_item = null;
  private JMenuItem mark_pattern_item = null;
  private JMenuItem mark_ambiguities_item = null;

  private BasePlotGroup base_plot_group;
}
