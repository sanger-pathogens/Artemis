/* SelectMenu.java
 *
 * created: Thu Jan 14 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/SelectMenu.java,v 1.13 2008-05-21 14:14:04 tjc Exp $
 **/

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.sequence.*;

import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.util.StringVector;
import uk.ac.sanger.artemis.components.genebuilder.GeneUtils;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.io.KeyVector;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.InvalidRelationException;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.Location;


import java.awt.Cursor;
import java.awt.Toolkit;
import java.awt.event.*;
import java.util.Vector;
import java.util.Enumeration;

import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.KeyStroke;


/**
 *  This menu has contains items such a "Select all", "Select none" and
 *  "Select by key".
 *
 *  @author Kim Rutherford
 **/

public class SelectMenu extends SelectionMenu 
{

  private static final long serialVersionUID = 1L;

  /**
   *  The EntryGroup object that was passed to the constructor.
   **/
  private EntryGroup entry_group;

  /**
   *  The GotoEventSource object that was passed to the constructor.
   **/
  private GotoEventSource goto_event_source;

  /**
   *  The shortcut for Select All
   **/
  final static KeyStroke SELECT_ALL_KEY =
    KeyStroke.getKeyStroke(KeyEvent.VK_A, 
                           Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()); //InputEvent.CTRL_MASK);

  /**
   *  The shortcut for Select None
   **/
  final static KeyStroke SELECT_NONE_KEY =
    KeyStroke.getKeyStroke (KeyEvent.VK_N, 
                           Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()); //InputEvent.CTRL_MASK);


  /**
   *  Warn if the user tries to select more than this number of features.
   **/
  final private int MAX_FEATURE_TO_SELECT_COUNT = 10000;

  /** busy cursor */
  private Cursor cbusy = new Cursor(Cursor.WAIT_CURSOR);
  /** done cursor */
  private Cursor cdone = new Cursor(Cursor.DEFAULT_CURSOR);


  /**
   *  Create a new SelectMenu object and all it's menu items.
   *  @param frame The JFrame that owns this JMenu.
   *  @param selection The Selection that the commands in the menu will
   *    operate on.
   *  @param goto_event_source The object the we will call makeBaseVisible ()
   *    on.
   *  @param entry_group The EntryGroup object where features and exons are
   *    selected from.
   *  @param base_plot_group The BasePlotGroup associated with this JMenu -
   *    needed to call getCodonUsageAlgorithm()
   *  @param menu_name The name of the new menu.
   **/
  public SelectMenu(final JFrame frame,
                    final Selection selection,
                    final GotoEventSource goto_event_source,
                    final EntryGroup entry_group,
                    final BasePlotGroup base_plot_group,
                    final String menu_name)
 {
   this(frame,selection,goto_event_source,entry_group,
        base_plot_group,null,null,menu_name);
 }

  /**
   *  Create a new SelectMenu object and all it's menu items.
   *  @param frame The JFrame that owns this JMenu.
   *  @param selection The Selection that the commands in the menu will
   *    operate on.
   *  @param goto_event_source The object the we will call makeBaseVisible ()
   *    on.
   *  @param entry_group The EntryGroup object where features and exons are
   *    selected from.
   *  @param base_plot_group The BasePlotGroup associated with this JMenu -
   *    needed to call getCodonUsageAlgorithm()
   *  @param menu_name The name of the new menu.
   **/
  public SelectMenu(final JFrame frame,
                    final Selection selection,
                    final GotoEventSource goto_event_source,
                    final EntryGroup entry_group,
                    final BasePlotGroup base_plot_group,
                    final AlignmentViewer alignQueryViewer,  
                    final AlignmentViewer alignSubjectViewer,
                    final String menu_name) 
 {
    super(frame, menu_name, selection);

    this.entry_group = entry_group;
    this.goto_event_source = goto_event_source;

    JMenuItem selector_item = new JMenuItem("Feature Selector ...");
    selector_item.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        new Selector(selection, goto_event_source, getEntryGroup(),
                     base_plot_group);
      }
    });
    add(selector_item);

    addSeparator();

    JMenuItem all_item = new JMenuItem("All");
    all_item.setAccelerator(SELECT_ALL_KEY);
    all_item.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        frame.setCursor(cbusy);
        selectAll();
        frame.setCursor(cdone);
      }
    });
    add(all_item);

    JMenuItem all_bases_item = new JMenuItem("All Bases");
    all_bases_item.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event) 
      {
        frame.setCursor(cbusy);
        selectAllBases();
        frame.setCursor(cdone);
      }
    });
    add(all_bases_item);


    if(alignQueryViewer != null || alignSubjectViewer != null)
    {
      JMenuItem all_diffs = new JMenuItem("All Features in Non-matching Regions");
      all_diffs.addActionListener(new ActionListener() 
      {
        public void actionPerformed(ActionEvent event) 
        {
          final JCheckBox cbox_strict = new JCheckBox("Strictly no overlapping regions",true);

          int n = JOptionPane.showConfirmDialog(null, cbox_strict,
                            "Selecting Features in Non-matching Regions",
                            JOptionPane.OK_CANCEL_OPTION,
                            JOptionPane.QUESTION_MESSAGE);

          if(n == JOptionPane.CANCEL_OPTION)
            return;

          frame.setCursor(cbusy);
          Vector diffs = null;
          if(alignQueryViewer == null || alignSubjectViewer == null)
          {
            if(alignQueryViewer != null)
              diffs = alignQueryViewer.getDifferenceCoords(false);
            else
              diffs = alignSubjectViewer.getDifferenceCoords(true);
          }
          else
          {
            final Vector diffs1 = alignQueryViewer.getDifferenceCoords(false);
            final Vector diffs2 = alignSubjectViewer.getDifferenceCoords(true);
            diffs = AddMenu.union(diffs1,diffs2);
          }

          if(diffs != null)
            selectAllFeatures(diffs,cbox_strict.isSelected());
          frame.setCursor(cdone);
        }
      });

      add(all_diffs);

    }

    JMenuItem none_item = new JMenuItem("None");
    none_item.setAccelerator(SELECT_NONE_KEY);
    none_item.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event) 
      {
        clearSelection();
      }
    });
    add(none_item);

    JMenuItem select_by_key_item = new JMenuItem("By Key");
    select_by_key_item.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event) 
      {
        selectByKeyPopup();
      }
    });

    add(select_by_key_item);

    JMenuItem select_non_pseudo_cds_item =
      new JMenuItem("CDS Features without /pseudogene");
    select_non_pseudo_cds_item.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)
      {
        // select all CDS features that do not have the /pseudo  or /pseudogene qualifier
        final FeaturePredicateConjunction predicate = new FeaturePredicateConjunction(
            new FeatureKeyQualifierPredicate(Key.CDS, "pseudo", false),
            new FeatureKeyQualifierPredicate(Key.CDS, "pseudogene", false),
            FeaturePredicateConjunction.AND);

        clearSelection();

        selectFeaturesByPredicate(predicate);
      }
    });
    add(select_non_pseudo_cds_item);

    
    final boolean isDatabaseGroup = GeneUtils.isDatabaseEntry( getEntryGroup() );
    final JMenuItem select_all_cds_item;
    if(isDatabaseGroup)
    {
      select_non_pseudo_cds_item.setEnabled(false);
      select_all_cds_item = new JMenuItem("All "+
          DatabaseDocument.EXONMODEL+" Features");
    }
    else
      select_all_cds_item = new JMenuItem("All CDS Features");
    select_all_cds_item.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event) 
      {
        final FeatureKeyPredicate predicate;
        if(isDatabaseGroup)
          predicate = new FeatureKeyPredicate(new Key(DatabaseDocument.EXONMODEL));
        else
          predicate = new FeatureKeyPredicate(Key.CDS);

        clearSelection();

        selectFeaturesByPredicate(predicate);
      }
    });
    add(select_all_cds_item);

    JMenuItem select_same_type_item = new JMenuItem("Same Key");
    select_same_type_item.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event) 
      {
        selectSameKey();
      }
    });
    add(select_same_type_item);

    JMenuItem select_matching_qualifiers_item =
      new JMenuItem("Features Matching Qualifier");
    select_matching_qualifiers_item.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event) 
      {
        selectMatchingQualifiers();
      }
    });
    add(select_matching_qualifiers_item);

    JMenuItem select_orf_item = new JMenuItem("Open Reading Frame");
    select_orf_item.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event) 
      {
        selectOrf();
      }
    });
    add(select_orf_item);

    final JMenuItem select_features_in_range_item =
      new JMenuItem("Features Overlapping Selection");
    select_features_in_range_item.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event) 
      {
        selectOverlappingFeatures();
      }
    });

    add(select_features_in_range_item);
    
    
    final JMenuItem select_features_contained_item =
        new JMenuItem("Features Within Selection");
    select_features_contained_item.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent event) 
        {
          selectContainedFeatures();
        }
      });

      add(select_features_contained_item);

    JMenuItem select_base_range_item = new JMenuItem("Base Range ...");
    select_base_range_item.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event) 
      {
        selectBaseRange();
      }
    });
    add(select_base_range_item);

    JMenuItem select_aa_range_in_feature_item = new JMenuItem("Feature AA Range ...");
    select_aa_range_in_feature_item.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        selectFeatureAARange();
      }
    });
    add(select_aa_range_in_feature_item);

    addSeparator();

    JMenuItem selection_toggle_item = new JMenuItem("Toggle Selection");
    selection_toggle_item.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        toggleSelection();
      }
    });
    add(selection_toggle_item);
  }


  /**
   *  Create a new SelectMenu object and all it's menu items.
   *  @param frame The JFrame that owns this JMenu.
   *  @param selection The Selection that the commands in the menu will
   *    operate on.
   *  @param goto_event_source The object the we will call makeBaseVisible ()
   *    on.
   *  @param entry_group The EntryGroup object where features and exons are
   *    selected from.
   *  @param base_plot_group The BasePlotGroup associated with this JMenu -
   *    needed to call getCodonUsageAlgorithm()
   **/
  public SelectMenu(final JFrame frame,
                    final Selection selection,
                    final GotoEventSource goto_event_source,
                    final EntryGroup entry_group,
                    final BasePlotGroup base_plot_group) 
  {
    this (frame, selection, goto_event_source, entry_group,
          base_plot_group, "Select");
  }


  /**
   *  Select all the features in the active entries.  If there are more than
   *  1000 features the user is asked for confirmation.
   **/
  private void selectAll() 
  {
    if(getEntryGroup().getAllFeaturesCount() >
       MAX_FEATURE_TO_SELECT_COUNT) 
    {
      final YesNoDialog dialog =
        new YesNoDialog(getParentFrame(),
                        "Are you sure you want to select " +
                         getEntryGroup().getAllFeaturesCount() +
                        " features?");

      if(!dialog.getResult ()) 
        return;
    }

    final FeatureVector all_features = getEntryGroup().getAllFeatures();

    if(getEntryGroup() instanceof SimpleEntryGroup) 
    {
      // special case for speed
      getSelection().set(all_features);
    } 
    else 
    {
      clearSelection();
      getSelection().add(all_features);
    }
  }

  /**
   *  Select all the bases.
   **/
  private void selectAllBases()
  {
    try 
    {
      final Strand strand = getEntryGroup().getBases().getForwardStrand();

      final MarkerRange new_range =
        strand.makeMarkerRangeFromPositions(1, strand.getSequenceLength());

      getSelection().setMarkerRange(new_range);
    }
    catch (OutOfRangeException e) 
    {
      throw new Error ("internal error - unexpected exception: " + e);
    }
  }


  /**
  *
  * Select all features in a given collection of regions.
  * @param regions   Vector of coordinates of non-matching region.
  * @param isStrict  true if not allowing any overlap with matching region.
  *
  */
  private void selectAllFeatures(final Vector regions, final boolean isStrict)
  {
    
    final FeaturePredicate predicate = new FeaturePredicate()
    {
      public boolean testPredicate(final Feature feature)
      {
        final Location loc = feature.getLocation();
        Key key = feature.getKey();
        if(key.equals("source"))
          return false;

        final int feat_start = loc.getFirstBase();
        final int feat_end   = loc.getLastBase();
//      final int feat_len   = feature.getBaseCount();

        int diff_start;
        int diff_end;
        Integer coords[];

        Enumeration eDiffs = regions.elements();
        while(eDiffs.hasMoreElements())
        {
          coords = (Integer[])eDiffs.nextElement();
          diff_start = coords[0].intValue();
          diff_end   = coords[1].intValue();

          if( isStrict && feat_start >= diff_start && feat_start <= diff_end &&
                           feat_end >= diff_start   && feat_end <= diff_end )
            return true;
          else if( !isStrict && 
                   ((feat_start >= diff_start && feat_start <= diff_end)  ||
                    (feat_end >= diff_start   && feat_end <= diff_end) ||
                    (diff_start >= feat_start && diff_start <= feat_end) ||
                    (diff_end >= feat_start   && diff_end <= feat_end )) )
            return true;
          else if(diff_start > feat_end)
            return false;
        }

        return false;
      }
    };

    selectFeaturesByPredicate(predicate);
  }

  /**
   *  Remove the all Features and FeatureSegments in the EntryGroup from the
   *  Selection.  If the current EntryGroup is a SimpleEntryGroup then this
   *  method just calls getSelection ().clear ();
   **/
  private void clearSelection () 
  {
    if (getEntryGroup () instanceof SimpleEntryGroup) {
      // special case for speed
      getSelection ().clear ();
    } else {
      getSelection ().setMarkerRange (null);

      final FeatureEnumeration test_enumerator = getEntryGroup ().features ();

      while (test_enumerator.hasMoreFeatures ()) {
        final Feature this_feature = test_enumerator.nextFeature ();

        getSelection ().remove (this_feature);
        getSelection ().removeSegmentsOf (this_feature);
      }
    }
  }

  /**
   *  Invert the selection - after calling this method the selection will
   *  contain only those features that were not in the selection before the
   *  call.
   **/
  private void toggleSelection () {
    final FeatureVector selected_features = getSelection ().getAllFeatures ();

    final FeatureVector all_features = getEntryGroup ().getAllFeatures ();

    final FeatureVector new_selection = new FeatureVector ();

    for (int i = 0 ; i < all_features.size () ; ++i) {
      final Feature this_feature = all_features.elementAt (i);

      if (!selected_features.contains (this_feature)) {
        new_selection.add (this_feature);
      }
    }

    clearSelection ();

    getSelection ().add (new_selection);
  }


  /**
   *  Popup a TextRequester component and then select all the features that
   *  have the same key as the user types.
   **/
  private void selectByKeyPopup () {
    final KeyChooser key_chooser;

    final EntryInformation default_entry_information =
      Options.getArtemisEntryInformation ();

    key_chooser = new KeyChooser (default_entry_information,
                                  new Key ("misc_feature"));

    key_chooser.getKeyChoice ().addItemListener (new ItemListener () {
      public void itemStateChanged (ItemEvent itemEvent) {
        selectByKey (key_chooser.getKeyChoice ().getSelectedItem ());
        key_chooser.setVisible (false);
        key_chooser.dispose ();
      }
    });

    key_chooser.getOKButton ().addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent actionEvent) {
        selectByKey (key_chooser.getKeyChoice ().getSelectedItem ());
        key_chooser.setVisible (false);
        key_chooser.dispose ();
      }
    });

    key_chooser.setVisible (true);
  }

  /**
   *  Select all the features that have the given key.
   **/
  private void selectByKey (final Key key) {
    final FeaturePredicate predicate = new FeatureKeyPredicate (key);

    selectFeaturesByPredicate (predicate);
  }

  /**
   *  Select all the features that have the given key and contains a qualifier
   *  with the given name and value.
   *  If key is null select any feature with the contains a qualifier with the
   *  given name and value.
   **/
  private void selectByQualifier (final Key key,
                                  final String name, final String value) {
    final FeaturePredicate predicate =
      new FeatureKeyQualifierPredicate (key, name, value);

    selectFeaturesByPredicate (predicate);
  }

  /**
   *  Select all the features that have the same key as the currently selected
   *  features.
   **/
  private void selectSameKey () {
    if (! checkForSelectionFeatures ()) {
      return;
    }

    final KeyVector seen_keys = new KeyVector ();

    final FeatureVector selected_features = getSelection ().getAllFeatures ();

    for (int i = 0 ; i < selected_features.size () ; ++i) {
      final Feature current_feature = selected_features.elementAt (i);

      if (!seen_keys.contains (current_feature.getKey ())) {
        seen_keys.add (current_feature.getKey ());
      }
    }

    clearSelection ();

    for (int i = 0 ; i < seen_keys.size () ; ++i) {
      selectByKey ((Key)seen_keys.get(i));
    }
  }

  /**
   *  Ask the user for a qualifier name, list all of the qualifier with that
   *  name in the currently selected feature then select all features that
   *  have that name and qualifier.
   **/
  private void selectMatchingQualifiers () {
    if (!checkForSelectionFeatures (getParentFrame (), getSelection ())) {
      return;
    }

    final FeatureVector selected_features = getSelection ().getAllFeatures ();

    if (selected_features.size () > 1) {
      new MessageDialog (getParentFrame (), "select only one feature");
      return;
    }

    final Feature selected_feature = selected_features.elementAt (0);

    final StringVector qualifier_names =
      Feature.getAllQualifierNames (selected_features);

    if (qualifier_names.size () == 0) {
      new MessageDialog (getParentFrame (), "feature has no qualifiers");
      return;
    }

    final ChoiceFrame name_choice_frame =
      new ChoiceFrame ("Select a qualifier name ...", qualifier_names);

    final JComboBox name_choice = name_choice_frame.getChoice ();

    class SelectListener implements ActionListener, ItemListener {
      public void doStuff () {
        selectMatchingQualifiersHelper (selected_feature,
                                        (String) name_choice.getSelectedItem ());
        name_choice_frame.setVisible (false);
        name_choice_frame.dispose ();
      }

      public void itemStateChanged (ItemEvent itemEvent) {
        doStuff ();
        // need to remove() because itemStateChanged() is called twice on the
        // Tru64 1.4.1 JVM
        name_choice.removeItemListener (this);
      }

      public void actionPerformed (ActionEvent actionEvent) {
        doStuff ();
        name_choice_frame.getOKButton ().removeActionListener (this);
      }
    };

    final SelectListener listener = new SelectListener ();

    name_choice.addItemListener (listener);
    name_choice_frame.getOKButton ().addActionListener (listener);

    name_choice_frame.setVisible (true);
  }

  /**
   *  Pop up a list of the values of the qualifiers with the given name in the
   *  given feature then select all features that have a qualifier with the
   *  same and value.
   **/
  private void selectMatchingQualifiersHelper (final Feature feature,
                                               final String name) {
    final Qualifier qualifier;

    try {
      qualifier = feature.getQualifierByName (name);
    } catch (InvalidRelationException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    }

    if (qualifier == null) {
      throw new Error ("internal error - unexpected null reference");
    }

    final StringVector qualifier_values = qualifier.getValues ();

    if (qualifier_values.size () == 0) {
      new MessageDialog (getParentFrame (),
                         "qualifier " + name + " has no values in ");
      return;
    }

    if (qualifier_values.size () == 1) {
      selectByQualifier (null, name, (String)qualifier_values.elementAt (0));
    } else {
      final ChoiceFrame value_choice_frame =
        new ChoiceFrame ("Select a qualifier value", qualifier_values);

      final JComboBox value_choice = value_choice_frame.getChoice ();

      value_choice.addItemListener (new ItemListener () {
        public void itemStateChanged (ItemEvent itemEvent) {
          selectByQualifier (null, name,
                             (String) value_choice.getSelectedItem ());
          value_choice_frame.setVisible (false);
          value_choice_frame.dispose ();
        }
      });

      value_choice_frame.getOKButton ().addActionListener (new ActionListener () {
        public void actionPerformed (ActionEvent actionEvent) {
          selectByQualifier (null, name,
                             (String) value_choice.getSelectedItem ());
          value_choice_frame.setVisible (false);
          value_choice_frame.dispose ();
        }
      });

      value_choice_frame.setVisible (true);
    }    
  }

  /**
   *  Select all the features that match the given predicate.
   **/
  private void selectFeaturesByPredicate (final FeaturePredicate predicate) {
    final FeatureVector new_selection_features = new FeatureVector ();

    final FeatureEnumeration feature_enum = getEntryGroup ().features ();

    while (feature_enum.hasMoreFeatures ()) {
      final Feature current_feature = feature_enum.nextFeature ();

      if (predicate.testPredicate (current_feature)) {
        new_selection_features.add (current_feature);
      }
    }

    clearSelection ();

    getSelection ().add (new_selection_features);
  }

  /**
   *  If a range of bases is selected, then select the ORF around those
   *  bases.
   **/
  private void selectOrf () {
    if (!checkForSelectionRange ()) {
      return;
    }

    final MarkerRange range = getSelection ().getMarkerRange ();
    final MarkerRange start_orf_range =
      Strand.getORFAroundMarker (range.getStart (), true);

    /*final*/ Marker end_marker;

    // get the marker of the first base of the last codon in the range (we
    // want to keep in the same frame)
    try {
      final int start_end_diff = range.getCount () - 1;

      final int mod_length = start_end_diff - start_end_diff % 3;

      end_marker = range.getStart ().moveBy (mod_length);
    } catch (OutOfRangeException e) {
      end_marker = range.getStart ();
    }

    final MarkerRange end_orf_range =
      Strand.getORFAroundMarker (end_marker, true);

    MarkerRange new_range = range;

    if (start_orf_range != null) {
      new_range = new_range.extendRange (start_orf_range);
    }

    if (end_orf_range != null) {
      new_range = new_range.extendRange (end_orf_range);
    }

    if (start_orf_range == null && end_orf_range == null &&
        range.getCount () <= 6) {   // 6 == two codons
      new MessageDialog (getParentFrame (),
                         "there is no open reading frame at the selected " +
                         "bases");
      return;
    }

    getSelection ().setMarkerRange (new_range);
  }

  /**
   *  This method will ask the user for a range of bases (such as 100..200),
   *  using a MarkerRangeRequester component, then selects that range.
   **/
  private void selectBaseRange () {
    final MarkerRangeRequester range_requester =
      new MarkerRangeRequester ("enter a range of bases (eg. 100..200):",
                                18, "");

    final MarkerRangeRequesterListener listener =
      new MarkerRangeRequesterListener () {
        public void actionPerformed (final MarkerRangeRequesterEvent event) {
          final MarkerRange marker_range =
            event.getMarkerRange (getEntryGroup ().getBases ());

          getSelection ().setMarkerRange (marker_range);

          makeBaseVisible (getSelection ().getStartBaseOfSelection ());
        }
      };

    range_requester.addMarkerRangeRequesterListener (listener);

    range_requester.setVisible (true);
  }

  /**
   *  This method will ask the user for a range of amino acids in a feature
   *  (such as 100..200), using a Requester component, then selects that
   *  range on the sequence.
   **/
  private void selectFeatureAARange () {
    if (!checkForSelectionFeatures (getParentFrame (), getSelection ())) {
      return;
    }

    final FeatureVector selected_features = getSelection ().getAllFeatures ();

    if (selected_features.size () > 1) {
      new MessageDialog (getParentFrame (), "select only one feature");
      return;
    }

    final Feature selected_feature = selected_features.elementAt (0);

    final MarkerRangeRequester range_requester =
      new MarkerRangeRequester ("enter a range of amino acids (eg. 100..200):",
                                18, "");

    final MarkerRangeRequesterListener listener =
      new MarkerRangeRequesterListener () {
        public void actionPerformed (final MarkerRangeRequesterEvent event) {
          final Range range = event.getRawRange ();

          if (range == null) {
            return;
          }

          final int start_position = range.getStart ();
          final int end_position = range.getEnd ();

          final int start_pos_in_feature = (start_position - 1) * 3 + 1;
          final int end_pos_in_feature = (end_position - 1) * 3 + 3;

          try {
            final Marker start_marker =
              selected_feature.getPositionInSequence (start_pos_in_feature);

            final Marker end_marker =
              selected_feature.getPositionInSequence (end_pos_in_feature);

            final MarkerRange marker_range =
              new MarkerRange (start_marker.getStrand (),
                               start_marker.getPosition (),
                               end_marker.getPosition ());

            getSelection ().setMarkerRange (marker_range);

            makeBaseVisible (getSelection ().getStartBaseOfSelection ());
          } catch (OutOfRangeException e) {
            new MessageDialog (getParentFrame (),
                               "amino acid range is out of range " +
                               "for this feature");
          }
        }
      };

    range_requester.addMarkerRangeRequesterListener (listener);

    range_requester.setVisible (true);
  }

  /**
   *  If there are some selected bases, select the feature in that range
   *  (replacing the orginal selection).  If some features are selected,
   *  select those features that overlap the selected features (and unselect
   *  the original features).
   **/
  private void selectOverlappingFeatures () {
    final MarkerRange selected_marker_range =
      getSelection ().getMarkerRange ();

    if (selected_marker_range == null) {
      final FeatureVector selected_features =
        getSelection ().getAllFeatures ();

      if (selected_features.size () == 0) {
        new MessageDialog (getParentFrame (), "nothing selected");
        return;
      } else {
        final FeatureVector new_selection = new FeatureVector ();

        for (int i = 0 ; i < selected_features.size () ; ++i) {
          final Feature this_feature = selected_features.elementAt (i);

          final Range this_feature_raw_range = this_feature.getMaxRawRange ();

          try {
            final FeatureVector overlapping_features =
              getEntryGroup ().getFeaturesInRange (this_feature_raw_range);

            for (int overlapping_feature_index = 0 ;
                 overlapping_feature_index < overlapping_features.size () ;
                 ++overlapping_feature_index) {
              final Feature overlapping_feature =
                overlapping_features.elementAt (overlapping_feature_index);

              if (!new_selection.contains (overlapping_feature)) {
                new_selection.add (overlapping_feature);
              }
            }
          } catch (OutOfRangeException e) {
            throw new Error ("internal error - unexpected exception: " + e);
          }
        }

        for (int i = 0 ; i < selected_features.size () ; ++i) {
          new_selection.remove (selected_features.elementAt (i));
        }

        getSelection ().set (new_selection);
      }
    } else {
      try {
        final Range raw_range = selected_marker_range.getRawRange ();
        getSelection ().set (getEntryGroup ().getFeaturesInRange (raw_range));
      } catch (OutOfRangeException e) {
        throw new Error ("internal error - unexpected exception: " + e);
      }
    }
  }

  /**
   *  If there are some selected bases or features select the features that are
   *  fully contained within that selection.
   **/
  private void selectContainedFeatures () {
    final MarkerRange selected_marker_range =
      getSelection ().getMarkerRange ();

    final FeatureVector featuresContainedByRange;
    try
    {
      if (selected_marker_range == null) {
        final FeatureVector selected_features = getSelection ().getAllFeatures ();
        if (selected_features.size () == 0) {
          new MessageDialog (getParentFrame (), "nothing selected");
          return;
        }

        featuresContainedByRange = new FeatureVector();
        for(int i=0; i<selected_features.size(); i++) {
          Range r = selected_features.elementAt(i).getMaxRawRange();
          final FeatureVector featuresContained = getFeaturesContainedByRange(r);
          for(int j=0; j<featuresContained.size(); j++) {
            final Feature f = featuresContained.elementAt(j);
            if(!featuresContainedByRange.contains(f) && !selected_features.contains(f))
              featuresContainedByRange.add(f);
          }
        }
      } else {
        featuresContainedByRange = getFeaturesContainedByRange(
            selected_marker_range.getRawRange());
      }
    } catch (OutOfRangeException e) {
      throw new Error ("internal error - unexpected exception: " + e);
    }
    getSelection ().set (featuresContainedByRange);
  }

  /**
   * Return the features that are contained by a given range
   * @param r
   * @return
   * @throws OutOfRangeException
   */
  private FeatureVector getFeaturesContainedByRange(Range r) throws OutOfRangeException
  {
    final FeatureVector featuresInRange =
        getEntryGroup ().getFeaturesInRange (r);
    final FeatureVector featuresContainedByRange = new FeatureVector();
    for(int i=0; i<featuresInRange.size(); i++)
    {
      Feature f = featuresInRange.elementAt(i);
      Range maxR =  f.getMaxRawRange(); 
      if(r.getStart() <= maxR.getStart() &&
         r.getEnd()   >= maxR.getEnd())
        featuresContainedByRange.add(f);
    }
    return featuresContainedByRange;
  }

  /**
   *  This method sends an GotoEvent to all the GotoEvent listeners that will
   *  make the given base visible.
   **/
  private void makeBaseVisible (final Marker base_marker) {
    goto_event_source.gotoBase (base_marker);
  }

  /**
   *  Return the EntryGroup that was passed to the constructor.
   **/
  private EntryGroup getEntryGroup () {
    return entry_group;
  }

}
