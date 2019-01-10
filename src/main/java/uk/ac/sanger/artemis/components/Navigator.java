/* Navigator.java
 *
 * created: Sun Jan 10 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/Navigator.java,v 1.4 2008-11-13 12:02:36 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.io.RangeVector;
import uk.ac.sanger.artemis.sequence.*;
import uk.ac.sanger.artemis.*;

import uk.ac.sanger.artemis.util.StringVector;

import java.awt.*;
import java.awt.event.*;

import javax.swing.*;

/**
 *  This component allows the user to navigate around the Entry.
 *
 *  @author Kim Rutherford
 *  @version $Id: Navigator.java,v 1.4 2008-11-13 12:02:36 tjc Exp $
 **/

public class Navigator extends JFrame
    implements EntryGroupChangeListener 
{
  private static final long serialVersionUID = 1L;

  /** selects the goto base function. */
  private final JRadioButton goto_base_button;

  /** selects the goto base pattern function. */
  private final JRadioButton goto_base_pattern_button =
      new JRadioButton ("Find Base Pattern:", true);;

  /** selects the find amino acid sequence function. */
  private final JRadioButton goto_aa_pattern_button =
      new JRadioButton ("Find Amino Acid String:", true);;

  /** selects the goto feature qualifier value function. */
  private final JRadioButton goto_qualifier_button;

  /** selects the goto gene name function. */
  private final JRadioButton goto_gene_name_button;

  /** selects the goto feature key function. */
  private final JRadioButton goto_key_button;

  /** pattern to search for if the user has selected the
      goto base function. */
  private final JTextField goto_base_text;

  /** pattern to search for if the user has selected the
      goto base pattern function. */
  private final JTextField goto_base_pattern_text;

  /** pattern to search for if the user has selected the
      goto amino acid function. */
  private final JTextField goto_aa_pattern_text;

  /** pattern to search for if the user has selected the
      goto qualifier value function. */
  private final JTextField goto_qualifier_textfield;

  /** pattern to search for if the user has selected the
      goto gene name function. */
  private final JTextField goto_gene_name_textfield;

  /** key to search for if the user has selected the
      goto key function. */
  private final JTextField goto_feature_key_textfield;

  /** selects if the search should start at first/last
      base or first/last feature (depending on the search type). */
  private JRadioButton start_at_an_end_button;

  /** selects if the search should start at the
      position of the current selection. */
  private final JRadioButton start_at_selection_button;

  /** If checked the search will go backwards. */
  private final JCheckBox search_backward_button;

  /** If checked the search will ignore the case of the query and subject. */
  private final JCheckBox ignore_case_button;

  /** If checked the search text is allowed to match a substring of a
      qualifier value. */
  private final JCheckBox partial_match_button;
  
  /** Check whether the match overlaps selected range */
  private final JCheckBox overlaps_with_selection;
  
  /** Search fwd strand */
  private final JCheckBox fwd_strand = new JCheckBox("Forward Strand ", true);;
  
  /** Search bwd strand */
  private final JCheckBox rev_strand = new JCheckBox("Reverse Strand", true);;

  /**
   *  The GotoEventSource object that was passed to the constructor.
   **/
  private final GotoEventSource goto_event_source;

  /**
   *  The EntryGroup object that was passed to the constructor.
   **/
  private final EntryGroup entry_group;

  /**
   *  This is the Selection that was passed to the constructor.
   **/
  private final Selection selection;
  
  private boolean disposeOnClose = false;
  
  /**
   *  Create a new Navigator component.
   *  @param selection The Selection that the commands will operate on.
   *  @param goto_event_source The object the we will call gotoBase () on.
   *  @param entry_group The EntryGroup object used when searching for
   *    qualifier text in features.
   **/
  public Navigator (final Selection selection,
                    final GotoEventSource goto_event_source,
                    final EntryGroup entry_group) 
  {
    super ("Artemis Navigator");

    this.selection = selection;
    this.entry_group = entry_group;
    this.goto_event_source = goto_event_source;

    final Font default_font = Options.getOptions ().getFont ();

    GridBagLayout gridbag = new GridBagLayout();
    getContentPane ().setLayout (gridbag);

    setFont (default_font);

    GridBagConstraints c = new GridBagConstraints();

    c.fill = GridBagConstraints.HORIZONTAL;
    c.anchor = GridBagConstraints.NORTH;
    c.weighty = 0;

    final int TEXT_FIELD_WIDTH = 25;

    final ButtonGroup button_group = new ButtonGroup ();

    goto_base_button = new JRadioButton ("Goto Base:", true);

    button_group.add (goto_base_button);

    final JPanel goto_base_panel = new JPanel (new FlowLayout (FlowLayout.LEFT));
    goto_base_panel.add (goto_base_button);
    c.gridwidth = 2;
    gridbag.setConstraints (goto_base_panel, c);
    getContentPane ().add (goto_base_panel);

    goto_base_text = new JTextField ("", TEXT_FIELD_WIDTH);
    c.gridwidth = GridBagConstraints.REMAINDER;
    gridbag.setConstraints (goto_base_text, c);
    getContentPane ().add (goto_base_text);

    goto_base_text.addKeyListener (new KeyAdapter () {
      public void keyTyped(final KeyEvent e) {
        goto_base_button.setSelected (true);
        if (e.getKeyChar () == '\n') {
          doGoto ();
        }
      }
    });


    goto_gene_name_button =
      new JRadioButton ("Goto Feature With Gene Name:", true);

    button_group.add (goto_gene_name_button);

    final JPanel goto_gene_name_panel = new JPanel (new FlowLayout (FlowLayout.LEFT));
    goto_gene_name_panel.add (goto_gene_name_button);
    c.gridwidth = 2;
    gridbag.setConstraints (goto_gene_name_panel, c);
    getContentPane ().add (goto_gene_name_panel);

    goto_gene_name_textfield = new JTextField ("", TEXT_FIELD_WIDTH);
    c.gridwidth = GridBagConstraints.REMAINDER;
    gridbag.setConstraints (goto_gene_name_textfield, c);
    getContentPane ().add (goto_gene_name_textfield);

    goto_gene_name_textfield.addKeyListener (new KeyAdapter () {
      public void keyTyped(final KeyEvent e) {
        goto_gene_name_button.setSelected (true);
        if (e.getKeyChar () == '\n') {
          doGoto ();
        }
      }
    });


    goto_qualifier_button =
      new JRadioButton ("Goto Feature With This Qualifier Value:", true);
    button_group.add (goto_qualifier_button);

    final JPanel goto_qualifier_panel = new JPanel (new FlowLayout (FlowLayout.LEFT));
    goto_qualifier_panel.add (goto_qualifier_button);
    c.gridwidth = 2;
    gridbag.setConstraints (goto_qualifier_panel, c);
    getContentPane ().add (goto_qualifier_panel);

    goto_qualifier_textfield = new JTextField ("", TEXT_FIELD_WIDTH);
    c.gridwidth = GridBagConstraints.REMAINDER;
    gridbag.setConstraints (goto_qualifier_textfield, c);
    getContentPane ().add (goto_qualifier_textfield);

    goto_qualifier_textfield.addKeyListener (new KeyAdapter () {
      public void keyTyped(final KeyEvent e) {
        goto_qualifier_button.setSelected (true);
        if (e.getKeyChar () == '\n') {
          doGoto ();
        }
      }
    });


    goto_key_button =
      new JRadioButton ("Goto Feature With This Key:", true);
    button_group.add (goto_key_button);

    final JPanel goto_key_panel = new JPanel (new FlowLayout (FlowLayout.LEFT));
    goto_key_panel.add (goto_key_button);
    c.gridwidth = 2;
    gridbag.setConstraints (goto_key_panel, c);
    getContentPane ().add (goto_key_panel);

    goto_feature_key_textfield = new JTextField ("", TEXT_FIELD_WIDTH);
    c.gridwidth = GridBagConstraints.REMAINDER;
    gridbag.setConstraints (goto_feature_key_textfield, c);
    getContentPane ().add (goto_feature_key_textfield);

    goto_feature_key_textfield.addKeyListener (new KeyAdapter () {
      public void keyTyped(final KeyEvent e) {
        goto_key_button.setSelected (true);
        if (e.getKeyChar () == '\n') {
          doGoto ();
        }
      }
    });


    goto_base_pattern_button.addItemListener(new ItemListener() {
      public void itemStateChanged(ItemEvent arg0) {
        rev_strand.setEnabled(goto_aa_pattern_button.isSelected() || 
            goto_base_pattern_button.isSelected());
        fwd_strand.setEnabled(goto_aa_pattern_button.isSelected() || 
            goto_base_pattern_button.isSelected());
      }
    });
    button_group.add (goto_base_pattern_button);

    final JPanel goto_base_pattern_panel = new JPanel (new FlowLayout (FlowLayout.LEFT));
    goto_base_pattern_panel.add (goto_base_pattern_button);
    c.gridwidth = 2;
    gridbag.setConstraints (goto_base_pattern_panel, c);
    getContentPane ().add (goto_base_pattern_panel);

    goto_base_pattern_text = new JTextField ("", TEXT_FIELD_WIDTH);
    c.gridwidth = GridBagConstraints.REMAINDER;
    gridbag.setConstraints (goto_base_pattern_text, c);
    getContentPane ().add (goto_base_pattern_text);

    goto_base_pattern_text.addKeyListener (new KeyAdapter () {
      public void keyTyped(final KeyEvent e) {
        goto_base_pattern_button.setSelected (true);
        if (e.getKeyChar () == '\n') {
          doGoto ();
        }
      }
    });


    goto_aa_pattern_button.addItemListener(new ItemListener() {
      public void itemStateChanged(ItemEvent arg0) {
        rev_strand.setEnabled(goto_aa_pattern_button.isSelected() || 
            goto_base_pattern_button.isSelected());
        fwd_strand.setEnabled(goto_aa_pattern_button.isSelected() || 
            goto_base_pattern_button.isSelected());
      }
    });
    button_group.add (goto_aa_pattern_button);

    final JPanel goto_aa_pattern_panel = new JPanel (new FlowLayout (FlowLayout.LEFT));
    goto_aa_pattern_panel.add (goto_aa_pattern_button);
    c.gridwidth = 2;
    gridbag.setConstraints (goto_aa_pattern_panel, c);
    getContentPane ().add (goto_aa_pattern_panel);

    goto_aa_pattern_text = new JTextField ("", TEXT_FIELD_WIDTH);
    c.gridwidth = GridBagConstraints.REMAINDER;
    gridbag.setConstraints (goto_aa_pattern_text, c);
    getContentPane ().add (goto_aa_pattern_text);

    goto_aa_pattern_text.addKeyListener (new KeyAdapter () {
      public void keyTyped(final KeyEvent e) {
        goto_aa_pattern_button.setSelected (true);
        if (e.getKeyChar () == '\n') {
          doGoto ();
        }
      }
    });

    goto_base_button.setSelected (true);

    final ButtonGroup start_position_button_group = new ButtonGroup ();
    final JPanel start_at_an_end_panel = new JPanel (new FlowLayout(FlowLayout.LEFT));

    start_at_an_end_button =
      new JRadioButton ("beginning (or end)", true);

    start_position_button_group.add (start_at_an_end_button);

    c.gridwidth = GridBagConstraints.REMAINDER;

    start_at_an_end_panel.add (new JLabel(" Start search at:"));
    start_at_an_end_panel.add (start_at_an_end_button);

    start_at_selection_button =
      new JRadioButton ("selection", false);
    start_at_selection_button.addItemListener(new ItemListener() {
      public void itemStateChanged(ItemEvent e) {
        overlaps_with_selection.setEnabled(start_at_selection_button.isSelected());
      }
    });

    start_position_button_group.add (start_at_selection_button);

    start_at_an_end_panel.add (start_at_selection_button);

    c.gridwidth = GridBagConstraints.REMAINDER;
    getContentPane ().add (Box.createVerticalStrut(5), c);
    getContentPane ().add (start_at_an_end_panel, c);


    final JPanel select_within_button_panel = new JPanel (new FlowLayout (FlowLayout.LEFT));
    overlaps_with_selection = new JCheckBox("Overlaps With Selection", false);
    overlaps_with_selection.setEnabled(start_at_selection_button.isSelected());
    select_within_button_panel.add(overlaps_with_selection);
    getContentPane ().add (select_within_button_panel, c);

    //
    // STRAND SELECTION
    final JPanel strand_panel = new JPanel (new FlowLayout (FlowLayout.LEFT));
    fwd_strand.setEnabled(goto_aa_pattern_button.isSelected() || 
        goto_base_pattern_button.isSelected());
    strand_panel.add(fwd_strand);

    rev_strand.setEnabled(goto_aa_pattern_button.isSelected() || 
        goto_base_pattern_button.isSelected());
    strand_panel.add(rev_strand);
    getContentPane ().add (strand_panel, c);


    //
    // OTHER OPTIONS
    final JPanel option_button_panel = new JPanel (new FlowLayout (FlowLayout.LEFT));
    search_backward_button = new JCheckBox ("Search Backward", false);
    ignore_case_button = new JCheckBox ("Ignore Case", true);
    partial_match_button = new JCheckBox ("Allow Substring Matches", true);

    option_button_panel.add (search_backward_button);
    option_button_panel.add (ignore_case_button);
    option_button_panel.add (partial_match_button);

    getContentPane ().add (option_button_panel, c);


    final JButton goto_button = new JButton ("Goto");
    goto_button.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent e) {
        doGoto ();
      }
    });


    final JButton clear_button = new JButton ("Clear");
    clear_button.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent e) {
        clear ();
      }
    });


    final JButton close_button = new JButton ("Close");
    close_button.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent e) {
        onClose();
      }
    });


    final JPanel close_and_goto_panel = new JPanel (
        new FlowLayout (FlowLayout.CENTER, 15, 5));

    close_and_goto_panel.add (goto_button);
    close_and_goto_panel.add (clear_button);
    close_and_goto_panel.add (close_button);

    gridbag.setConstraints (close_and_goto_panel, c);
    getContentPane ().add (close_and_goto_panel);


    addWindowListener (new WindowAdapter () {
      public void windowClosing (WindowEvent event) {
        onClose();
      }
    });

    getEntryGroup ().addEntryGroupChangeListener (this);

    pack ();
    Utilities.centreFrame(this);
    setVisible (true);
  }

  /**
   *  Implementation of the EntryGroupChangeListener interface.  We listen to
   *  EntryGroupChange events so that we can get rid of the Navigator when the
   *  EntryGroup is no longer in use (for example when the EntryEdit is
   *  closed).
   **/
  public void entryGroupChanged (final EntryGroupChangeEvent event) {
    switch (event.getType ()) {
    case EntryGroupChangeEvent.DONE_GONE:
      getEntryGroup ().removeEntryGroupChangeListener (this);
      dispose ();
      break;
    }
  }
  
  private void onClose()
  {
    if(disposeOnClose)
    {
      getEntryGroup().removeEntryGroupChangeListener(Navigator.this);
      dispose();
    }
    else
      setVisible(false);
  }
  
  protected void setDisposeOnClose(boolean disposeOnClose)
  {
    this.disposeOnClose = disposeOnClose;
  }

  /**
   *  This method finds the given pattern in the Bases of the given EntryGroup
   *  and returns a MarkerRange for it.
   *  @param pattern This is the pattern to search for.
   *  @param entry_group The EntryGroup to search.
   *  @param selection The current selection.
   *  @param start_at_end If true or if there is nothing in the Selection then
   *    the search will start at the first or last base (depending on value of
   *    the next parameter), otherwise the search will start at the Selection.
   *  @param search_backwards If true the search will move from last base to
   *    first base, otherwise first to last.
   *  @return The range that matches the pattern or null if there is no match.
   **/
  private static MarkerRange findBasePattern (final BasePattern pattern,
                                             final EntryGroup entry_group,
                                             final Selection selection,
                                             final boolean start_at_end,
                                             final boolean search_backwards,
                                             final boolean search_fwd_strand,
                                             final boolean search_bwd_strand) {
    // if start_at_end is false we want to start the search at the selection
    final Marker selection_base = selection.getLowestBaseOfSelection ();

    final Marker start_position;

    if (start_at_end || selection_base == null) {
      // null means start the search from the beginning
      start_position = null;
    } else {
      start_position = selection_base;
    }

    final MarkerRange match_range =
      pattern.findMatch (entry_group.getBases (),
                         start_position,
                         entry_group.getSequenceLength (),
                         search_backwards,
                         search_fwd_strand, search_bwd_strand);

    if (match_range == null) {
      return null;
    } else {
      return match_range;
    }
  }

  /**
   *  This method finds the given amino acid seqeunce, sets the selection to
   *  the matching bases and then sends a GotoEvent to all the GotoEvent
   *  listeners that will make the first base of the match visible.
   *  @param sequence This is the pattern to search for.
   *  @param entry_group The EntryGroup to search.
   *  @param selection The current selection.
   *  @param start_at_end If true or if there is nothing in the Selection then
   *    the search will start at first or last base (depending on value of the
   *    next parameter), otherwise the search will start at the Selection.
   *  @param search_backwards If true the search will move from last base to
   *    first base, otherwise first to last.
   *  @return The range that matches the pattern or null if there is no match.
   **/
  private static MarkerRange findAminoAcidSequence (final AminoAcidSequence sequence,
                                                   final EntryGroup entry_group,
                                                   final Selection selection,
                                                   final boolean start_at_end,
                                                   final boolean search_backwards,
                                                   final boolean search_fwd_strand,
                                                   final boolean search_bwd_strand) {
    // if start_at_end is false we want to start the search at the selection
    final Marker selection_base = selection.getLowestBaseOfSelection ();

    final Marker start_position;

    if (start_at_end || selection_base == null) {
      // null means start the search from the beginning
      start_position = null;
    } else {
      start_position = selection_base;
    }

    return sequence.findMatch (entry_group.getBases (),
                               start_position, search_backwards, 
                               search_fwd_strand, search_bwd_strand);
  }

  /**
   *  This method will perform the selected goto function.
   **/
  private void doGoto () {
    if (goto_base_button.isSelected ()) {
        doGotoBase ();
      return;
    }

    if (goto_base_pattern_button.isSelected ()) {
      if(!rev_strand.isSelected() && !fwd_strand.isSelected())
        new MessageDialog (this, "Select a strand to search.");
      else
        doGotoBasePattern ();
      return;
    }

    if (goto_aa_pattern_button.isSelected ()) {
      if(!rev_strand.isSelected() && !fwd_strand.isSelected())
        new MessageDialog (this, "Select a strand to search.");
      else
        doGotoAAPattern ();
      return;
    }

    if (goto_gene_name_button.isSelected ()) {
      final StringVector qualifiers_to_search = new StringVector ();
      qualifiers_to_search.add (Options.getOptions ().getAllGeneNames ());

      doGotoQualifierValue (qualifiers_to_search,
                            goto_gene_name_textfield.getText ().trim ());
      return;
    }

    if (goto_qualifier_button.isSelected ()) {
      doGotoQualifierValue (null, goto_qualifier_textfield.getText ().trim ());
      return;
    }

    if (goto_key_button.isSelected ()) {
      doGotoKey ();
      return;
    }
  }

  /**
   *  Clear all the JTextField components.
   **/
  private void clear () {
    goto_base_text.setText ("");
    goto_base_pattern_text.setText ("");
    goto_aa_pattern_text.setText ("");
    goto_qualifier_textfield.setText ("");
    goto_gene_name_textfield.setText ("");
    goto_feature_key_textfield.setText ("");
  }

  /**
   *  Show a message dialog saying that the search text field is empty.
   **/
  private void searchTextEmptyError () {
    new MessageDialog (this, "You have not entered a value to search for");
  }

  /**
   *  Go to the base that the user typed in to the goto_base_text TextArea.
   **/
  private void doGotoBase () {
    final String number_string = goto_base_text.getText ().trim ();

    if (number_string.length () == 0) {
      new MessageDialog (this, "you have not entered a number to go to");
      return;
    }

    try {
      final int destination_base =
        Integer.valueOf (number_string).intValue ();

      final MarkerRange destination_range =
        goto_event_source.gotoBase (destination_base);

      if (destination_range == null) {
        new MessageDialog (this,
                           "the base is out of range for this sequence");
      } else {
        // success select that base
        getSelection ().setMarkerRange (destination_range);
      }
    } catch (NumberFormatException e) {
      new MessageDialog (this, "Cannot understand this number: " +
                         goto_base_text.getText ());
    }
  }

  /**
   *  Go to the base pattern that the user typed in to the
   *  goto_base_pattern_text TextArea.
   **/
  private void doGotoBasePattern () {
    try {
      final String pattern_string =
        goto_base_pattern_text.getText ().trim ();

      if (pattern_string.length () == 0) {
        new MessageDialog (this, "you have not entered a pattern to go to");
        return;
      }

      final BasePattern pattern = new BasePattern (pattern_string);

      final boolean start_at_an_end = start_at_an_end_button.isSelected ();

      final MarkerRange match_range =
        findBasePattern (pattern, getEntryGroup (),
                         getSelection (), start_at_an_end,
                         search_backward_button.isSelected (),
                         fwd_strand.isSelected(), rev_strand.isSelected());

      if (match_range == null) {
        new MessageDialog (this, "reached the end of sequence");
      } else {
        
        if(start_at_selection_button.isSelected() && overlaps_with_selection.isSelected()) {
          if(!overlapsRange(getSelection().getSelectionRanges(), match_range)) {
            new MessageDialog (this, "Next match ("+ 
               match_range.getRawRange().toString()+") outside selection.");
            return;
          }
        }
        
        getSelection ().setMarkerRange (match_range);

        final Marker first_selected_base =
          getSelection ().getLowestBaseOfSelection ();

        goto_event_source.gotoBase (first_selected_base);
        start_at_selection_button.setSelected (true);
      }
    } catch (BasePatternFormatException e) {
      new MessageDialog (this,
                         "Illegal base pattern: " +
                         goto_base_pattern_text.getText ());
    }
  }

  /**
   *  Find the amino acid pattern that the user typed in to the
   *  goto_aa_pattern_text TextArea.
   **/
  private void doGotoAAPattern () {
    final String pattern_string =
      goto_aa_pattern_text.getText ().trim ().toUpperCase();

    if (pattern_string.length () == 0) {
      new MessageDialog (this, "you have not entered a pattern to go to");
      return;
    }

    final AminoAcidSequence pattern = new AminoAcidSequence (pattern_string);

    final boolean start_at_an_end = start_at_an_end_button.isSelected ();
    final boolean search_backwards = search_backward_button.isSelected ();

    final MarkerRange match_range =
      findAminoAcidSequence (pattern, getEntryGroup (),
                             getSelection (), start_at_an_end,
                             search_backwards,
                             fwd_strand.isSelected(), rev_strand.isSelected());

    if (match_range == null)
      new MessageDialog (this, "reached the end of sequence");
    else {
      
      if(start_at_selection_button.isSelected() && overlaps_with_selection.isSelected()) {
        if(!overlapsRange(getSelection().getSelectionRanges(), match_range)) {
          new MessageDialog (this, "Next match ("+ 
             match_range.getRawRange().toString()+") outside selection.");
          return;
        }
      }
      
      start_at_selection_button.setSelected (true);
      getSelection ().setMarkerRange (match_range);
      goto_event_source.gotoBase (getSelection ().getLowestBaseOfSelection ());
    }
  }

  /**
   * Check if a MarkerRange is contained in a RangeVector
   * @param ranges
   * @param match_range
   * @return
   */
  private boolean overlapsRange(final RangeVector ranges,
                                final MarkerRange match_range)
  {
    for(int i=0; i<ranges.size(); i++) {
      Range r = (Range)ranges.get(i);

      if( ( match_range.getRawStart().getRawPosition() >=  r.getStart() && 
            match_range.getRawStart().getRawPosition() <= r.getEnd() ) ||
          ( match_range.getRawEnd().getRawPosition() >=  r.getStart() && 
            match_range.getRawEnd().getRawPosition() <= r.getEnd() ))
        return true;
    }  
    return false;
  }
  
  /**
   *  Select the next feature that contains the text given by the user in the
   *  goto_qualifier_textfield TextArea.
   **/
  private void doGotoQualifierValue (final StringVector qualifiers_to_search,
                                     final String search_text) {
    if (search_text.equals ("")) {
      searchTextEmptyError ();
      return;
    }

    final FeatureVector selected_features =
      getSelection ().getAllFeatures ();

    final int index;

    if (selected_features.size () == 0 ||
        start_at_an_end_button.isSelected ()) {
      index = -1;
    } else {
      index = getEntryGroup ().indexOf (selected_features.elementAt (0));
    }

    final int first_search_feature_index;

    Feature found_feature = null;

    if (search_backward_button.isSelected ()) {
      if (index == -1) {
        // nothing was selected so start the search at the first feature
        first_search_feature_index =
          getEntryGroup ().getAllFeaturesCount () - 1;
      } else {
        first_search_feature_index = index - 1;
      }

      for (int i = first_search_feature_index ; i >= 0 ; --i) {
        final Feature this_feature = getEntryGroup ().featureAt (i);

        if (this_feature.containsText (search_text,
                                       ignore_case_button.isSelected (),
                                       partial_match_button.isSelected (),
                                       qualifiers_to_search)) {
          found_feature = this_feature;
          break;
        }
      }
    } else {
      if (index == -1) {
        // nothing was selected so start the search at the first feature
        first_search_feature_index = 0;
      } else {
        first_search_feature_index = index + 1;
      }

      for (int i = first_search_feature_index ;
           i < getEntryGroup ().getAllFeaturesCount () ;
           ++i) {
        final Feature this_feature = getEntryGroup ().featureAt (i);

        if (this_feature.containsText (search_text,
                                       ignore_case_button.isSelected (),
                                       partial_match_button.isSelected (),
                                       qualifiers_to_search)) {

          found_feature = this_feature;
          break;
        }
      }
    }

    if (found_feature == null) {
      getSelection ().clear ();
      start_at_an_end_button.setSelected (true);

      new MessageDialog (this, "text not found");
    } else {
      getSelection ().set (found_feature);
      goto_event_source.gotoBase (getSelection ().getLowestBaseOfSelection ());
      start_at_selection_button.setSelected (true);
    }
  }

  /**
   *  Select the next feature that has the key given by the user in the
   *  goto_feature_key_textfield TextArea.
   **/
  private void doGotoKey () {
    final FeatureVector selected_features =
      getSelection ().getAllFeatures ();

    final int index;

    if (selected_features.size () == 0 ||
        start_at_an_end_button.isSelected ()) {
      index = -1;
    } else {
      index = getEntryGroup ().indexOf (selected_features.elementAt (0));
    }

    final int first_search_feature_index;

    final String search_key_string =
      goto_feature_key_textfield.getText ().trim ();

    if (search_key_string.equals ("")) {
      searchTextEmptyError ();
      return;
    }

    Feature found_feature = null;

    if (search_backward_button.isSelected ()) {
      if (index == -1) {
        // nothing was selected so start the search at the first feature
        first_search_feature_index =
          getEntryGroup ().getAllFeaturesCount () - 1;
      } else {
        first_search_feature_index = index - 1;
      }

      for (int i = first_search_feature_index ; i >= 0 ; --i) {
        final Feature this_feature = getEntryGroup ().featureAt (i);

        if (keyMatches (this_feature, search_key_string)) {
          found_feature = this_feature;
          break;
        }
      }
    } else {
      if (index == -1) {
        // nothing was selected so start the search at the first feature
        first_search_feature_index = 0;
      } else {
        first_search_feature_index = index + 1;
      }

      for (int i = first_search_feature_index ;
           i < getEntryGroup ().getAllFeaturesCount () ;
           ++i) {
        final Feature this_feature = getEntryGroup ().featureAt (i);

        if (keyMatches (this_feature, search_key_string)) {
          found_feature = this_feature;
          break;
        }
      }
    }

    if (found_feature == null) {
      getSelection ().clear ();
      start_at_an_end_button.setSelected (true);

      new MessageDialog (this, "key not found");
    } else {
      getSelection ().set (found_feature);
      goto_event_source.gotoBase (getSelection ().getLowestBaseOfSelection ());
      start_at_selection_button.setSelected (true);
    }
  }

  /**
   *  Returns true if and only if the given feature has search_key_string as
   *  it's Key.
   **/
  private boolean keyMatches (final Feature test_feature,
                             final String search_key_string) {
    final String feature_key_string;

    if (ignore_case_button.isSelected ()) {
      feature_key_string = test_feature.getKey ().toString ().toLowerCase ();
    } else {
      feature_key_string = test_feature.getKey ().toString ();
    }

    if (feature_key_string.equals (search_key_string.toLowerCase ())) {
      return true;
    } else {
      return false;
    }
  }

  /**
   *  Return the Selection object that was passed to the constructor.
   **/
  private Selection getSelection () {
    return selection;
  }

  /**
   *  Return the EntryGroup object that was passed to the constructor.
   **/
  private EntryGroup getEntryGroup () {
    return entry_group;
  }

}
