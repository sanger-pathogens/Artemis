/* Selector.java
 *
 * created: Tue Apr 11 2000
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2000,2001,2002  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/Selector.java,v 1.11 2008-01-17 16:00:31 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.sequence.AminoAcidSequence;
import uk.ac.sanger.artemis.sequence.Strand;
import uk.ac.sanger.artemis.sequence.Bases;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.GFFDocumentEntry;
import uk.ac.sanger.artemis.io.Location;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.io.RangeVector;

import uk.ac.sanger.artemis.util.StringVector;

import java.awt.*;
import java.awt.event.*;
import java.util.Collections;
import java.util.Comparator;
import java.util.StringTokenizer;

import javax.swing.*;

/**
 *  This component allows the user to set the selection by filtering the
 *  features in an EntryGroup on key and contents.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: Selector.java,v 1.11 2008-01-17 16:00:31 tjc Exp $
 **/

public class Selector extends JFrame
    implements EntryGroupChangeListener {
  
  private static final long serialVersionUID = 1L;
  private JCheckBox by_key_button;
  private JCheckBox by_qualifier_button;
  private JCheckBox by_motif_button;
  private JCheckBox less_than_bases_button;
  private JCheckBox greater_than_bases_button;
  private JCheckBox less_than_exons_button;
  private JCheckBox greater_than_exons_button;
  private JCheckBox ignore_case_checkbox;
  private JCheckBox match_any_word_checkbox;
  private JCheckBox forward_strand_checkbox;
  private JCheckBox reverse_strand_checkbox;
  private JCheckBox intron_pattern_button;
  private KeyChoice key_selector;
  private QualifierChoice qualifier_selector;
  private JTextField qualifier_text;
  private JTextField motif_text;
  private JTextField less_than_bases_text;
  private JTextField greater_than_bases_text;
  private JTextField less_than_exons_text;
  private JTextField greater_than_exons_text;

  /**
   *  If checked the search text is allowed to match a substring of a
   *  qualifier value.
   **/
  final private JCheckBox partial_match_checkbox;

  /**
   *  The EntryGroup object that was passed to the constructor.
   **/
  final private EntryGroup entry_group;

  /**
   *  This is the Selection that was passed to the constructor.
   **/
  final private Selection selection;
  
  public static org.apache.log4j.Logger logger4j = 
    org.apache.log4j.Logger.getLogger(Selector.class);
  
  /**
   *  Create a new Selector that van set the given Selection.
   *  @param selection The Selection that the commands in the menu will
   *    operate on.
   *  @param entry_group The component will choose features from this
   *    EntryGroup.
   *  @param goto_event_source The object the we will call gotoBase() on.
   *  @param base_plot_group The BasePlotGroup associated with this JMenu -
   *    needed to call getCodonUsageAlgorithm()
   **/
  public Selector(final Selection selection,
                   final GotoEventSource goto_event_source,
                   final EntryGroup entry_group,
                   final BasePlotGroup base_plot_group) 
  {
    super("Artemis Feature Selector");

    this.selection = selection;
    this.entry_group = entry_group;

    final Font default_font = Options.getOptions().getFont();

    setFont(default_font);

    GridBagLayout gridbag = new GridBagLayout();
    getContentPane().setLayout(gridbag);

    GridBagConstraints c = new GridBagConstraints();

    c.fill = GridBagConstraints.HORIZONTAL;
    c.anchor = GridBagConstraints.WEST;
    c.weighty = 0;
    c.gridwidth = GridBagConstraints.REMAINDER;

    final JLabel top_label = new JLabel("Select by:");
    gridbag.setConstraints(top_label, c);
    getContentPane().add(top_label);


    by_key_button = new JCheckBox("Key: ", false);
    final JPanel by_key_panel = new JPanel();
    by_key_panel.setLayout(new FlowLayout(FlowLayout.LEFT));
    by_key_panel.add(by_key_button);
    c.gridwidth = 2;
    gridbag.setConstraints(by_key_panel, c);
    getContentPane().add(by_key_panel);

//  final EntryInformation default_entry_information =
//    Options.getArtemisEntryInformation();

    Entry default_entry = getEntryGroup().getDefaultEntry();
    if(default_entry == null)
      default_entry = getEntryGroup().elementAt(0);
  
    final EntryInformation default_entry_information =
                        default_entry.getEntryInformation();

    key_selector = new KeyChoice(default_entry_information);

    c.gridwidth = GridBagConstraints.REMAINDER;
    gridbag.setConstraints(key_selector, c);
    getContentPane().add(key_selector);

    key_selector.addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent _) 
      {
        by_key_button.setSelected(true);
      }
    });

    by_key_button.addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent e) 
      {
        if(!by_key_button.isSelected()) 
          by_qualifier_button.setSelected(false);
      }
    });


    by_qualifier_button = new JCheckBox("Qualifier: ", false);
    final JPanel by_qualifier_panel = new JPanel();
    by_qualifier_panel.setLayout(new FlowLayout(FlowLayout.LEFT));
    by_qualifier_panel.add(by_qualifier_button);
    c.gridwidth = 2;
    gridbag.setConstraints(by_qualifier_panel, c);
    getContentPane().add(by_qualifier_panel);

    boolean isGFF = false;
    if(getEntryGroup().getDefaultEntry() != null &&
       getEntryGroup().getDefaultEntry().getEMBLEntry() instanceof GFFDocumentEntry)
      isGFF = true;
    qualifier_selector = new QualifierChoice(default_entry_information,
                                              key_selector.getSelectedItem(),
                                              null, isGFF);
    
    c.gridwidth = GridBagConstraints.REMAINDER;
    gridbag.setConstraints(qualifier_selector, c);
    getContentPane().add(qualifier_selector);

    qualifier_selector.addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent _) 
      {
        by_qualifier_button.setSelected(true);
        by_key_button.setSelected(true);
      }
    });

    key_selector.addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent _) 
      {
        qualifier_selector.setKey(key_selector.getSelectedItem());
      }
    });

    by_qualifier_button.addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent e) 
      {
        if(by_qualifier_button.isSelected())
        {
          if(!by_key_button.isSelected()) 
            by_key_button.setSelected(true);
        }
      }
    });


    final JLabel qualifier_text_label =
      new JLabel("   Containing this text: ");
    c.gridwidth = 2;
    gridbag.setConstraints(qualifier_text_label, c);
    getContentPane().add(qualifier_text_label);

    qualifier_text = new JTextField("", 18);

    c.gridwidth = GridBagConstraints.REMAINDER;
    gridbag.setConstraints(qualifier_text, c);
    getContentPane().add(qualifier_text);


    ignore_case_checkbox = new JCheckBox("Ignore Case", true);

    final JPanel ignore_case_panel = new JPanel();

    ignore_case_panel.setLayout(new FlowLayout(FlowLayout.LEFT));
    ignore_case_panel.add(ignore_case_checkbox);

    c.gridwidth = 2;
    gridbag.setConstraints(ignore_case_panel, c);
    getContentPane().add(ignore_case_panel);


    partial_match_checkbox = new JCheckBox("Allow Partial Match", true);

    final JPanel partial_match_panel = new JPanel();

    partial_match_panel.setLayout(new FlowLayout(FlowLayout.LEFT));
    partial_match_panel.add(partial_match_checkbox);

    c.gridwidth = GridBagConstraints.REMAINDER;
    gridbag.setConstraints(partial_match_panel, c);
    getContentPane().add(partial_match_panel);


    match_any_word_checkbox = new JCheckBox("Match Any Word", false);

    final JPanel match_any_word_panel = new JPanel();

    match_any_word_panel.setLayout(new FlowLayout(FlowLayout.LEFT));
    match_any_word_panel.add(match_any_word_checkbox);

    c.gridwidth = GridBagConstraints.REMAINDER;
    gridbag.setConstraints(match_any_word_panel, c);
    getContentPane().add(match_any_word_panel);

    andSeparator(gridbag, c);

    less_than_bases_button = new JCheckBox("Up to: ", false);
    final JPanel less_than_bases_panel = new JPanel();
    less_than_bases_panel.setLayout(new FlowLayout(FlowLayout.LEFT));
    less_than_bases_panel.add(less_than_bases_button);
    c.gridwidth = 3;
    gridbag.setConstraints(less_than_bases_panel, c);
    getContentPane().add(less_than_bases_panel);

    less_than_bases_text = new JTextField("", 18);

    gridbag.setConstraints(less_than_bases_text, c);
    getContentPane().add(less_than_bases_text);

    final JLabel less_than_bases_label = new JLabel(" bases long");

    c.gridwidth = GridBagConstraints.REMAINDER;
    gridbag.setConstraints(less_than_bases_label, c);
    getContentPane().add(less_than_bases_label);

    andSeparator(gridbag, c);

    greater_than_bases_button = new JCheckBox("At least: ", false);
    final JPanel greater_than_bases_panel = new JPanel();
    greater_than_bases_panel.setLayout(new FlowLayout(FlowLayout.LEFT));
    greater_than_bases_panel.add(greater_than_bases_button);
    c.gridwidth = 2;
    gridbag.setConstraints(greater_than_bases_panel, c);
    getContentPane().add(greater_than_bases_panel);

    greater_than_bases_text = new JTextField("", 18);

    gridbag.setConstraints(greater_than_bases_text, c);
    getContentPane().add(greater_than_bases_text);

    final JLabel greater_than_bases_label = new JLabel(" bases long");

    c.gridwidth = GridBagConstraints.REMAINDER;
    gridbag.setConstraints(greater_than_bases_label, c);
    getContentPane().add(greater_than_bases_label);

    andSeparator(gridbag, c);
    
    less_than_exons_button = new JCheckBox("Up to: ", false);
    final JPanel less_than_exons_panel = new JPanel();
    less_than_exons_panel.setLayout(new FlowLayout(FlowLayout.LEFT));
    less_than_exons_panel.add(less_than_exons_button);
    c.gridwidth = 3;
    gridbag.setConstraints(less_than_exons_panel, c);
    getContentPane().add(less_than_exons_panel);

    less_than_exons_text = new JTextField("", 18);

    gridbag.setConstraints(less_than_exons_text, c);
    getContentPane().add(less_than_exons_text);

    final JLabel less_than_exons_label = new JLabel(" exons long");

    c.gridwidth = GridBagConstraints.REMAINDER;
    gridbag.setConstraints(less_than_exons_label, c);
    
    
    getContentPane().add(less_than_exons_label);

    andSeparator(gridbag, c);

    greater_than_exons_button = new JCheckBox("At least: ", false);
    final JPanel greater_than_exons_panel = new JPanel();
    greater_than_exons_panel.setLayout(new FlowLayout(FlowLayout.LEFT));
    greater_than_exons_panel.add(greater_than_exons_button);
    c.gridwidth = 2;
    gridbag.setConstraints(greater_than_exons_panel, c);
    getContentPane().add(greater_than_exons_panel);

    greater_than_exons_text = new JTextField("", 18);

    gridbag.setConstraints(greater_than_exons_text, c);
    getContentPane().add(greater_than_exons_text);

    final JLabel greater_than_exons_label = new JLabel(" exons long");
    
    c.gridwidth = GridBagConstraints.REMAINDER;
    gridbag.setConstraints(greater_than_exons_label, c);
    getContentPane().add(greater_than_exons_label);

    andSeparator(gridbag, c);
    
    intron_pattern_button = new JCheckBox("Contains introns without GT/GC start and AG end", false);
    final JPanel intron_pattern_panel = new JPanel();
    intron_pattern_panel.setLayout(new FlowLayout(FlowLayout.LEFT));
    intron_pattern_panel.add(intron_pattern_button);
    c.gridwidth = GridBagConstraints.REMAINDER;
    gridbag.setConstraints(intron_pattern_panel, c);
    getContentPane().add(intron_pattern_panel);
    
        
    {
      // leave a blank line to give the user a visual clue
      final JLabel blank_label = new JLabel("");

      gridbag.setConstraints(blank_label, c);
      getContentPane().add(blank_label);

      final JLabel and_label = new JLabel("And by:");

      gridbag.setConstraints(and_label, c);
      getContentPane().add(and_label);
    }


    by_motif_button = new JCheckBox("Amino acid motif: ", false);
    final JPanel by_motif_panel = new JPanel();
    by_motif_panel.setLayout(new FlowLayout(FlowLayout.LEFT));
    by_motif_panel.add(by_motif_button);
    c.gridwidth = 2;
    gridbag.setConstraints(by_motif_panel, c);
    getContentPane().add(by_motif_panel);

    motif_text = new JTextField("", 18);

    c.gridwidth = GridBagConstraints.REMAINDER;
    gridbag.setConstraints(motif_text, c);
    getContentPane().add(motif_text);

    {
      // leave a blank line to give the user a visual clue
      final JLabel blank_label = new JLabel("");

      gridbag.setConstraints(blank_label, c);
      getContentPane().add(blank_label);
    }

    final JPanel strand_panel = new JPanel();

    strand_panel.setLayout(new FlowLayout(FlowLayout.LEFT));

    forward_strand_checkbox = new JCheckBox("Forward Strand Features", true);
    strand_panel.add(forward_strand_checkbox);
    forward_strand_checkbox.addItemListener(new ItemListener()
    {
      public void itemStateChanged(ItemEvent _)
      {
        if(!forward_strand_checkbox.isSelected() &&
            !reverse_strand_checkbox.isSelected()) 
          // make sure one of the strand is always selected
          reverse_strand_checkbox.setSelected(true);
      }
    });


    reverse_strand_checkbox = new JCheckBox("Reverse Strand Features", true);
    strand_panel.add(reverse_strand_checkbox);
    reverse_strand_checkbox.addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent _) 
      {
        if(!reverse_strand_checkbox.isSelected() &&
            !forward_strand_checkbox.isSelected()) 
          // make sure one of the strand is always selected
          forward_strand_checkbox.setSelected(true);
      }
    });


    c.gridwidth = GridBagConstraints.REMAINDER;
    gridbag.setConstraints(strand_panel, c);
    getContentPane().add(strand_panel);





    final JButton select_button = new JButton("Select");

    select_button.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent e) 
      {
        selection.set(getSelected());
      }
    });

    final JButton view_button = new JButton("View");

    view_button.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent e)
      {
        final FeaturePredicate predicate =
          new FeatureFromVectorPredicate(getSelected());

        String title = "All features";

        if(by_key_button.isSelected()) 
          title += " with key \"" + key_selector.getSelectedItem() + "\"";

        if(by_qualifier_button.isSelected()) 
        {
          if(qualifier_text.getText().trim().length() > 0)
            title += " with qualifier \"" + qualifier_selector.getSelectedItem() +
                     "\" containing text \"" + qualifier_text.getText() + "\"";
          else
            title += " with qualifier \"" +
                     qualifier_selector.getSelectedItem() + "\"";
        }

        if(forward_strand_checkbox.isSelected() &&
           !reverse_strand_checkbox.isSelected()) 
          title += " on the forward strand";

        if(!forward_strand_checkbox.isSelected() &&
           reverse_strand_checkbox.isSelected()) 
          title += " on the reverse strand";

        if(by_motif_button.isSelected()) 
          title += " with motif: " + motif_text.getText().trim().toUpperCase();
        
        
        if(less_than_bases_button.isSelected())
          title += " at most " + less_than_bases_text.getText().trim() +
                   " bases long";
        
        if(greater_than_bases_button.isSelected())
          title += " at least " + greater_than_bases_text.getText().trim() +
                   " bases long";
        
        if(less_than_exons_button.isSelected()) 
          title += " at most " + less_than_exons_text.getText().trim() +
                   " exons long";
   
        
        if(greater_than_exons_button.isSelected()) 
          title += " at least " + greater_than_exons_text.getText().trim() +
                   " exons long";
        
        if(intron_pattern_button.isSelected())
          title += " containing introns without GT/GC start and AG end";
        
        final FilteredEntryGroup filtered_entry_group =
          new FilteredEntryGroup(entry_group, predicate, title);

        final FeatureListFrame feature_list_frame =
          new FeatureListFrame(title,
                               getSelection(),
                               goto_event_source, filtered_entry_group,
                               base_plot_group);

        feature_list_frame.setVisible(true);
      }
    });


    final JButton close_button = new JButton("Close");

    close_button.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent e) 
      {
        Selector.this.dispose();
      }
    });

    final FlowLayout flow_layout =
      new FlowLayout(FlowLayout.CENTER, 15, 5);

    final JPanel bottom_button_panel = new JPanel(flow_layout);

    bottom_button_panel.add(select_button);
    bottom_button_panel.add(view_button);
    bottom_button_panel.add(close_button);

    gridbag.setConstraints(bottom_button_panel, c);
    getContentPane().add(bottom_button_panel);

    addWindowListener(new WindowAdapter() 
    {
      public void windowClosing(WindowEvent event) 
      {
        getEntryGroup().removeEntryGroupChangeListener(Selector.this);
        Selector.this.dispose();
      }
    });

    getEntryGroup().addEntryGroupChangeListener(this);
    pack();

    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
    setLocation(new Point((screen.width  - getSize().width) / 2,
                          (screen.height - getSize().height) / 2));

    setVisible(true);
  }

  private void andSeparator(GridBagLayout gridbag, 
                            GridBagConstraints c)
  {
    // leave a blank line to give the user a visual clue
    final JLabel blank_label = new JLabel("");

    gridbag.setConstraints(blank_label, c);
    getContentPane().add(blank_label);

    final JLabel and_label = new JLabel("And:");
    
    gridbag.setConstraints(and_label, c);
    getContentPane().add(and_label);
  }
  
  /**
   *  Implementation of the EntryGroupChangeListener interface.  We listen to
   *  EntryGroupChange events so that we can get rid of the Navigator when the
   *  EntryGroup is no longer in use(for example when the EntryEdit is
   *  closed).
   **/
  public void entryGroupChanged(final EntryGroupChangeEvent event) 
  {
    switch(event.getType())
    {
      case EntryGroupChangeEvent.DONE_GONE:
       getEntryGroup().removeEntryGroupChangeListener(this);
       dispose();
       break;
    }
  }

  /**
   *  Return those features that match the current setting of the Selector.
   **/
  private FeatureVector getSelected() 
  {
    final FeaturePredicate key_and_qualifier_predicate;
    final FeaturePredicate motif_predicate;

    if(by_key_button.isSelected()) 
    {
      if(by_qualifier_button.isSelected()) 
      {
        final String search_text = qualifier_text.getText().trim();
        final String qualifier_name =
         (String) qualifier_selector.getSelectedItem();

        if(search_text.length() == 0) 
        {
          key_and_qualifier_predicate =
            new FeatureKeyQualifierPredicate(key_selector.getSelectedItem(),
                                             qualifier_name,
                                             true);
        }
        else
        {
          if(match_any_word_checkbox.isSelected()) 
          {
            final FeaturePredicateVector temp_predicates =
              new FeaturePredicateVector();

            //final StringVector words =
            //  StringVector.getStrings(search_text, " ");
            
            final StringTokenizer tok = new StringTokenizer(search_text, " \n");

            while(tok.hasMoreTokens()) 
            {
              final String this_word = tok.nextToken().trim();
              final FeaturePredicate new_predicate =
                new FeatureKeyQualifierPredicate(key_selector.getSelectedItem(),
                                                 qualifier_name,
                                                 this_word,
                                                 partial_match_checkbox.isSelected(),
                                                 ignore_case_checkbox.isSelected());

              temp_predicates.add(new_predicate);
            }

            key_and_qualifier_predicate =
              new FeaturePredicateConjunction(temp_predicates,
                                              FeaturePredicateConjunction.OR);
          } 
          else 
          {
            key_and_qualifier_predicate =
              new FeatureKeyQualifierPredicate(key_selector.getSelectedItem(),
                                               qualifier_name,
                                               search_text,
                                               partial_match_checkbox.isSelected(),
                                               ignore_case_checkbox.isSelected());
            
          }
        }
      } 
      else 
      {
        final String search_text = qualifier_text.getText().trim();
        if(search_text.length() == 0) 
        {
          key_and_qualifier_predicate =
            new FeatureKeyPredicate(key_selector.getSelectedItem());
        } 
        else 
        {
          if(match_any_word_checkbox.isSelected())
          {
            final FeaturePredicateVector temp_predicates =
              new FeaturePredicateVector();
            
            final StringVector words = 
              StringVector.getStrings(search_text, " ");
            

            for(int i = 0 ; i < words.size() ; ++i) 
            {
              final String this_word =(String)words.elementAt(i);
              final FeaturePredicate new_predicate =
                new FeatureKeyQualifierPredicate(key_selector.getSelectedItem(),
                                                 null, // match any qialifier
                                                 this_word,
                                                 partial_match_checkbox.isSelected(),
                                                 ignore_case_checkbox.isSelected());
              temp_predicates.add(new_predicate);
            }

            key_and_qualifier_predicate =
              new FeaturePredicateConjunction(temp_predicates,
                                              FeaturePredicateConjunction.OR);
          } 
          else 
          {
            key_and_qualifier_predicate =
              new FeatureKeyQualifierPredicate(key_selector.getSelectedItem(),
                                               null, // match any qialifier
                                               search_text,
                                               partial_match_checkbox.isSelected(),
                                               ignore_case_checkbox.isSelected());
          }
        }
      }
    } 
    else 
      key_and_qualifier_predicate = null;

    if(by_motif_button.isSelected()) 
    {
      final AminoAcidSequence amino_acid_sequence =
        new AminoAcidSequence(motif_text.getText().trim());

      motif_predicate =
        new FeaturePatternPredicate(amino_acid_sequence);
    } 
    else 
      motif_predicate = null;


    FeaturePredicate less_than_bases_predicate = null;
    if(less_than_bases_button.isSelected() &&
       less_than_bases_text.getText().trim().length() > 0)
    {
      try 
      {
        final int less_than_bases_int =
          Integer.valueOf(less_than_bases_text.getText().trim()).intValue();

        less_than_bases_predicate = new FeaturePredicate() 
        {
          public boolean testPredicate(final Feature feature) 
          {
            if(feature.getBaseCount() <= less_than_bases_int) 
              return true;
            else 
              return false;
          }
        };
      } 
      catch(NumberFormatException e) 
      {
        new MessageDialog(this,
                           "warning this is not a number: " +
                           less_than_bases_text.getText());
        less_than_bases_predicate = null;
      }
    } 


    FeaturePredicate greater_than_bases_predicate = null;
    if(greater_than_bases_button.isSelected() &&
       greater_than_bases_text.getText().trim().length() > 0) 
    {
      try 
      {
        final int greater_than_bases_int =
          Integer.valueOf(greater_than_bases_text.getText().trim()).intValue();

        greater_than_bases_predicate = new FeaturePredicate() 
        {
          public boolean testPredicate(final Feature feature) 
          {
            if(feature.getBaseCount() >= greater_than_bases_int) 
              return true;
            else 
              return false;
          }
        };
      } 
      catch(NumberFormatException e) 
      {
        new MessageDialog(this,
                          "warning this is not a number: " +
                          greater_than_bases_text.getText());
        greater_than_bases_predicate = null;
      }
    } 



    FeaturePredicate less_than_exons_predicate = null;
    if(less_than_exons_button.isSelected() &&
       less_than_exons_text.getText().trim().length() > 0)
    {
      try
      {
        final int less_than_exons_int =
          Integer.valueOf(less_than_exons_text.getText().trim()).intValue();

        less_than_exons_predicate = new FeaturePredicate()
        {
          public boolean testPredicate(final Feature feature)
          {
            if(feature.getSegments().size() <= less_than_exons_int) 
              return true;
            else 
              return false;
          }
        };
      } 
      catch(NumberFormatException e) 
      {
        new MessageDialog(this,
                          "warning this is not a number: " +
                          less_than_exons_text.getText());
        less_than_exons_predicate = null;
      }
    } 


    FeaturePredicate greater_than_exons_predicate = null;
    if(greater_than_exons_button.isSelected() &&
       greater_than_exons_text.getText().trim().length() > 0)
    {
      try 
      {
        final int greater_than_exons_int =
          Integer.valueOf(greater_than_exons_text.getText().trim()).intValue();

        greater_than_exons_predicate = new FeaturePredicate()
        {
          public boolean testPredicate(final Feature feature)
          {
            if(feature.getSegments().size() >= greater_than_exons_int) 
              return true;
            else 
              return false;
          }
        };
      } 
      catch(NumberFormatException e) 
      {
        new MessageDialog(this,
                          "warning this is not a number: " +
                          greater_than_exons_text.getText());
        greater_than_exons_predicate = null;
      }
    } 
    
    
    FeaturePredicate intron_pattern_predicate = null;
    if(intron_pattern_button.isSelected())
      intron_pattern_predicate = getIntronPredicate();
    
    FeaturePredicate predicate = null;
    if(!by_key_button.isSelected() && !by_motif_button.isSelected())
    {
      // default to selecting all features
      predicate = new FeaturePredicate() 
      {
        public boolean testPredicate(final Feature feature)
        {
          return true;
        }
      };
    } 
    else
    {
      if(motif_predicate != null && key_and_qualifier_predicate != null)
      {
        predicate =
          new FeaturePredicateConjunction(key_and_qualifier_predicate,
                                          motif_predicate,
                                          FeaturePredicateConjunction.AND);
      } 
      else 
      {
        if(motif_predicate != null) 
          predicate = motif_predicate;
        else
          predicate = key_and_qualifier_predicate;
      }
    }

    if(less_than_bases_predicate != null) 
    {
      predicate =
        new FeaturePredicateConjunction(predicate,
                                        less_than_bases_predicate,
                                        FeaturePredicateConjunction.AND);
    }

    if(greater_than_bases_predicate != null)
    {
      predicate =
        new FeaturePredicateConjunction(predicate,
                                         greater_than_bases_predicate,
                                         FeaturePredicateConjunction.AND);
    }

    if(less_than_exons_predicate != null)
    {
      predicate =
        new FeaturePredicateConjunction(predicate,
                                        less_than_exons_predicate,
                                        FeaturePredicateConjunction.AND);
    }

    if(greater_than_exons_predicate != null)
    {
      predicate =
        new FeaturePredicateConjunction(predicate,
                                        greater_than_exons_predicate,
                                        FeaturePredicateConjunction.AND);
    }

    if(intron_pattern_predicate !=null)
    {
      predicate =
        new FeaturePredicateConjunction(predicate,
                                        intron_pattern_predicate, 
                                        FeaturePredicateConjunction.AND);
    }
    
    final FeatureEnumeration test_enumerator = entry_group.features();
    final FeatureVector return_features = new FeatureVector();

    while(test_enumerator.hasMoreFeatures())
    {
      final Feature this_feature = test_enumerator.nextFeature();

      if(predicate.testPredicate(this_feature))
      {
        if(this_feature.isForwardFeature())
        {
          if(forward_strand_checkbox.isSelected()) 
            return_features.add(this_feature);
        } 
        else
        {
          if(reverse_strand_checkbox.isSelected()) 
            return_features.add(this_feature);
        }
      }
    }

    return return_features;
  }
  
  protected static FeaturePredicate getIntronPredicate()
  {
    return new FeaturePredicate()
    {
      public boolean testPredicate(final Feature feature)
      { 
        if(feature.getSegments().size() < 2)
          return false;
        
        final Location location = feature.getLocation().copy();
        final RangeVector ranges = location.getRanges();
        Collections.sort(ranges, new Comparator<Range>(){
          public int compare(final Range r1, final Range r2)
          {
            return r1.getStart()-r2.getStart();
          }
        });
        
        final Strand strand = feature.getStrand();
        for(int i = 0; i < ranges.size () -1; ++i) 
        {
          final int end_of_range_1 =
            ((Range)ranges.elementAt(i)).getEnd ();
          final int start_of_range_2 =
            ((Range)ranges.elementAt(i+1)).getStart ();
          
          // ignore - the exons overlap so there is no room for an intron
          if(end_of_range_1 > start_of_range_2) 
            continue;

          try
          {
            Range feature_range = new Range(end_of_range_1 + 1,
                                            start_of_range_2 - 1);

            final char bases[];
            if(location.isComplement())
            {
              final char tmp_bases[] = strand.getRawSubSequenceC(feature_range);
              final char tmp2_bases[] = new char[feature_range.getCount()];
              
              for(int j=0; j<tmp2_bases.length; j++)
                tmp2_bases[j] = tmp_bases[j];
              bases = Bases.reverseComplement(  tmp2_bases );
            }
            else
              bases = strand.getRawSubSequenceC(feature_range);
                         
            if(bases.length < 3)
              return true;

            
            int length = feature_range.getCount();
            if( bases[0] != 'g' ||
               (bases[1] != 't' && bases[1] != 'c') ||
                bases[length-1] != 'g' ||
                bases[length-2] != 'a' )
            {
              if(bases[length-1] != 'g' || bases[length-2] != 'a')
                if(location.isComplement())
                  logger4j.info("INTRON SPLICE SITE MISSING AG: "+end_of_range_1);
                else
                  logger4j.info("INTRON SPLICE SITE MISSING AG: "+start_of_range_2);
              
              else
                if(location.isComplement())
                  logger4j.info("INTRON SPLICE SITE MISSING GT/GC: "+start_of_range_2);
                else
                  logger4j.info("INTRON SPLICE SITE MISSING GT/GC: "+end_of_range_1);
              return true;
            }
            
          }
          catch(uk.ac.sanger.artemis.util.OutOfRangeException oore)
          {           
          }
        }

        return false;
      }
    };
  }

  /**
   *  Return the Selection object that was passed to the constructor.
   **/
  private Selection getSelection() 
  {
    return selection;
  }

  /**
   *  Return the EntryGroup object that was passed to the constructor.
   **/
  private EntryGroup getEntryGroup() 
  {
    return entry_group;
  }

}
