/* FeatureInfo.java
 *
 * created: Sat Dec 19 1998
 *
 * This file is part of Artemis
 *
 * Copyright (C) 1998,1999,2000,2002  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/FeatureInfo.java,v 1.1 2004-06-09 09:46:39 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.sequence.*;
import uk.ac.sanger.artemis.plot.*;

import java.awt.*;
import java.awt.event.*;

import javax.swing.*;

/**
 *  This component displays a summary of the statistical information about one
 *  feature (such as correlation score, codon usage and gc content).
 *
 *  @author Kim Rutherford
 *  @version $Id: FeatureInfo.java,v 1.1 2004-06-09 09:46:39 tjc Exp $
 **/
public class FeatureInfo extends JFrame
    implements EntryChangeListener, FeatureChangeListener {
  /**
   *  Create a new FeatureInfo component.
   *  @param feature The Feature that this component is showing information
   *    about.
   **/
  public FeatureInfo (final Feature feature,
                      final CodonUsageAlgorithm codon_usage_algorithm) {
    super ("Feature infomation: " + feature.getIDString ());

    this.feature = feature;
    this.entry = feature.getEntry ();
    this.codon_usage_algorithm = codon_usage_algorithm;

    setBackground (new Color (210, 210, 210));
    getContentPane ().setBackground (new Color (210, 210, 210));

    codon_info_areas = new JTextArea [4][4];
    base_count_info_areas = new JTextArea [4];

    final JPanel centre_panel = new JPanel ();
    centre_panel.setLayout (new BorderLayout ());

    makeCountList ();
    getContentPane ().add (aa_count_panel, "West");

    makeMiscInfogrid ();
    centre_panel.add (misc_info_panel, "South");

    makeCodonInfogrid ();
    centre_panel.add (codon_info_panel, "Center");

    getContentPane ().add (centre_panel, "Center");

    button_panel = new JPanel ();
    final JButton close_button = new JButton ("Close");

    close_button.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent e) {
        stopListening ();
        FeatureInfo.this.dispose ();
      }
    });

    button_panel.add (close_button);

    getContentPane ().add (button_panel, "South");

    updateComponents ();

    getFeature ().addFeatureChangeListener (this);
    getFeature ().getEntry ().addEntryChangeListener (this);

    addWindowListener (new WindowAdapter () {
      public void windowClosing (WindowEvent event) {
        stopListening ();
        FeatureInfo.this.dispose ();
      }
    });

    pack ();

    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
    setLocation (new Point ((screen.width - getSize ().width) / 2,
                            (screen.height - getSize ().height) / 2));

    setVisible (true);
  }

  /**
   *  Remove this object as a feature change and entry change listener.
   **/
  private void stopListening () {
    getEntry ().removeEntryChangeListener (this);
    getFeature ().removeFeatureChangeListener (this);
  }

  /**
   *  Implementation of the EntryChangeListener interface.  We listen to
   *  EntryChange events so we can notify the user if of this component if the
   *  feature gets deleted.
   **/
  public void entryChanged (EntryChangeEvent event) {
    switch (event.getType ()) {
    case EntryChangeEvent.FEATURE_DELETED:
      if (event.getFeature () == getFeature ()) {
        stopListening ();
        dispose ();
      }
      break;
    default:
      // do nothing;
      break;
    }
  }

  /**
   *  Implementation of the FeatureChangeListener interface.  We need to
   *  listen to feature change events from the Features in this object.
   *  @param event The change event.
   **/
  public void featureChanged (FeatureChangeEvent event) {
    updateComponents ();
  }

  /**
   *  Create all the labels in the misc_info_panel.
   **/
  private void makeMiscInfogrid () {
    misc_info_panel = new JPanel ();

    GridLayout grid = new GridLayout (0, 1);
    misc_info_panel.setLayout (grid);

    molecular_weight_label = new JLabel ();
    misc_info_panel.add (molecular_weight_label);

    correlation_scores_label = new JLabel ();
    misc_info_panel.add (correlation_scores_label);

    if (codon_usage_algorithm != null) {
      usage_scores_label = new JLabel ();
      misc_info_panel.add (usage_scores_label);
    }
  }

  /**
   *  Create a List containing a count of the occurrences of each codon in the
   *  feature.
   **/
  private void makeCountList () {
    aa_count_panel = new JPanel ();

    final JPanel aa_count_sub_panel1 = new JPanel ();
    aa_count_sub_panel1.setLayout (new GridLayout (0, 1));
    aa_count_panel.add (aa_count_sub_panel1);

    final JPanel aa_count_sub_panel2 = new JPanel ();
    aa_count_sub_panel2.setLayout (new GridLayout (0, 1));
    aa_count_panel.add (aa_count_sub_panel2);

    aa_count_list = new JLabel [AminoAcidSequence.symbol_count];

    for (int i = 0 ; i < AminoAcidSequence.symbol_count / 2 ; ++i) {
      aa_count_list[i] = new JLabel ();
      aa_count_sub_panel1.add (aa_count_list[i]);
    }

    // add a blank label so that the columns balance
    aa_count_sub_panel1.add (new JLabel (""));

    for (int i = AminoAcidSequence.symbol_count / 2 ;
         i < AminoAcidSequence.symbol_count ;
         ++i) {
      aa_count_list[i] = new JLabel ();
      aa_count_sub_panel2.add (aa_count_list[i]);
    }
  }

  /**
   *  Create all the labels in the codon_info_panel.
   **/
  private void makeCodonInfogrid () {
    codon_info_panel = new JPanel ();

    GridBagLayout gridbag = new GridBagLayout();

    codon_info_panel.setLayout (gridbag);

    GridBagConstraints c = new GridBagConstraints();

    final int GRID_WIDTH = 7;

    c.gridwidth = GRID_WIDTH;
    c.anchor = GridBagConstraints.NORTH;
    c.insets = new Insets (1, 1, 1, 1);
    c.fill = GridBagConstraints.NONE;


    // create the first row with some dummy labels to pad it out

    // the first row has the bases T,C,A,G as column titles
    final JLabel dummy_label1 = new JLabel ("");
    gridbag.setConstraints (dummy_label1, c);
    codon_info_panel.add (dummy_label1);

    final JLabel dummy_label2 = new JLabel ("");
    gridbag.setConstraints (dummy_label2, c);
    codon_info_panel.add (dummy_label2);

    for (int i = 0 ; i < Bases.letter_index.length ; ++i) {
      final JLabel new_label =
        new JLabel (String.valueOf (Bases.letter_index[i]).toUpperCase ());
      gridbag.setConstraints (new_label, c);
      codon_info_panel.add (new_label);
    }

    final JLabel dummy_label3 = new JLabel ("");
    c.gridwidth = GridBagConstraints.REMAINDER;
    gridbag.setConstraints (dummy_label3, c);
    codon_info_panel.add (dummy_label3);

    c.anchor = GridBagConstraints.WEST;
    c.gridwidth = GRID_WIDTH;

    final int ROWS = 4;

    // create the main rows

    for (int row = 0 ; row < ROWS ; ++row) {
      {
        final JLabel new_label =
          new JLabel (String.valueOf (Bases.letter_index[row]).toUpperCase ());

        new_label.setBorder (new javax.swing.border.EmptyBorder (1,1,1,1));
        gridbag.setConstraints (new_label, c);
        codon_info_panel.add (new_label);
      }

      {
        base_count_info_areas[row] = new JTextArea (5, 15);
        base_count_info_areas[row].setEditable (false);

        final JScrollPane scroll_pane =
          new JScrollPane (base_count_info_areas[row]);
        gridbag.setConstraints (scroll_pane, c);
        codon_info_panel.add (scroll_pane);

        final Color parent_background = getBackground ();
        base_count_info_areas[row].setBackground (parent_background);
      }

      for (int column = 2 ; column < GRID_WIDTH - 1 ; ++column) {

        final JTextArea new_text = new JTextArea (5, 12);
        new_text.setText ("r: " + row + " c: " + column);

        new_text.setEditable (false);

        new_text.setBackground (Color.white);

        final JScrollPane scroll_pane = new JScrollPane (new_text);
        gridbag.setConstraints (scroll_pane, c);
        codon_info_panel.add (scroll_pane);

        codon_info_areas[row][column - 2] = new_text;
      }

      c.gridwidth = GridBagConstraints.REMAINDER;

      {
        final JTextArea new_text = new JTextArea (5, 2);
        new_text.setText (Bases.letter_index[0] + "\n" +
                          Bases.letter_index[1] + "\n" +
                          Bases.letter_index[2] + "\n" +
                          Bases.letter_index[3]);

        new_text.setEditable (false);

        final JScrollPane scroll_pane = new JScrollPane (new_text);
        gridbag.setConstraints (scroll_pane, c);
        codon_info_panel.add (scroll_pane);
      }

      c.gridwidth = GRID_WIDTH;
    }
  }

  /**
   *  Set the List components with the information from the feature.
   **/
  private void updateComponents () {
    updateAACountList ();

    updateCodonInfoAreas ();

    updateMolecularWeightLabel ();

    updateCorrelationScoresLabel ();

    updateUsageScoresLabel ();

    updateBaseCounts ();

    validate ();
  }

  /**
   *  Update the aa_count_list JTextArea.
   **/
  private void updateAACountList () {
    for (int i = 0 ; i < AminoAcidSequence.symbol_count ; ++i) {
      if (i == AminoAcidSequence.symbol_count - 1 &&
          getFeature ().getResidueCount (i) == 0) {
        // don't include Selenocysteine in the list in the list if there are
        // none
        aa_count_list[i].setText ("");
        continue;
      }

      final String three_letter_abbreviation =
        AminoAcidSequence.getThreeLetterAbbreviation (i);

      final char one_letter_abbreviation =
        Character.toUpperCase (AminoAcidSequence.getSymbolFromIndex (i)) ;

      final String label_string =
        three_letter_abbreviation + " (" + one_letter_abbreviation +"): " +
        getFeature ().getResidueCount (i);

      aa_count_list[i].setText (label_string);
    }
  }

  /**
   *  Update all the JTextArea objects in codon_info_areas.
   **/
  private void updateCodonInfoAreas () {
    for (int first = 0 ; first < 4 ; ++first) {

      for (int second = 0 ; second < 4 ; ++second) {

        codon_info_areas[first][second].setText ("");

        for (int third = 0 ; third < 4 ; ++third) {

          final char this_codon_symbol =
            AminoAcidSequence.getCodonTranslation (Bases.letter_index[first],
                                                   Bases.letter_index[second],
                                                   Bases.letter_index[third]);

          final String this_codon_string =
            AminoAcidSequence.getThreeLetterAbbreviation (this_codon_symbol);

          final int this_codon_count =
            getFeature ().getCodonCount (first,second,third);

          final int this_codon_index =
            AminoAcidSequence.getSymbolIndex (this_codon_symbol);

          final int this_amino_acid_count =
            getFeature ().getResidueCount (this_codon_index);

          final String percent_string;

          if (this_amino_acid_count == this_codon_count) {
            percent_string = "ALL";
          } else {
            if (this_amino_acid_count < this_codon_count) {
              // some of the amino acids aren't coded for by the standard
              // codon - probably because of the /transl_except qualifier
              percent_string = "---";
            } else {
              final int percentage =
                100 * this_codon_count / this_amino_acid_count;

              percent_string = String.valueOf (percentage) + "%";
            }
          }

          final String first_part_of_line =
            this_codon_string + " " + this_codon_count;

          codon_info_areas[first][second].append (first_part_of_line);

          final int line_length_without_spaces =
            first_part_of_line.length () + percent_string.length ();

          // add some spaces so that everything lines up nicely
          for (int i = 0 ; i < 11 - line_length_without_spaces ; ++i) {
            codon_info_areas[first][second].append (" ");
          }

          codon_info_areas[first][second].append (percent_string);

          if (third != 3) {
            // insert a newline on all but the last line
            codon_info_areas[first][second].append ("\n");
          }
        }
      }
    }
  }

  /**
   *  Update molecular_weight_label.
   **/
  private void updateMolecularWeightLabel () {
    final String new_text =
      "Mol weight: " +
      getFeature ().getTranslation ().getMolecularWeight () + "  Start: " +
      getFeature ().getFirstBase () + "  End: " +
      getFeature ().getLastBase () + "  Bases: " +
      getFeature ().getBaseCount () + "  AA length: " +
      getFeature ().getTranslation ().length ();
    molecular_weight_label.setText (new_text);
  }


  /**
   *  The method updates the components in base_count_info_areas.
   **/
  private void updateBaseCounts () {
    final int translation_base_count =
      getFeature ().getTranslationBases ().length ();

    final int amino_acid_count = getFeature ().getTranslation ().length ();

    for (int i = 0 ; i < 4 ; ++i) {
      base_count_info_areas[i].setText ("");
      {
        final String string =
          "ALL:" + updateBaseCountsFormatter (getFeature ().getBaseCount (i),
                                              translation_base_count);

        base_count_info_areas[i].append (string);
      }

      for (int codon_base_index = 0;
           codon_base_index < 3 ;
           ++codon_base_index) {
        String label = null;
        switch (codon_base_index) {
        case 0: label = "\n1st:"; break;
        case 1: label = "\n2nd:"; break;
        case 2: label = "\n3rd:"; break;
        }

        final int base_count =
          getFeature ().getPositionalBaseCount (codon_base_index, i);

        final String string =
          label + updateBaseCountsFormatter (base_count, amino_acid_count);

        base_count_info_areas[i].append (string);
      }
    }
  }

  /**
   *  A helper method for updateBaseCounts ().
   **/
  private String updateBaseCountsFormatter (int base_count,
                                            int total) {
    final String count_string = "     " + base_count;

    final String percent_string;

    if (base_count < total) {
      percent_string = " " + 100 * base_count / total + "%";
    } else {
      percent_string = "ALL";
    }

    final int count_width = 5;

    final int percent_width = 3;

    return
      count_string.substring (count_string.length () - count_width) + " " +
      percent_string.substring (percent_string.length () - percent_width);
  }

  /**
   *  Update correlation_scores_label.
   **/
  private void updateCorrelationScoresLabel () {
    final int c_total =
      getFeature ().getBaseCount (Bases.getIndexOfBase ('c'));
    final int g_total =
      getFeature ().getBaseCount (Bases.getIndexOfBase ('g'));

    final String c3_score;
    final String g1_score;
    final String g3_score;

    final int c3_count =
      getFeature ().getPositionalBaseCount (2, Bases.getIndexOfBase ('c'));
    final int g1_count =
      getFeature ().getPositionalBaseCount (0, Bases.getIndexOfBase ('g'));
    final int g3_count =
      getFeature ().getPositionalBaseCount (2, Bases.getIndexOfBase ('g'));

    if (c_total == 0) {
      c3_score = "ALL";
    } else {
      c3_score =
        String.valueOf (1000 * (3 * c3_count - c_total) / c_total / 10.0);
    }

    if (g_total == 0) {
      g1_score = "ALL";
    } else {
      g1_score =
        String.valueOf (1000 * (3 * g1_count - g_total) / g_total / 10.0);
    }

    if (g_total == 0) {
      g3_score = "ALL";
    } else {
      g3_score =
        String.valueOf (1000 * (3 * g3_count - g_total) / g_total / 10.0);
    }

    final double cor1_2_score =
      ((int) getFeature ().get12CorrelationScore () * 10) / 10.0;

    final double gc_percent =
      ((int) (getFeature ().getPercentGC () * 100.0)) / 100.0;

    final String new_label =
      "position 1/2 score = " + cor1_2_score + "  " +
      "C3/G1/G3 (o-e)/e = " + c3_score + " " + g1_score + " " + g3_score +
      "  " + gc_percent + "% GC";

    correlation_scores_label.setText (new_label);
  }

  /**
   *  Update usage_scores_label.
   **/
  private void updateUsageScoresLabel () {
    if (codon_usage_algorithm != null) {
      final String new_label =
        "usage score = " + codon_usage_algorithm.getFeatureScore (feature);

      usage_scores_label.setText (new_label);
    }
  }

  /**
   *  Return the feature this component is showing information about.
   **/
  private Feature getFeature () {
    return feature;
  }

  /**
   *  Return the Entry that contains the Feature this object is displaying.
   **/
  private Entry getEntry () {
    return entry;
  }

  /**
   *  The Feature that this component is showing information about.
   **/
  private Feature feature = null;

  /**
   *  The Entry that contains the Feature this object is displaying.
   **/
  private Entry entry;

  /**
   *  A panel containing the labels with the amino acid counts for the
   *  feature.
   **/
  private JPanel aa_count_panel = null;

  /**
   *  An array Labels of showing the amino acid counts for the feature.  There
   *  is one Label for each symbol (24).
   **/
  private JLabel [] aa_count_list = null;

  /**
   *  A panel containing general statistics about the whole sequence.
   **/
  private JPanel misc_info_panel = null;

  /**
   *  A Label containing molecular weight of the protein, start and end
   *  positions and length.
   **/
  private JLabel molecular_weight_label = null;

  /**
   *  A Label containing the correlation scores for this feature.
   **/
  private JLabel correlation_scores_label = null;

  /**
   *  A Label containing the usage scores for this feature.
   **/
  private JLabel usage_scores_label = null;

  /**
   *  Statistics about the codon frequency.
   **/
  private JPanel codon_info_panel = null;

  /**
   *  A panel to hold the close button.
   **/
  private JPanel button_panel = null;

  /**
   *  These components each hold information about 4 codons each.
   **/
  private JTextArea [] [] codon_info_areas;

  /**
   *  These components display position information about the base counts.
   *  See positional_base_counts and base_counts.
   **/
  private JTextArea [] base_count_info_areas;

  /**
   *  The CodonUsageAlgorithm reference that was passed to the constructor.
   *  (Used by updateUsageScoresLabel ()).
   **/
  private final CodonUsageAlgorithm codon_usage_algorithm;
}
