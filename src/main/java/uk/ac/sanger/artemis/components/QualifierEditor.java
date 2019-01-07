/* QualifierEditor.java
 *
 * created: Tue Oct 23 2001
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2000,2001,2002  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/QualifierEditor.java,v 1.5 2007-07-09 13:07:38 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;

import uk.ac.sanger.artemis.util.ReadOnlyException;
import uk.ac.sanger.artemis.io.GFFDocumentEntry;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierInfo;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.io.QualifierParseException;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.EntryInformationException;

import java.awt.*;
import java.awt.event.*;
import java.util.Vector;

import javax.swing.*;

/**
 *  This component allows qualifiers to be added to or replaced in several
 *  features at once.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: QualifierEditor.java,v 1.5 2007-07-09 13:07:38 tjc Exp $
 **/

public class QualifierEditor extends JFrame {
  /**
   *  Create a new QualifierEditor for the given features.
   **/
  public QualifierEditor (final FeatureVector features,
                          final EntryGroup entry_group) {
    super (getFrameTitle (features));

    this.features = features;
    this.entry_group = entry_group;

    final Font font = Options.getOptions ().getFont ();

    final Feature first_feature = features.elementAt (0);

    final EntryInformation entry_information =
      first_feature.getEntry ().getEntryInformation ();

    setFont (font);

    boolean isGFF = false;
    if(entry_group.getDefaultEntry().getEMBLEntry() instanceof GFFDocumentEntry)
      isGFF = true;
    final QualifierChoice qualifier_choice =
      new QualifierChoice (entry_information, first_feature.getKey (), null, isGFF);

    final JPanel outer_qualifier_choice_panel = new JPanel ();
    final JPanel qualifier_choice_panel = new JPanel ();
    outer_qualifier_choice_panel.setLayout (new BorderLayout ());

    outer_qualifier_choice_panel.add (qualifier_choice_panel, "West");

    final JButton qualifier_button = new JButton ("Insert qualifier:");

    qualifier_button.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent e) {
        final String qualifier_name =
          (String) qualifier_choice.getSelectedItem ();
        final QualifierInfo qualifier_info =
          entry_information.getQualifierInfo (qualifier_name);
        
        if (qualifier_info == null) {
          new MessageDialog (QualifierEditor.this, "internal error: no " +
                             "qualifier info for " + qualifier_name);
        } else {
          qualifier_text_area.append ("/" + qualifier_name);
          
          switch (qualifier_info.getType ()) {
          case QualifierInfo.QUOTED_TEXT:
            qualifier_text_area.append ("=\"\"");
            break;
            
          case QualifierInfo.NO_VALUE:
          case QualifierInfo.OPTIONAL_QUOTED_TEXT:
            break;
            
          default:
            qualifier_text_area.append ("=");
          }
          
          qualifier_text_area.append ("\n");
        }
      }
      
    });

    qualifier_choice_panel.add (qualifier_button);

    qualifier_choice_panel.add (qualifier_choice);

    getContentPane ().add (outer_qualifier_choice_panel, "North");

    qualifier_text_area = new QualifierTextArea ();

    add_button.setFont (getFont ());
    replace_button.setFont (getFont ());
    close_button.setFont (getFont ());

    add_button.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent e) {
        addOrReplace (false);
      }
    });

    replace_button.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent e) {
        addOrReplace (true);
      }
    });

    close_button.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent e) {
        dispose ();
      }
    });

    button_panel.setFont (getFont ());

    button_panel.add (replace_button);
    button_panel.add (add_button);
    button_panel.add (close_button);

    getContentPane ().add (qualifier_text_area, "Center");
    getContentPane ().add (button_panel, "South");

    addWindowListener (new WindowAdapter () {
      public void windowClosing (WindowEvent event) {
        dispose ();
      }
    });

    pack ();

    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
    setLocation (new Point ((screen.width - getSize ().width) / 2,
                            (screen.height - getSize ().height) / 2));
  }

  /**
   *  Add to or replace the qualifiers in all the features that were passed to
   *  the constructor with the qualifiers from qualifier_text_area.
   *  @param replace If false any existing qualifiers of the same name in the
   *    features are left unchanged.  If true existing qualifiers of the same
   *    name in the features will be deleted.
   **/
  private void addOrReplace (final boolean replace) {
    try 
    {
      entry_group.getActionController().startAction();

      // this will contain one QualifierVector object for each Feature in the
      // features vector (in the same order)
      final Vector qualifier_vector_vector = new Vector();
      final int features_size = features.size();

      for(int i = 0; i<features_size; ++i) 
      {
        final Feature this_feature = features.elementAt(i);

        final Entry this_feature_entry = this_feature.getEntry();

        if(this_feature_entry == null)
        {
          // feature has already been deleted
          qualifier_vector_vector.addElement(null);
        } 
        else 
        {
          final EntryInformation entry_information =
            this_feature_entry.getEntryInformation();

          try 
          {
            final QualifierVector qualifiers =
              qualifier_text_area.getParsedQualifiers(entry_information);

            qualifier_vector_vector.addElement(qualifiers);

          }
          catch(QualifierParseException e) 
          {
            new MessageDialog(this,
                              "error while parsing: " + e.getMessage());
            return;
          }
        }
      }

      if(qualifier_vector_vector.size() != features_size) 
        throw new Error("Internal error in QualifierEditor.add() - " +
                        "mismatched array sizes");

      for(int feature_index = 0; feature_index < features_size;
          ++feature_index) 
      {
        final Feature this_feature = features.elementAt(feature_index);
        this_feature.resetColour();

        if(qualifier_vector_vector.elementAt(feature_index) == null)
          continue;

        final QualifierVector qualifier_vector =
          (QualifierVector)qualifier_vector_vector.elementAt(feature_index);

        final int qualifier_vector_size = qualifier_vector.size();

        for(int qualifier_index = 0; qualifier_index < qualifier_vector_size;
            ++qualifier_index) 
        {
          final Qualifier this_qualifier =
            (Qualifier)qualifier_vector.elementAt(qualifier_index);

          if(replace) 
          {
            try 
            {
              this_feature.setQualifier(this_qualifier);
            } 
            catch(EntryInformationException e) 
            {
              new MessageDialog(this,
                                "failed to add qualifiers to: " +
                                this_feature.getIDString() + ": " +
                                e.getMessage());
            } 
            catch(ReadOnlyException e) 
            {
              new MessageDialog(this,
                                "failed to add qualifiers to read-only " +
                                "feature: " + this_feature.getIDString() +
                                ": " + e.getMessage());
            }
          } 
          else 
          {
            try
            {
              this_feature.addQualifierValues(this_qualifier);
            } 
            catch(EntryInformationException e) 
            {
              new MessageDialog(this,
                                "failed to add qualifiers to: " +
                                this_feature.getIDString() + ": " +
                                e.getMessage());
            } 
            catch(ReadOnlyException e) 
            {
              new MessageDialog(this,
                                "failed to add qualifiers to read-only " +
                                "feature: " + this_feature.getIDString() +
                                ": " + e.getMessage());
            }

          }
        }
      }
    } 
    finally 
    {
      entry_group.getActionController().endAction();
    }
  }

  /**
   *  Return an appropriate String to use for the title of this JFrame.
   **/
  static private String getFrameTitle (final FeatureVector features) {
    boolean etc_flag = false;

    final StringBuffer buffer = new StringBuffer ();

    final int MAX_LENGTH = 80;

    for (int i = 0 ; i < features.size () ; ++i) {
      final Feature this_feature = features.elementAt (i);

      final String feature_name = this_feature.getIDString ();

      if (feature_name == null) {
        etc_flag = true;
        continue;
      }

      if (buffer.length () + feature_name.length () < MAX_LENGTH) {
        if (buffer.length () == 0) {
          buffer.append ("Add or replace qualifiers of: ");
          buffer.append (feature_name);
        } else {
          buffer.append (", ").append (feature_name);
        }
      } else {
        etc_flag = true;
        break;
      }
    }

    if (buffer.length () == 0) {
      buffer.append ("Add or replace qualifiers");
    } else {
      if (etc_flag) {
        buffer.append (", ...");
      }
    }

    return buffer.toString ();
  }

  private QualifierTextArea qualifier_text_area;

  private JButton add_button = new JButton ("Add");
  private JButton replace_button = new JButton ("Replace");
  private JButton close_button = new JButton ("Close");

  private FlowLayout flow_layout =
    new FlowLayout (FlowLayout.CENTER, 25, 5);

  private JPanel button_panel = new JPanel (flow_layout);

  /**
   *  The Feature objects that were passed to the constructor.
   **/
  private FeatureVector features = new FeatureVector ();

  /**
   *  The EntryGroup that contains the Features (passed to the constructor).
   **/
  private EntryGroup entry_group;
}
