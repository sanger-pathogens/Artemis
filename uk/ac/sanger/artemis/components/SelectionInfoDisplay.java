/* SelectionInfoDisplay.java
 *
 * created: Tue Dec 15 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/SelectionInfoDisplay.java,v 1.7 2005-01-06 11:21:06 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.sequence.*;

import uk.ac.sanger.artemis.io.StreamQualifier;
import uk.ac.sanger.artemis.io.QualifierInfo;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.util.OutOfRangeException;

import java.awt.*;
import java.awt.event.*;
import java.util.Vector;

import javax.swing.*;

/**
 *  This class displays information about the selection in a Label.
 *
 *  @author Kim Rutherford
 *  @version $Id: SelectionInfoDisplay.java,v 1.7 2005-01-06 11:21:06 tjc Exp $
 **/
public class SelectionInfoDisplay extends CanvasPanel
    implements SelectionChangeListener 
{
  /**
   *  Create a new SelectionInfoDisplay component.
   *  @param entry_group The EntryGroup that this component will display.
   *  @param selection The Selection that this component will display.
   **/
  public SelectionInfoDisplay(final EntryGroup entry_group,
                              final Selection selection) 
  {
    this.entry_group = entry_group;
    this.selection = selection;

    getSelection().addSelectionChangeListener(this);

    addComponentListener(new ComponentAdapter() 
    {
      public void componentResized(ComponentEvent e) 
      {
        setSize(getSize().width,
                getFontHeight()+1);
        repaint();
      }
      public void componentShown(ComponentEvent e) 
      {
        setSize(getSize().width,
                getFontHeight()+1);
        repaint();
      }
    });

    setBackground(new Color(230,230,230));
    setSize(80, getFontHeight());
  }

  /**
   *
   **/
  public Dimension getPreferredSize() 
  {
    return new Dimension(10, getFontHeight());
  }

  /**
   *
   **/
  public Dimension getMinimumSize() 
  {
    return new Dimension(10, getFontHeight());
  }

  /**
   *  Implementation of the SelectionChangeListener interface.  We listen to
   *  SelectionChange events so that we can update the list to reflect the
   *  current selection.
   **/
  public void selectionChanged(SelectionChangeEvent event) 
  {
    repaint();
  }

  /**
   *  Draw the label.
   **/
  public void paintComponent(final Graphics g) 
  {
    super.paintComponent(g);
    if(!isVisible ()) 
      return;

    final FeatureVector features = getSelection().getAllFeatures();
    final StringBuffer new_text  = new StringBuffer();

    new_text.append(markerRangeText(getSelection(), entry_group));

    int base_total = 0;
    int aa_total = 0;
    boolean saw_a_non_cds = false;

    if(features.size () > 0) 
    {
      final StringBuffer feature_names = new StringBuffer ();

      if(features.size () < 100) 
      {
        if(features.size () > 1)
          feature_names.append ("  (");

        // show up to 10 features
        for(int i = 0 ; i < features.size () ; ++i) 
        {
          final Feature current_feature = features.elementAt (i);

          base_total += current_feature.getBaseCount ();
          aa_total   += current_feature.getAACount ();

          if(!current_feature.getKey().equals ("CDS")) 
            saw_a_non_cds = true;

          if(i < 10) 
          {
            if(i != 0) 
              feature_names.append(' ');
            
            feature_names.append(current_feature.getIDString ());
          }
        }

        if(features.size() > 10) 
          feature_names.append ("...");

        if(features.size() == 1) 
        {
          // only one feature so append some qualifiers
          feature_names.append("  (");
          feature_names.append(getQualifierString(features.elementAt (0)));
        }
        feature_names.append(")");
      }

      if(features.size() == 1) 
        new_text.append ("Selected feature:  ");
      else 
        new_text.append (features.size () + " selected features  ");

      if(features.size() < 100) 
      {
        // show a count of the number of bases and amino acids in the selected
        // feature
        if(features.size() > 1) 
        {
          new_text.append ("total bases " + base_total);

          if(!saw_a_non_cds) 
            new_text.append ("  total amino acids " + aa_total);
        }
        else
        {
          if(features.size() == 1) 
          {
            new_text.append ("bases " + base_total);
            if(!saw_a_non_cds) 
              new_text.append ("  amino acids " + aa_total);
          }
        }
      }

      new_text.append("  ");
      new_text.append(feature_names.toString());
    }


    String text = new_text.toString();

    if(text.length () > 0) 
    {
      if (text.length () > 150) 
      {
        // call substring () just to keep the String a sensible length
        text = text.substring (0, 150);
      }
    } 
    else 
      text = "Nothing selected";

//  g.setColor(new Color (230, 230, 230));
//  g.fillRect(0, 0, getCanvasWidth(), getCanvasHeight());

    g.setColor(Color.black);
    g.drawString(text, 2, getFontMaxAscent() + 1);
  }

  /**
   *  Return a String containing the qualifiers (except the /note) of given
   *  feature.
   **/
  private String getQualifierString(final Feature feature) 
  {
    final QualifierVector qualifiers = feature.getQualifiers();

    // if we see a /note or /fasta_file it gets saved for later
    final Vector saved_qualifiers = new Vector();

    final StringBuffer string_buffer = new StringBuffer();

    final EntryInformation entry_information =
      feature.getEntry().getEntryInformation();

    for (int i = 0 ; i < qualifiers.size () ; ++i) 
    {
      final Qualifier this_qualifier = (Qualifier)qualifiers.elementAt (i);

      if (this_qualifier.getName ().equals ("note") ||
          this_qualifier.getName ().endsWith ("_file")) 
      {
        // save the qualifier and put it last
        saved_qualifiers.addElement (this_qualifier);
        continue;
      }

      if (string_buffer.length () > 0)
      {
        // put spaces between the qualifiers
        string_buffer.append (' ');
      }

      final QualifierInfo qualifier_info =
        entry_information.getQualifierInfo (this_qualifier.getName ());

      string_buffer.append (StreamQualifier.toString (qualifier_info,
                                                      this_qualifier));
    }

    for(int i = 0 ; i < saved_qualifiers.size () ; ++i) 
    {
      if(string_buffer.length () > 0) 
      {
        // put spaces between the qualifiers
        string_buffer.append (' ');
      }

      final Qualifier this_qualifier =
        (Qualifier) saved_qualifiers.elementAt (i);

      final QualifierInfo qualifier_info =
        entry_information.getQualifierInfo (this_qualifier.getName ());

      string_buffer.append (StreamQualifier.toString (qualifier_info,
                                                      this_qualifier));
    }

    return string_buffer.toString ();
  }

  /**
   *  If the selection contains a MarkerRange (a range of bases) then return a
   *  String in this form: "Selected Bases on forward strand: <start>..<end> ",
   *  otherwise return "".  This method is also used by the SelectionViewer
   *  class.
   *  @param current_selection The selection to summarize.
   *  @param entry_group The EntryGroup that the selection refers to.
   **/
  static String markerRangeText(final Selection current_selection,
                                final EntryGroup entry_group) 
  {
    final MarkerRange marker_range = current_selection.getMarkerRange ();

    if (marker_range == null) 
      return "";
    else
    {
      final StringBuffer buffer = new StringBuffer ();

      final int start_pos = marker_range.getStart ().getPosition ();
      final int end_pos = marker_range.getEnd ().getPosition ();

      if(marker_range.getStrand ().isForwardStrand ()) 
      {
        if(start_pos == end_pos) 
        {
          buffer.append ("One selected base on forward strand: " +
                         start_pos + "  ");
        }
        else 
        {
          buffer.append ((end_pos - start_pos + 1) +
                         " selected bases on forward strand: " +
                         start_pos + ".." + end_pos + " ");
        }
      }
      else 
      {
        if(start_pos == end_pos) 
        {
          buffer.append ("One selected base on reverse strand: " +
                         marker_range.getStart ().getPosition () +
                         "  = complement (" +
                         marker_range.getEnd ().getRawPosition () + ") ");
        }
        else 
        {
          buffer.append ((end_pos - start_pos + 1) +
                         " selected bases on reverse strand: " +
                         marker_range.getStart ().getPosition () + ".." +
                         marker_range.getEnd ().getPosition () +
                         "  = complement (" +
                         marker_range.getEnd ().getRawPosition () + ".." +
                         marker_range.getStart ().getRawPosition () + ") ");
        }
      }

      if(marker_range.getCount () >= 3) 
      {
        // check if this codon is in a feature

        final Range raw_range = marker_range.getRawRange ();

        final FeatureVector features;

        try 
        {
          features = entry_group.getFeaturesInRange (raw_range);
        } 
        catch (OutOfRangeException e) 
        {
          return "illegal marker range";
        }

        for(int i = 0 ; i < features.size () ; ++i)
        {
          final Feature this_feature = features.elementAt (i);

          if(!this_feature.isProteinFeature ()) 
          {
            // only display protein features
            continue;
          }

          if(marker_range.isForwardMarker () !=
             this_feature.isForwardFeature ()) 
          {
            // only display feature positions in features on the same strand
            // as the selected bases
            continue;
          }

          final String this_feature_id = this_feature.getIDString ();

          final Marker start_marker = marker_range.getStart ();

          final Marker end_marker = marker_range.getEnd ();

          final int start_position_in_feature =
            this_feature.getFeaturePositionFromMarker (start_marker);

          final int start_codon_position =
            start_position_in_feature - this_feature.getCodonStart () + 1;

          String start_string = null;

          if(start_position_in_feature != -1) 
          {
            if (start_codon_position % 3 == 0) 
            {
              start_string = "codon " + (start_codon_position / 3 + 1);
            }
          }

          String end_string = null;

          final int end_position_in_feature =
            this_feature.getFeaturePositionFromMarker (end_marker);

          final int end_codon_position =
            end_position_in_feature - this_feature.getCodonStart () - 1;

          if (end_position_in_feature != -1) 
          {
            if (start_codon_position != end_codon_position &&
                end_codon_position % 3 == 0)
            {
              end_string = "codon " + (end_codon_position / 3 + 1);
            }
          }

          if(!(start_string == null && end_string == null)) 
          {
            buffer.append (" (");

            if(start_string != null) 
              buffer.append (start_string);

            if(start_string != null && end_string != null) 
              buffer.append (" to ");

            if(end_string != null) 
              buffer.append (end_string);

            buffer.append(" in feature " + this_feature.getIDString () +
                          ") ");
          }
        }
      }

      return buffer.append(" ").toString();
    }
  }

  /**
   *  Return the Selection reference that was passed to the constructor.
   **/
  private Selection getSelection () 
  {
    return selection;
  }

  /**
   *  Return the Selection object that was passed to the constructor.
   **/
  private EntryGroup getEntryGroup () 
  {
    return entry_group;
  }

  /**
   *  The reference of the EntryGroup object that was passed to the
   *  constructor.
   **/
  private EntryGroup entry_group;

  /**
   *  This is a reference to the Selection object that created this
   *  component.
   **/
  private final Selection selection;
}

