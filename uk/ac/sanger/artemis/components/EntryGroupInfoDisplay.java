/* EntryGroupInfoDisplay.java
 *
 * created: Fri Mar 12 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/EntryGroupInfoDisplay.java,v 1.2 2004-12-14 15:53:11 tjc Exp $
 **/

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.util.StringVector;
import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.Location;
import uk.ac.sanger.artemis.io.FuzzyRange;
import uk.ac.sanger.artemis.io.RangeVector;
import uk.ac.sanger.artemis.io.InvalidRelationException;

import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.sequence.*;

import java.util.Hashtable;
import java.util.Enumeration;
import java.awt.event.*;

import javax.swing.JFrame;

/**
 *  This component will show general information about an EntryGroup.  It will
 *  show the sequence length, GC content and a summary of the active entries.
 *
 *  @author Kim Rutherford
 *  @version $Id: EntryGroupInfoDisplay.java,v 1.2 2004-12-14 15:53:11 tjc Exp $
 **/

public class EntryGroupInfoDisplay
    implements FeatureChangeListener, EntryChangeListener,
               EntryGroupChangeListener, SequenceChangeListener 
{
  /**
   *  Summarize both strands.
   **/
  final static public int BOTH = 3;

  /**
   *  Create a new EntryGroupInfoDisplay object to summarize all features on
   *  both strands of the given EntryGroup.
   *  @param parent_frame The reference of the parent JFrame.
   *  @param entry_group The EntryGroup to view.
   **/
  public EntryGroupInfoDisplay(final JFrame parent_frame,
                               final EntryGroup entry_group)
  {
    this(parent_frame, entry_group, BOTH);
  }

  /**
   *  Create a new EntryGroupInfoDisplay object to summarize the given
   *  EntryGroup.
   *  @param parent_frame The reference of the parent JFrame.
   *  @param entry_group The EntryGroup to view.
   *  @param strand_flag Indicates which strand to summarize.  BOTH means
   *    summarize all features on both strands.  FORWARD means summarize
   *    forward strand only.  REVERSE means summarize reverse strand only.
   **/
  public EntryGroupInfoDisplay(final JFrame parent_frame,
                               EntryGroup entry_group,
                               final int strand_flag) 
  {
    this.strand_flag = strand_flag;
    this.parent_frame = parent_frame;
    this.entry_group = entry_group;

    if(strand_flag == Bases.FORWARD) 
    {
      final FeaturePredicate forward_feature_predicate =
        new FeaturePredicate() 
      {
        public boolean testPredicate(final Feature test_feature) 
        {
          return test_feature.isForwardFeature();
        }
      };

      final FilteredEntryGroup filtered_entry_group =
        new FilteredEntryGroup(entry_group, forward_feature_predicate,
                               "forward strand features");

      entry_group = filtered_entry_group;
    }
    else
    {
      if(strand_flag == Bases.REVERSE)
      {
        final FeaturePredicate reverse_feature_predicate =
          new FeaturePredicate()
        {
          public boolean testPredicate(final Feature test_feature) 
          {
            return !test_feature.isForwardFeature();
          }
        };

        final FilteredEntryGroup filtered_entry_group =
          new FilteredEntryGroup(entry_group, reverse_feature_predicate,
                                 "reverse strand features");

        entry_group = filtered_entry_group;
      } 
      else
      {
        if(strand_flag != BOTH) 
          throw new Error("internal error - illegal argument");
      }
    }

    file_viewer = new FileViewer(getFrameName());
    updateView();

    entry_group.addEntryChangeListener(this);
    entry_group.addFeatureChangeListener(this);
    entry_group.addEntryGroupChangeListener(this);
    entry_group.getBases().addSequenceChangeListener(this,
                                                     Bases.MAX_PRIORITY);

    file_viewer.addWindowListener(new WindowAdapter() 
    {
      public void windowClosed(WindowEvent event) 
      {
        stopListening();
      }
    });
  }

  /**
   *  Return the String that should be used for the title of the Overview Frame
   **/
  private String getFrameName() 
  {
    if(entry_group instanceof FilteredEntryGroup) 
    {
      final FilteredEntryGroup filtered_entry_group =
                                 (FilteredEntryGroup)entry_group;

      final String filter_name = filtered_entry_group.getFilterName();

      if(filter_name == null)
        return "Artemis overview of: " + parent_frame.getTitle();
      else
        return "Artemis overview of: " + parent_frame.getTitle() +
          " (" + filter_name + ")";
    } 
    else 
    {
      final StringBuffer buffer = new StringBuffer("Artemis Overview");

      if(entry_group.size() > 0) 
      {
        buffer.append(" of:");

        for(int i = 0 ; i < entry_group.size() ; ++i) 
        {
          final Entry this_entry = entry_group.elementAt(i); 

          if(this_entry.getName() == null) 
            buffer.append(" [no name]");
          else
            buffer.append(" " + this_entry.getName());
        }
      }

      return buffer.toString();
    }
  }
  
  /**
   *  Remove this object as a selection change listener.
   **/
  private void stopListening() 
  {
    entry_group.removeEntryChangeListener(this);
    entry_group.removeEntryGroupChangeListener(this);
    entry_group.removeFeatureChangeListener(this);
    entry_group.getBases().removeSequenceChangeListener(this);
  }

  /**
   *  Implementation of the FeatureChangeListener interface.   We listen to
   *  FeatureChange events so that we can update the display if qualifiers
   *  change.
   **/
  public void featureChanged(final FeatureChangeEvent event) 
  {
    updateView();
  }

  /**
   *  Implementation of the EntryGroupChangeListener interface.  We listen to
   *  EntryGroupChange events so that we can update the display if entries
   *  are added, deleted, activated or deactivated.
   **/
  public void entryGroupChanged(final EntryGroupChangeEvent event) 
  {
    if(event.getType() == EntryGroupChangeEvent.DONE_GONE) 
    {
      stopListening();
      file_viewer.dispose();
    }
    else
      updateView();
  }

  /**
   *  Implementation of the EntryChangeListener interface.  We listen to
   *  EntryChange events so that we can update the display if features are
   *  added or deleted.
   **/
  public void entryChanged(final EntryChangeEvent event) 
  {
    updateView();
  }

  /**
   *  Implementation of the SequenceChangeListener interface.  We listen to
   *  SequenceChange events so that we can update the display if the sequence
   *  changes.
   **/
  public void sequenceChanged(final SequenceChangeEvent event) 
  {
    updateView();
  }

  /**
   *  Update the FileViewer component to reflect the current state of the
   *  EntryGroup.
   **/
  private void updateView() 
  {
    final String frame_name = getFrameName();
    file_viewer.setTitle(frame_name);
    
    final StringBuffer buffer = new StringBuffer();

    buffer.append(frame_name + "\n\n");
    buffer.append("Number of bases: " +
                  entry_group.getSequenceLength() + "\n");
    buffer.append("Number of features in the active entries: " +
                  entry_group.getAllFeaturesCount() + "\n\n");

    final Bases entry_group_bases = entry_group.getBases();

    // cds_count - gene_count = pseudo gene count
    int gene_count = 0;
    int cds_count  = 0;
    int spliced_gene_count = 0;
    int spliced_gene_bases_count = 0;

    StringBuffer gene_bases_buffer = new StringBuffer();
    int gene_bases_count_with_introns = 0;
    int cds_bases_count = 0;
    int exon_count      = 0;
    int intron_count    = 0;
    int partial_count   = 0;

    final Hashtable table = new Hashtable();
    final FeatureEnumeration feature_enumerator = entry_group.features();

    while(feature_enumerator.hasMoreFeatures())
    {
      final Feature this_feature = feature_enumerator.nextFeature();

      final Key key = this_feature.getKey();
      final String key_string = key.toString();

      final Location location  = this_feature.getLocation();
      final RangeVector ranges = location.getRanges();

      for(int i = 0; i < ranges.size(); ++i) 
        if(ranges.elementAt(i) instanceof FuzzyRange)
          ++partial_count;

      try 
      {
        String colour = this_feature.getValueOfQualifier("colour");

        if(colour == null || colour.length() == 0)
          colour = "no colour";

        if(table.containsKey(key_string)) 
        {
          final Hashtable colour_table = (Hashtable)table.get(key_string);
          final Integer colour_value   = (Integer)colour_table.get(colour);

          if(colour_value == null) 
              colour_table.put(colour, new Integer(1));
          else 
          {
            final int old_value = ((Integer)colour_value).intValue();
            colour_table.put(colour, new Integer(old_value + 1));
          }
        } 
        else
        {
          final Hashtable colour_table = new Hashtable();
          colour_table.put(colour, new Integer(1));
          table.put(key_string, colour_table);
        }
      } 
      catch(InvalidRelationException e) 
      {
        throw new Error("internal error - unexpected exception: " + e);
      }

      if(this_feature.isCDS()) 
      {
        cds_count++;
        cds_bases_count += this_feature.getBaseCount();

        final Qualifier pseudo_qualifier;

        try 
        {
          pseudo_qualifier = this_feature.getQualifierByName("pseudo");
        } 
        catch(InvalidRelationException e) 
        {
          throw new Error("internal error - unexpected exception: " + e);
        }

        if(pseudo_qualifier == null) 
        {
          ++gene_count;

          final FeatureSegmentVector segments = this_feature.getSegments();

          if(segments.size() > 1)
          {
            ++spliced_gene_count;
            spliced_gene_bases_count += this_feature.getBaseCount();
            intron_count += segments.size() - 1;
          }

          exon_count += segments.size();

          gene_bases_buffer.append(this_feature.getBases());

          gene_bases_count_with_introns += this_feature.getRawLastBase() -
                                           this_feature.getRawFirstBase() + 1;
        }
      }
    }

    final String gene_bases    = gene_bases_buffer.toString();
    final int gene_bases_count = gene_bases.length();

    final int pseudo_gene_count = cds_count - gene_count;
    final int pseudo_gene_bases_count = cds_bases_count - gene_bases_count;

    if(cds_count > 0) 
    {
      if(gene_count > 0) 
      {
        buffer.append("Genes (CDS features without a /pseudo qualifier):\n");

        final int non_spliced_gene_count = gene_count - spliced_gene_count;
        final int non_spliced_gene_bases_count =
          gene_bases_count - spliced_gene_bases_count;

        if(spliced_gene_count > 0) 
        {
          buffer.append("   spliced:\n");
          buffer.append("      count: " + spliced_gene_count + "\n");
          buffer.append("      bases: " + spliced_gene_bases_count + "\n");
          buffer.append("      introns: " + intron_count + "\n");
        }

        if(non_spliced_gene_count > 0) 
        {
          buffer.append("   non-spliced:\n");
          buffer.append("      count: " + non_spliced_gene_count + "\n");
          buffer.append("      bases: " + non_spliced_gene_bases_count +
                        "\n");
        }

        buffer.append("   all:\n");
        buffer.append("      count: " + gene_count + "\n");
        buffer.append("      partials: " + partial_count + "\n");
        if(exon_count == gene_count) 
        {
          buffer.append("      bases: "
                        + gene_bases_count + "\n");
        } 
        else
        {
          buffer.append("      bases(excluding introns): "
                        + gene_bases_count + "\n");
          buffer.append("      bases(including introns): "
                        + gene_bases_count_with_introns + "\n");
          buffer.append("      exons: " + exon_count + "\n");
          buffer.append("      average exon length: " +
                        10L * gene_bases_count / exon_count / 10.0 + "\n");
          buffer.append("      average intron length: " +
                        10L * (gene_bases_count_with_introns -
                               gene_bases_count) /
                        (exon_count - gene_count) / 10.0 + "\n");
          buffer.append("      average number of exons per gene: " +
                        100L * exon_count / gene_count / 100.00 + "\n");
        }
        buffer.append("      density: " +
                      1000000L * gene_count /
                      entry_group.getSequenceLength() / 1000.0 +
                      " genes per kb   (" +
                      entry_group.getSequenceLength() / gene_count +
                       " bases per gene)\n");
        buffer.append("      average length: " +
                      gene_bases_count / gene_count + "\n");
        buffer.append("      average length (including introns): " +
                      gene_bases_count_with_introns / gene_count + "\n");
        buffer.append("      coding percentage: " +
                      1000L * gene_bases_count /
                      entry_group.getSequenceLength() / 10.0 + "\n");
        buffer.append("      coding percentage (including introns): " +
                      1000L * gene_bases_count_with_introns /
                      entry_group.getSequenceLength() / 10.0 + "\n\n");

        buffer.append("      gene sequence composition:\n\n");

        final StringVector gene_base_summary =
          SelectionViewer.getBaseSummary(gene_bases);

        for(int i = 0 ; i < gene_base_summary.size() ; ++i)
        {
          buffer.append("         ");
          buffer.append(gene_base_summary.elementAt(i)).append("\n");
        }

        buffer.append("\n");
      }

      if(pseudo_gene_count > 0)
      {
        buffer.append("Pseudo genes (CDS features with a /pseudo " +
                      "qualifier):\n");

        buffer.append("   count: " + pseudo_gene_count + "\n");
        buffer.append("   bases: " + pseudo_gene_bases_count + "\n");
        buffer.append("   average length: " +
                      pseudo_gene_bases_count / pseudo_gene_count + "\n\n");
      }

      if(pseudo_gene_count > 0) 
      {
        buffer.append("All CDS features:\n");
        buffer.append("   count: " + cds_count + "\n");
        buffer.append("   bases: " + cds_bases_count + "\n");
        buffer.append("   average length: " +
                      cds_bases_count / cds_count + "\n\n");
      }
    }

    final Strand strand;

    if(strand_flag == Bases.FORWARD || strand_flag == BOTH) 
      strand = entry_group.getBases().getForwardStrand();
    else
      strand = entry_group.getBases().getReverseStrand();

    final StringVector base_summary =
      SelectionViewer.getBaseSummary(strand.getStrandBases());

    buffer.append("\nOverall sequence composition:\n\n");

    for(int i = 0 ; i < base_summary.size() ; ++i) 
      buffer.append(base_summary.elementAt(i)).append("\n");

    buffer.append("\nSummary of the active entries:\n");

    final Enumeration e = table.keys();

    while(e.hasMoreElements()) 
    {
      final String this_key = (String)e.nextElement();
      final Hashtable colour_table = (Hashtable)table.get(this_key);

      buffer.append(this_key + ": ");

      final StringBuffer colour_string = new StringBuffer();

      final Enumeration colour_enum = colour_table.keys();

      int total = 0;

      while(colour_enum.hasMoreElements()) 
      {
        final String this_colour = (String) colour_enum.nextElement();

        final int colour_count =
          ((Integer) colour_table.get(this_colour)).intValue();

        total += colour_count;

        final String end_string;

        if(this_colour.equals("no colour")) 
          end_string = "no colour";
        else
          end_string = "colour: " + this_colour;

        if(colour_count == 1) 
          colour_string.append("  one has " + end_string + "\n");
        else
          colour_string.append("  " + colour_count + " have " +
                               end_string + "\n");
      }

      buffer.append(total + "\n");
      buffer.append(colour_string);
    }

    file_viewer.setText(buffer.toString());
  }

  /**
   *  This is the EntryGroup object that we are viewing.
   **/
  private EntryGroup entry_group;

  /**
   *  The strand indicator that was passed to the constructor.
   **/
  private int strand_flag;

  /**
   *  The FileViewer object that is displaying the EntryGroup.
   **/
  private FileViewer file_viewer;

  /**
   *  The Frame that was passed to the constructor.
   **/
  private JFrame parent_frame;
}
