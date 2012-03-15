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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/AddMenu.java,v 1.42 2009-06-01 09:49:07 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;

import uk.ac.sanger.artemis.plot.CodonUsageAlgorithm;

import uk.ac.sanger.artemis.sequence.BasePattern;
import uk.ac.sanger.artemis.sequence.BasePatternFormatException;
import uk.ac.sanger.artemis.sequence.Bases;
import uk.ac.sanger.artemis.sequence.MarkerRange;
import uk.ac.sanger.artemis.sequence.MarkerRangeVector;
import uk.ac.sanger.artemis.sequence.Strand;
import uk.ac.sanger.artemis.util.*;
import uk.ac.sanger.artemis.components.genebuilder.GeneUtils;
import uk.ac.sanger.artemis.io.ChadoCanonicalGene;
import uk.ac.sanger.artemis.io.DatabaseDocumentEntry;
import uk.ac.sanger.artemis.io.GFFDocumentEntry;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.InvalidRelationException;
import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.io.RangeVector;
import uk.ac.sanger.artemis.io.Location;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.LocationParseException;
import uk.ac.sanger.artemis.io.EntryInformationException;

import java.awt.Cursor;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.Comparator;
import java.util.Vector;
import java.util.Enumeration;
import java.util.regex.Pattern;

import javax.swing.Box;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JTextField;
import javax.swing.KeyStroke;


/**
 *  A Menu with commands that add new features/entries to an EntryGroup.  This
 *  should have been called CreateMenu.
 *
 *  @author Kim Rutherford
 *  @version $Id: AddMenu.java,v 1.42 2009-06-01 09:49:07 tjc Exp $
 **/
public class AddMenu extends SelectionMenu 
{

  private static final long serialVersionUID = 1L;

  /** The GotoEventSource object that was passed to the constructor. */
  private GotoEventSource goto_event_source = null;

  /** The EntryGroup object that was passed to the constructor. */
  private EntryGroup entry_group;

  private BasePlotGroup base_plot_group;
  
  /**
   *  The shortcut for "Create From Base Range".
   **/
  final static KeyStroke CREATE_FROM_BASE_RANGE_KEY =
    KeyStroke.getKeyStroke (KeyEvent.VK_C,
                            Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()); //InputEvent.CTRL_MASK);

  final static public int CREATE_FROM_BASE_RANGE_KEY_CODE = KeyEvent.VK_C;

  /** busy cursor */
  private Cursor cbusy = new Cursor(Cursor.WAIT_CURSOR);
  /** done cursor */
  private Cursor cdone = new Cursor(Cursor.DEFAULT_CURSOR);

  /*private AlignmentViewer alignQueryViewer;
  private AlignmentViewer alignSubjectViewer;*/

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
   *  @param base_plot_group The AlignmentViewer associated with this JMenu
   *
   *  @param menu_name The name of the new menu.
   **/
  public AddMenu (final JFrame frame,
                  final Selection selection,
                  final EntryGroup entry_group,
                  final GotoEventSource goto_event_source,
                  final BasePlotGroup base_plot_group,
                  final String menu_name) 
  {
    this(frame,selection,entry_group,
         goto_event_source,base_plot_group,null,null,menu_name);
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
   *  @param menu_name The name of the new menu.
   **/
  public AddMenu(final JFrame frame,
                 final Selection selection,
                 final EntryGroup entry_group,
                 final GotoEventSource goto_event_source,
                 final BasePlotGroup base_plot_group,
                 final AlignmentViewer alignQueryViewer, 
                 final AlignmentViewer alignSubjectViewer,
                 final String menu_name)
  {
    super (frame, menu_name, selection);

    /*this.alignQueryViewer   = alignQueryViewer;
    this.alignSubjectViewer = alignSubjectViewer;*/
    this.entry_group = entry_group;
    this.base_plot_group = base_plot_group;

    final JMenuItem new_feature_item = new JMenuItem ("New Feature");
    new_feature_item.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent event) {
        makeNewFeature ();
      }
    });

    add (new_feature_item);

   
    final JMenuItem create_feature_from_range_item =
      new JMenuItem ("Feature From Base Range");
    
    if(!GeneUtils.isDatabaseEntry(entry_group))
      create_feature_from_range_item.setAccelerator (CREATE_FROM_BASE_RANGE_KEY);
    create_feature_from_range_item.addActionListener(new ActionListener () {
      public void actionPerformed (ActionEvent event) {
        createFeatureFromBaseRange (getParentFrame (), getSelection (),
                                    entry_group, getGotoEventSource ());
      }
    });
    add (create_feature_from_range_item);

    final JMenuItem create_gene_model_from_range_item = new JMenuItem(
          "Gene Model From Base Range");  
    if(GeneUtils.isDatabaseEntry(entry_group))
      create_gene_model_from_range_item.setAccelerator(CREATE_FROM_BASE_RANGE_KEY);
    create_gene_model_from_range_item.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)
      {
        final Entry default_entry = entry_group.getDefaultEntry();
        if(default_entry == null)
        {
          new MessageDialog(frame, "There is no default entry.");
          return;
        }
        final uk.ac.sanger.artemis.io.Entry entry = default_entry.getEMBLEntry();
        if( !(entry instanceof GFFDocumentEntry || entry instanceof DatabaseDocumentEntry) )
        {
          new MessageDialog(frame, 
              "Expecting a GFF entry. The default entry "+
              default_entry.getName()+" does not support this.");
          return;
        }
        entry_group.getActionController ().startAction ();
        GeneUtils.createGeneModel(getParentFrame(), getSelection(),
            entry_group, getGotoEventSource());
        entry_group.getActionController ().endAction ();
      }
    });
    add (create_gene_model_from_range_item);   

    if(alignQueryViewer != null || alignSubjectViewer != null)
    {
      JMenuItem create_difference_feature  =
        new JMenuItem("Features From Non-matching Regions");
      create_difference_feature.addActionListener(new ActionListener()
      {
        public void actionPerformed (ActionEvent event) 
        {
          frame.setCursor(cbusy);
          Vector diffs = null;
          String comparisonNote = "";
          if(alignQueryViewer == null || alignSubjectViewer == null)
          {
            final Entry sequence_entry;
            if(alignQueryViewer != null)
            {
              diffs = alignQueryViewer.getDifferenceCoords(false);
              sequence_entry = alignQueryViewer.getSubjectEntryGroup().getSequenceEntry();
            }
            else
            {
              diffs = alignSubjectViewer.getDifferenceCoords(true);
              sequence_entry = alignSubjectViewer.getQueryEntryGroup().getSequenceEntry();
            }

            comparisonNote = comparisonNote + sequence_entry.getName();
          }
          else    // multi-comparison
          {
            Vector diffs1 = alignQueryViewer.getDifferenceCoords(false);
            Vector diffs2 = alignSubjectViewer.getDifferenceCoords(true);
          
            Entry sequence_entry;
            sequence_entry = alignQueryViewer.getSubjectEntryGroup().getSequenceEntry();
            comparisonNote = comparisonNote + sequence_entry.getName();

            sequence_entry = alignSubjectViewer.getQueryEntryGroup().getSequenceEntry();
            comparisonNote = comparisonNote + " and " + sequence_entry.getName();
  
            diffs = union(diffs1,diffs2);
          }
          
          createFeatures(diffs, frame, comparisonNote);
          frame.setCursor(cdone);
        }
      });

      add (create_difference_feature);
    }

    final JMenuItem create_intron_features_item =
      new JMenuItem ("Intron Features");
    create_intron_features_item.addActionListener(new ActionListener () {
      public void actionPerformed (ActionEvent event) {
        createIntronFeatures (getParentFrame (), getSelection (),
                              entry_group);
      }
    });
    
    add (create_intron_features_item);
    if(GeneUtils.isDatabaseEntry(entry_group))
      create_intron_features_item.setEnabled(false);
    
    final JMenuItem create_intergenic_features_item =
      new JMenuItem ("Intergenic Features");
    create_intergenic_features_item.addActionListener(new ActionListener () {
      public void actionPerformed (ActionEvent event) {
        createIntergenicFeatures(getParentFrame (), entry_group);
      }
    });
    add (create_intergenic_features_item);

    final JMenuItem create_exon_features_item =
      new JMenuItem ("Exon Features");
    create_exon_features_item.addActionListener(new ActionListener () {
      public void actionPerformed (ActionEvent event) {
        createExonFeatures (getParentFrame (), getSelection (),
                            entry_group);
      }
    });

    add (create_exon_features_item);
    if(GeneUtils.isDatabaseEntry(entry_group))
      create_exon_features_item.setEnabled(false);

    final JMenuItem create_gene_features_item =
      new JMenuItem ("Gene Features");
    create_gene_features_item.addActionListener(new ActionListener () {
      public void actionPerformed (ActionEvent event) {
        createGeneFeatures (getParentFrame (), getSelection (),
                            entry_group);
      }
    });

    add (create_gene_features_item);
    if(GeneUtils.isDatabaseEntry(entry_group))
      create_gene_features_item.setEnabled(false);

    addSeparator ();

    final JMenuItem new_entry_item = new JMenuItem ("New Entry");
    new_entry_item.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent event) {
        makeNewEntry ();
      }
    });

    add (new_entry_item);

    addSeparator ();

    final JMenuItem mark_orfs_with_size_item = new JMenuItem ("Mark Open Reading Frames ...");

    mark_orfs_with_size_item.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent event) {
        markORFSWithSize (false, frame);
      }
    });

    add (mark_orfs_with_size_item);


    final JMenuItem mark_empty_orfs_with_size_item = new JMenuItem ("Mark Empty ORFs ...");

    mark_empty_orfs_with_size_item.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent event) {
        markORFSWithSize (true, frame);
      }
    });

    add (mark_empty_orfs_with_size_item);


    final JMenuItem mark_orfs_range_item = new JMenuItem ("Mark ORFs In Range ...");
    mark_orfs_range_item.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent event) {
        markOpenReadingFramesInRange ();
      }
    });

    add (mark_orfs_range_item);


    final JMenuItem mark_pattern_item = new JMenuItem ("Mark From Pattern ...");
    mark_pattern_item.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent event) {
        makeFeaturesFromPattern ();
      }
    });

    add (mark_pattern_item);

    final JMenuItem mark_ambiguities_item = new JMenuItem ("Mark Ambiguities");
    mark_ambiguities_item.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent event) {
        markAmbiguities ();
      }
    });

    add (mark_ambiguities_item);
  }

  /**
  *
  * Find the union of coordinates in two Vecor objects.
  * @param v1 Vector of Integer coordinates
  * @param v2 Vector of Integer coordinates
  *  
  */
  protected static Vector union(final Vector v1, final Vector v2)
  {
    final Vector union = new Vector();

    for(int i=0; i<v1.size(); i++)
    {
      Integer[] imatch = (Integer[])v1.get(i);
      int istart = imatch[0].intValue();
      int iend   = imatch[1].intValue();

      for(int j=0; j<v2.size(); j++)
      {
        Integer[] jmatch = (Integer[])v2.get(j);
        int jstart = jmatch[0].intValue();
        int jend   = jmatch[1].intValue();

        if( (istart >= jstart && istart <= jend) ||
            (iend >= jstart && iend <= jend) ||
            (jstart > istart && jend < iend) )
        {
          if(jstart > istart)
            istart = jstart;
          
          if(iend < jend)
            jend = iend;

          final Integer coords[] = new Integer[2];
          coords[0] = new Integer(istart);
          coords[1] = new Integer(jend);
          union.add(coords);
        }
      }
    }

    return union;
  }

  /**
  *
  * Create features from Vector of coordinates.
  * @param diffs Vector of coordinates to create feature from.
  * @param frame The JFrame that owns this JMenu.
  *
  */
  private void createFeatures(Vector diffs, JFrame frame, String name) 
  {
    Enumeration eDiffs = diffs.elements();
    while(eDiffs.hasMoreElements())
    {
      Integer coords[] = (Integer[])eDiffs.nextElement();
      int start = coords[0].intValue();
      int end   = coords[1].intValue();
//    System.out.println(start+" "+end);
         
      final Entry default_entry = entry_group.getDefaultEntry();
      if(default_entry == null) 
      {
        new MessageDialog(frame, "There is no default entry");
        return;
      }

      Location loc = null;
      Feature temp_feature;
      try 
      {
        loc = new Location(start+".."+end);
        Key misc_feature = new Key("misc_feature");
        temp_feature = default_entry.createFeature(misc_feature, loc);
        Qualifier note = new Qualifier("note",
                                       "Automatically generated region of difference with "+
                                       name);
        temp_feature.setQualifier(note);
      } 
      catch(EntryInformationException e) 
      {
        // use the default key instead
        final Key default_key =
          default_entry.getEntryInformation().getDefaultKey();

        try
        {
          temp_feature =
              default_entry.createFeature(default_key, loc);
        }
        catch(EntryInformationException einfo)
        {
          throw new Error("internal error - unexpected exception: " + einfo);
        }
        catch(ReadOnlyException eRead)
        {
          new MessageDialog(frame, "feature not created: " +
                          "the default entry is read only");
        }
        catch(OutOfRangeException eout)
        {
          throw new Error("internal error - unexpected exception: " + eout);
        }

      }
      catch(ReadOnlyException e) 
      {
        new MessageDialog(frame, "feature not created: " +
                        "the default entry is read only");
      } 
      catch(OutOfRangeException e) 
      {
        throw new Error("internal error - unexpected exception: " + e);
      }
      catch(LocationParseException  lpe)
      {
        throw new Error("internal error - unexpected exception: " + lpe);
      }
    }
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

          
          final JFrame edit_frame = new JFrame("Artemis Feature Edit: " + 
              new_feature.getIDString() +
              (new_feature.isReadOnly() ?
                  "  -  (read only)" :
                  ""));
          
          final FeatureEdit feature_edit = new FeatureEdit(new_feature, entry_group,
              getSelection(), getGotoEventSource(), edit_frame);
          
          edit_frame.addWindowListener(new WindowAdapter() 
          {
            public void windowClosing(WindowEvent event) 
            {
              feature_edit.stopListening();
              edit_frame.dispose();
            }
          });
          
          edit_frame.getContentPane().add(feature_edit);
          edit_frame.pack();

          //final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
          //edit_frame.setLocation(new Point((screen.width - edit_frame.getSize().width)/2,
          //                            (screen.height - edit_frame.getSize().height)/2));
          Utilities.centreFrame(edit_frame);
          
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

          edit_frame.setVisible(true);
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
      } else if(default_entry.isReadOnly()) {
        new MessageDialog (frame, "The default entry is read only");
        return;
      }
        

      try {
        final Location new_location = range.createLocation ();

        /*final*/ Feature temp_feature;
        QualifierVector qualifiers = null;
        final boolean isDatabaseEntry = GeneUtils.isDatabaseEntry(entry_group);
        
        if(isDatabaseEntry)
        {
          String uniquename = GeneUtils.promptForUniquename(entry_group, 
                                     range.isForwardMarker(), range.getRawRange());
          Qualifier qualifier = new Qualifier("ID", uniquename);
          qualifiers = new QualifierVector();
          qualifiers.add(qualifier);
        }
        
        try 
        {
          final Key key;
          if(isDatabaseEntry)
            key = new Key("region");
          else
            key = Key.CDS;
          
          if(qualifiers == null)
            temp_feature = default_entry.createFeature (key, new_location);
          else
            temp_feature = default_entry.createFeature (key, new_location, qualifiers);
          
        /*if(isDatabaseEntry)
          {
            final ChadoCanonicalGene chado_gene = new ChadoCanonicalGene();
            chado_gene.setGene(temp_feature.getEmblFeature());
            ((uk.ac.sanger.artemis.io.GFFStreamFeature)
              (temp_feature.getEmblFeature())).setChadoGene(chado_gene);
          }*/
        }
        catch (EntryInformationException e) 
        {
          // use the default key instead

          final Key default_key =
            default_entry.getEntryInformation ().getDefaultKey ();

          try 
          {
            if(qualifiers == null)
              temp_feature =
                default_entry.createFeature (default_key, new_location);
            else
              temp_feature =
                default_entry.createFeature (default_key, new_location, qualifiers);
          } 
          catch (EntryInformationException ex) 
          {
            throw new Error ("internal error - unexpected exception: " + ex);
          }
        }

        final Feature new_feature = temp_feature;

        selection.setMarkerRange (null);
        selection.set (new_feature);
        
        final ActionListener cancel_listener = new ActionListener() 
        {
          public void actionPerformed(ActionEvent e)
          {
            try 
            {
              new_feature.removeFromEntry ();
              selection.setMarkerRange(range);
            }
            catch (ReadOnlyException exception) 
            {
              throw new Error("internal error - unexpected exception: " +
                              exception);
            }
          }
        };
          
        EditMenu.editSelectedFeatures(entry_group, selection, goto_event_source,
            new_feature, cancel_listener, null);

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

        if (!(selection_feature.isProteinFeature () ||
              selection_feature.getKey().equals("5'UTR") ||
              selection_feature.getKey().equals("3'UTR"))) {
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

        int select = 0;
        for (int range_index = 0 ;
             range_index < cds_ranges.size () - 1 ;
             ++range_index) {
          final int end_of_range_1 =
            ((Range)cds_ranges.elementAt(range_index)).getEnd ();
          final int start_of_range_2 =
            ((Range)cds_ranges.elementAt(range_index + 1)).getStart ();

          if (end_of_range_1 > start_of_range_2) {
            // ignore - the exons overlap so there is no room for an intron
            continue;
          }

          Range new_range = null;

          try {
            new_range = new Range (end_of_range_1 + 1,
                                   start_of_range_2 - 1);
          } catch (OutOfRangeException e) {
            Object[] options = { "CANCEL", "IGNORE", "IGNORE ALL"};

            if(select != 2)
            {
              select = JOptionPane.showOptionDialog(null,
                          "Found overlapping CDS\n"+e,
                          "Out of Range",
                           JOptionPane.YES_NO_CANCEL_OPTION,
                           JOptionPane.WARNING_MESSAGE,
                           null,
                           options,
                           options[0]);
              if(select == 0)
                throw new Error ("internal error - unexpected exception: " + e);
             
              continue;
            }
          }

          final RangeVector intron_ranges = new RangeVector ();

          intron_ranges.add (new_range);

          final Key intron_key = new Key ("intron");
          final Location intron_location =
            new Location (intron_ranges, cds_location.isComplement ());
          final QualifierVector qualifiers = new QualifierVector ();

          try {
            StringVector sysNames = Options.getOptions().getSystematicQualifierNames();
            for(int i=0; i<sysNames.size(); i++) {
              Qualifier qual = selection_feature.getQualifierByName((String) sysNames.get(i));
              if(qual != null && qual.getValues() != null && qual.getValues().size() > 0) {
                qualifiers.addQualifierValues(qual);
                break;
              }
            }
          } catch (InvalidRelationException e) {
            throw new Error ("internal error - unexpected exception: " + e);
          }
          
          try {
            selection_feature.getEntry ().createFeature (intron_key,
                                                         intron_location,
                                                         qualifiers);
          } catch (ReadOnlyException e) {
            throw new Error ("internal error - unexpected exception: " + e);
          } catch (EntryInformationException e) {
            throw new Error ("internal error - unexpected exception: " + e);
          } catch (OutOfRangeException e) {
            Object[] options = { "CANCEL", "IGNORE", "IGNORE ALL"};

            if(select != 2)
            {
              select = JOptionPane.showOptionDialog(null, 
                          "Found overlapping CDS\n"+e,
                          "Out of Range",
                           JOptionPane.YES_NO_CANCEL_OPTION,
                           JOptionPane.WARNING_MESSAGE,
                           null,
                           options,
                           options[0]);
              if(select == 0)
                throw new Error ("internal error - unexpected exception: " + e);
            }
          }
        }
      }
    } finally {
      entry_group.getActionController ().endAction ();
    }
  }

  
  /**
   *  Create intergenic regions between CDS.
   *  @param frame The Frame to use for MessageDialog components.
   *  @param entry_group The EntryGroup to create the features in.
   **/
  static void createIntergenicFeatures (final JFrame frame,
                                    final EntryGroup entry_group) 
  {
    try 
    {
      entry_group.getActionController ().startAction ();

      final FeaturePredicate predicate;
      
      if(GeneUtils.isDatabaseEntry(entry_group))
        predicate = new FeatureKeyPredicate(new Key("gene"));
      else
      {
        final FeaturePredicateVector temp_predicates =
          new FeaturePredicateVector();
        
        temp_predicates.add(new FeatureKeyPredicate(Key.CDS));
        temp_predicates.add(new FeatureKeyPredicate(new Key("tRNA")));
        
        predicate =
          new FeaturePredicateConjunction(temp_predicates,
                                          FeaturePredicateConjunction.OR);
      }
      
      FeatureVector cdsFeatures = new FeatureVector ();
      final FeatureEnumeration feature_enum = entry_group.features ();
      while (feature_enum.hasMoreFeatures ()) 
      {
        final Feature current_feature = feature_enum.nextFeature ();
        if(predicate.testPredicate (current_feature)) 
          cdsFeatures.add (current_feature);
      }

      cdsFeatures = cdsFeatures.sort(feature_comparator);
      
      int prevEnd = 0;
      Entry newEntry = null;
      boolean prevForward = true;
      
      // search for intergenic regions
      for(int i=0; i < cdsFeatures.size (); i++)
      {
        final Feature this_feature = cdsFeatures.elementAt(i);
        final Location cds_location = this_feature.getLocation();
        final Range r = cds_location.getTotalRange();
        
        int currentStart = r.getStart();
        
        if(i==0 && r.getStart()==1)
        {
          prevEnd = r.getEnd();
          prevForward = this_feature.isForwardFeature();
          // check for overlapping CDS
          if(i<cdsFeatures.size()-1)
          {
        	  int next = i+1; 	
            while(getTotalRange((Feature)cdsFeatures.elementAt(next)).getStart() <= prevEnd+1)
            {
              Feature f = (Feature)cdsFeatures.elementAt(next);
              
              if(prevEnd < getTotalRange(f).getEnd())
              {
                prevEnd = getTotalRange(f).getEnd();
                prevForward = f.isForwardFeature();
              }
              
              i = next;
              next++;
            }
          }
          continue;
        }
        
        try 
        {
          Range new_range = new Range(prevEnd + 1,
                                      currentStart - 1);
          Location location = new Location(new_range);
          final Key key;
          
          if(GeneUtils.isDatabaseEntry(entry_group))
            key = new Key ("region");            
          else
            key = new Key ("misc_feature");
          
          
          // intergenic regions (IGR) - flanking CDS 4 possible:
          // IGR-F (forward):  cds> IGR cds>
          // IGR-R (reverse): <cds IGR <cds
          // IGR-B (both): <cds IGR cds>
          // IGR-X: cds> IGR <cds
          String note;
          
          if(this_feature.isForwardFeature())
          {
            if(prevForward)
              note = "IGR-F";
            else
              note = "IGR-B";
          }
          else
          {
            if(prevForward)
              note = "IGR-X";
            else
              note = "IGR-R";
          }
          final QualifierVector qualifiers = new QualifierVector ();
          final Qualifier qualifier = new Qualifier("note", note);
          qualifiers.add(qualifier);
          
          if(newEntry == null)
            newEntry = entry_group.createEntry("intergenic");
          
          newEntry.createFeature(key, location, qualifiers);
          prevEnd = r.getEnd();
          prevForward = this_feature.isForwardFeature();
          
          // check for overlapping CDS
          if(i<cdsFeatures.size()-1)
          {
        	  int next = i+1;
            while( next < cdsFeatures.size() &&
                  getTotalRange((Feature)cdsFeatures.elementAt(next)).getStart() <= prevEnd+1)
            {
              Feature f = (Feature)cdsFeatures.elementAt(next);
              
              if(prevEnd < getTotalRange(f).getEnd())
              {
                prevEnd = getTotalRange(f).getEnd();
                prevForward = f.isForwardFeature();
              }
              i = next;
              next++;
            }
          }
          if(i==cdsFeatures.size()-1)
          {
            if(entry_group.getSequenceLength() > r.getEnd())
            {
              new_range = new Range(prevEnd + 1,
                  entry_group.getSequenceLength());
              location = new Location(new_range);
              newEntry.createFeature(key,
                                        location, qualifiers);
            }
          }
        } 
        catch (OutOfRangeException e) {}
        catch(ReadOnlyException e)  {}
        catch(EntryInformationException e) {}
        
      }
    } 
    finally 
    {
      entry_group.getActionController ().endAction ();
    }
  }

  private static Range getTotalRange(Feature f)
  {
    return ((Range) f.getLocation().getTotalRange() );
  }
  
  /**
   *  This is used by getSortedFeaturesInRange().
   **/
  final private static Comparator feature_comparator = new Comparator()
  {
    /**
     *  Compare two Objects with respect to ordering.
     *  @return a negative number if feature1_object is less than
     *    feature2_object ; a positive number if feature1_object is greater
     *    than feature2_object; else 0
     **/
    public int compare(final Object feature1_object,
                       final Object feature2_object)
    {
      final Feature feature1 =(Feature) feature1_object;
      final Feature feature2 =(Feature) feature2_object;

      final int feature1_start = feature1.getLocation().getTotalRange().getStart();
      final int feature2_start = feature2.getLocation().getTotalRange().getStart();

      if(feature1_start < feature2_start)
        return -1;
      else if(feature1_start > feature2_start)
        return 1;

      final int feature1_end = feature1.getLocation().getTotalRange().getEnd();
      final int feature2_end = feature2.getLocation().getTotalRange().getEnd();
      
      if(feature1_end < feature2_end)
        return -1;
      else
        return 1;
    }
  };

  
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
          final Range this_range = (Range)cds_ranges.elementAt (range_index);

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

        try 
        {
          
          if(selection_feature.getEmblFeature() instanceof GFFStreamFeature)
          {
            String uniquename = GeneUtils.promptForUniquename(entry_group, 
                                 selection_feature.isForwardFeature());
          
            Qualifier qualifier = new Qualifier("ID", uniquename);
            qualifiers.setQualifier(qualifier);
          }
          
          Feature feature = 
            selection_feature.getEntry ().createFeature (gene_key,
                                                       gene_location,
                                                       qualifiers);
          
          
          if(feature.getEmblFeature() instanceof GFFStreamFeature)
          {
            ChadoCanonicalGene chado_gene = new ChadoCanonicalGene();
            chado_gene.setGene(feature.getEmblFeature());
            ((uk.ac.sanger.artemis.io.GFFStreamFeature)
                (feature.getEmblFeature())).setChadoGene(chado_gene);
          }
          
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
  private void markORFSWithSize (final boolean empty_only, final JFrame frame) {
    final int default_minimum_orf_size =
      Options.getOptions ().getMinimumORFSize ();

    Box inputBox = Box.createVerticalBox();
    inputBox.add(new JLabel("Minimum open reading frame size:"));
    JTextField minSize = new JTextField(String.valueOf (default_minimum_orf_size));
    inputBox.add(minSize);
    JCheckBox useFastaBoundary = new JCheckBox("break at contig boundaries (multiple fasta only)", false);
    inputBox.add(useFastaBoundary);
    
    int select = JOptionPane.showConfirmDialog(getParentFrame(), 
                                inputBox, "ORF options", 
                                JOptionPane.OK_CANCEL_OPTION,
                                JOptionPane.QUESTION_MESSAGE);
    
    if(select == JOptionPane.CANCEL_OPTION)
      return;
    
    final String requester_text = minSize.getText().trim();

    if (requester_text.length () == 0) 
      return;
   
    try
    {
      final int minimum_orf_size =
            Integer.valueOf (requester_text).intValue ();

      markOpenReadingFrames(minimum_orf_size, empty_only, 
                            useFastaBoundary.isSelected(), frame);
    } 
    catch (NumberFormatException e) 
    {
      new MessageDialog(getParentFrame(),
                        "this is not a number: " + requester_text);
    }
  }

  /**
   *  Create a new Feature for each open reading frame.
   *  @param minimum_orf_size All the returned ORFs will be at least this many
   *    amino acids long.
   *  @param empty_only If true only those ORFS that don't already contain a
   *    segment will be marked.
   **/
  private void markOpenReadingFrames (final int minimum_orf_size,
                                      final boolean empty_only,
                                      final boolean isMultiFasta,
                                      final JFrame frame) {
    frame.setCursor(cbusy);
    try {
      final Entry new_entry =
        entry_group.createEntry ("ORFS_" + minimum_orf_size + '+');

      if(isMultiFasta)   // ensure ORF's do not cross fasta boundaries
      {
        final FeatureKeyPredicate predicate =
          new FeatureKeyPredicate(new Key("fasta_record"));
        final FeatureVector fasta_features = new FeatureVector ();

        final FeatureEnumeration feature_enum = entry_group.features ();

        while (feature_enum.hasMoreFeatures ()) 
        {
          final Feature current_feature = feature_enum.nextFeature ();
          if (predicate.testPredicate (current_feature))
            fasta_features.add (current_feature);
        }
        
        for(int i=0;i<fasta_features.size();i++)
        {
          final int start = fasta_features.elementAt(i).getFirstBase();
          final int last  = fasta_features.elementAt(i).getLastBase();

          final MarkerRange forward_range = 
            entry_group.getBases().getForwardStrand().
            makeMarkerRangeFromPositions(start, last);

          markOpenReadingFrames(new_entry, forward_range, minimum_orf_size,
                                empty_only, last, start);
          
          int length = entry_group.getBases().getLength();
          final MarkerRange backward_range = 
            entry_group.getBases().getReverseStrand().
            makeMarkerRangeFromPositions((length-last+1), (length-start+1));

          markOpenReadingFrames(new_entry, backward_range, minimum_orf_size,
                                empty_only, length-start+1, length-last+1);
        }
      }
      else
      {
        final int sequence_length = entry_group.getSequenceLength();

        final Strand forward_strand = entry_group.getBases().getForwardStrand();

        final MarkerRange forward_range = forward_strand
            .makeMarkerRangeFromPositions(1, sequence_length);

        markOpenReadingFrames(new_entry, forward_range, minimum_orf_size,
            empty_only, sequence_length, 1);

        final Strand backward_strand = entry_group.getBases()
            .getReverseStrand();

        final MarkerRange backward_range = backward_strand
            .makeMarkerRangeFromPositions(1, sequence_length);

        markOpenReadingFrames(new_entry, backward_range, minimum_orf_size,
            empty_only, sequence_length, 1);
      }
    } catch (OutOfRangeException e) {
      frame.setCursor(cdone);
      throw new Error ("internal error - unexpected OutOfRangeException");
    }
    frame.setCursor(cdone);
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
                                 false, entry_group.getSequenceLength(), 1);


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
                                      final boolean empty_only,
                                      final int sequence_end,
                                      final int sequence_start) {
    final MarkerRange [] forward_orf_ranges =
      Strand.getOpenReadingFrameRanges (search_range, minimum_orf_size, sequence_end,
          sequence_start);

    String uniquename = GeneUtils.promptForUniquename(entry_group, search_range.isForwardMarker());
    
    final Key key;
    if(GeneUtils.isDatabaseEntry(entry_group))
      key = new Key("region");
    else
      key = Key.CDS;
    
    for(int i = 0 ; i < forward_orf_ranges.length ; ++i) 
    {
      final MarkerRange this_range = forward_orf_ranges[i];
      final Feature new_feature;

      try 
      {
        QualifierVector qualifiers = null;
        if(entry.getEMBLEntry() instanceof 
            uk.ac.sanger.artemis.io.DatabaseDocumentEntry)
        {
          Qualifier qualifier = new Qualifier("ID", uniquename+Integer.toString(i+1));
          qualifiers = new QualifierVector();
          qualifiers.setQualifier(qualifier);
        }
        new_feature = makeFeatureFromMarkerRange (entry, this_range, key, qualifiers);
      } 
      catch (EntryInformationException e) 
      {
        e.printStackTrace();
        new MessageDialog (getParentFrame (), "cannot continue: " +
                           "the default entry does not support CDS features");
        return;
      } 
      catch (ReadOnlyException e) 
      {
        new MessageDialog (getParentFrame (), "cannot continue: " +
                           "the default entry is read only");
        return;
      }

      if(empty_only && overlapsAnActiveSegment (new_feature)) 
      {
        try 
        {
          new_feature.removeFromEntry ();
        } 
        catch (ReadOnlyException exception) 
        {
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
                                              final Key key,
                                              QualifierVector qualifiers)
      throws EntryInformationException, ReadOnlyException 
  {
    try 
    {
      final Location new_location = range.createLocation ();

      if(qualifiers == null)
      {
        qualifiers = new QualifierVector ();
        qualifiers.setQualifier (new Qualifier ("note", "none"));
      }

      final Feature new_feature =
        entry.createFeature (key, new_location, qualifiers);

      final CodonUsageAlgorithm codon_usage_algorithm =
        base_plot_group.getCodonUsageAlgorithm ();

      if(codon_usage_algorithm != null) 
      {
        int score =
          (int) (codon_usage_algorithm.getFeatureScore (new_feature) * 50);

        if(score < 0) 
          score = 0;

        if(score > 100) 
          score = 100;

        final String score_string = String.valueOf (score);
        new_feature.addQualifierValues (new Qualifier ("score",
                                                       score_string));

        final int var_colour = 255 - score * 5 / 2;
        final String colour_string = var_colour + " " + var_colour + " 255";
        new_feature.addQualifierValues (new Qualifier ("colour",
                                                       colour_string));
      }

      return new_feature;
    } 
    catch (OutOfRangeException e) 
    {
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
  private void makeFeaturesFromPattern (final BasePattern pattern) 
  {
    final MarkerRangeVector matches =
      pattern.findMatches (entry_group.getBases (),
                           null,        // search from start
                           entry_group.getSequenceLength ());

    if (matches.size () == 0) 
    {
      new MessageDialog (getParentFrame (),
                         "no matches found for: " + pattern);
      return;
    }

    final int TOO_MANY_MATCHES = 100;

    if (matches.size () > TOO_MANY_MATCHES) 
    {
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

    String uniquename = null;
    
    if(entry_group.getDefaultEntry().getEMBLEntry() instanceof 
        uk.ac.sanger.artemis.io.DatabaseDocumentEntry)
      uniquename = GeneUtils.promptForUniquename(entry_group, true); 
    
    for (int i = 0 ; i < matches.size () ; ++i) 
    {
      try 
      {
        QualifierVector qualifiers = null;
        
        if(uniquename != null)
        {
          Qualifier qualifier = new Qualifier("ID", uniquename+Integer.toString(i+1));
          qualifiers = new QualifierVector();
          qualifiers.setQualifier(qualifier);
        }
        
        final Feature new_feature =
          makeFeatureFromMarkerRange (new_entry, matches.elementAt (i), key, qualifiers);
        new_feature.setQualifier (new Qualifier ("note", pattern.toString ()));
      } 
      catch (EntryInformationException e) 
      {
        new MessageDialog (getParentFrame (), "cannot continue: " +
                           e.getMessage ());
        return;
      } 
      catch (ReadOnlyException e) 
      {
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
    final Key unsure;
    if(entry_group.getSequenceEntry() != null &&
       entry_group.getSequenceEntry().getEMBLEntry() instanceof DatabaseDocumentEntry)
      unsure = new Key ("region");
    else
      unsure = new Key ("unsure");
    
    Pattern p = Pattern.compile("^[nN]+$"); // pattern match for all n's
    
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

          if(p.matcher(unsure_bases).matches())
          {
            final Feature feature =
              new_entry.createFeature (new Key("gap"), location, qualifiers);
          }
          else
          {
            qualifiers.setQualifier (new Qualifier ("note", unsure_bases));
            final Feature feature =
              new_entry.createFeature (unsure, location, qualifiers);
          }
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

}
