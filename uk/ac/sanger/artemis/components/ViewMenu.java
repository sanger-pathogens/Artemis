/* ViewMenu.java
 *
 * created: Tue Dec 29 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/ViewMenu.java,v 1.15 2009-04-22 08:50:13 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.sequence.*;

import uk.ac.sanger.artemis.util.*;
import uk.ac.sanger.artemis.components.filetree.FileList;
import uk.ac.sanger.artemis.components.filetree.RemoteFileNode;
import uk.ac.sanger.artemis.components.genebuilder.GeneUtils;

import uk.ac.sanger.artemis.io.ChadoCanonicalGene;
import uk.ac.sanger.artemis.io.DocumentEntry;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.ValidateFeature;
import uk.ac.sanger.artemis.io.InvalidRelationException;
import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.j2ssh.FTProgress;
import uk.ac.sanger.artemis.j2ssh.FileTransferProgressMonitor;

import java.io.*;
import java.awt.*;
import java.awt.event.*;

import javax.swing.*;

import com.sshtools.j2ssh.sftp.FileAttributes;

/**
 *  A popup menu with viewing commands.
 *
 *  @author Kim Rutherford
 *  @version $Id: ViewMenu.java,v 1.15 2009-04-22 08:50:13 tjc Exp $
 **/

public class ViewMenu extends SelectionMenu 
{
  /** */
  private static final long serialVersionUID = 1L;

  /**
   *  The EntryGroup that was passed to the constructor.
   **/
  private EntryGroup entry_group = null;

  /**
   *  The Selection that was passed to the constructor.
   **/
  private Selection selection = null;

  /**
   *  The GotoEventSource that was passed to the constructor.
   **/
  //private GotoEventSource goto_event_source = null;

  private BasePlotGroup base_plot_group;
  
  /**
   *  Create a new ViewMenu object.
   *  @param frame The JFrame that owns this JMenu.
   *  @param selection The Selection that the commands in the menu will
   *    operate on.
   *  @param entry_group The EntryGroup object where new features/entries will
   *    be added.
   *  @param goto_event_source The object that we will call makeBaseVisible()
   *    on.
   *  @param base_plot_group The BasePlotGroup associated with this JMenu -
   *    needed to call getCodonUsageAlgorithm()
   *  @param menu_name The name of the new menu.
   **/
  public ViewMenu(final JFrame frame,
                  final Selection selection,
                  final GotoEventSource goto_event_source,
                  final EntryGroup entry_group,
                  final BasePlotGroup base_plot_group,
                  final String menu_name) 
  {
    super(frame, menu_name, selection);

    this.entry_group = entry_group;
    this.selection = selection;

    this.base_plot_group = base_plot_group;

    final JMenuItem plot_features_item = new JMenuItem("Feature Plots");
    plot_features_item.setAccelerator(PLOT_FEATURES_KEY);
    plot_features_item.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event) 
      {
        plotSelectedFeatures(getParentFrame(), getSelection());
      }
    });

    final JMenuItem view_feature_item = new JMenuItem("Selected Features");
    view_feature_item.setAccelerator(VIEW_FEATURES_KEY);
    view_feature_item.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        viewSelectedFeatures(getParentFrame(), getSelection());
      }
    });

    final JMenuItem view_selection_item = new JMenuItem("Selection");
    view_selection_item.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent event) {
        new SelectionViewer(getSelection(), entry_group);
      }
    });

    final JMenuItem feature_info_item = new JMenuItem("Feature Statistics");
    feature_info_item.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        viewSelectedFeatureInfo();
      }
    });

    final SelectionSubMenu view_bases = new SelectionSubMenu(this, "Bases");
    final JMenuItem view_bases_item = new JMenuItem("Bases Of Selection");
    view_bases_item.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event) 
      {
        viewSelectedBases(true, false);
      }
    });
    view_bases.add(view_bases_item);

    final JMenuItem view_bases_as_fasta_item =
      new JMenuItem("Bases Of Selection As FASTA");
    view_bases_as_fasta_item.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        viewSelectedBases(false, false);
      }
    });
    view_bases.add(view_bases_as_fasta_item);
    view_bases.addSeparator();

    final JMenuItem view_exon_bases =
        new JMenuItem("Bases Of Selected Exons As FASTA");
    view_exon_bases.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        viewSelectedBases(false, true);
      }
    });
    view_bases.add(view_exon_bases);

    final SelectionSubMenu view_aa = new SelectionSubMenu(this, "Amino Acids");
    final JMenuItem view_aa_item = new JMenuItem("Amino Acids Of Selection");
    view_aa_item.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        viewSelectedAminoAcids(true);
      }
    });
    view_aa.add(view_aa_item);

    final JMenuItem view_aa_as_fasta_item =
      new JMenuItem("Amino Acids Of Selection As FASTA");
    view_aa_as_fasta_item.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        viewSelectedAminoAcids(false);
      }
    });
    view_aa.add(view_aa_as_fasta_item);

    final JMenuItem overview_item = new JMenuItem("Overview");
    overview_item.setAccelerator(OVERVIEW_KEY);
    overview_item.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        new EntryGroupInfoDisplay(getParentFrame(), entry_group);
      }
    });

    final JMenuItem forward_overview_item = new JMenuItem("Forward Strand Overview");
    forward_overview_item.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)
      {
        new EntryGroupInfoDisplay(getParentFrame(), entry_group,
                                  Bases.FORWARD);
      }
    });

    final JMenuItem reverse_overview_item = new JMenuItem("Reverse Strand Overview");
    reverse_overview_item.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event)
      {
        new EntryGroupInfoDisplay(getParentFrame(), entry_group,
                                  Bases.REVERSE);
      }
    });

    final JMenuItem view_cds_item = new JMenuItem("CDS Genes And Products");
    if(GeneUtils.isDatabaseEntry(entry_group))
      view_cds_item.setEnabled(false);
    view_cds_item.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        final FeaturePredicate feature_predicate =
          new FeatureKeyPredicate(Key.CDS);

        final String filter_name =
          "CDS features (filtered from: " +
          getParentFrame().getTitle() + ")";

        final FilteredEntryGroup filtered_entry_group =
          new FilteredEntryGroup(entry_group, feature_predicate, filter_name);

        final FeatureListFrame feature_list_frame =
          new FeatureListFrame(filter_name,
                                selection, goto_event_source,
                                filtered_entry_group,
                                base_plot_group);

        feature_list_frame.getFeatureList().setShowGenes(true);
        feature_list_frame.getFeatureList().setShowProducts(true);

        feature_list_frame.setVisible(true);
      }
    });

    JMenu search_results_menu = null;

    search_results_menu = new SelectionSubMenu(this, "Search Results");

    final boolean sanger_options =
      Options.getOptions().getPropertyTruthValue("sanger_options");

    final ExternalProgramVector external_programs =
      Options.getOptions().getExternalPrograms();

    final StringVector external_program_names = new StringVector();

    for(int i = 0 ; i < external_programs.size() ; ++i) 
    {
      final ExternalProgram external_program =
        external_programs.elementAt(i);

      final String new_name = external_program.getName();

      if(!external_program_names.contains(new_name)) 
        external_program_names.add(new_name);
    }

    for(int i = 0 ; i < external_program_names.size() ; ++i) 
    {
      final String external_program_name =
        (String)external_program_names.elementAt(i);

      final JMenuItem new_menu =
        makeSearchResultsMenu(external_program_name, false, sanger_options);
      search_results_menu.add(new_menu);
    }

    if(sanger_options) 
    {
      search_results_menu.addSeparator();

      for(int i = 0 ; i < external_program_names.size() ; ++i)
      {
        final String external_program_name =
          (String)external_program_names.elementAt(i);

        final JMenuItem new_menu =
          makeSearchResultsMenu(external_program_name, true, sanger_options);
        search_results_menu.add(new_menu);
      }
    }

    final int MAX_FILTER_FEATURE_COUNT = 10000;

    final JMenu feature_filters_menu = new SelectionSubMenu(this, "Feature Filters");

    final JMenuItem bad_start_codons_item =
      new JMenuItem("Suspicious Start Codons ...");
    bad_start_codons_item.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        if(checkEntryGroupSize(MAX_FILTER_FEATURE_COUNT))
          showBadStartCodons(getParentFrame(), selection, 
                             entry_group, goto_event_source,
                             base_plot_group);
      }
    });

    final JMenuItem bad_stop_codons_item =
      new JMenuItem("Suspicious Stop Codons ...");
    bad_stop_codons_item.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        if(checkEntryGroupSize(MAX_FILTER_FEATURE_COUNT)) 
          showBadStopCodons(getParentFrame(), selection,
                            entry_group, goto_event_source,
                            base_plot_group);
      }
    });

    final JMenuItem stop_codons_in_translation =
      new JMenuItem("Stop Codons In Translation ...");
    stop_codons_in_translation.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        if(checkEntryGroupSize(MAX_FILTER_FEATURE_COUNT)) 
          showStopsInTranslation(getParentFrame(), selection,
                                 entry_group, goto_event_source,
                                 base_plot_group);
      }
    });
    
    
    final JMenuItem intronsSpliceSite =
        new JMenuItem("Introns without GT/GC start and AG end ...");
    intronsSpliceSite.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        if(checkEntryGroupSize(MAX_FILTER_FEATURE_COUNT)) 
          showIntrons(selection, entry_group, goto_event_source,
              base_plot_group);
      }
    });
    
    
    final JMenuItem geneModelCheck =
        new JMenuItem("Gene model boundary check ...");
    geneModelCheck.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        if(checkEntryGroupSize(MAX_FILTER_FEATURE_COUNT)) 
          geneBoundaryCheck(selection, entry_group, goto_event_source,
              base_plot_group);
      }
    });
    
    
    final JMenuItem validate =
        new JMenuItem("Validation checks ...");
    validate.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        if(checkEntryGroupSize(MAX_FILTER_FEATURE_COUNT))
        {
          ValidateFeature gffTest = new ValidateFeature(getEntryGroup());
          gffTest.featureListErrors(entry_group, selection, goto_event_source, base_plot_group);
        }
      }
    });
    

    final JMenuItem bad_feature_keys_item =
      new JMenuItem("Non EMBL Keys ...");
    bad_feature_keys_item.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        if(checkEntryGroupSize(MAX_FILTER_FEATURE_COUNT))
          showNonEMBLKeys(getParentFrame(), selection,
                          entry_group, goto_event_source,
                          base_plot_group);
      }
    });

    final JMenuItem duplicated_keys_item =
      new JMenuItem("Duplicated Features ...");
    duplicated_keys_item.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        if(checkEntryGroupSize(MAX_FILTER_FEATURE_COUNT)) 
          showDuplicatedFeatures(getParentFrame(), selection,
                                 entry_group, goto_event_source,
                                 base_plot_group);
      }
    });

    final JMenuItem overlapping_cds_features_item;
    
    if(GeneUtils.isDatabaseEntry(entry_group))
      overlapping_cds_features_item =
        new JMenuItem("Overlapping "+DatabaseDocument.EXONMODEL+" Features ...");
    else
      overlapping_cds_features_item =   
        new JMenuItem("Overlapping CDS Features ...");
    overlapping_cds_features_item.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        if(checkEntryGroupSize(MAX_FILTER_FEATURE_COUNT)) 
          showOverlappingCDSs(getParentFrame(), selection,
                              entry_group, goto_event_source,
                              base_plot_group);
      }
    });

    final JMenuItem same_stop_cds_features_item;
    if(GeneUtils.isDatabaseEntry(entry_group))
      same_stop_cds_features_item = new JMenuItem(
          DatabaseDocument.EXONMODEL+"s Sharing Stop Codons ...");
    else
      same_stop_cds_features_item = new JMenuItem("CDSs Sharing Stop Codons ...");
    same_stop_cds_features_item.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        if(checkEntryGroupSize(MAX_FILTER_FEATURE_COUNT)) 
          showFeaturesWithSameStopCodons(getParentFrame(), 
                                         selection, entry_group,
                                         goto_event_source,
                                         base_plot_group);
      }
    });

    final JMenuItem missing_qualifier_features_item =
      new JMenuItem("Features Missing Required Qualifiers ...");
    missing_qualifier_features_item.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event) 
      {
        if(checkEntryGroupSize(MAX_FILTER_FEATURE_COUNT)) 
          showMissingQualifierFeatures(getParentFrame(), selection,
                                       entry_group, goto_event_source,
                                       base_plot_group);
      }
    });
    
    final JMenuItem all_filters_item =
        new JMenuItem("Apply All Filters Above ...");
    all_filters_item.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent event) 
        {
          if(!checkEntryGroupSize(MAX_FILTER_FEATURE_COUNT))
            return;
          
          showBadStartCodons(getParentFrame(), selection, 
              entry_group, goto_event_source,
              base_plot_group);
          
          showBadStopCodons(getParentFrame(), selection,
              entry_group, goto_event_source,
              base_plot_group);
          
          showStopsInTranslation(getParentFrame(), selection,
              entry_group, goto_event_source,
              base_plot_group);
          
          showNonEMBLKeys(getParentFrame(), selection,
              entry_group, goto_event_source,
              base_plot_group);
          
          showDuplicatedFeatures(getParentFrame(), selection,
              entry_group, goto_event_source,
              base_plot_group);
          
          showOverlappingCDSs(getParentFrame(), selection,
              entry_group, goto_event_source,
              base_plot_group);
          
          showFeaturesWithSameStopCodons(getParentFrame(), 
              selection, entry_group,
              goto_event_source,
              base_plot_group);
          
          showMissingQualifierFeatures(getParentFrame(), selection,
              entry_group, goto_event_source,
              base_plot_group);
          
          showFilterByMultipleID(getParentFrame(), selection,
              entry_group, goto_event_source,
              base_plot_group);
        }
      });

    final JMenuItem filter_by_key_item =
      new JMenuItem("Filter By Key ...");
    filter_by_key_item.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        if(checkEntryGroupSize(MAX_FILTER_FEATURE_COUNT)) 
          showFilterByKey(getParentFrame(), selection,
                          entry_group, goto_event_source,
                          base_plot_group);
      }
    });


    final JMenuItem filter_by_multiple_sys_id =
      new JMenuItem("Duplicate Systematic Name Qualifier ...");
    filter_by_multiple_sys_id.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        showFilterByMultipleID(getParentFrame(), selection,
                               entry_group, goto_event_source,
                               base_plot_group);
      }
    });
    
    
    final JMenuItem filter_by_selection_item =
      new JMenuItem("Selected Features ...");
    filter_by_selection_item.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        showFilterBySelection(getParentFrame(),
                              selection, entry_group, goto_event_source,
                              base_plot_group);
      }
    });

    feature_filters_menu.add(bad_start_codons_item);
    feature_filters_menu.add(bad_stop_codons_item);
    feature_filters_menu.add(stop_codons_in_translation);
    feature_filters_menu.add(intronsSpliceSite);
    if(GeneUtils.isGFFEntry( getEntryGroup() ))
      feature_filters_menu.add(geneModelCheck);

    feature_filters_menu.add(bad_feature_keys_item);
    feature_filters_menu.add(duplicated_keys_item);
    feature_filters_menu.add(overlapping_cds_features_item);
    feature_filters_menu.add(same_stop_cds_features_item);
    feature_filters_menu.add(missing_qualifier_features_item);
    feature_filters_menu.add(filter_by_multiple_sys_id);
    feature_filters_menu.add(validate);
    
    feature_filters_menu.addSeparator();
    feature_filters_menu.add(all_filters_item);
    feature_filters_menu.addSeparator();
    feature_filters_menu.add(filter_by_key_item);
    feature_filters_menu.add(filter_by_selection_item);

    add(view_feature_item);
    add(view_selection_item);
    addSeparator();
    if(search_results_menu != null)
      add(search_results_menu);
    
    add(view_cds_item);
    add(feature_filters_menu);
    addSeparator();
    add(overview_item);
    add(forward_overview_item);
    add(reverse_overview_item);
    addSeparator();
    add(view_bases);
    add(view_aa);
    addSeparator();
    add(feature_info_item);
    add(plot_features_item);
    
    //If "view Custom Annotation" is selected in menu or
    //view_Custom_Annotation is set up in option file,add submenus to the View menu. 
    //@Author Luj 19/08/2013 (start)
        if (
                //Options.getOptions().getPropertyTruthValue("view_Custom_Annotation") ||
               (System.getProperty("viewCustomAnnotation")!=null &&
                System.getProperty("viewCustomAnnotation").equals("true"))) {
            JMenu view_Customized_Annotation_menu = null;
            view_Customized_Annotation_menu = new JMenu("View Customized Annotation Results");
            final JMenuItem view_GAMOLA_BLAST = new JMenuItem("View GAMOLA BLAST Annotation for Selection");
            view_GAMOLA_BLAST.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent event) {
                    new SelectionCustomViewer(getSelection(), entry_group, "Blast_database");
                }
            });
            view_Customized_Annotation_menu.add(view_GAMOLA_BLAST);

            final JMenuItem view_GAMOLA_COG = new JMenuItem("View GAMOLA COG Annotation for Selection");
            view_GAMOLA_COG.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent event) {
                    new SelectionCustomViewer(getSelection(), entry_group, "COG_database");
                }
            });
            view_Customized_Annotation_menu.add(view_GAMOLA_COG);

            final JMenuItem view_GAMOLA_PFAM = new JMenuItem("View GAMOLA PFAM Annotation for Selection");
            view_GAMOLA_PFAM.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent event) {
                    new SelectionCustomViewer(getSelection(), entry_group, "PFam_database");
                }
            });
            view_Customized_Annotation_menu.add(view_GAMOLA_PFAM);

            final JMenuItem view_GAMOLA_TIGR = new JMenuItem("View GAMOLA TIGR Annotation for Selection");
            view_GAMOLA_TIGR.addActionListener(new ActionListener() {
                public void actionPerformed(ActionEvent event) {
                    new SelectionCustomViewer(getSelection(), entry_group, "TIGRfam_database");
                }
            });
            view_Customized_Annotation_menu.add(view_GAMOLA_TIGR);
            addSeparator();
            add(view_Customized_Annotation_menu);
            view_Customized_Annotation_menu.add(view_GAMOLA_BLAST);
            view_Customized_Annotation_menu.add(view_GAMOLA_COG);
            view_Customized_Annotation_menu.add(view_GAMOLA_PFAM);
            view_Customized_Annotation_menu.add(view_GAMOLA_TIGR);

        }
  }
  //@Author Luj 19/08/2013 (end)
  

  
  /**
   *  Create a new ViewMenu object.
   *  @param entry_edit The EntryEdit that owns this JMenu.
   *  @param selection The Selection that the commands in the menu will
   *    operate on.
   *  @param entry_group The EntryGroup object where new features/entries will
   *    be added.
   *  @param goto_event_source The object that we will call makeBaseVisible()
   *    on.
   *  @param base_plot_group The BasePlotGroup associated with this JMenu -
   *    needed to call getCodonUsageAlgorithm()
   **/
  public ViewMenu(final JFrame frame,
                  final Selection selection,
                  final GotoEventSource goto_event_source,
                  final EntryGroup entry_group,
                  final BasePlotGroup base_plot_group) 
  {
    this(frame, selection, goto_event_source, entry_group,
          base_plot_group, "View");
  }

  /**
   *  The shortcut for Show Feature Plots.
   **/
  final static KeyStroke PLOT_FEATURES_KEY =
    KeyStroke.getKeyStroke (KeyEvent.VK_W, 
                            Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()); //InputEvent.CTRL_MASK);

  final static public int PLOT_FEATURES_KEY_CODE = KeyEvent.VK_W;

  /**
   *  The shortcut for View Selected Features.
   **/
  final static KeyStroke VIEW_FEATURES_KEY =
    KeyStroke.getKeyStroke (KeyEvent.VK_V, 
                            Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()); //InputEvent.CTRL_MASK);

  final static public int VIEW_FEATURES_KEY_CODE = KeyEvent.VK_V;

  /**
   *  The shortcut for Show Overview.
   **/
  final static KeyStroke OVERVIEW_KEY =
    KeyStroke.getKeyStroke (KeyEvent.VK_O, 
                            Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()); //InputEvent.CTRL_MASK);

  final static public int OVERVIEW_KEY_CODE = KeyEvent.VK_O;

  /**
   *  The shortcut for View FASTA in browser.
   **/
  final static KeyStroke FASTA_IN_BROWSER_KEY =
    KeyStroke.getKeyStroke (KeyEvent.VK_F, 
                            Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()); //InputEvent.CTRL_MASK);

  final static public int FASTA_IN_BROWSER_KEY_CODE = KeyEvent.VK_F;

  /**
   *  The shortcut for View FASTA.
   **/
  final static KeyStroke VIEW_FASTA_KEY =
    KeyStroke.getKeyStroke (KeyEvent.VK_R, 
                            Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()); //InputEvent.CTRL_MASK);

  final static public int VIEW_FASTA_KEY_CODE = KeyEvent.VK_R;

  /**
   *  The shortcut for View BLASTP in browser.
   **/
  final static KeyStroke BLASTP_IN_BROWSER_KEY =
    KeyStroke.getKeyStroke (KeyEvent.VK_B, 
                            Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()); //InputEvent.CTRL_MASK);

  final static public int BLASTP_IN_BROWSER_KEY_CODE = KeyEvent.VK_B;

  /**
   *  The shortcut for View BLASTP.
   **/
  final static KeyStroke VIEW_BLASTP_KEY =
    KeyStroke.getKeyStroke (KeyEvent.VK_BACK_QUOTE , 
                            Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()); //InputEvent.CTRL_MASK);

  final static public int VIEW_BLASTP_KEY_CODE = KeyEvent.VK_BACK_QUOTE;

  /**
   *  The shortcut for View HTH.
   **/
  final static KeyStroke VIEW_HTH_KEY =
    KeyStroke.getKeyStroke (KeyEvent.VK_H, 
                            Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()); //InputEvent.CTRL_MASK);

  final static public int VIEW_HTH_KEY_CODE = KeyEvent.VK_H;

  /**
   *  Make a JMenuItem for viewing the results of running the given program.
   *  @param send_to_browser if true the results should be sent straight to
   *    the web browser rather than using a SearchResultViewer object.
   *  @param sanger_options true if the sanger_options is set to true in the
   *    options file.
   **/
  private JMenuItem makeSearchResultsMenu(final String program_name,
                                          final boolean send_to_browser,
                                          final boolean sanger_options) 
  {
    final String suffix;
    if(send_to_browser)
      suffix = new String(" results (in browser)");
    else
      suffix = new String(" results");

    final JMenuItem new_menu = new JMenuItem(program_name + suffix);

    if ((sanger_options && send_to_browser || !sanger_options)
        && program_name.equals("fasta"))
      new_menu.setAccelerator (FASTA_IN_BROWSER_KEY);
    else 
    {
      if ((sanger_options && send_to_browser || !sanger_options)
          && program_name.equals("blastp")) 
        new_menu.setAccelerator (BLASTP_IN_BROWSER_KEY);
      else 
      {
        if(program_name.equals("fasta"))
          new_menu.setAccelerator (VIEW_FASTA_KEY);
        else if(program_name.equals("blastp")) 
          new_menu.setAccelerator (VIEW_BLASTP_KEY);
        else if(program_name.equals("hth")) 
          new_menu.setAccelerator (VIEW_HTH_KEY);
      }
    }

    new_menu.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event) 
      {
        viewExternalResults(getParentFrame(), getSelection(),
                            program_name, send_to_browser);
      }
    });

    return new_menu;
  }

  /**
   *  Popup a FeatureListFrame containing the non-pseudo CDS features that
   *  have invalid start codons.
   *  @param parent_frame The parent JFrame.
   *  @param selection The Selection to pass to the FeatureList.
   *  @param entry_group The EntryGroup to pass to the FilteredEntryGroup.
   *  @param goto_source The GotoEventSource to pass to the FeatureList.
   *  @param base_plot_group The BasePlotGroup associated with this JMenu -
   *    needed to call getCodonUsageAlgorithm()
   **/
  protected static void showBadStartCodons (final JFrame parent_frame,
                                         final Selection selection,
                                         final EntryGroup entry_group,
                                         final GotoEventSource goto_source,
                                         final BasePlotGroup base_plot_group) 
  {
    final FeaturePredicate cds_predicate;
    
    if(GeneUtils.isDatabaseEntry(entry_group))
      cds_predicate = new FeatureKeyPredicate(new Key(DatabaseDocument.EXONMODEL));
    else
      cds_predicate =
          new FeaturePredicateConjunction(
              new FeatureKeyQualifierPredicate(Key.CDS, "pseudo", false),
              new FeatureKeyQualifierPredicate(Key.CDS, "pseudogene", false),
              FeaturePredicateConjunction.AND);

    final FeaturePredicate feature_predicate = new FeaturePredicate ()
    {
      public boolean testPredicate (final Feature feature) 
      {
        if(!cds_predicate.testPredicate (feature))
          return false;

        if(feature.hasValidStartCodon (true)) 
          return false;
        else 
          return true;
      }
    };

    final String filter_name;
    
    if(GeneUtils.isDatabaseEntry(entry_group))
      filter_name =
        DatabaseDocument.EXONMODEL+" features with suspicious start codons (filtered from: " +
        parent_frame.getTitle () + ")";
    else
      filter_name =
        "CDS features with suspicious start codons (filtered from: " +
        parent_frame.getTitle () + ")";

    final FilteredEntryGroup filtered_entry_group =
      new FilteredEntryGroup (entry_group, feature_predicate, filter_name);

    final FeatureListFrame feature_list_frame =
      new FeatureListFrame (filter_name,
                            selection, goto_source, filtered_entry_group,
                            base_plot_group);

    feature_list_frame.setVisible (true);
  }

  /**
   *  Popup a FeatureListFrame containing the non-pseudo CDS features that
   *  have invalid stop codons.
   *  @param parent_frame The parent JFrame.
   *  @param selection The Selection to pass to the FeatureList.
   *  @param entry_group The EntryGroup to pass to the FilteredEntryGroup.
   *  @param goto_source The GotoEventSource to pass to the FeatureList.
   *  @param base_plot_group The BasePlotGroup associated with this JMenu -
   *    needed to call getCodonUsageAlgorithm()
   **/
  protected static void showBadStopCodons (final JFrame parent_frame,
                                        final Selection selection,
                                        final EntryGroup entry_group,
                                        final GotoEventSource goto_source,
                                        final BasePlotGroup
                                          base_plot_group) 
  {
    final FeaturePredicate cds_predicate;

    if(GeneUtils.isDatabaseEntry(entry_group))
      cds_predicate = new FeatureKeyPredicate(new Key(DatabaseDocument.EXONMODEL));
    else
      cds_predicate =
          new FeaturePredicateConjunction(
              new FeatureKeyQualifierPredicate(Key.CDS, "pseudo", false),
              new FeatureKeyQualifierPredicate(Key.CDS, "pseudogene", false),
              FeaturePredicateConjunction.AND);

    final FeaturePredicate feature_predicate =  new FeaturePredicate () 
    {
      public boolean testPredicate (final Feature feature) 
      {
        if(!cds_predicate.testPredicate (feature)) 
          return false;

        if(feature.hasValidStopCodon (true)) 
          return false;
        else 
          return true;
      }
    };

    final String filter_name;
    
    if(GeneUtils.isDatabaseEntry(entry_group))
      filter_name=
        DatabaseDocument.EXONMODEL+" features with suspicious stop codons (filtered from: " +
        parent_frame.getTitle () + ")";
    else
      filter_name =
        "CDS features with suspicious stop codons (filtered from: " +
        parent_frame.getTitle () + ")";

    final FilteredEntryGroup filtered_entry_group =
      new FilteredEntryGroup (entry_group, feature_predicate, filter_name);

    final FeatureListFrame feature_list_frame =
      new FeatureListFrame (filter_name,
                            selection, goto_source, filtered_entry_group,
                            base_plot_group);

    feature_list_frame.setVisible (true);
  }

  /**
   *  Popup a FeatureListFrame containing the non-pseudo CDS features that
   *  contain a stop codon in the translation.
   *  @param parent_frame The parent Frame.
   *  @param selection The Selection to pass to the FeatureList.
   *  @param entry_group The EntryGroup to pass to the FilteredEntryGroup.
   *  @param goto_source The GotoEventSource to pass to the FeatureList.
   *  @param base_plot_group The BasePlotGroup associated with this Menu -
   *    needed to call getCodonUsageAlgorithm()
   **/
  private static void showStopsInTranslation(final Frame parent_frame,
                                             final Selection selection,
                                             final EntryGroup entry_group,
                                             final GotoEventSource goto_source,
                                             final BasePlotGroup
                                               base_plot_group) 
  {
    final FeaturePredicate cds_predicate;

    if(GeneUtils.isDatabaseEntry(entry_group))
      cds_predicate = new FeatureKeyPredicate(new Key(DatabaseDocument.EXONMODEL));
    else
      cds_predicate =
          new FeaturePredicateConjunction(
              new FeatureKeyQualifierPredicate(Key.CDS, "pseudo", false),
              new FeatureKeyQualifierPredicate(Key.CDS, "pseudogene", false),
              FeaturePredicateConjunction.AND);

    final FeaturePredicate feature_predicate = new FeaturePredicate () 
    {
      public boolean testPredicate (final Feature feature) 
      {
        if(!cds_predicate.testPredicate (feature)) 
          return false;

        final AminoAcidSequence amino_acids = feature.getTranslation ();
        if(amino_acids.containsStopCodon ())
          return true;
        else
          return false;
      }
    };

    final String filter_name;
    if(GeneUtils.isDatabaseEntry(entry_group))
      filter_name =
        DatabaseDocument.EXONMODEL+" features with stop codon(s) in translation (filtered from: " +
        parent_frame.getTitle () + ")";
    else
      filter_name =
        "CDS features with stop codon(s) in translation (filtered from: " +
        parent_frame.getTitle () + ")";

    final FilteredEntryGroup filtered_entry_group =
      new FilteredEntryGroup (entry_group, feature_predicate, filter_name);

    final FeatureListFrame feature_list_frame =
      new FeatureListFrame (filter_name,
                            selection, goto_source, filtered_entry_group,
                            base_plot_group);

    feature_list_frame.setVisible (true);
  }

  /**
   *  Popup a FeatureListFrame containing the features that contains 
   *  introns without GT/GC start and AG end.
   *  @param selection The Selection to pass to the FeatureList.
   *  @param entry_group The EntryGroup to pass to the FilteredEntryGroup.
   *  @param goto_source The GotoEventSource to pass to the FeatureList.
   *  @param base_plot_group The BasePlotGroup associated with this JMenu -
   *    needed to call getCodonUsageAlgorithm()
   **/
  protected static void showIntrons(final Selection selection,
                                    final EntryGroup entry_group,
                                    final GotoEventSource goto_source,
                                    final BasePlotGroup base_plot_group)
  {
    final FeaturePredicate feature_predicate = Selector.getIntronPredicate();
    final String filter_name = "Contains introns without GT/GC start and AG end";
    final FilteredEntryGroup filtered_entry_group =
        new FilteredEntryGroup (entry_group, feature_predicate, filter_name);

      final FeatureListFrame feature_list_frame =
        new FeatureListFrame (filter_name,
                              selection, goto_source, filtered_entry_group,
                              base_plot_group);

    feature_list_frame.setVisible (true);
  }
  

  /**
   *  Popup a FeatureListFrame containing the gene models with boundaries that
   *  need fixing.
   *  @param selection The Selection to pass to the FeatureList.
   *  @param entry_group The EntryGroup to pass to the FilteredEntryGroup.
   *  @param goto_source The GotoEventSource to pass to the FeatureList.
   *  @param base_plot_group The BasePlotGroup associated with this JMenu -
   *    needed to call getCodonUsageAlgorithm()
   **/
  protected static void geneBoundaryCheck(
                                    final Selection selection,
                                    final EntryGroup entry_group,
                                    final GotoEventSource goto_source,
                                    final BasePlotGroup base_plot_group)
  {
    final FeaturePredicate feature_predicate = new FeaturePredicate() {
      // return true if boundary need fixing
      public boolean testPredicate(Feature feature)
      {
        if(!GeneUtils.isGFFEntry(entry_group) ||
           !feature.getKey().equals("gene"))
          return false;
        ChadoCanonicalGene chadoGene =
            ((GFFStreamFeature)feature.getEmblFeature()).getChadoGene();
        
        if(chadoGene == null)
          return false;
        if(GeneUtils.isBoundaryOK(chadoGene) > 0)
          return true;
        
        return !GeneUtils.isStrandOK(chadoGene);
      }
    };
    
    final String filter_name = "Gene Model Boundary";
    final FilteredEntryGroup filtered_entry_group =
        new FilteredEntryGroup (entry_group, feature_predicate, filter_name);

    final FeatureListFrame feature_list_frame =
        new FeatureListFrame (filter_name,
                              selection, goto_source, filtered_entry_group,
                              base_plot_group);

    feature_list_frame.setVisible (true);
  }
  
  /**
   *  Popup a FeatureListFrame containing the features that have non-EMBL keys.
   *  @param parent_frame The parent JFrame.
   *  @param selection The Selection to pass to the FeatureList.
   *  @param entry_group The EntryGroup to pass to the FilteredEntryGroup.
   *  @param goto_source The GotoEventSource to pass to the FeatureList.
   *  @param base_plot_group The BasePlotGroup associated with this JMenu -
   *    needed to call getCodonUsageAlgorithm()
   **/
  protected static void showNonEMBLKeys(final JFrame parent_frame,
                                        final Selection selection,
                                        final EntryGroup entry_group,
                                        final GotoEventSource goto_source,
                                        final BasePlotGroup
                                        base_plot_group) 
  {
    final FeaturePredicate feature_predicate =
      new FeaturePredicate () {
        public boolean testPredicate (final Feature feature) {
          if (feature.hasValidEMBLKey ()) {
            return false;
          } else {
            return true;
          }
        }
      };

    final String filter_name =
      "features with a non-EMBL key (filtered from: " +
      parent_frame.getTitle () + ")";

    final FilteredEntryGroup filtered_entry_group =
      new FilteredEntryGroup (entry_group, feature_predicate, filter_name);

    final FeatureListFrame feature_list_frame =
      new FeatureListFrame (filter_name,
                            selection, goto_source, filtered_entry_group,
                            base_plot_group);

    feature_list_frame.setVisible (true);
  }


  /**
   *  Popup a FeatureListFrame containing the features that have the same key
   *  and location as another features (ie. duplicates).
   *  @param parent_frame The parent JFrame.
   *  @param selection The Selection to pass to the FeatureList.
   *  @param entry_group The EntryGroup to pass to the FilteredEntryGroup.
   *  @param goto_source The GotoEventSource to pass to the FeatureList.
   *  @param base_plot_group The BasePlotGroup associated with this JMenu -
   *    needed to call getCodonUsageAlgorithm()
   **/
  protected static void showDuplicatedFeatures (final JFrame parent_frame,
                                             final Selection selection,
                                             final EntryGroup entry_group,
                                             final GotoEventSource goto_source,
                                             final BasePlotGroup
                                               base_plot_group)
  {
    final FeaturePredicate feature_predicate =
      new FeaturePredicate () {
        public boolean testPredicate (final Feature feature) {
          final Entry feature_entry = feature.getEntry ();

          final int feature_index = feature_entry.indexOf (feature);

          if (feature_index + 1 == feature_entry.getFeatureCount ()) {
            // last in the Entry
            return false;
          }

          final Feature next_feature =
            feature_entry.getFeature (feature_index + 1);

          if (feature.getKey ().equals (next_feature.getKey ()) &&
              feature.getLocation ().equals (next_feature.getLocation ())) {
            return true;
          } else {
            return false;
          }
        }
      };

    final String filter_name =
      "duplicated Features (filtered from: " +
      parent_frame.getTitle () + ")";

    final FilteredEntryGroup filtered_entry_group =
      new FilteredEntryGroup (entry_group, feature_predicate, filter_name);

    final FeatureListFrame feature_list_frame =
      new FeatureListFrame (filter_name,
                            selection, goto_source, filtered_entry_group,
                            base_plot_group);

    feature_list_frame.setVisible (true);
  }

  /**
   *  Popup a FeatureListFrame containing those CDS features that overlap with
   *  the next feature.
   *  @param parent_frame The parent JFrame.
   *  @param selection The Selection to pass to the FeatureList.
   *  @param entry_group The EntryGroup to pass to the FilteredEntryGroup.
   *  @param goto_source The GotoEventSource to pass to the FeatureList.
   *  @param base_plot_group The BasePlotGroup associated with this JMenu -
   *    needed to call getCodonUsageAlgorithm()
   **/
  protected static void showOverlappingCDSs(final JFrame parent_frame,
                                            final Selection selection,
                                            final EntryGroup entry_group,
                                            final GotoEventSource goto_source,
                                            final BasePlotGroup base_plot_group) 
  {
    final Key key;
    if(GeneUtils.isDatabaseEntry(entry_group))
      key = new Key(DatabaseDocument.EXONMODEL);
    else
      key = Key.CDS;
    
    final FeatureKeyPredicate cds_predicate = new FeatureKeyPredicate (key);

    final FeaturePredicate feature_predicate =
      new FeaturePredicate () {
        public boolean testPredicate (final Feature test_feature) {
          if (!cds_predicate.testPredicate (test_feature)) {
            return false;
          }

          final Range feature_range = test_feature.getMaxRawRange ();

          final FeatureVector overlapping_features;

          try {
            overlapping_features =
              entry_group.getFeaturesInRange (feature_range);
          } catch (OutOfRangeException e) {
            throw new Error ("internal error - unexpected exception: " + e);
          }

          for (int i = 0 ; i < overlapping_features.size () ; ++i) {
            final Feature current_feature = overlapping_features.elementAt (i);

            if (current_feature != test_feature &&
                cds_predicate.testPredicate (current_feature)) {
              return true;
            }
          }

          return false;
        }
      };

    final String filter_name =
      "overlapping "+key.getKeyString()+" features (filtered from: " +
      parent_frame.getTitle () + ")";

    final FilteredEntryGroup filtered_entry_group =
      new FilteredEntryGroup (entry_group, feature_predicate, filter_name);

    final FeatureListFrame feature_list_frame =
      new FeatureListFrame (filter_name,
                            selection, goto_source, filtered_entry_group,
                            base_plot_group);

    feature_list_frame.setVisible (true);
  }

  /**
   *  Popup a FeatureListFrame containing those CDS features that overlap with
   *  the next feature.
   *  @param parent_frame The parent JFrame.
   *  @param selection The Selection to pass to the FeatureList.
   *  @param entry_group The EntryGroup to pass to the FilteredEntryGroup.
   *  @param goto_source The GotoEventSource to pass to the FeatureList.
   *  @param base_plot_group The BasePlotGroup associated with this JMenu -
   *    needed to call getCodonUsageAlgorithm()
   **/
  private static void
    showFeaturesWithSameStopCodons(final JFrame parent_frame,
                                   final Selection selection,
                                   final EntryGroup entry_group,
                                   final GotoEventSource goto_source,
                                   final BasePlotGroup base_plot_group) 
  {
    final Key key;
    if(GeneUtils.isDatabaseEntry(entry_group))
      key = new Key(DatabaseDocument.EXONMODEL);
    else
      key = Key.CDS;
    final FeatureKeyPredicate cds_predicate = new FeatureKeyPredicate (key);

    final FeaturePredicate feature_predicate =
      new FeaturePredicate () {
        public boolean testPredicate (final Feature test_feature) {
          if (!cds_predicate.testPredicate (test_feature)) {
            return false;
          }

          final Range feature_range = test_feature.getMaxRawRange ();

          final FeatureVector overlapping_features;

          try {
            overlapping_features =
              entry_group.getFeaturesInRange (feature_range);
          } catch (OutOfRangeException e) {
            throw new Error ("internal error - unexpected exception: " + e);
          }

          for (int i = 0 ; i < overlapping_features.size () ; ++i) {
            final Feature current_feature = overlapping_features.elementAt (i);

            if (current_feature == test_feature) {
              continue;
            }

            if (current_feature != test_feature &&
                cds_predicate.testPredicate (current_feature) &&
                (current_feature.getLastBase () ==
                 test_feature.getLastBase ()) &&
                (current_feature.getSegments ().lastElement ().getFrameID () ==
                 test_feature.getSegments ().lastElement ().getFrameID ())) {
              return true;
            }
          }
          return false;
        }
      };

    final String filter_name =
      key.getKeyString()+" features with the same stop codon as another (filtered from: " +
      parent_frame.getTitle () + ")";

    final FilteredEntryGroup filtered_entry_group =
      new FilteredEntryGroup (entry_group, feature_predicate, filter_name);

    final FeatureListFrame feature_list_frame =
      new FeatureListFrame (filter_name,
                            selection, goto_source, filtered_entry_group,
                            base_plot_group);

    feature_list_frame.setVisible (true);
  }

  /**
   *  Popup a FeatureListFrame containing the features that are missing
   *  required EMBL qualifiers.
   *  @param parent_frame The parent JFrame.
   *  @param selection The Selection to pass to the FeatureList.
   *  @param entry_group The EntryGroup to pass to the FilteredEntryGroup.
   *  @param goto_source The GotoEventSource to pass to the FeatureList.
   *  @param base_plot_group The BasePlotGroup associated with this JMenu -
   *    needed to call getCodonUsageAlgorithm()
   **/
  protected static void showMissingQualifierFeatures(final JFrame parent_frame,
                      final Selection selection, final EntryGroup entry_group,
                      final GotoEventSource goto_source,
                      final BasePlotGroup base_plot_group)
  {
    final FeaturePredicate feature_predicate =
      new FeaturePredicate () {
        public boolean testPredicate (final Feature feature) {
          if (feature.hasRequiredQualifiers ()) {
            return false;
          } else {
            return true;
          }
        }
      };

    final String filter_name =
      "features that are missing a required EMBL " +
      "qualifier (filtered from: " +
      parent_frame.getTitle () + ")";

    final FilteredEntryGroup filtered_entry_group =
      new FilteredEntryGroup (entry_group, feature_predicate, filter_name);

    final FeatureListFrame feature_list_frame =
      new FeatureListFrame (filter_name,
                            selection, goto_source, filtered_entry_group,
                            base_plot_group);

    feature_list_frame.setVisible (true);
  }

  /**
   *  Popup a FeatureListFrame containing only those features that have the
   *  key choosen by the user.
   *  @param parent_frame The parent JFrame.
   *  @param selection The Selection to pass to the FeatureList.
   *  @param entry_group The EntryGroup to pass to the FilteredEntryGroup.
   *  @param goto_source The GotoEventSource to pass to the FeatureList.
   *  @param base_plot_group The BasePlotGroup associated with this JMenu -
   *    needed to call getCodonUsageAlgorithm()
   **/
  private static void showFilterByKey (final JFrame parent_frame,
                                      final Selection selection,
                                      final EntryGroup entry_group,
                                      final GotoEventSource goto_source,
                                      final BasePlotGroup base_plot_group) 
  {
    final KeyChooser key_chooser = 
      new KeyChooser (Options.getArtemisEntryInformation (),
                                  new Key ("misc_feature"));

    key_chooser.getKeyChoice ().addItemListener (new ItemListener () {
      public void itemStateChanged (ItemEvent e) {
        if (e.getStateChange () == ItemEvent.SELECTED) {
          showFilterByKeyHelper (parent_frame,
                                 key_chooser.getKeyChoice ().getSelectedItem (),
                                 selection, entry_group, goto_source,
                                 base_plot_group);
          key_chooser.setVisible (false);
          key_chooser.dispose ();
        }
      }
    });

    key_chooser.getOKButton ().addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent _) {
        showFilterByKeyHelper (parent_frame,
                               key_chooser.getKeyChoice ().getSelectedItem (),
                               selection, entry_group, goto_source,
                               base_plot_group);
        key_chooser.setVisible (false);
        key_chooser.dispose ();
      }
    });

    key_chooser.setVisible (true);
  }

  /**
   *  Popup a FeatureListFrame containing only those features that have the
   *  key choosen by the user.
   *  @param parent_frame The parent JFrame.
   *  @param key The key to use in the filter.
   *  @param selection The Selection to pass to the FeatureList.
   *  @param entry_group The EntryGroup to pass to the FilteredEntryGroup.
   *  @param goto_source The GotoEventSource to pass to the FeatureList.
   *  @param base_plot_group The BasePlotGroup associated with this JMenu -
   *    needed to call getCodonUsageAlgorithm()
   **/
  private static void showFilterByKeyHelper (final JFrame parent_frame,
                      final Key key, final Selection selection,
                      final EntryGroup entry_group, 
                      final GotoEventSource goto_source,
                      final BasePlotGroup base_plot_group)
  {
    final String filter_name =
      "features with key: " + key + " (filtered from: " +
      parent_frame.getTitle () + ")";

    final FilteredEntryGroup filtered_entry_group =
      new FilteredEntryGroup (entry_group,
                              new FeatureKeyPredicate (key), filter_name);

    final FeatureListFrame feature_list_frame =
      new FeatureListFrame (filter_name,
                            selection, goto_source, filtered_entry_group,
                            base_plot_group);

    feature_list_frame.setVisible (true);
  }
  
  
  /**
   *  Popup a FeatureListFrame containing the features that have multiple
   *  systematic ID qualifiers (e.g. two or more systematic_id or locus_tag).
   *  @param parent_frame The parent JFrame.
   *  @param selection The Selection to pass to the FeatureList.
   *  @param entry_group The EntryGroup to pass to the FilteredEntryGroup.
   *  @param goto_source The GotoEventSource to pass to the FeatureList.
   *  @param base_plot_group The BasePlotGroup associated with this JMenu -
   *    needed to call getCodonUsageAlgorithm()
   **/
  protected static void showFilterByMultipleID(final JFrame parent_frame,
                      final Selection selection, final EntryGroup entry_group,
                      final GotoEventSource goto_source,
                      final BasePlotGroup base_plot_group)
  {
    final FeaturePredicate feature_predicate = new FeaturePredicate () 
    {
      public boolean testPredicate (final Feature feature) 
      {
        StringVector names = Options.getOptions().getSystematicQualifierNames();
        
        for(int i=0; i<names.size(); i++)
        {
          try
          {
            Qualifier qualifier = feature.getQualifierByName((String) names.get(i));
            if(qualifier != null && qualifier.getValues().size() > 1)
              return true;
          }
          catch (InvalidRelationException e){}
        }
        return false;
      }
    };

    final String filter_name =
      "features that have multiple systematic name qualifiers " +
      "(filtered from: " + parent_frame.getTitle () + ")";

    final FilteredEntryGroup filtered_entry_group =
      new FilteredEntryGroup (entry_group, feature_predicate, filter_name);

    final FeatureListFrame feature_list_frame =
      new FeatureListFrame (filter_name,
                            selection, goto_source, filtered_entry_group,
                            base_plot_group);

    feature_list_frame.setVisible (true);
  }

  

  /**
   *  Popup a FeatureListFrame containing only those features that are
   *  currently in the Selection.
   *  @param parent_frame The parent JFrame.
   *  @param selection The Selection to pass to the FeatureList and to use to
   *    filter the entry_group.
   *  @param entry_group The EntryGroup to pass to the FilteredEntryGroup.
   *  @param goto_source The GotoEventSource to pass to the FeatureList.
   *  @param base_plot_group The BasePlotGroup associated with this JMenu -
   *    needed to call getCodonUsageAlgorithm()
   **/
  private static void showFilterBySelection (final JFrame parent_frame,
                      final Selection selection, final EntryGroup entry_group,
                      final GotoEventSource goto_source,
                      final BasePlotGroup base_plot_group)
  {
    final FeaturePredicate predicate =
      new FeatureFromVectorPredicate (selection.getAllFeatures ());

    final String filter_name =
      "features from the selection (filtered from: " +
      parent_frame.getTitle () + ")";

    final FilteredEntryGroup filtered_entry_group =
      new FilteredEntryGroup (entry_group, predicate, filter_name);

    final FeatureListFrame feature_list_frame =
      new FeatureListFrame (filter_name,
                            selection, goto_source, filtered_entry_group,
                            base_plot_group);

    feature_list_frame.setVisible (true);
  }

  /**
   *  viewSelectedFeatures(), plotSelectedFeatures() etc. will only show this
   *  many features.
   **/
  private static final int MAXIMUM_SELECTED_FEATURES = 25;

  /**
   *  Open a view window for each of the selected features.  The viewer will
   *  listen for feature change events and update itself.
   *  @param frame The JFrame to use for MessageDialog components.
   *  @param selection The Selection containing the features to merge.
   **/
  protected static void viewSelectedFeatures (final JFrame frame,
                                    final Selection selection) 
  {
    final FeatureVector features_to_view = selection.getAllFeatures ();

    if (features_to_view.size () > MAXIMUM_SELECTED_FEATURES) {
      new MessageDialog (frame, "warning: only viewing the first " +
                         MAXIMUM_SELECTED_FEATURES + " selected features");
    }

    for (int i = 0 ;
         i < features_to_view.size () && i < MAXIMUM_SELECTED_FEATURES ;
         ++i) {
      final Feature selection_feature = features_to_view.elementAt (i);

      new FeatureViewer (selection_feature);
    }
  }

  /**
   *  Open an plot viewing window for each of the selected features.  The plot
   *  viewer will listen for feature change events and update itself.
   *  @param frame The JFrame to use for MessageDialog components.
   *  @param selection The Selection containing the features to plot.
   **/
  protected static void plotSelectedFeatures (final JFrame frame,
                                    final Selection selection) 
  {
    final FeatureVector features_to_plot = selection.getAllFeatures ();

    if (features_to_plot.size () > MAXIMUM_SELECTED_FEATURES) {
      new MessageDialog (frame, "warning: only showing plots for the first " +
                         MAXIMUM_SELECTED_FEATURES + " selected features");
    }

    for(int i = 0;
        i < features_to_plot.size () && i < MAXIMUM_SELECTED_FEATURES;
        ++i) 
    {
      final Feature selection_feature = features_to_plot.elementAt (i);

      new FeaturePlotGroup (selection_feature);
    }
  }

  /**
   *  Show the output file from an external program (like fasta) for the
   *  selected Feature objects.  The name of the file to read is stored in a
   *  feature qualifier.  The qualifier used is the program name plus "_file".
   *  @param frame The JFrame to use for MessageDialog components.
   *  @param selection The Selection containing the features to merge.
   *  @param send_to_browser if true the results should be sent straight to
   *    the web browser rather than using a SearchResultViewer object.
   **/
  protected static void viewExternalResults (final JFrame frame,
                                   final Selection selection,
                                   final String program_name,
                                   final boolean send_to_browser) 
  {
    final FeatureVector features_to_view = selection.getAllFeatures ();

    if (features_to_view.size () > MAXIMUM_SELECTED_FEATURES) {
      new MessageDialog (frame, "warning: only viewing results from " +
                         "the first " + MAXIMUM_SELECTED_FEATURES +
                         " selected features");
    }

    for (int i = 0 ;
         i < features_to_view.size () && i < MAXIMUM_SELECTED_FEATURES ;
         ++i) {
      final Feature this_feature = features_to_view.elementAt (i);

      Qualifier qualifier = null;
      try
      {
        qualifier = this_feature.getQualifierByName(program_name + "_file");
      }
      catch(InvalidRelationException e1)
      {
        // TODO Auto-generated catch block
        e1.printStackTrace();
      }
      
      if(qualifier == null)
        continue;
      StringVector qualifier_values = qualifier.getValues();
      String qualifier_value;
      if(qualifier_values == null || qualifier_values.size() == 0) {
        new MessageDialog (frame,
                           "Message",
                           "No " + program_name + " results for " +
                           this_feature.getIDString ());
        continue;
      }
      
      for(int j=0; j<qualifier_values.size(); j++)
      {
        qualifier_value = (String) qualifier_values.get(j);

        // strip off database name
        int ind;
        if((ind = qualifier_value.indexOf(':')) > -1)
          qualifier_value = qualifier_value.substring(ind + 1);

        String file_name = qualifier_value;

        // remove the program name string (if any) from the qualifier, so that
        // we can try to read the file with and without the prefix.
        if(file_name.startsWith(program_name + File.separatorChar)||
           file_name.startsWith(program_name + "/")) 
        {
          file_name = file_name.substring (program_name.length () + 1);
        }

       
        try
        {
          Document document = getSearchDocument(this_feature, program_name, file_name);

          if(document == null)
          {
            final String message_string = "No " + program_name
                + " results for " + this_feature.getIDString()
                + " (file not found: " + qualifier_value + ")";

            new MessageDialog(frame, message_string);

            continue;
          }

          if(send_to_browser)
          {
            String fileName = document.toString();
            if(document instanceof ZipFileDocument)
                fileName = ((ZipFileDocument)document).writeTmpFile(
                    getDocumentContents(document));
            
            SearchResultViewer.sendToBrowser(fileName);
          }
          else
          {
            new SearchResultViewer(program_name + " results for "
                + this_feature.getIDString() + " from " + document.getName(), document);
          }
        }
        catch(ExternalProgramException e)
        {
          new MessageDialog(frame, "error while open results file: " + e);
        }
        catch(IOException e)
        {
          new MessageDialog(frame, "error while open results file: " + e);
        }
      }
    }
  }
  
  /**
   * Get the corresponding fasta/blast search result. This looks initially
   * in program zip files, and then local file documents and finally, if a 
   * remote file connection has been made, it looks in the remote file system.
   * @param this_feature
   * @param program
   * @param fileName
   * @return
   * @throws IOException
   */
  public static Document getSearchDocument(Feature this_feature, String program, String fileName) 
          throws IOException
  {      
    Entry entry = this_feature.getEntry();
    Document root_document = entry.getRootDocument();

    if(root_document == null)
      root_document = new FileDocument(new File("."));
    Document document = null;

    final File dir_name = new File(program);
    File rootFile;
    if(root_document.getLocation() instanceof File)
      rootFile = (File)root_document.getLocation();
    else
      rootFile = new File(".");
    
    final Document[] possible_documents = new Document[] {
        new ZipFileDocument(new File(rootFile, program+File.separatorChar+program+".zip"), fileName),
        new ZipFileDocument(new File(rootFile, program+".zip"), fileName),
        new ZipFileDocument(new File(program+".zip"), fileName),
        new ZipFileDocument(new File(program, program+".zip"), fileName),
        root_document.append(program).append(fileName),
        root_document.append(fileName),
        new FileDocument(new File(fileName)),
        new FileDocument(dir_name).append(fileName),
        new FileDocument(new File(System.getProperty("user.dir")))
            .append(program).append(fileName),
        };

    for(int k = 0; k < possible_documents.length; ++k)
    {
      final Document this_document = possible_documents[k];
      if(this_document.readable())
      {
        document = this_document;
        break;
      }
      else
      {
        final File gzip_file = new File(this_document.toString() + ".gz");
        final Document gzip_document = new FileDocument(gzip_file);

        if(gzip_document.readable())
        {
          document = gzip_document;
          break;
        }
      }
    }

    // if not found locally try SSH remote site
    if(document == null && 
        ((DocumentEntry) (entry.getEMBLEntry())).getDocument() instanceof RemoteFileDocument)
    {
      File fdata = new File(program, fileName);
      // check on the remote side and scp the file over
      document = checkRemoteNode(this_feature, program, fileName, fdata.getParentFile());
    }
    return document;
  }

  private static String getDocumentContents(Document document) throws IOException
  {
    final BufferedReader buffered_reader = new BufferedReader(document.getReader());
    String line;
    final StringBuffer line_buffer = new StringBuffer();
    while((line = buffered_reader.readLine()) != null) 
      line_buffer.append(line).append('\n');

    buffered_reader.close();  
    return line_buffer.toString();
  }
  
  /**
   *  Open a FeatureInfo component for each of the selected features.  The
   *  new component will listen for feature change events and update itself.
   **/
  private void viewSelectedFeatureInfo () 
  {
    final FeatureVector features_to_view = getSelection ().getAllFeatures ();

    if (features_to_view.size () > MAXIMUM_SELECTED_FEATURES) {
      new MessageDialog (getParentFrame (),
                         "warning: only viewing the statistics for " +
                         "the first " + MAXIMUM_SELECTED_FEATURES +
                         " selected features");
    }

    for (int i = 0 ;
         i < features_to_view.size () && i < MAXIMUM_SELECTED_FEATURES ;
         ++i) {
      final Feature selection_feature = features_to_view.elementAt (i);

      new FeatureInfo (selection_feature,
                       base_plot_group.getCodonUsageAlgorithm ());
    }
  }

  /**
   *  View the bases of the selected features.  This creates a
   *  FeatureBaseViewer for each feature.
   *  @param include_numbers If true then the sequence will be numbered
   *    (every second line of the display will be numbers rather than
   *    sequence).
   **/
  private void viewSelectedBases (final boolean include_numbers, final boolean selectedExonsOnly) {
    if (getSelection ().isEmpty ()) {
      new MessageDialog (getParentFrame (), "Nothing selected");
      return;
    }

    if (selection.getMarkerRange () == null) {
      final FeatureVector fs = getSelection ().getAllFeatures ();

      if(selectedExonsOnly)
      {
        if (fs.size () > 1) 
          new MessageDialog (getParentFrame (),
               "warning: only viewing bases for the selected exons of one feature");
        new FeatureBaseViewer (fs.elementAt(0), include_numbers, 
            getSelection ().getSelectedSegments());
      }
      else
      {
        if (fs.size () > MAXIMUM_SELECTED_FEATURES) 
          new MessageDialog (getParentFrame (),
               "warning: only viewing bases for " +
               "the first " + MAXIMUM_SELECTED_FEATURES +
               " selected features");
        for (int i = 0; i < fs.size () && i < MAXIMUM_SELECTED_FEATURES; ++i)
          new FeatureBaseViewer (fs.elementAt (i), include_numbers, null);
      }
    } 
    else 
    {
      final SequenceViewer sequence_viewer =
        new SequenceViewer ("Selected bases", include_numbers);

      final String bases = getSelection ().getSelectionText ();
      sequence_viewer.setSequence (null, bases);
    }
  }

  /**
   *  View the bases of the selected CDS features.  This creates a
   *  FeatureBaseViewer for each CDS feature.
   *  @param include_numbers If true then the amino acids will be numbered
   *    (every second line of the display will be numbers rather than
   *    sequence).
   **/
  private void viewSelectedAminoAcids (final boolean include_numbers) {
    if (getSelection ().isEmpty ()) {
      new MessageDialog (getParentFrame (), "Nothing selected");
      return;
    }

    final MarkerRange range = selection.getMarkerRange ();

    if (range == null) {
      final FeatureVector features_to_view = getSelection ().getAllFeatures ();

      if (features_to_view.size () > MAXIMUM_SELECTED_FEATURES) {
        new MessageDialog (getParentFrame (),
                           "warning: only viewing amino acids for " +
                           "the first " + MAXIMUM_SELECTED_FEATURES +
                           " selected features");
      }

      for (int i = 0 ;
           i < features_to_view.size () && i < MAXIMUM_SELECTED_FEATURES ;
           ++i) {
        final Feature this_feature = features_to_view.elementAt (i);

        new FeatureAminoAcidViewer (this_feature, include_numbers);
      }
    } else {
      final SequenceViewer sequence_viewer =
        new SequenceViewer ("Selected bases (translated)", include_numbers);

      final String bases = getSelection ().getSelectionText ();

      final AminoAcidSequence amino_acids =
        AminoAcidSequence.getTranslation (bases, true);

      sequence_viewer.setSequence (null, amino_acids.toString ());
    }
  }

  /**
   *  Check that the number of features in the EntryGroup is less than the
   *  given number.  If it is less return true.  If not then popup a
   *  YesNoDialog asking for confirmation and return true if and only if the
   *  user chooses yes.
   **/
  private boolean checkEntryGroupSize (final int max_size) {
    final int feature_count = getEntryGroup ().getAllFeaturesCount ();

    if (feature_count < max_size) {
      return true;
    } else {
      final YesNoDialog dialog =
        new YesNoDialog (getParentFrame (),
                         "there are " + feature_count + " features in the " +
                         "active entries - continue?");

      return dialog.getResult ();
    }
  }

  /**
   *  Return the EntryGroup that was passed to the constructor.
   **/
  private EntryGroup getEntryGroup () 
  {
    return entry_group;
  }
  
  
  /**
   * Look for the file on the remote side and copy and transfer is
   * locally.
   * @param this_feature
   * @param program_name
   * @param file_name
   * @param dir_name
   * @return
   */
  private static Document checkRemoteNode(final Feature this_feature,
                               final String program,
                               String file_name,
                               final File dir_name)
  {
    RemoteFileDocument doc = (RemoteFileDocument) (((DocumentEntry) this_feature
        .getEntry().getEMBLEntry()).getDocument());

    RemoteFileNode node = doc.getRemoteFileNode();
    Document document = null;
    FileList flist = new FileList();
    String path = node.getPathName().concat(
        "/" + program + "/" + program + ".zip");
    FileAttributes attr = flist.stat(path);
    
    if(attr != null && attr.isFile())
    {
      File fn = flist.getZipEntryContents(path, file_name, dir_name);
      if(fn == null)
        return null;
      return new FileDocument(fn);
    }

    path = node.getPathName().concat(
          "/" + program + "/" + file_name);
    attr = flist.stat(path);

    if(attr == null || !attr.isFile())
    {
      file_name = file_name + ".gz";
      path = node.getPathName().concat(
        "/" + program + "/" + file_name);
      attr = flist.stat(path);
    }

    
    if(attr != null && attr.isFile())
    {
      FileTransferProgressMonitor monitor = new FileTransferProgressMonitor(
          null);
      FTProgress progress = monitor.add(node.getFile());

      RemoteFileNode fileNode = new RemoteFileNode("", file_name, null,
          node.getPathName().concat("/" + program), false);
      document = new RemoteFileDocument(fileNode);
      byte contents[] = fileNode.getFileContents(progress);
      File fn = new File(dir_name.getAbsoluteFile(), file_name);
      
      writeByteFile(contents, fn);
      ((RemoteFileDocument)document).setString(fn.getAbsolutePath());
      monitor.close();
    }
    return document;
  }
  
  private static boolean writeByteFile(byte[] contents, File fn)
  {
    if(fn.exists())
    {
      int n = JOptionPane.showConfirmDialog(null,
                                 "Overwrite \n"+fn.getName()+"?",
                                 "Overwrite File",
                                 JOptionPane.YES_NO_OPTION);
      if(n == JOptionPane.NO_OPTION)
        return false;
    }
    else if(!fn.getParentFile().canWrite())
    {
      if(!fn.getParentFile().mkdir())
        JOptionPane.showMessageDialog(null,"Cannot write "+fn.getName()+" to "+
                        fn.getParentFile().getAbsolutePath(),
                        "Write Permission Denied",
                        JOptionPane.WARNING_MESSAGE);
    }

    try
    {
      FileOutputStream out = new FileOutputStream(fn);
      out.write(contents);
      out.close(); 
    }
    catch(FileNotFoundException fnfe) { return false;}
    catch(IOException ioe) { return false;}
  
    return true;
  }


}
