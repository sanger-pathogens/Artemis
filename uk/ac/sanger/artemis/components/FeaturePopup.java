/* FeaturePopup.java
 *
 * created: Wed Oct 21 1998
 *
 * This file is part of Artemis
 *
 * Copyright(C) 1998,1999,2000,2001,2002  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/FeaturePopup.java,v 1.1 2004-06-09 09:46:45 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.util.StringVector;

import java.io.*;
import java.awt.MenuItem;
import java.awt.CheckboxMenuItem;
import java.awt.event.*;
import javax.swing.*;

/**
 *  FeaturePopup class
 *
 *  @author Kim Rutherford
 *  @version $Id: FeaturePopup.java,v 1.1 2004-06-09 09:46:45 tjc Exp $
 *
 **/

public class FeaturePopup extends JPopupMenu 
{

  /**
   *  The reference of the EntryGroup object that was passed to the
   *  constructor.
   **/
  private EntryGroup entry_group;

  /**
   *  This is the Selection object that was passed to the constructor.
   **/
  final private Selection selection;

  /**
   *  This is a reference to the GotoEventSource object that was passed to the
   *  constructor.
   **/
  private GotoEventSource goto_event_source;

  /**
   *  The reference of the object that created this popup.
   **/
  private DisplayComponent owner;

  /**
   *  If the parent component of this popup is a FeatureDisplay then this will
   *  contain it's reference.
   **/
  private FeatureDisplay feature_display = null;

  /**
   *  If the parent component of this popup is a FeatureList then this will
   *  contain it's reference.
   **/
  private FeatureList feature_list = null;

  /**
   *  Set by the constructor to be the(possibly) empty vector of selected
   *  features.
   **/
  private FeatureVector selection_features;

  /**
   *  Set by the constructor to be the(possibly) empty vector of selected
   *  features.
   **/
  private FeatureSegmentVector selection_segments;

  private JCheckBoxMenuItem show_labels_item = null;
  private JCheckBoxMenuItem one_line_per_entry_item = null;
  private JCheckBoxMenuItem show_forward_frame_lines_item = null;
  private JCheckBoxMenuItem show_reverse_frame_lines_item = null;
  private JCheckBoxMenuItem show_start_codons_item = null;
  private JCheckBoxMenuItem show_stop_codons_item = null;
  private JCheckBoxMenuItem show_feature_arrows_item = null;
  private JCheckBoxMenuItem show_feature_borders_item = null;
  private JCheckBoxMenuItem frame_features_item = null;
  private JCheckBoxMenuItem source_features_item = null;
  private JCheckBoxMenuItem rev_comp_display_item = null;
  private JCheckBoxMenuItem base_colours_item = null;
  private JCheckBoxMenuItem correlation_scores_item = null;
  private JCheckBoxMenuItem show_genes_item = null;
  private JCheckBoxMenuItem show_products_item = null;
  private JCheckBoxMenuItem show_qualifiers_item = null;
  private JMenuItem entry_group_menu_item = null;
  private JMenuItem select_menu_item = null;
  private JMenuItem add_menu_item = null;
  private JMenuItem view_menu_item = null;
  private JMenuItem edit_menu_item = null;
  private JMenuItem goto_menu_item = null;
  private JMenuItem write_menu_item = null;
  private JMenuItem run_menu_item = null;
  private JMenuItem broadcast_item = null;
  private JMenuItem raise_feature_item = null;
  private JMenuItem lower_feature_item = null;
  private JMenuItem smallest_to_front_item = null;
  private JMenuItem zoom_to_selection_item = null;
  private JMenuItem score_cutoffs_item = null;
  private JMenuItem select_visible_range = null;
  private JMenuItem save_feature_list_item = null;
  private JMenuItem select_visible_features = null;

  private BasePlotGroup base_plot_group = null;

  /**
   *  Create a new FeaturePopup object.
   *  @param owner The component where this popup was popped up from.
   *  @param selection The selection to use for this popup.
   *  @param base_plot_group The BasePlotGroup associated with this JMenu -
   *    needed to call getCodonUsageAlgorithm()
   **/
  public FeaturePopup(final DisplayComponent owner,
                      final EntryGroup entry_group,
                      final Selection selection,
                      final GotoEventSource goto_event_source,
                      final BasePlotGroup base_plot_group) 
  {
    super(getMenuName(owner));

    this.owner = owner;
    this.entry_group = entry_group;
    this.selection = selection;
    this.goto_event_source = goto_event_source;
    this.base_plot_group = base_plot_group;

    selection_features = selection.getSelectedFeatures();
    selection_segments = selection.getSelectedSegments();

    makeSubMenus();

    addGenericItems();

    if(owner instanceof FeatureDisplay) 
    {
      feature_display =(FeatureDisplay) owner;
      addFeatureDisplayItems();
    } 
    else 
    {
      // must be a FeatureList
      feature_list =(FeatureList) owner;
      addFeatureListItems();
    }

    maybeAdd(raise_feature_item);
    maybeAdd(lower_feature_item);
    maybeAdd(smallest_to_front_item);
    maybeAdd(zoom_to_selection_item);
    maybeAdd(select_visible_range);
    maybeAdd(select_visible_features);
    maybeAdd(score_cutoffs_item);
    maybeAdd(save_feature_list_item);
    addSeparator();
    maybeAdd(entry_group_menu_item);
    maybeAdd(select_menu_item);
    maybeAdd(goto_menu_item);
    maybeAdd(view_menu_item);
    maybeAdd(edit_menu_item);
    maybeAdd(add_menu_item);
    maybeAdd(write_menu_item);
    maybeAdd(run_menu_item);
    addSeparator();
    maybeAdd(show_labels_item);
    maybeAdd(one_line_per_entry_item);
    maybeAdd(show_forward_frame_lines_item);
    maybeAdd(show_reverse_frame_lines_item);
    maybeAdd(show_start_codons_item);
    maybeAdd(show_stop_codons_item);
    maybeAdd(show_feature_arrows_item);
    maybeAdd(show_feature_borders_item);
    maybeAdd(frame_features_item);
    maybeAdd(source_features_item);
    maybeAdd(rev_comp_display_item);
    maybeAdd(base_colours_item);
    maybeAdd(correlation_scores_item);
    maybeAdd(show_genes_item);
    maybeAdd(show_qualifiers_item);
    maybeAdd(show_products_item);
    addSeparator();
    maybeAdd(broadcast_item);
  }

  /**
   *  Rename the name String to use for this JMenu.
   **/
  private static String getMenuName(final DisplayComponent owner) 
  {
    if(owner instanceof FeatureDisplay) 
      return "Feature Viewer JMenu";
    else
      return "Feature List JMenu";
  }

  /**
   *  Add an item only if it isn't null.
   **/
  private void maybeAdd(JMenuItem item) 
  {
    if(item != null) 
      add(item);
  }

  /**
   *  Create those menu items that are relevant to all components.
   **/
  private void addGenericItems() 
  {
    if(selection_features.size() > 0 || selection_segments.size() > 0) {

    }
  }

  /**
   *  Create the Edit, Add and View sub menus.
   **/
  public void makeSubMenus() 
  {
    final JFrame frame = owner.getParentFrame();
    entry_group_menu_item = new EntryGroupMenu(frame, getEntryGroup());

    select_menu_item = new SelectMenu(frame, selection,
                                      getGotoEventSource(),
                                      getEntryGroup(),
                                      base_plot_group);

    view_menu_item = new ViewMenu(frame, selection,
                                  getGotoEventSource(),
                                  getEntryGroup(),
                                  base_plot_group);

    goto_menu_item = new GotoMenu(frame, selection,
                                  getGotoEventSource(),
                                  getEntryGroup());

    if(Options.readWritePossible()) 
    {
      edit_menu_item = new EditMenu(frame, selection,
                                    getGotoEventSource(),
                                    getEntryGroup(),
                                    base_plot_group);
      if(entry_group instanceof SimpleEntryGroup) 
        add_menu_item = new AddMenu(frame, selection,
                                    getEntryGroup(),
                                    getGotoEventSource(),
                                    base_plot_group);

      write_menu_item = new WriteMenu(frame, selection,
                                      getEntryGroup());
      if(Options.isUnixHost()) 
        run_menu_item = new RunMenu(frame, selection);
    }
  }

  /**
   *  Create those menu items that are relevant only to FeatureDisplay objects.
   **/
  private void addFeatureDisplayItems() 
  {
    show_start_codons_item = new JCheckBoxMenuItem("Start Codons");
    show_start_codons_item.setState(feature_display.getShowStartCodons());
    show_start_codons_item.addItemListener(new ItemListener()
    {
      public void itemStateChanged(ItemEvent event) 
      {
        feature_display.setShowStartCodons(show_start_codons_item.getState());
      }
    });

    show_stop_codons_item = new JCheckBoxMenuItem("Stop Codons");
    show_stop_codons_item.setState(feature_display.getShowStopCodons());
    show_stop_codons_item.addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent event) 
      {
        feature_display.setShowStopCodons(show_stop_codons_item.getState());
      }
    });

    show_feature_arrows_item = new JCheckBoxMenuItem("Feature Arrows");
    show_feature_arrows_item.setState(feature_display.getShowFeatureArrows());
    show_feature_arrows_item.addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent event) 
      {
        feature_display.setShowFeatureArrows(show_feature_arrows_item.getState());
      }
    });

    show_feature_borders_item = new JCheckBoxMenuItem("Feature Borders");
    show_feature_borders_item.setState(feature_display.getShowFeatureBorders());
    show_feature_borders_item.addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent event) 
      {
        feature_display.setShowFeatureBorders(show_feature_borders_item.getState());
      }
    });

    show_labels_item = new JCheckBoxMenuItem("Feature Labels");
    show_labels_item.setState(feature_display.getShowLabels());
    show_labels_item.addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent e) 
      {
        feature_display.setShowLabels(show_labels_item.getState());
      }
    });

    one_line_per_entry_item = new JCheckBoxMenuItem("One Line Per Entry");
    one_line_per_entry_item.setState(feature_display.getOneLinePerEntryFlag());
    one_line_per_entry_item.addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent e) 
      {
        final boolean new_state = one_line_per_entry_item.getState();
        if(new_state && getEntryGroup().size() > 8) 
          feature_display.setShowLabels(false);
        feature_display.setOneLinePerEntry(new_state);
      }
    });

    show_forward_frame_lines_item =
      new JCheckBoxMenuItem("Forward Frame Lines");
    show_forward_frame_lines_item.setState(feature_display.getShowForwardFrameLines());
    show_forward_frame_lines_item.addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent e) 
      {
        feature_display.setShowForwardFrameLines(show_forward_frame_lines_item.getState());
      }
    });

    show_reverse_frame_lines_item =
      new JCheckBoxMenuItem("Reverse Frame Lines");
    show_reverse_frame_lines_item.setState(feature_display.getShowReverseFrameLines());
    show_reverse_frame_lines_item.addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent e) 
      {
        feature_display.setShowReverseFrameLines(show_reverse_frame_lines_item.getState());
      }
    });

    frame_features_item = new JCheckBoxMenuItem("All Features On Frame Lines");
    frame_features_item.setState(feature_display.getFrameFeaturesFlag());
    frame_features_item.addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent e) 
      {
        feature_display.setFrameFeaturesFlag(frame_features_item.getState());
      }
    });

    source_features_item = new JCheckBoxMenuItem("Show Source Features");
    source_features_item.setState(feature_display.getShowSourceFeatures());
    source_features_item.addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent e) 
      {
        feature_display.setShowSourceFeatures(source_features_item.getState());
      }
    });

    rev_comp_display_item = new JCheckBoxMenuItem("Flip Display");
 
    rev_comp_display_item.setState(feature_display.isRevCompDisplay());
    rev_comp_display_item.addItemListener(new ItemListener()
    {
      public void itemStateChanged(ItemEvent e) 
      {
        feature_display.setRevCompDisplay(rev_comp_display_item.getState());
      }
    });

    base_colours_item = new JCheckBoxMenuItem("Colourise Bases");
    base_colours_item.setState(feature_display.getShowBaseColours());
    base_colours_item.addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent e) 
      {
        feature_display.setShowBaseColours(base_colours_item.getState());
      }
    });

    smallest_to_front_item =
      new JMenuItem("Smallest Features In Front");
    smallest_to_front_item.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent e) 
      {
        // clear the selection because selected features will always be on
        // top - which is not usually what is wanted
        selection.clear();
        feature_display.smallestToFront();
      }
    });

    score_cutoffs_item = new JMenuItem("Set Score Cutoffs ...");

    score_cutoffs_item.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent e) 
      {
        final ScoreChangeListener minimum_listener =
          new ScoreChangeListener() 
          {
            public void scoreChanged(final ScoreChangeEvent event) 
            {
              feature_display.setMinimumScore(event.getValue());
            }
          };

        final ScoreChangeListener maximum_listener =
          new ScoreChangeListener() 
          {
            public void scoreChanged(final ScoreChangeEvent event) 
            {
              feature_display.setMaximumScore(event.getValue());
            }
          };

        final ScoreChanger score_changer =
          new ScoreChanger("Score Cutoffs",
                            minimum_listener, maximum_listener,
                            0, 100);

        score_changer.setVisible(true);
      }
    });

    if(selection_features.size() > 0 || selection_segments.size() > 0) 
    {
      raise_feature_item = new JMenuItem("Raise Selected Features");
      raise_feature_item.addActionListener(new ActionListener() 
      {
        public void actionPerformed(ActionEvent event) 
        {
          raiseSelection();
        }
      });

      lower_feature_item = new JMenuItem("Lower Selected Features");
      lower_feature_item.addActionListener(new ActionListener() 
      {
        public void actionPerformed(ActionEvent event) 
        {
          lowerSelection();
        }
      });
    }

    if(!selection.isEmpty()) 
    {
      zoom_to_selection_item = new JMenuItem("Zoom to Selection");
      zoom_to_selection_item.addActionListener(new ActionListener() 
      {
        public void actionPerformed(ActionEvent event) 
        {
          zoomToSelection((FeatureDisplay) owner);
        }
      });
    }


    select_visible_range =
      new JMenuItem("Select Visible Range");
    select_visible_range.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent e) 
      {
        selection.setMarkerRange(feature_display.getVisibleMarkerRange());
      }
    });

    select_visible_features =
      new JMenuItem("Select Visible Features");
    select_visible_features.addActionListener(new ActionListener() {
      public void actionPerformed(ActionEvent e) {
        selection.set(feature_display.getCurrentVisibleFeatures());
      }
    });
  }

  /**
   *  Create those menu items that are relevant only to FeatureList objects.
   **/
  private void addFeatureListItems() 
  {
    if(Options.getOptions().readWritePossible()) 
    {
      save_feature_list_item = new JMenuItem("Save List To File ...");
      save_feature_list_item.addActionListener(new ActionListener() 
      {
        public void actionPerformed(ActionEvent e) 
        {
          saveFeatureList();
        }
      });
    }

    correlation_scores_item = new JCheckBoxMenuItem("Show Correlation Scores");
    correlation_scores_item.setState(feature_list.getCorrelationScores());
    correlation_scores_item.addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent e) 
      {
        feature_list.setCorrelationScores(correlation_scores_item.getState());
      }
    });

    show_genes_item = new JCheckBoxMenuItem("Show Gene Names");
    show_genes_item.setState(feature_list.getShowGenes());
    show_genes_item.addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent e) 
      {
        feature_list.setShowGenes(show_genes_item.getState());
      }
    });

    show_products_item = new JCheckBoxMenuItem("Show Products");
    show_products_item.setState(feature_list.getShowProducts());
    show_products_item.addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent e) 
      {
        if(show_products_item.getState()) 
          feature_list.setShowQualifiers(false);
        
        feature_list.setShowProducts(show_products_item.getState());
      }
    });
    
    show_qualifiers_item = new JCheckBoxMenuItem("Show Qualifiers");
    show_qualifiers_item.setState(feature_list.getShowQualifiers());
    show_qualifiers_item.addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent e) 
      {
        feature_list.setShowQualifiers(show_qualifiers_item.getState());
        if(show_qualifiers_item.getState()) 
          feature_list.setShowProducts(false);
      }
    });
  }

  /**
   *  Save the text of the feature list to a file.
   **/
  private void saveFeatureList() 
  {
    final JFrame frame = owner.getParentFrame();
    final StickyFileChooser file_dialog = new StickyFileChooser();

    file_dialog.setDialogTitle("Choose save file ...");
    file_dialog.setDialogType(JFileChooser.SAVE_DIALOG);
    final int status = file_dialog.showOpenDialog(frame);

    if(status != JFileChooser.APPROVE_OPTION ||
       file_dialog.getSelectedFile() == null) 
      return;

    final File write_file =
      new File(file_dialog.getCurrentDirectory(),
                file_dialog.getSelectedFile().getName());
    
    if(write_file.exists()) 
    {
      final YesNoDialog yes_no_dialog =
        new YesNoDialog(frame,
                         "this file exists: " + write_file +
                         " overwrite it?");
      if(yes_no_dialog.getResult()) 
      {
        // yes - continue
      }
      else 
        return;
    }
    
    try 
    {
      final PrintWriter writer =
        new PrintWriter(new FileWriter(write_file));
      
      final StringVector list_strings = feature_list.getListStrings();
      
      for(int i = 0 ; i < list_strings.size() ; ++i) 
        writer.println(list_strings.elementAt(i));
      
      writer.close();
    } 
    catch(IOException e) 
    {
      new MessageDialog(frame, "error while writing: " + e.getMessage());
    }
  }

  /**
   *  Raise the selected features. (FeatureDisplay only.)
   **/
  private void raiseSelection() 
  {
    final FeatureVector features_to_raise = selection.getAllFeatures();

    for(int i = 0 ; i < features_to_raise.size() ; ++i) 
    {
      final Feature selection_feature = features_to_raise.elementAt(i);
      feature_display.raiseFeature(selection_feature);
    }
  }

  /**
   *  Lower the selected features. (FeatureDisplay only.)
   **/
  private void lowerSelection() 
  {
    final FeatureVector features_to_lower = selection.getAllFeatures();

    for(int i = 0 ; i < features_to_lower.size() ; ++i) 
    {
      final Feature selection_feature = features_to_lower.elementAt(i);
      feature_display.lowerFeature(selection_feature);
    }
  }

  /**
   *  Zoom the FeatureDisplay to the selection.
   **/
  static void zoomToSelection(final FeatureDisplay feature_display) 
  {
    final Selection selection = feature_display.getSelection();

    if(selection.isEmpty()) 
      return;

    // why bother in this case?
    if(feature_display.getEntryGroup().getSequenceLength() < 1000) 
      return;

    int first_base;
    int last_base;

    final FeatureSegmentVector segments = selection.getSelectedSegments();

    if(segments.size() == 1) 
    {
      // special case - zoom to the feature instead
      first_base = segments.elementAt(0).getFeature().getRawFirstBase();
      last_base  = segments.elementAt(0).getFeature().getRawLastBase();
    }
    else
    {
      first_base = selection.getLowestBaseOfSelection().getRawPosition();
      last_base  = selection.getHighestBaseOfSelection().getRawPosition();
    }

    if(first_base < 250) 
      first_base = 250;
    else 
      first_base -= 250;

    last_base += 250;

    feature_display.setFirstAndLastBase(first_base, last_base);
  }

  /**
   *  Return the EntryGroup object that this FeatureDisplay is displaying.
   **/
  private EntryGroup getEntryGroup() 
  {
    return entry_group;
  }

  /**
   *  Return an object that implements the GotoEventSource interface and is
   *  for the sequence that this DisplayComponent is displaying.
   **/
  public GotoEventSource getGotoEventSource() 
  {
    return goto_event_source;
  }

  /**
   *  Return a new CheckboxMenuItem unless the VM is 1.2 or worse(ie. 1.3 or
   *  1.4) on GNU/Linux, in which case return a new ArtemisCheckboxMenuItem
   **/
  private MenuItem makeCheckboxMenuItem(final String item_name) 
  {
    if(Options.getOptions().isBuggyLinuxVM() &&
        Options.getOptions().getPropertyTruthValue("buggy_linux_vm_fix")) 
      return new ArtemisCheckboxMenuItem(item_name);
    else 
      return new CheckboxMenuItem(item_name);
  }

  /**
   *  Add a ItemListener to a ArtemisCheckboxMenuItem or a CheckboxMenuItem.
   **/
  private void addMenuItemListener(final MenuItem menu_item,
                                    final ItemListener listener)
  {
    if(menu_item instanceof ArtemisCheckboxMenuItem) 
     ((ArtemisCheckboxMenuItem)menu_item).addItemListener(listener);
    else 
     ((CheckboxMenuItem)menu_item).addItemListener(listener);
  }

  /**
   *  Set the state of a ArtemisCheckboxMenuItem or a CheckboxMenuItem.
   **/
  private void setCheckboxMenuItemState(final MenuItem menu_item,
                                         final boolean state) 
  {
    if(menu_item instanceof ArtemisCheckboxMenuItem) 
     ((ArtemisCheckboxMenuItem)menu_item).setState(state);
    else 
     ((CheckboxMenuItem)menu_item).setState(state);
  }

  /**
   *  Get the state of a ArtemisCheckboxMenuItem or a CheckboxMenuItem.
   **/
  private boolean getCheckboxMenuItemState(final MenuItem menu_item) 
  {
    if(menu_item instanceof ArtemisCheckboxMenuItem) 
      return((ArtemisCheckboxMenuItem)menu_item).getState();
    else 
      return((CheckboxMenuItem)menu_item).getState();
  }

}
