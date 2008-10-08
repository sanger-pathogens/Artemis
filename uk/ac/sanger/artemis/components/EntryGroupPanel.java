/* EntryGroupPanel.java
 *
 * created: Wed Jun 21 2000
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/EntryGroupPanel.java,v 1.9 2008-10-08 15:31:48 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.components.genebuilder.GeneUtils;
import uk.ac.sanger.artemis.sequence.*;

import java.awt.Component;
import java.awt.event.*;
import javax.swing.*;

/**
 *  A JPanel that can show an EntryGroup(in some way).
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: EntryGroupPanel.java,v 1.9 2008-10-08 15:31:48 tjc Exp $
 **/

abstract public class EntryGroupPanel extends CanvasPanel 
{
  /** shortcut for View FASTA in browser.  */
  final int fasta_in_browser_key = KeyEvent.VK_F;

  /** shortcut for View FASTA.  */
  final int view_fasta_key = KeyEvent.VK_R;

  /** shortcut for View BLASTP in browser. */
  final int blastp_in_browser_key = KeyEvent.VK_B;

  /** shortcut for View BLASTP. */
  final int view_blastp_key = KeyEvent.VK_TAB;

  /** EntryGroup this component is displaying */
  private EntryGroup entry_group;

  /**
   *  This is a reference to the Selection object that was passed to the
   *  constructor.
   **/
  private Selection selection;

  /**
   *  This is a reference to the GotoEventSource object that was passed to the
   *  constructor.
   **/
  private GotoEventSource goto_event_source;

  /**
   *  The BasePlotGroup object that was passed to the constructor.
   **/
  private BasePlotGroup base_plot_group;

  /**
   *  Create a new EntryGroupPanel for the given EntryGroup.
   *  @param entry_group The EntryGroup that this component will display.
   *  @param selection The Selection object for this component.  Selected
   *    objects will be highlighted.
   *  @param base_plot_group The BasePlotGroup associated with this JMenu -
   *    needed to call getCodonUsageAlgorithm()
   **/
  public EntryGroupPanel(final EntryGroup entry_group,
                         final Selection selection,
                         final GotoEventSource goto_event_source,
                         final BasePlotGroup base_plot_group) 
  {
    this.entry_group       = entry_group;
    this.selection         = selection;
    this.goto_event_source = goto_event_source;
    this.base_plot_group   = base_plot_group;

    addKeyListener(new KeyAdapter() 
    {
      public void keyPressed(final KeyEvent event) 
      {
        handleKeyPress(event);
      }
    });
  }

  /**
   *  Return an object that implements the GotoEventSource interface and which
   *  is relative to the sequence that this DisplayComponent is displaying.
   **/
  public GotoEventSource getGotoEventSource() 
  {
    return goto_event_source;
  }

  /**
   *  Walk up through parents of this JComponent and return the JFrame at the
   *  top.
   **/
  public JFrame getParentFrame() 
  {
    Component parent = this.getParent();

    while(parent != null && !(parent instanceof JFrame)) 
      parent = parent.getParent();

    return (JFrame)parent;
  }

  /**
   *  Return the EntryGroup object that this FeatureDisplay is displaying.
   **/
  public EntryGroup getEntryGroup() 
  {
    return entry_group;
  }

  /**
   *  Returns the Selection that was passed to the constructor.
   **/
  public Selection getSelection() 
  {
    return selection;
  }

  /**
   *  Returns the BasePlotGroup that was passed to the constructor.
   **/
  public BasePlotGroup getBasePlotGroup()
  {
    return base_plot_group;
  }

  /**
   *  This method sends an GotoEvent to all the GotoEvent listeners that will
   *  make the start of the first selected Feature or FeatureSegment visible.
   **/
  protected void makeSelectionVisible() 
  {
    final Marker first_base = getSelection().getStartBaseOfSelection();

    final GotoEvent new_event = new GotoEvent(this, first_base);

    getGotoEventSource().sendGotoEvent(new_event);
  }

  /**
   *  Return true if and only if the given MouseEvent(a mouse press) should
   *  pop up a JPopupMenu.
   **/
  protected boolean isMenuTrigger(final MouseEvent event)
  {
    if(event.isPopupTrigger() ||
       (event.getModifiers() & InputEvent.BUTTON3_MASK) != 0) 
      return true;
    else 
      return false;
  }
  
  /**
   *  Handle key press events.
   **/
  private void handleKeyPress(final KeyEvent event) 
  {
    // this is done so that menu shortcuts don't cause each action to be
    // performed twice   
    if(!event.isShiftDown() && event.getModifiers() != 0)
    {
      if(getParentFrame() instanceof MultiComparator &&
         event.getModifiers() == InputEvent.ALT_MASK)
      {
        int event_key_code = event.getKeyCode();
        switch(event_key_code)
        {
        case EditMenu.TRIM_FEATURES_KEY_CODE:
          EditMenu.trimSelected(getParentFrame(), getSelection(),
              getEntryGroup(), false, true);
          break;
        case EditMenu.TRIM_FEATURES_TO_NEXT_ANY_KEY_CODE:
          EditMenu.trimSelected(getParentFrame(), getSelection(),
              getEntryGroup(), true, true);
          break;
        case EditMenu.EXTEND_TO_PREVIOUS_STOP_CODON_KEY_CODE:
          EditMenu.extendToORF(getParentFrame(), getSelection(),
              getEntryGroup(), false);
          break;
        case EditMenu.UNDO_KEY_CODE:
          EditMenu.undo(getParentFrame(), getSelection(), getEntryGroup());
          break;
        }
      }
      return;
    }
    
    final FeatureVector selected_features = getSelection().getAllFeatures();

    if(selected_features.size() == 0) 
    {
      final MarkerRange marker_range = getSelection().getMarkerRange();

      if(marker_range != null) 
      {
        switch(event.getKeyCode()) 
        {
          case AddMenu.CREATE_FROM_BASE_RANGE_KEY_CODE:
            if(!GeneUtils.isDatabaseEntry(entry_group))
              AddMenu.createFeatureFromBaseRange(getParentFrame(),
                                                 getSelection(),
                                                 entry_group,
                                                 getGotoEventSource());
            else
            {
              entry_group.getActionController ().startAction ();
              GeneUtils.createGeneModel(getParentFrame(), getSelection(),
                  entry_group, getGotoEventSource());
              entry_group.getActionController ().endAction ();
            }
            break;
        }
      }

    } 
    else
    {
      final Feature first_selected_feature = selected_features.elementAt(0);

      final int feature_index =
        getEntryGroup().indexOf(first_selected_feature);

      final int event_key_code = event.getKeyCode();

      switch(event_key_code) 
      {
        case KeyEvent.VK_UP:
        case KeyEvent.VK_DOWN:
        case KeyEvent.VK_LEFT:
        case KeyEvent.VK_RIGHT:
          {
            if((event_key_code == KeyEvent.VK_LEFT ||
                 event_key_code == KeyEvent.VK_RIGHT) &&
                !event.isShiftDown()) 
              break;
          
            final FeatureSegmentVector selected_segments =
                        getSelection().getSelectedSegments();

            if(event.isShiftDown()) 
            {
              if(selected_segments.size() > 0) 
              {
                final FeatureSegment selected_segment =
                               selected_segments.elementAt(0);

                final Feature selected_segment_feature =
                                selected_segment.getFeature();

                final FeatureSegmentVector feature_segments =
                       selected_segment_feature.getSegments();

                final int segment_index =
                  feature_segments.indexOf(selected_segment);

                if(event_key_code == KeyEvent.VK_UP ||
                   event_key_code == KeyEvent.VK_LEFT &&
                   !selected_segment_feature.isForwardFeature() ||
                   event_key_code == KeyEvent.VK_RIGHT &&
                   selected_segment_feature.isForwardFeature()) 
                {
                  final int segment_count = feature_segments.size();
                  final int new_index = segment_index + 1;
                  if(segment_index < segment_count - 1) 
                    selection.set(feature_segments.elementAt(new_index));
                } 
                else
                {
                  if(segment_index > 0) 
                  {
                    final int new_index = segment_index - 1;
                    selection.set(feature_segments.elementAt(new_index));
                  }
                }
              }
            } 
            else
            {
              if(event_key_code == KeyEvent.VK_DOWN) 
              {
                if(feature_index <
                    getEntryGroup().getAllFeaturesCount() - 1) 
                {
                  selection.set(getEntryGroup().featureAt(feature_index + 1));
                  makeSelectionVisible();
                }
              }
              else
              {
                if(feature_index > 0) 
                {
                  selection.set(getEntryGroup().featureAt(feature_index - 1));
                  makeSelectionVisible();
                }
              }
            }
          }
          break;
        case EditMenu.EDIT_FEATURES_KEY_CODE:
          EditMenu.editSelectedFeatures(getParentFrame(),
                                        getEntryGroup(),
                                        getSelection(),
                                        getGotoEventSource());
          break;
//      case EditMenu.UNDO_KEY_CODE:
//        EditMenu.undo(getParentFrame(), getSelection(), getEntryGroup());
//        break;
        case EditMenu.MERGE_FEATURES_KEY_CODE:
          EditMenu.mergeFeatures(getParentFrame(), getSelection(),
                                  getEntryGroup());
          break;
        case EditMenu.DUPLICATE_KEY_CODE:
          EditMenu.duplicateFeatures(getParentFrame(), getSelection(),
                                       getEntryGroup());
          break;
        case EditMenu.DELETE_FEATURES_KEY_CODE:
          EditMenu.deleteSelectedFeatures(getParentFrame(), getSelection(),
                                           getEntryGroup());
          break;
//      case EditMenu.TRIM_FEATURES_KEY_CODE:
//        EditMenu.trimSelected(getParentFrame(), getSelection(),
//                              getEntryGroup(), false, true);
//       break;
//      case EditMenu.TRIM_FEATURES_TO_NEXT_ANY_KEY_CODE:
//        EditMenu.trimSelected(getParentFrame(), getSelection(),
//                               getEntryGroup(), true, true);
//        break;
//      case EditMenu.EXTEND_TO_PREVIOUS_STOP_CODON_KEY_CODE:
//        EditMenu.extendToORF(getParentFrame(), getSelection(),
//                              getEntryGroup(), false);
//        break;
        case AddMenu.CREATE_FROM_BASE_RANGE_KEY_CODE:
          AddMenu.createFeatureFromBaseRange(getParentFrame(), getSelection(),
                                              entry_group,
                                              getGotoEventSource());
          break;
        case ViewMenu.VIEW_FEATURES_KEY_CODE:
          ViewMenu.viewSelectedFeatures(getParentFrame(), getSelection());
          break;
        case ViewMenu.PLOT_FEATURES_KEY_CODE:
          ViewMenu.plotSelectedFeatures(getParentFrame(), getSelection());
          break;
        case ViewMenu.OVERVIEW_KEY_CODE:
          new EntryGroupInfoDisplay(getParentFrame(), entry_group);
          break;
        case ViewMenu.FASTA_IN_BROWSER_KEY_CODE:
          viewResults("fasta", true);
          break;
        case ViewMenu.BLASTP_IN_BROWSER_KEY_CODE:
          viewResults("blastp", true);
          break;
        case ViewMenu.VIEW_FASTA_KEY_CODE:
          viewResults("fasta", false);
          break;
        case ViewMenu.VIEW_BLASTP_KEY_CODE:
          viewResults("blastp", false);
          break;
        case ViewMenu.VIEW_HTH_KEY_CODE:
          viewResults("hth", false);
          break;
      }
    }
  }

  /**
   *  Show the output file from an external program (like fasta) for the
   *  selected Feature objects.  The name of the file to read is stored in a
   *  feature qualifier.  The qualifier used is the program name plus "_file".
   *  @param send_to_browser if true the results should be sent straight to
   *    the web browser rather than using a SearchResultViewer object.
   **/
  private void viewResults(final String program_name,
                           final boolean send_to_browser) 
  {
    final boolean sanger_options =
      Options.getOptions().getPropertyTruthValue("sanger_options");

    if(sanger_options && send_to_browser) 
      ViewMenu.viewExternalResults(getParentFrame(), getSelection(),
                                   program_name, true);
    else 
      ViewMenu.viewExternalResults(getParentFrame(), getSelection(),
                                   program_name, false);
  }

}
