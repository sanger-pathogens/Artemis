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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/FeaturePopup.java,v 1.22 2008-01-17 09:57:06 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.components.genebuilder.GeneUtils;
import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.util.StringVector;

import java.io.*;
import java.awt.Component;
import java.awt.BorderLayout;
import java.awt.event.*;

import javax.swing.*;

/**
 *  FeaturePopup class
 *
 *  @author Kim Rutherford
 *  @version $Id: FeaturePopup.java,v 1.22 2008-01-17 09:57:06 tjc Exp $
 *
 **/

public class FeaturePopup extends JPopupMenu 
{
  /** */
  private static final long serialVersionUID = 1L;

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
  private BasePlotGroup base_plot_group = null;
  private JMenuItem feature_display_menus[] = null;
  private int SCORE_MINIMUM = Integer.MAX_VALUE;
  private int SCORE_MAXIMUM = Integer.MIN_VALUE;

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

    final JMenuItem action_menus[] = makeSubMenus();
    JMenuItem feature_list_menus[] = null;

    if(owner instanceof FeatureDisplay)
    {
      feature_display = (FeatureDisplay)owner;
      feature_display_menus = addFeatureDisplayItems();
      for(int i = 0; i<feature_display_menus.length; i++)
        if(!(feature_display_menus[i] instanceof JCheckBoxMenuItem))
          maybeAdd(feature_display_menus[i]);
    }
    else // must be a FeatureList
    {
      feature_list = (FeatureList)owner;
      feature_list_menus = addFeatureListItems();
      for(int i=0; i<feature_list_menus.length; i++)
        if(!(feature_list_menus[i] instanceof JCheckBoxMenuItem))
          maybeAdd(feature_list_menus[i]);
    }
    addSeparator();
    

    final JMenuItem miValidate = new JMenuItem("Validation report...");
    miValidate.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent arg0)
      {
        FeatureVector features = null;
        if(selection.getAllFeatures().size() < 1)
        {
          int status = JOptionPane.showConfirmDialog(owner.getParentFrame(), 
              "No features selected. Validate all features.", "Select Features", 
              JOptionPane.OK_CANCEL_OPTION);
          if(status == JOptionPane.CANCEL_OPTION)
            return;
        }
        else
          features = selection.getAllFeatures();
        new ValidateViewer(getEntryGroup(), features);
      }     
    });
    maybeAdd(miValidate);
    addSeparator();
    

    for(int i = 0; i<action_menus.length; i++)
      maybeAdd(action_menus[i]);

    addSeparator();
    if(owner instanceof FeatureDisplay) 
    {
      for(int i = 0; i<feature_display_menus.length; i++)
        if((feature_display_menus[i] instanceof JCheckBoxMenuItem))
          maybeAdd(feature_display_menus[i]);
    }
    else 
    {
      for(int i=0; i<feature_list_menus.length; i++)
        if((feature_list_menus[i] instanceof JCheckBoxMenuItem))
          maybeAdd(feature_list_menus[i]);
    }

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
   *  Create the Edit, Add and View sub menus.
   **/
  private JMenuItem[] makeSubMenus()
  {
    final JMenuItem[] action_menus = new JMenuItem[8];

    final JFrame frame = owner.getParentFrame();
    action_menus[0] = new EntryGroupMenu(frame, getEntryGroup());

    AlignmentViewer alignQueryViewer   = null;
    AlignmentViewer alignSubjectViewer = null;

    if(frame instanceof MultiComparator)
    {
      MultiComparator mc = (MultiComparator)frame;

      // determine which FeatureDisplay
      int ifeature_display;
      for(ifeature_display = 0; ifeature_display < mc.getEntryGroupArray().length;
          ++ifeature_display)
      {
        final EntryGroup this_entry_group = mc.getEntryGroupArray()[ifeature_display];
        if(this_entry_group == getEntryGroup())
          break;
      }

      if(ifeature_display==0)
        alignQueryViewer = null;
      else
        alignQueryViewer = mc.getAlignmentViewerArray()[ifeature_display-1];

      if(ifeature_display == mc.getEntryGroupArray().length-1)
        alignSubjectViewer = null;
      else
        alignSubjectViewer = mc.getAlignmentViewerArray()[ifeature_display];

      action_menus[1] =
          new SelectMenu(frame, selection,
                         getGotoEventSource(),
                         getEntryGroup(),
                         base_plot_group,
                         alignQueryViewer, alignSubjectViewer,
                         "Select");
    }
    else
      action_menus[1] = new SelectMenu(frame, selection,
                                       getGotoEventSource(),
                                       getEntryGroup(),
                                       base_plot_group);

    action_menus[2] = new ViewMenu(frame, selection,
                                  getGotoEventSource(),
                                  getEntryGroup(),
                                  base_plot_group);

    action_menus[3] = new GotoMenu(frame, selection,
                                  getGotoEventSource(),
                                  getEntryGroup());

    if(Options.readWritePossible()) 
    {
      action_menus[4] = new EditMenu(frame, selection,
                                    getGotoEventSource(),
                                    getEntryGroup(),
                                    base_plot_group, owner);
      if(entry_group instanceof SimpleEntryGroup) 
      {
        if(frame instanceof MultiComparator)
          action_menus[5] = new AddMenu(frame, selection,
                                    getEntryGroup(),
                                    getGotoEventSource(),
                                    base_plot_group,
                                    alignQueryViewer, alignSubjectViewer,
                                    "Create");
        else
          action_menus[5] = new AddMenu(frame, selection,
                                    getEntryGroup(),
                                    getGotoEventSource(),
                                    base_plot_group);
      }

      action_menus[6] = new WriteMenu(frame, selection,
                                      getEntryGroup());
      if(Options.isUnixHost()) 
        action_menus[7] = new RunMenu(frame, selection);
    }
    return action_menus;
  }

  /**
   *  Create those menu items that are relevant only to FeatureDisplay objects.
   **/
  private JMenuItem[] addFeatureDisplayItems() 
  {
    final JMenuItem[] feature_display_menus;
    
    final boolean isDatabaseGroup = GeneUtils.isDatabaseEntry( getEntryGroup() );
    final boolean isGFFGroup      = GeneUtils.isGFFEntry( getEntryGroup() );
    
    if(isDatabaseGroup || isGFFGroup)
      feature_display_menus = new JMenuItem[22];
    else
      feature_display_menus = new JMenuItem[21];
    
    feature_display_menus[0] = new JCheckBoxMenuItem("Start Codons");
    ((JCheckBoxMenuItem)feature_display_menus[0]).setState(
                               feature_display.getShowStartCodons());
    feature_display_menus[0].addItemListener(new ItemListener()
    {
      public void itemStateChanged(ItemEvent event) 
      {
        feature_display.setShowStartCodons(
                   ((JCheckBoxMenuItem)feature_display_menus[0]).getState());
      }
    });

    feature_display_menus[1] = new JCheckBoxMenuItem("Stop Codons");
    ((JCheckBoxMenuItem)feature_display_menus[1]).setState(
                               feature_display.getShowStopCodons());
    feature_display_menus[1].addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent event) 
      {
        feature_display.setShowStopCodons(
                  ((JCheckBoxMenuItem)feature_display_menus[1]).getState());
      }
    });

    feature_display_menus[2] = new JCheckBoxMenuItem("Feature Arrows");
    ((JCheckBoxMenuItem)feature_display_menus[2]).setState(
                               feature_display.getShowFeatureArrows());
    feature_display_menus[2].addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent event) 
      {
        feature_display.setShowFeatureArrows(
                  ((JCheckBoxMenuItem)feature_display_menus[2]).getState());
      }
    });

    feature_display_menus[3] = new JCheckBoxMenuItem("Feature Borders");
    ((JCheckBoxMenuItem)feature_display_menus[3]).setState(
                                   feature_display.getShowFeatureBorders());
    feature_display_menus[3].addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent event) 
      {
        feature_display.setShowFeatureBorders(
                  ((JCheckBoxMenuItem)feature_display_menus[3]).getState());
      }
    });

    feature_display_menus[4] = new JCheckBoxMenuItem("Feature Labels");
    ((JCheckBoxMenuItem)feature_display_menus[4]).setState(feature_display.getShowLabels());
    feature_display_menus[4].addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent e) 
      {
        feature_display.setShowLabels(
                  ((JCheckBoxMenuItem)feature_display_menus[4]).getState());
      }
    });

    feature_display_menus[5] = new JCheckBoxMenuItem("One Line Per Entry");
    ((JCheckBoxMenuItem)feature_display_menus[5]).setState(
                                 feature_display.getOneLinePerEntryFlag());
    feature_display_menus[5].addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent e) 
      {
        final boolean new_state =
                  ((JCheckBoxMenuItem)feature_display_menus[5]).getState();
        if(new_state && getEntryGroup().size() > 8) 
          feature_display.setShowLabels(false);
        feature_display.setOneLinePerEntry(new_state);
        if(new_state)
          ((JCheckBoxMenuItem)feature_display_menus[6]).setState(false);
      }
    });
    
    feature_display_menus[6] = new JCheckBoxMenuItem("Feature Stack View");
    ((JCheckBoxMenuItem)feature_display_menus[6]).setState(
                                 feature_display.getFeatureStackViewFlag());
    feature_display_menus[6].addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent e) 
      {
        final boolean new_state = 
                  ((JCheckBoxMenuItem)feature_display_menus[6]).getState();
        if(new_state && getEntryGroup().size() > 8) 
          feature_display.setShowLabels(false);
        feature_display.setFeatureStackViewFlag(new_state);
        if(new_state)
          ((JCheckBoxMenuItem)feature_display_menus[5]).setState(false);
      }
    });

    feature_display_menus[7] = new JCheckBoxMenuItem("Forward Frame Lines");
    ((JCheckBoxMenuItem)feature_display_menus[7]).setState(
                                feature_display.getShowForwardFrameLines());
    feature_display_menus[7].addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent e) 
      {
        feature_display.setShowForwardFrameLines(
                  ((JCheckBoxMenuItem)feature_display_menus[7]).getState());
      }
    });

    feature_display_menus[8] = new JCheckBoxMenuItem("Reverse Frame Lines");
    ((JCheckBoxMenuItem)feature_display_menus[8]).setState(
                               feature_display.getShowReverseFrameLines());
    feature_display_menus[8].addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent e) 
      {
        feature_display.setShowReverseFrameLines(
                  ((JCheckBoxMenuItem)feature_display_menus[8]).getState());
      }
    });

    feature_display_menus[9] = new JCheckBoxMenuItem("All Features On Frame Lines");
    ((JCheckBoxMenuItem)feature_display_menus[9]).setState(
                                            feature_display.getFrameFeaturesFlag());
    feature_display_menus[9].addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent e) 
      {
        feature_display.setFrameFeaturesFlag(
                  ((JCheckBoxMenuItem)feature_display_menus[9]).getState());
      }
    });

    feature_display_menus[10] = new JCheckBoxMenuItem("Show Source Features");
    ((JCheckBoxMenuItem)feature_display_menus[10]).setState(
                                    feature_display.getShowSourceFeatures());
    feature_display_menus[10].addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent e) 
      {
        feature_display.setShowSourceFeatures(
                  ((JCheckBoxMenuItem)feature_display_menus[10]).getState());
      }
    });

    feature_display_menus[11] = new JCheckBoxMenuItem("Flip Display");
    ((JCheckBoxMenuItem)feature_display_menus[11]).setState(
                                  feature_display.isRevCompDisplay());
    feature_display_menus[11].addItemListener(new ItemListener()
    {
      public void itemStateChanged(ItemEvent e) 
      {
        feature_display.setRevCompDisplay(
                  ((JCheckBoxMenuItem)feature_display_menus[11]).getState());
      }
    });

    feature_display_menus[12] = new JCheckBoxMenuItem("Colourise Bases");
    ((JCheckBoxMenuItem)feature_display_menus[12]).setState(
                                   feature_display.getShowBaseColours());
    feature_display_menus[12].addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent e) 
      {
        feature_display.setShowBaseColours(
                  ((JCheckBoxMenuItem)feature_display_menus[12]).getState());
      }
    });

    feature_display_menus[13] = new JMenuItem("Smallest Features In Front");
    feature_display_menus[13].addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent e) 
      {
        // clear the selection because selected features will always be on
        // top - which is not usually what is wanted
        selection.clear();
        feature_display.smallestToFront();
      }
    });

    feature_display_menus[14] = new JMenuItem("Set Score Cutoffs ...");
    feature_display_menus[14].addActionListener(new ActionListener() 
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

        FeatureVector features = getEntryGroup().getAllFeatures();
        int score;
        for(int i=0; i<features.size(); i++)
        {
          score = features.elementAt(i).getScore();
          if(score > -1)
          {
            if(score < SCORE_MINIMUM)
              SCORE_MINIMUM = score;
            else if(score > SCORE_MAXIMUM)
              SCORE_MAXIMUM = score;
          }
        }
        
        if(SCORE_MINIMUM == Integer.MAX_VALUE)
          SCORE_MINIMUM = 0;
        if(SCORE_MAXIMUM == Integer.MIN_VALUE)
          SCORE_MAXIMUM = 100;
        
        final ScoreChanger score_changer =
          new ScoreChanger("Score Cutoffs",
                            minimum_listener, maximum_listener,
                            SCORE_MINIMUM, SCORE_MAXIMUM);

        score_changer.setVisible(true);
      }
    });

    feature_display_menus[15] = new JMenuItem("Raise Selected Features");
    feature_display_menus[15].addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        raiseSelection();
      }
    });

    feature_display_menus[16] = new JMenuItem("Lower Selected Features");
    feature_display_menus[16].addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        lowerSelection();
      }
    });

    feature_display_menus[17] = new JMenuItem("Zoom to Selection");
    feature_display_menus[17].addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        zoomToSelection((FeatureDisplay) owner);
      }
    });

    feature_display_menus[18] = new JMenuItem("Select Visible Range");
    feature_display_menus[18].addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent e) 
      {
        selection.setMarkerRange(feature_display.getVisibleMarkerRange());
      }
    });

    feature_display_menus[19] = new JMenuItem("Select Visible Features");
    feature_display_menus[19].addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent e) 
      {
        selection.set(feature_display.getCurrentVisibleFeatures());
      }
    });

    feature_display_menus[20] = new JMenuItem("Frame Line Features ...");
    feature_display_menus[20].addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        final DefaultListModel listModel = new DefaultListModel();
        final Object protein_keys[] = ((FeatureDisplay)owner).getProteinKeys();

        for(int i=0; i<protein_keys.length; i++)
          listModel.addElement(protein_keys[i]);

        JPanel frame_keys = new JPanel(new BorderLayout());
        
        final JLabel label = new JLabel("Features Displayed on the Frame Lines:");
        final JList protein_list = new JList(listModel);
        final JScrollPane jsp    = new JScrollPane(protein_list);

        JButton remove_butt = new JButton("REMOVE");
        remove_butt.addActionListener(new ActionListener()
        {
          public void actionPerformed(ActionEvent e)
          {
            while(!protein_list.isSelectionEmpty())
              listModel.remove(protein_list.getSelectedIndex());
          }
        });

        Box bdown = Box.createVerticalBox();
        bdown.add(label);
        bdown.add(jsp);
        bdown.add(remove_butt);
        frame_keys.add(bdown, BorderLayout.CENTER);


        final KeyChoice key_choice =
           new KeyChoice(entry_group.elementAt(0).getEntryInformation(),
                         new uk.ac.sanger.artemis.io.Key("CDS"));

        JButton add_butt = new JButton("ADD");
        add_butt.addActionListener(new ActionListener()
        {
          public void actionPerformed(ActionEvent e)
          {
            listModel.addElement(key_choice.getSelectedItem().toString());
          }
        });

        bdown = Box.createVerticalBox();
        bdown.add(Box.createVerticalGlue());
        bdown.add(key_choice);
        bdown.add(add_butt);
        frame_keys.add(bdown, BorderLayout.EAST);

        int select = JOptionPane.showConfirmDialog(null, frame_keys,
                                "Frame Line Features ...",
                                 JOptionPane.OK_CANCEL_OPTION,
                                 JOptionPane.QUESTION_MESSAGE);

        if(select == JOptionPane.CANCEL_OPTION)
          return;
        
        ((FeatureDisplay)owner).setProteinKeys(listModel.toArray());
      }
    });
    
    if(isDatabaseGroup || isGFFGroup)
    {
      feature_display_menus[21] = new JMenuItem("Show/Hide Features ...");
      feature_display_menus[21].addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          GeneUtils.defineShowHideGeneFeatures(feature_display.getEntryGroup()
              .getAllFeatures());
        }
      });
    }
    
    return feature_display_menus;
  }


  /**
  *
  * Override show() 
  *
  */
  public void show(Component invoker, int x, int y)
  {
    if(owner instanceof FeatureDisplay)
    {
      if(selection_features.size() > 0 || selection_segments.size() > 0)
      {
        feature_display_menus[14].setEnabled(true);
        feature_display_menus[15].setEnabled(true);
//      feature_display_menus[16].setEnabled(true);
      }
      else
      {
        feature_display_menus[14].setEnabled(false);
        feature_display_menus[15].setEnabled(false);
//      feature_display_menus[16].setEnabled(false);
      }
    }

    super.show(invoker, x, y);
  }

  /**
   *  Create those menu items that are relevant only to FeatureList objects.
   **/
  private JMenuItem[] addFeatureListItems() 
  {
    final JMenuItem feature_list_menus[] = new JMenuItem[7];
    if(Options.readWritePossible()) 
    {
      feature_list_menus[0] = new JMenuItem("Save List To File ...");
      feature_list_menus[0].addActionListener(new ActionListener() 
      {
        public void actionPerformed(ActionEvent e) 
        {
          saveFeatureList();
        }
      });
    }

    feature_list_menus[1] = new JCheckBoxMenuItem("Show Correlation Scores");
    ((JCheckBoxMenuItem)feature_list_menus[1]).setState(
                                        feature_list.getCorrelationScores());
    feature_list_menus[1].addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent e) 
      {
        feature_list.setCorrelationScores(
                      ((JCheckBoxMenuItem)feature_list_menus[1]).getState());
      }
    });

    feature_list_menus[2] = new JCheckBoxMenuItem("Show Gene Names");
    ((JCheckBoxMenuItem)feature_list_menus[2]).setState(feature_list.getShowGenes());
    feature_list_menus[2].addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent e) 
      {
        boolean show_sysid = ((JCheckBoxMenuItem)feature_list_menus[3]).getState();
        boolean show_genes = ((JCheckBoxMenuItem)feature_list_menus[2]).getState();
        if(show_sysid && show_genes)
        {
          ((JCheckBoxMenuItem)feature_list_menus[3]).setState(false);
          feature_list.setShowSystematicID(false);
        }

        feature_list.setShowGenes(show_genes);
      }
    });

    feature_list_menus[3] = new JCheckBoxMenuItem("Show Systematic ID");
    ((JCheckBoxMenuItem)feature_list_menus[3]).setState(feature_list.getShowSysID());
    feature_list_menus[3].addItemListener(new ItemListener()
    {
      public void itemStateChanged(ItemEvent e)
      {
        boolean show_sysid = ((JCheckBoxMenuItem)feature_list_menus[3]).getState();
        boolean show_genes = ((JCheckBoxMenuItem)feature_list_menus[2]).getState();
        if(show_genes && show_sysid)
        {
          ((JCheckBoxMenuItem)feature_list_menus[2]).setState(false);
          feature_list.setShowGenes(false);
        }

        feature_list.setShowSystematicID(show_sysid);
      }
    });

    feature_list_menus[4] = new JCheckBoxMenuItem("Show Products");
    ((JCheckBoxMenuItem)feature_list_menus[4]).setState(feature_list.getShowProducts());
    feature_list_menus[4].addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent e) 
      {
        boolean show_products = ((JCheckBoxMenuItem)feature_list_menus[4]).getState();
        if(show_products) 
        {
          ((JCheckBoxMenuItem)feature_list_menus[5]).setState(false);
          feature_list.setShowQualifiers(false);
        }
        feature_list.setShowProducts(show_products);
      }
    });
    
    feature_list_menus[5] = new JCheckBoxMenuItem("Show Qualifiers");
    ((JCheckBoxMenuItem)feature_list_menus[5]).setState(feature_list.getShowQualifiers());
    feature_list_menus[5].addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent e) 
      {
        boolean show_qualifiers = ((JCheckBoxMenuItem)feature_list_menus[5]).getState();
        feature_list.setShowQualifiers(show_qualifiers);
        if(show_qualifiers) 
        { 
          ((JCheckBoxMenuItem)feature_list_menus[4]).setState(false);
          feature_list.setShowProducts(false);
        }
      }
    });
    
    feature_list_menus[6] = new JMenuItem("Show Selected Qualifiers ...");
    feature_list_menus[6].addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        final StringVector sv = feature_list.getShowUserDefinedQualifier();
        final Object display_names[];
        
        if(sv != null && sv.size() > 0)
          display_names = sv.toArray();
        else
        {
          display_names = new Object[1];
          display_names[0] = "note";
        }
        
        final String[] description = { 
        		"Qualifiers to display in feature list displays."};
        ListSelectionPanel displayListSelectionPanel =
          new ListSelectionPanel(entry_group, display_names,
        		  description, false, new Key("gene"), null);
        int select = JOptionPane.showConfirmDialog(null, 
            displayListSelectionPanel,
            "Select Qualifier to Display",
             JOptionPane.OK_CANCEL_OPTION,
             JOptionPane.QUESTION_MESSAGE);
        
        if(select == JOptionPane.CANCEL_OPTION)
          return;

        feature_list.setShowUserDefinedQualifier(
        		displayListSelectionPanel.getResultString());
      }
    });

    return feature_list_menus;
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
    final int status = file_dialog.showSaveDialog(frame);

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

}
