/* MultiComparator.java
 *
 * created: Tue Sep 11 2001
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2001,2002  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/MultiComparator.java,v 1.22 2008-09-03 10:58:33 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.chado.ChadoTransactionManager;
import uk.ac.sanger.artemis.components.alignment.BamView;
import uk.ac.sanger.artemis.components.alignment.FileSelectionDialog;
import uk.ac.sanger.artemis.components.filetree.FileManager;
import uk.ac.sanger.artemis.components.genebuilder.GeneUtils;
import uk.ac.sanger.artemis.components.variant.VCFview;
import uk.ac.sanger.artemis.util.RemoteFileDocument;
import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.util.InputStreamProgressListener;
import uk.ac.sanger.artemis.io.EntryInformationException;
import uk.ac.sanger.artemis.io.DocumentEntryFactory;
import uk.ac.sanger.artemis.io.DocumentEntry;

import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Image;
import java.awt.MediaTracker;
import java.awt.Point;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.IOException;
import java.io.File;
import java.net.URL;
import java.util.List;

import javax.swing.JFrame;
import javax.swing.JMenuBar;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;


/**
 *  This JFrame component contains an arbitrary number of AlignmentViewer
 *  components and FeatureDisplay components along with ComparatorGlue objects
 *  to keep them synchronized.
 *  @author Kim Rutherford
 **/

public class MultiComparator extends JFrame 
{

  private static final long serialVersionUID = 1L;

  /**
   *  An array of EntryGroup objects to display(set by the constructor).  The
   *  array will contain at least one EntryGroup and will be exactly one
   *  element longer than comparison_data_array.
   **/
  private EntryGroup[] entry_group_array = null;

  /**
   *  An array of ComparisonData objects to display(set by the constructor).
   *  The array will be exactly one element shorter than entry_group_array.
   **/
  private ComparisonData[] comparison_data_array = null;

  /**
   *  An array of Selection objects - one per EntryGroup.  Created by
   *  makeSelectionArray()
   **/
  private Selection[] selection_array = null;

  /**
   *  An array of GotoEventSource objects - one per EntryGroup.  Created by
   *  makeGotoEventSourceArray()
   **/
  private GotoEventSource[] goto_event_source_array = null;

  /**
   *  An array of FeatureDisplay objects - one per EntryGroup.  Created in the
   *  constructor.
   **/
  private FeatureDisplay[] feature_display_array = null;

  private JPanel bamPanel[] = null;
  private JPanel vcfPanel[] = null;
  private Dimension dimensionAlignViewer = null;
  
  /**
   *  An array of AlignmentViewer objects - one per ComparisonData.  Created
   *  in the constructor.
   **/
  private AlignmentViewer[] alignment_viewer_array = null;

  /**
   *  An array of AlignmentViewer objects - one per ComparisonData.  Created
   *  in the constructor.
   **/
  private ComparatorGlue[] comparator_glue_array = null;

  /**
   *  An array of BasePlotGroup objects - one per FeatureDisplay/EntryGroup.
   *  Created in the constructor.
   **/
  private BasePlotGroup[] base_plot_group_array = null;

  /** This menu is populated by makeFileMenu(). */
  private final JMenu file_menu = new JMenu("File");

  /**
   *  The EntrySourceVector reference that is created in the constructor.
   **/
  private EntrySourceVector entry_sources;

  /** Used to show the progress of loading file. */
  private InputStreamProgressListener progress_listener;
  
  private GridBagLayout layout = new GridBagLayout();

  /**
   *  Initialise entry_group_array and comparison_data_array and create all
   *  the FeatureDisplay and AlignmentViewer components.
   *  @param entry_group_array The EntryGroup components to display.
   *  @param comparison_data_array The ComparisonData to display.  This array
   *    must have one element less than entry_group_array.
   *  @param progress_listener The object to which InputStreamProgressEvents
   *    will be send while reading.  Currently unused.
   **/
  public MultiComparator(final EntryGroup[] entry_group_array,
                         final ComparisonData[] comparison_data_array,
                         final InputStreamProgressListener
                            progress_listener) 
  {
    super();

    setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);

    if(entry_group_array.length != comparison_data_array.length + 1) 
      throw new Error("internal error - " +
                      "MultiComparator got illegal arguments");

    this.entry_group_array = entry_group_array;
    this.comparison_data_array = comparison_data_array;
    this.progress_listener = progress_listener;
    this.entry_sources = Utilities.getEntrySources(this, null);

    final StringBuffer title_buffer = new StringBuffer("ACT: ");
    for(int i = 0; i < entry_group_array.length - 1; ++i) 
      title_buffer.append(
          entry_group_array[i].getDefaultEntry().getName()).append(" vs ");

    final EntryGroup last_entry_group =
      entry_group_array[entry_group_array.length - 1];

    title_buffer.append(last_entry_group.getDefaultEntry().getName());

    setTitle(title_buffer.toString());

    makeSelectionArray();
    makeGotoEventSourceArray();

    feature_display_array =
      new FeatureDisplay[entry_group_array.length];
    base_plot_group_array =
      new BasePlotGroup[entry_group_array.length];
    bamPanel = 
      new JPanel[entry_group_array.length];
    vcfPanel =
      new JPanel[entry_group_array.length];
    alignment_viewer_array =
      new AlignmentViewer[comparison_data_array.length];
    comparator_glue_array =
      new ComparatorGlue[comparison_data_array.length];

    for(int i = 0; i < getEntryGroupArray().length; ++i) 
    {
      final EntryGroup this_entry_group = getEntryGroupArray()[i];

      final BasePlotGroup base_plot_group =
        new BasePlotGroup(this_entry_group, this,
                          getSelectionArray()[i],
                          getGotoEventSourceArray()[i]);

      final int scroll_bar_position;
      if(i == getEntryGroupArray().length - 1 &&
          getEntryGroupArray().length == 2) 
      {
        // put scrollbar at the top normally, but at the bottom of the lower
        // sequence if there's two entries on display
        scroll_bar_position = FeatureDisplay.SCROLLBAR_AT_BOTTOM;
      }
      else 
        scroll_bar_position = FeatureDisplay.SCROLLBAR_AT_TOP;

      feature_display_array[i] =
        new FeatureDisplay(this_entry_group,
                            getSelectionArray()[i],
                            getGotoEventSourceArray()[i],
                            base_plot_group,
                            scroll_bar_position);

      feature_display_array[i].setShowLabels(false);
      feature_display_array[i].setHardLeftEdge(false);

      if(getEntryGroupArray().length > 2) 
      {
        feature_display_array[i].setShowForwardFrameLines(false);
        feature_display_array[i].setShowReverseFrameLines(false);
      }

      feature_display_array[i].addDisplayAdjustmentListener(base_plot_group);

      base_plot_group_array[i] = base_plot_group;
      bamPanel[i] = new JPanel();
      vcfPanel[i] = new JPanel();
      this_entry_group.ref();
    }

    for(int i = 0 ; i < getComparisonDataArray().length ; ++i) 
    {
      alignment_viewer_array[i] =
        new AlignmentViewer(getFeatureDisplayArray()[i] ,
                             getFeatureDisplayArray()[i + 1],
                             getComparisonDataArray()[i]);

      comparator_glue_array[i] =
        new ComparatorGlue(this,
                            getFeatureDisplayArray()[i],
                            getFeatureDisplayArray()[i + 1],
                            getAlignmentViewerArray()[i]);
    }

    setFont(getDefaultFont());
    makeMenus();
    getContentPane().setLayout(layout);

    GridBagConstraints c = new GridBagConstraints();
    c.gridwidth  = GridBagConstraints.REMAINDER;
    c.fill       = GridBagConstraints.BOTH;
    c.anchor     = GridBagConstraints.NORTH;
    c.gridheight = 1;
    c.weightx    = 1;
    c.weighty    = 0;

    for(int i = 0; i < getFeatureDisplayArray().length; ++i) 
    {
      if(!(i == getEntryGroupArray().length - 1 &&
            getEntryGroupArray().length == 2)) 
      {
        // put graph above the sequence in this case
        c.weighty = 0;
        getContentPane().add(base_plot_group_array[i], c);
        c.weighty = 1;
        bamPanel[i].setVisible(false);
        vcfPanel[i].setVisible(false);
        getContentPane().add(bamPanel[i], c);
        getContentPane().add(vcfPanel[i], c);
      }

      c.weighty = 0;
      getContentPane().add(feature_display_array[i], c);

      if(i == getEntryGroupArray().length - 1 &&
            getEntryGroupArray().length == 2) 
      {
        // put graph below the sequence in this case
        c.weighty = 0;
        getContentPane().add(base_plot_group_array[i], c);
        c.weighty = 1;
        bamPanel[i].setVisible(false);
        vcfPanel[i].setVisible(false);
        getContentPane().add(bamPanel[i], c);
        getContentPane().add(vcfPanel[i], c);
      }

      if(i < getAlignmentViewerArray().length) 
      {
        c.fill = GridBagConstraints.BOTH;
        c.weighty = 0.5;
        getContentPane().add(getAlignmentViewerArray()[i], c);
      }
    }

    for(int i = 0 ; i < getEntryGroupArray().length ; ++i) 
    {
      final EntryGroupChangeListener change_listener =
                                 new EntryGroupChangeListener() 
      {
          public void entryGroupChanged(final EntryGroupChangeEvent event) 
          {
            switch(event.getType())
            {
              case EntryGroupChangeEvent.ENTRY_ADDED:
              case EntryGroupChangeEvent.ENTRY_DELETED:
                makeFileMenu();
                break;
            }
          }
        };

      getEntryGroupArray()[i].addEntryGroupChangeListener(change_listener);

      final EntryChangeListener entry_change_listener =
                                 new EntryChangeListener() 
      {
          public void entryChanged(final EntryChangeEvent event)
          {
            switch(event.getType()) 
            {
              case EntryChangeEvent.NAME_CHANGED:
                makeFileMenu();
                break;
            }
          }
        };

      getEntryGroupArray()[i].addEntryChangeListener(entry_change_listener);
    }

    addWindowListener(new WindowAdapter() 
    {
      public void windowClosing(WindowEvent event) 
      {
        closeComparator();
      }
    });

    final URL icon_url = Splash.class.getResource("/images/act.gif");
    if(icon_url != null) 
    {
      Toolkit toolkit = Toolkit.getDefaultToolkit();
      final Image icon_image = toolkit.getImage(icon_url);
      MediaTracker tracker = new MediaTracker(this);
      tracker.addImage(icon_image, 0);

      try 
      {
        tracker.waitForAll();
        setIconImage(icon_image);
      } 
      catch(InterruptedException e) 
      {
        // ignore and continue
      }
    }

    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();

    int screen_height = screen.height;
    int screen_width = screen.width;

    if(getAlignmentViewerArray().length > 0) 
    {
      if(screen_width <= 900 || screen_height <= 700) 
        setSize(screen_width * 9 / 10, screen_height * 9 / 10);
      else 
        setSize(900, 700);
    }

    setLocation(new Point((screen.width - getSize().width) / 2,
                          (screen.height - getSize().height) / 2));
    
    
    for(int i=0; i<getEntryGroupArray().length; i++)
      loadFromProperty(i);
  }
  
  /**
   * Load BAM/VCF files using the system properties flag -Dbam1, -Dbam2
   * for sequences 1, 2...
   * @param index the feature display to associate the files with
   */
  private void loadFromProperty(final int index)
  {
    final String bamProperty = "bam"+(index+1);
    if(System.getProperty(bamProperty) != null)
    {
      SwingWorker worker = new SwingWorker()
      {
        public Object construct()
        {
          MultiComparator.this.setCursor(new Cursor(Cursor.WAIT_CURSOR));
          String ngs[] = System.getProperty(bamProperty).split("[\\s,]");
          FileSelectionDialog fileChooser = new FileSelectionDialog(ngs);
          List<String> listBams = fileChooser.getFiles(".*\\.bam$");
          List<String> vcfFiles = fileChooser.getFiles(VCFview.VCFFILE_SUFFIX);
          loadBamAndVcf(listBams, vcfFiles, index);
          MultiComparator.this.setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
          return null;
        }
      };
      worker.start();
    }
  }

  /**
   *  Move all the sequences(and comparisons) up one place.  The first and
   *  last EntryGroups should be the same(reference).
   **/
  private void rotateSequences() 
  {
    final EntryGroup[] new_entry_groups =
      new EntryGroup[getEntryGroupArray().length];

    final ComparisonData[] new_comparison_data =
      new ComparisonData[getComparisonDataArray().length];

    for(int i = 1 ; i < new_entry_groups.length ; ++i) 
      new_entry_groups[i - 1] = getEntryGroupArray()[i];

    new_entry_groups[new_entry_groups.length - 1] = new_entry_groups[0];

    for(int i = 1; i < new_comparison_data.length; ++i) 
      new_comparison_data[i - 1] = getComparisonDataArray()[i];

    new_comparison_data[new_comparison_data.length - 1] =
      getComparisonDataArray()[0];

    final MultiComparator new_comparator =
      new MultiComparator(new_entry_groups, new_comparison_data,
                           progress_listener);

    setVisible(false);
    new_comparator.setVisible(true);
    closeComparator();
  }

  /**
   *  If there are no unsaved changes, close this EntryEdit.  Otherwise ask
   *  the user first.
   **/
  private void closeComparator() 
  {
    for(int i = 0 ; i < getEntryGroupArray().length ; ++i)
    {
      final EntryGroup entry_group = getEntryGroupArray()[i];

      if(entry_group.hasUnsavedChanges() &&
          entry_group.refCount() == 1) 
      {
        final YesNoDialog yes_no_dialog =
          new YesNoDialog(this,
                           "there are unsaved changes - really close?");

        if(!yes_no_dialog.getResult()) 
          return;
      }
    }

    setVisible(false);
    
    for(int i = 0 ; i < getEntryGroupArray().length ; ++i) 
    {
      final EntryGroup entry_group = getEntryGroupArray()[i];
      if(!GeneUtils.isDatabaseEntry(entry_group))
        entry_group.unref();
    }

    dispose();
  }

  /**
   *  Make one Selection object each EntryGroup.
   **/
  private void makeSelectionArray() 
  {
    final int length = getEntryGroupArray().length;
    selection_array = new Selection[length];
    for(int i = 0 ; i < length ; ++i) 
      selection_array[i] = new Selection(null);
  }

  /**
   *  Make subject_goto_event_source and query_goto_event_source and wire them
   *  togeather so that goto events do the right thing.
   **/
  private void makeGotoEventSourceArray() 
  {
    final int length = getEntryGroupArray().length;
    goto_event_source_array = new GotoEventSource[length];

    for(int i = 0 ; i < length ; ++i) 
    {
      goto_event_source_array[i] =
        new SimpleGotoEventSource(getEntryGroupArray()[i]) 
      {
          /**
           *  Send the given event to all the GotoListeners.
           **/
          public void sendGotoEvent(final GotoEvent goto_event) 
          {
            // temporarily disable moving so that the other FeatureDisplay
            // doesn't start moving(because of the listeners set up in
            // addDisplayListeners()
            super.sendGotoEvent(goto_event);
          }
        };
    }
  }

  /**
   *  Save the given entry, prompting for a file name if necessary.
   **/
  private void saveEntry(final Entry entry) 
  {
    if(entry.getName() == null) 
    {
      final EntryFileDialog file_dialog = new EntryFileDialog(this, true);

      file_dialog.saveEntry(entry, true, true,
                            true, DocumentEntryFactory.ANY_FORMAT);
    }
    else 
    {
      try
      {
        entry.save(DocumentEntryFactory.ANY_FORMAT);

        // save it back to ssh server
        if(((DocumentEntry)entry.getEMBLEntry()).getDocument()
                                   instanceof RemoteFileDocument)
        {
           RemoteFileDocument node =
               (RemoteFileDocument)(((DocumentEntry)entry.getEMBLEntry()).getDocument());

           File file = new File( ((DocumentEntry)entry.getEMBLEntry()).getDocument().getLocation().toString() );
           node.saveEntry(file);
        }

      }
      catch(IOException e)
      {
        new MessageDialog(this, "error while saving: " + e);
        return;
      } 
      catch(EntryInformationException e) 
      {
        new MessageDialog(this, "error while saving: " + e);
        return;
      }
    }
  }

  /**
   *  Return the name to use for the sub JMenu for the given EntryGroup.
   **/
  private String makeNewSubMenuName(final EntryGroup entry_group) 
  {
    final Entry sequence_entry = entry_group.getSequenceEntry();

    if(sequence_entry == null) 
      return "(no name)";
    else
    {
      final String sequence_name = sequence_entry.getName();
      if(sequence_name == null) 
        return "(no name)";
      else 
        return sequence_name;
    }

  }

  /**
   *  Make and add the menus for this component.
   **/
  private void makeMenus() 
  {
    //final Font default_font = getDefaultFont();
    final JMenuBar menu_bar = new JMenuBar();

    setJMenuBar(menu_bar);

    makeFileMenu();
    menu_bar.add(file_menu);

    final JMenu entries_menu = new JMenu("Entries");
    entries_menu.setMnemonic(KeyEvent.VK_N);
    menu_bar.add(entries_menu);
    final JMenu select_menu = new JMenu("Select");
    select_menu.setMnemonic(KeyEvent.VK_S);
    menu_bar.add(select_menu);
    final JMenu view_menu = new JMenu("View");
    view_menu.setMnemonic(KeyEvent.VK_V);
    menu_bar.add(view_menu);
    final JMenu goto_menu = new JMenu("Goto");
    goto_menu.setMnemonic(KeyEvent.VK_O);
    menu_bar.add(goto_menu);
    final JMenu edit_menu = new JMenu("Edit");
    edit_menu.setMnemonic(KeyEvent.VK_E);
    menu_bar.add(edit_menu);
    final JMenu create_menu = new JMenu("Create");
    create_menu.setMnemonic(KeyEvent.VK_C);
    menu_bar.add(create_menu);
    /*final JMenu write_menu = new JMenu("Write");
    write_menu.setMnemonic(KeyEvent.VK_W);
    menu_bar.add(write_menu);*/
    JMenu run_menu = new JMenu("Run");
    run_menu.setMnemonic(KeyEvent.VK_R);
    menu_bar.add(run_menu);
   
    final JMenu graph_menu = new JMenu("Graph");
    graph_menu.setMnemonic(KeyEvent.VK_G);
    menu_bar.add(graph_menu);

    for(int i = 0 ; i < getEntryGroupArray().length ; ++i)
    {
      final EntryGroup entry_group = getEntryGroupArray()[i];
      final String sub_menu_name = makeNewSubMenuName(entry_group);

      AlignmentViewer alignQueryViewer;
      if(i==0)
        alignQueryViewer = null;
      else
        alignQueryViewer = getAlignmentViewerArray()[i-1];

      AlignmentViewer alignSubjectViewer;
      if(i == getEntryGroupArray().length-1)
        alignSubjectViewer = null;
      else
        alignSubjectViewer = getAlignmentViewerArray()[i];

      final EntryGroupMenu this_entries_menu =
        new EntryGroupMenu(this,
                            getEntryGroupArray()[i],
                            sub_menu_name);
      entries_menu.add(this_entries_menu);

      final SelectMenu this_select_menu =
        new SelectMenu(this,
                        getSelectionArray()[i],
                        getGotoEventSourceArray()[i],
                        getEntryGroupArray()[i],
                        getBasePlotGroupArray()[i],
                        alignQueryViewer, alignSubjectViewer,
                        sub_menu_name);
      select_menu.add(this_select_menu);

      final ViewMenu this_view_menu =
        new ViewMenu(this,
                      getSelectionArray()[i],
                      getGotoEventSourceArray()[i],
                      getEntryGroupArray()[i],
                      getBasePlotGroupArray()[i],
                      sub_menu_name);
      view_menu.add(this_view_menu);

      final GotoMenu this_goto_menu =
        new GotoMenu(this,
                      getSelectionArray()[i],
                      getGotoEventSourceArray()[i],
                      getEntryGroupArray()[i],
                      sub_menu_name);
      goto_menu.add(this_goto_menu);

      AddMenu this_create_menu = null;

      if(Options.readWritePossible()) 
      {
        final EditMenu this_edit_menu =
          new EditMenu(this,
                        getSelectionArray()[i],
                        getGotoEventSourceArray()[i],
                        getEntryGroupArray()[i],
                        getBasePlotGroupArray()[i],
                        sub_menu_name,
                        getFeatureDisplayArray()[i]);
        edit_menu.add(this_edit_menu);

        this_create_menu =
          new AddMenu(this, getSelectionArray()[i],
                      getEntryGroupArray()[i],
                      getGotoEventSourceArray()[i],
                      getBasePlotGroupArray()[i],
                      alignQueryViewer, alignSubjectViewer,
                      sub_menu_name);
        create_menu.add(this_create_menu);

        final RunMenu this_run_menu =
            new RunMenu(this, getSelectionArray()[i],
                         sub_menu_name);
        run_menu.add(this_run_menu);
      }

      final GraphMenu this_graph_menu =
        new GraphMenu(this,
                       getEntryGroupArray()[i],
                       getBasePlotGroupArray()[i],
                       getFeatureDisplayArray()[i],
                       sub_menu_name, null, i+1);
      graph_menu.add(this_graph_menu);
    }

    final JMenuItem resize = new JMenuItem("Adjust panel heights ...");
    resize.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent arg0)
      {
        new ActPanelResizer(MultiComparator.this, layout);
      }
    });
    view_menu.add(resize);

    final JMenu display_menu = new JMenu("Display");
    display_menu.setMnemonic(KeyEvent.VK_D);

    final JMenuItem hide_on_frame_lines_item =
      new JMenuItem("Hide Frame Lines");
    hide_on_frame_lines_item.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        for(int i = 0 ; i < getFeatureDisplayArray().length ; ++i) 
        {
          getFeatureDisplayArray()[i].setShowForwardFrameLines(false);
          getFeatureDisplayArray()[i].setShowReverseFrameLines(false);
        }
      }
    });

    display_menu.add(hide_on_frame_lines_item);

    final JMenuItem show_on_frame_lines_item =
      new JMenuItem("Show Frame Lines");
    show_on_frame_lines_item.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        for(int i = 0 ; i < getFeatureDisplayArray().length ; ++i) 
        {
          getFeatureDisplayArray()[i].setShowForwardFrameLines(true);
          getFeatureDisplayArray()[i].setShowReverseFrameLines(true);
        }
      }
    });

    display_menu.add(show_on_frame_lines_item);

    menu_bar.add(display_menu);
  }

 
  /**
   *  Print menu 
   **/ 
  private void printMenu()
  {
    JMenuItem printImage = new JMenuItem("Save As Image Files (png/svg)...");
    printImage.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        PrintACT pact = new PrintACT(MultiComparator.this);
        pact.print();
      }
    });
    file_menu.add(printImage);
    
    JMenuItem printPS = new JMenuItem("Print...");
    printPS.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        PrintACT pact = new PrintACT(MultiComparator.this);
        pact.doPrintActions();
      }
    });
    file_menu.add(printPS);

// print preview
    JMenuItem printPreview = new JMenuItem("Print Preview");
    file_menu.add(printPreview);
    printPreview.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        PrintACT pact = new PrintACT(MultiComparator.this);
        pact.printPreview();
      }
    });
  }


  /**
   *  Make a new File menu replacing the current one(if any).
   **/
  private void makeFileMenu() 
  {
    file_menu.removeAll();
    file_menu.setMnemonic(KeyEvent.VK_F);
    final EntryGroup[] entry_group_array = getEntryGroupArray();

    JMenuItem popFileManager = new JMenuItem("Show File Manager ...");
    popFileManager.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)
      {
        if(ActMain.filemanager == null)
          ActMain.filemanager = new FileManager(MultiComparator.this,null);
        else
          ActMain.filemanager.setVisible(true);
      }
    });
    file_menu.add(popFileManager);

    if(entry_group_array[0] ==
       entry_group_array[entry_group_array.length - 1]) 
    {
      final JMenuItem rotate_button = new JMenuItem("Rotate Sequences");
      rotate_button.addActionListener(new ActionListener() 
      {
        public void actionPerformed(ActionEvent event) 
        {
          rotateSequences();
        }
      });

      file_menu.add(rotate_button);
      file_menu.addSeparator();
    }

    for(int i = 0; i < getEntryGroupArray().length; ++i) 
    {
      final EntryGroup entry_group = getEntryGroupArray()[i];
      final String new_menu_name = makeNewSubMenuName(entry_group);
      final JMenu entry_group_menu = new JMenu(new_menu_name);
      file_menu.add(entry_group_menu);

      for(int source_index = 0; source_index < entry_sources.size();
          ++source_index) 
      {
        final EntrySource this_source =
          entry_sources.elementAt(source_index);

        if(this_source.isFullEntrySource()) 
          continue;

        String entry_source_name = this_source.getSourceName();
        String menu_item_name = null;

        if(entry_source_name.equals("Filesystem")) 
          menu_item_name = "Read An Entry ...";
        else 
          menu_item_name = "Read An Entry From " + entry_source_name + " ...";

        final JMenuItem read_entry = new JMenuItem(menu_item_name);

        read_entry.addActionListener(new ActionListener() 
        {
          public void actionPerformed(ActionEvent event) 
          {
            readAnEntry(this_source,entry_group);
          }
        });

        entry_group_menu.add(read_entry);
      }

      entry_group_menu.addSeparator();

      if(entry_group == null || entry_group.size() == 0) 
      {
        // don't create a menu
      }
      else 
      {
        final JMenu save_entry_menu = new JMenu("Save Entry");

        for(int entry_index = 0; entry_index < entry_group.size();
            ++entry_index) 
        {
          final Entry this_entry = entry_group.elementAt(entry_index);
          String entry_name = this_entry.getName();

          if(entry_name == null) 
            entry_name = "no name";

          final ActionListener save_entry_listener =
            new ActionListener() 
          {
            public void actionPerformed(final ActionEvent event) 
            {
              MultiComparator.this.saveEntry(this_entry);
            }
          };

          final JMenuItem save_entry_item = new JMenuItem(entry_name);
          save_entry_item.addActionListener(save_entry_listener);
          save_entry_menu.add(save_entry_item);
        }

        entry_group_menu.add(save_entry_menu);

        if(GeneUtils.isDatabaseEntry(entry_group))
        {   
          final ChadoTransactionManager ctm = new ChadoTransactionManager();
          entry_group.addFeatureChangeListener(ctm);
          entry_group.addEntryChangeListener(ctm);
          entry_group.getBases().addSequenceChangeListener(ctm, 0);
          ctm.setEntryGroup(entry_group);

          final Selection selection = getSelectionArray()[i];
          final GotoEventSource goto_event_source = getGotoEventSourceArray()[i];
          final BasePlotGroup base_plot_group = getBasePlotGroupArray()[i];
            
          final JMenuItem commit_entry_item = new JMenuItem("Commit ...");
          commit_entry_item.addActionListener(new ActionListener() 
          {
            public void actionPerformed(final ActionEvent event) 
            {
              EntryEdit.commitToDatabase(entry_group, ctm, MultiComparator.this, 
                  selection, goto_event_source,base_plot_group);
            }
          });
          entry_group_menu.add(commit_entry_item);
        }
        
        final JMenuItem save_all_menu = new JMenuItem("Save All");

        save_all_menu.addActionListener(new ActionListener() 
        {
          public void actionPerformed(ActionEvent event) 
          {
            for(int entry_index = 0; entry_index < entry_group.size();
                ++entry_index) 
            {
              final Entry this_entry = entry_group.elementAt(entry_index);
              MultiComparator.this.saveEntry(this_entry);
            }
          }
        });

        entry_group_menu.add(save_all_menu);
        entry_group_menu.addSeparator();
        
        
        final WriteMenu this_write_menu =
          new WriteMenu(this,
                         getSelectionArray()[i],
                         getEntryGroupArray()[i]);
        entry_group_menu.add(this_write_menu);
        entry_group_menu.addSeparator();
      }

      JMenuItem read_bam_file = new JMenuItem("Read BAM / VCF ...");
      final int index = i;
      read_bam_file.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          FileSelectionDialog fileChooser = new FileSelectionDialog(
              null, false, "BAM / VCF View", "BAM / VCF");
          List<String> listBams = fileChooser.getFiles(".*\\.bam$");
          List<String> vcfFiles = fileChooser.getFiles(VCFview.VCFFILE_SUFFIX);
          
          loadBamAndVcf(listBams, vcfFiles, index);
        }
      });
      entry_group_menu.add(read_bam_file);
      
      final JMenuItem edit_subject = new JMenuItem("Edit In Artemis");
      edit_subject.addActionListener(new ActionListener() 
      {
        public void actionPerformed(ActionEvent event) 
        {
          final EntryEdit entry_edit = new EntryEdit(entry_group);
          entry_edit.setVisible(true);
        }
      });

      entry_group_menu.add(edit_subject);
    }

    file_menu.addSeparator();
    printMenu();
    file_menu.addSeparator();

    final JMenuItem close_button = new JMenuItem("Close");
    close_button.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        closeComparator();
      }
    });

    file_menu.add(close_button);
  }

  private void loadBamAndVcf(List<String> listBams, 
                             List<String> vcfFiles, 
                             int index)
  {
    final JPanel thisBamPanel = bamPanel[index];
    final JPanel thisVCFPanel = vcfPanel[index];
    final FeatureDisplay feature_display = feature_display_array[index];
    final EntryGroup entry_group = getEntryGroupArray()[index];
    
    if(listBams.size() > 0)
    {
      thisBamPanel.removeAll();
      thisBamPanel.setVisible(true);
    
      BamView bamView;
      try
      {
        bamView = new BamView(listBams, null, 2000, feature_display, 
            entry_group.getBases(), thisBamPanel, null);
      }
      catch(Exception ex)
      {
        JOptionPane.showMessageDialog(null,
          ex.getMessage(),
          "Error",
          JOptionPane.ERROR_MESSAGE);
        return;
      }

      if(dimensionAlignViewer == null)
        dimensionAlignViewer = alignment_viewer_array[0].getSize();
      thisBamPanel.setPreferredSize(new Dimension(500, dimensionAlignViewer.height/2));
      
      bamView.setDisplay(feature_display.getFirstVisibleForwardBase(), 
                       feature_display.getLastVisibleForwardBase(), null);
   
      bamView.getJspView().getVerticalScrollBar().setValue(
        bamView.getJspView().getVerticalScrollBar().getMaximum());
    
      feature_display.addDisplayAdjustmentListener(bamView);
      feature_display.getSelection().addSelectionChangeListener(bamView);
      MultiComparator.this.validate();
    }

    if (vcfFiles.size() > 0)
    {
      thisVCFPanel.removeAll();
      thisVCFPanel.setVisible(true);
      
      VCFview vcfView = new VCFview(null, thisVCFPanel, vcfFiles,
          feature_display.getMaxVisibleBases(), 1, null, null,
          null, feature_display);
      
      feature_display.addDisplayAdjustmentListener(vcfView);
      feature_display.getSelection().addSelectionChangeListener(vcfView);
      MultiComparator.this.validate();
    }  
  }
  
  /**
   *  Read an entry
   **/
  private void readAnEntry(final EntrySource this_source,
                           final EntryGroup entry_group)
  {
    try
    {
      final Entry new_entry = this_source.getEntry(entry_group.getBases(),
                                                   true);
      if(new_entry != null)
        entry_group.add(new_entry);
    }
    catch(final OutOfRangeException e)
    {
      new MessageDialog(MultiComparator.this,
             "read failed: one of the features " +
             "in the entry has an out of " +
             "range location: " +
             e.getMessage());
    }
    catch(final IOException e)
    {
      new MessageDialog(MultiComparator.this,
             "read failed due to an IO error: " +
             e.getMessage());
    }
  }

  /**
   *  Return the current default font(from Diana.options).
   **/
  private Font getDefaultFont() 
  {
    return Options.getOptions().getFont();
  }

  /**
   *  Return the EntryGroup objects that were passed to the constructor.
   **/
  protected EntryGroup[] getEntryGroupArray() 
  {
    return entry_group_array;
  }

  /**
   *  Return the ComparisonData objects that were passed to the constructor.
   **/
  private ComparisonData[] getComparisonDataArray() 
  {
    return comparison_data_array;
  }

  /**
   *  Return the AlignmentViewer objects that were created in the constructor.
   **/
  protected AlignmentViewer[] getAlignmentViewerArray() 
  {
    return alignment_viewer_array;
  }

  /**
   *  Return the FeatureDisplay objects that were created in the constructor.
   **/
  protected FeatureDisplay[] getFeatureDisplayArray() 
  {
    return feature_display_array;
  }

  /**
   *  Return the Selection objects that were created in the constructor.
   **/
  private Selection[] getSelectionArray() 
  {
    return selection_array;
  }

  /**
   *  Return the GotoEventSource objects that were created in the constructor.
   **/
  private GotoEventSource[] getGotoEventSourceArray() 
  {
    return goto_event_source_array;
  }

  /**
   *  Return the BasePlotGroup objects that were created in the constructor.
   **/
  protected BasePlotGroup[] getBasePlotGroupArray() 
  {
    return base_plot_group_array;
  }
  
  /**
   *  Return the Bam JPanel objects that were created in the constructor.
   **/
  protected JPanel[] getBamPanelArray() 
  {
    return bamPanel;
  }
  
  /**
   *  Return the VCF JPanel objects that were created in the constructor.
   **/
  protected JPanel[] getVcfPanelArray()
  {
    return vcfPanel;
  }

}
