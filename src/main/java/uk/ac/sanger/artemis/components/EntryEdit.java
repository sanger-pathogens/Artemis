/* EntryEdit.java
 *
 * created: Fri Oct  9 1998
 *
 * This file is part of Artemis
 *
 * Copyright(C) 1998,1999,2000,2001  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/EntryEdit.java,v 1.82 2009-09-24 12:42:16 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.chado.ChadoTransactionManager;
import uk.ac.sanger.artemis.chado.CommitFrame;
import uk.ac.sanger.artemis.circular.DNADraw;
import uk.ac.sanger.artemis.components.alignment.BamView;
import uk.ac.sanger.artemis.components.alignment.FileSelectionDialog;
import uk.ac.sanger.artemis.components.alignment.LookSeqPanel;
import uk.ac.sanger.artemis.components.filetree.FileList;
import uk.ac.sanger.artemis.components.filetree.FileManager;
import uk.ac.sanger.artemis.components.variant.VCFview;
import uk.ac.sanger.artemis.editor.BigPane;
import uk.ac.sanger.artemis.editor.FastaTextPane;
import uk.ac.sanger.artemis.editor.HitInfo;
import uk.ac.sanger.artemis.sequence.SequenceChangeEvent;
import uk.ac.sanger.artemis.sequence.SequenceChangeListener;

import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.util.ReadOnlyException;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.Document;
import uk.ac.sanger.artemis.util.IconManager;
import uk.ac.sanger.artemis.io.DatabaseInferredFeature;
import uk.ac.sanger.artemis.io.DocumentEntry;
import uk.ac.sanger.artemis.io.DocumentEntryFactory;
import uk.ac.sanger.artemis.io.EntryInformationException;
import uk.ac.sanger.artemis.io.DatabaseDocumentEntry;
import uk.ac.sanger.artemis.io.GenbankTblOutputStream;
import uk.ac.sanger.artemis.io.IndexFastaStream;
import uk.ac.sanger.artemis.io.InvalidRelationException;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Image;
import java.awt.Insets;
import java.awt.MediaTracker;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import java.io.IOException;

import javax.swing.*;

import com.sshtools.j2ssh.sftp.FileAttributes;

import java.util.Arrays;
import java.util.List;
import java.util.Vector;


/**
 *  Each object of this class is used to edit an EntryGroup object.
 *  @author Kim Rutherford
 *  @version $Id: EntryEdit.java,v 1.82 2009-09-24 12:42:16 tjc Exp $
 */
public class EntryEdit extends JFrame
    implements EntryGroupChangeListener, EntryChangeListener
{

  /** */
  private static final long serialVersionUID = 1L;

  /** The shortcut for Delete Selected Features. */
  final static KeyStroke SAVE_DEFAULT_KEY =
    KeyStroke.getKeyStroke(KeyEvent.VK_S,
                           Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()); // InputEvent.CTRL_MASK);

  /**
   *  A vector containing the Entry objects that this 
   *  EntryEdit object knows about.
   **/
  private EntryGroup entry_group;

  /**
   *  Created by the constructor to pass to those objects
   *  that are interested in GotoEvents.
   **/
  private GotoEventSource goto_event_source;

  private final JMenuBar menu_bar  = new JMenuBar();
  private final JMenu    file_menu = new JMenu("File");

  private EntryGroupDisplay group_display;
  private FeatureDisplay one_line_per_entry_display;
  private FeatureDisplay feature_display;
  private FeatureDisplay base_display;
  private BasePlotGroup  base_plot_group;
  private FeatureList feature_list;

  /** This Object contains the current selection. */
  private Selection selection = null;

  /** Alignment panel */
  private BamView bamView;
  private JPanel bamPanel;
  private VCFview vcfView;
  private JPanel vcfPanel;
  private JSplitPane lowerSplitPane;
  private JSplitPane ngSplitPane;
  
 /**
  *  The EntrySourceVector reference that is created in the constructor.
  **/
  private EntrySourceVector entry_sources;

  private SelectionInfoDisplay selection_info;
 
  private CommitButton commitButton;
  private static org.apache.log4j.Logger logger4j = 
    org.apache.log4j.Logger.getLogger(EntryEdit.class);
  
  /**
   *  Create a new EntryEdit object and JFrame.
   *  @param entry_group The EntryGroup object that this component is editing.
   */
  public EntryEdit(final EntryGroup entry_group) 
  {
    super("Artemis Entry Edit");
    setName("EntryEdit");

    setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
    entry_group.ref();
    this.entry_group = entry_group;

    // XXX add a InputStreamProgressListener
    this.entry_sources     = Utilities.getEntrySources(this, null);
    this.goto_event_source = new SimpleGotoEventSource(getEntryGroup());

    selection = new Selection(null);

    getEntryGroup().addFeatureChangeListener(selection);
    getEntryGroup().addEntryChangeListener(selection);
    getEntryGroup().addEntryGroupChangeListener(this);
    getEntryGroup().addEntryChangeListener(this);

    final Box box_panel = Box.createVerticalBox();
    final Box xBox = Box.createHorizontalBox();
    group_display = new EntryGroupDisplay(this);
    xBox.add(group_display);
    box_panel.add(xBox);
    
    if(getEntryGroup().getDefaultEntry() != null) 
    {
      final String name = getEntryGroup().getDefaultEntry().getName();
      if(name != null) 
        setTitle("Artemis Entry Edit: " + name);
      
      final JButton validate = new JButton("\u2713");
      validate.setToolTipText("Validate selected entries");
      validate.setFont(validate.getFont().deriveFont(18f));
      validate.setPreferredSize(new Dimension(25,25));
      validate.setMargin(new Insets(0, 0, 0, 0));
      validate.addActionListener(new ActionListener(){
        public void actionPerformed(ActionEvent arg0)
        {
          SwingUtilities.invokeLater(new Runnable() 
          {
            public void run() 
            {
              setCursor(new Cursor(Cursor.WAIT_CURSOR));
              new ValidateViewer(getEntryGroup(), 
                    getEntryGroup().getAllFeatures());
              setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
            }
          });
        }
      });
      xBox.add(validate);
      
      if(getEntryGroup().getDefaultEntry().getEMBLEntry() instanceof DatabaseDocumentEntry)
      {
        ChadoTransactionManager ctm = getDatabaseDocumentEntry().getChadoTransactionManager();
        
        getEntryGroup().addFeatureChangeListener(ctm);
        getEntryGroup().addEntryChangeListener(ctm);
        getEntryGroup().getBases().addSequenceChangeListener(ctm, 0);
        ctm.setEntryGroup(getEntryGroup());
        
        if(!getEntryGroup().getDefaultEntry().getEMBLEntry().isReadOnly())
        {
          commitButton = new CommitButton();
          getEntryGroup().addFeatureChangeListener(commitButton);
          getEntryGroup().addEntryChangeListener(commitButton);
          getEntryGroup().getBases().addSequenceChangeListener(commitButton, 0);
          xBox.add(commitButton);
        }

        if(DatabaseDocument.CHADO_INFER_CDS)
          DatabaseInferredFeature.addListenersToEntryGroup(getEntryGroup());
      }
    }

    final Font default_font = getDefaultFont();

    addWindowListener(new WindowAdapter() 
    {
      public void windowClosing(WindowEvent event) 
      {
        closeEntryEdit();
      }
    });

    getContentPane().setLayout(new BorderLayout());
    getContentPane().add(box_panel, "North");
    
    menu_bar.setFont(default_font);

    selection_info =
      new SelectionInfoDisplay(getEntryGroup(), getSelection());
    box_panel.add(selection_info);


    final boolean entry_buttons_option =
      Options.getOptions().getPropertyTruthValue("show_entry_buttons");

    group_display.setVisible(entry_buttons_option);


    final JPanel mainPanel = new JPanel(new BorderLayout());
    // set minimum size - this is then the smallest the split
    // pane divider can make this panel
    mainPanel.setMinimumSize(new Dimension(100,100));
    final Box main_box_panel = Box.createVerticalBox();
    mainPanel.add(main_box_panel, BorderLayout.NORTH);
    
    base_plot_group =
      new BasePlotGroup(getEntryGroup(), this, getSelection(),
                        getGotoEventSource());

    bamPanel = new JPanel();
    vcfPanel = new JPanel();
    ngSplitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT, 
        bamPanel, vcfPanel);
    ngSplitPane.setBorder(null);

    Dimension minimumSize = new Dimension(0, 0);
    bamPanel.setMinimumSize(minimumSize);
    vcfPanel.setMinimumSize(minimumSize);

    
    lowerSplitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT, 
        ngSplitPane, mainPanel);
    lowerSplitPane.setResizeWeight(0.);

    setNGDivider();
    
    final JSplitPane splitPane = new JSplitPane(JSplitPane.VERTICAL_SPLIT, 
                                             base_plot_group, lowerSplitPane);
    splitPane.setDividerSize(0);
    splitPane.setResizeWeight(0.);
    splitPane.setBorder(null);
    base_plot_group.setVisible(true);

    // lookseq read alignment
    LookSeqPanel lookseqPanel = null;
    JScrollPane jspLookSeq = null;
    if(Options.getOptions().getProperty("lookseq") != null)
    {
      lookseqPanel = new LookSeqPanel();
      jspLookSeq = new JScrollPane(lookseqPanel,
          JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
          JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
      jspLookSeq.setPreferredSize(new Dimension(jspLookSeq.getPreferredSize().width, 200));
      jspLookSeq.setVisible(false);
      main_box_panel.add(jspLookSeq);
    }


    // one line per entry display
    one_line_per_entry_display =
      new FeatureDisplay(getEntryGroup(), getSelection(),
                         getGotoEventSource(), base_plot_group);
    
    // read alignment panel
    //main_box_panel.add(bamPanel);
    
    one_line_per_entry_display.setShowLabels(false);
    one_line_per_entry_display.setOneLinePerEntry(true);
    one_line_per_entry_display.setVisible(false);

    // one line per entry expander button
    final JButton one_line_display_button = new JButton(">>");
    final Box one_line_button_box_across = setExpanderButton(one_line_display_button);
    one_line_display_button.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)
      {
        if(one_line_per_entry_display.isVisible())
          one_line_display_button.setText(">>");
        else
          one_line_display_button.setText("<<");
        one_line_per_entry_display.setVisible(!one_line_per_entry_display.isVisible());
        validate();
      }
    });
    main_box_panel.add(one_line_button_box_across);
    main_box_panel.add(one_line_per_entry_display);

    // feature display
    feature_display =
      new FeatureDisplay(getEntryGroup(), getSelection(),
                         getGotoEventSource(), base_plot_group);

    final Options options = Options.getOptions();

    if(options.getProperty("overview_feature_labels") != null) 
    {
      final boolean option_value =
        options.getPropertyTruthValue("overview_feature_labels");
      feature_display.setShowLabels(option_value);
    }

    if(options.getProperty("overview_one_line_per_entry") != null) 
    {
      final boolean option_value =
        options.getPropertyTruthValue("overview_one_line_per_entry");
      feature_display.setOneLinePerEntry(option_value);
    }
    
    if(options.getProperty("overview_feature_stack_view") != null) 
    {
      final boolean option_value =
        options.getPropertyTruthValue("overview_feature_stack_view");
      feature_display.setFeatureStackViewFlag(option_value);
    }

    feature_display.addDisplayAdjustmentListener(base_plot_group);
    feature_display.addDisplayAdjustmentListener(one_line_per_entry_display);

    one_line_per_entry_display.addDisplayAdjustmentListener(feature_display);

    // feature display expander button
    final JButton feature_display_button = new JButton("<<");
    final Box feature_display_button_box_across = setExpanderButton(feature_display_button);
    feature_display_button.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)
      {
        if(feature_display.isVisible())
          feature_display_button.setText(">>");
        else
          feature_display_button.setText("<<");
        feature_display.setVisible(!feature_display.isVisible());
      }
    });
    main_box_panel.add(feature_display_button_box_across);
    main_box_panel.add(feature_display);
    feature_display.setVisible(true);

    // base display
    base_display =
      new FeatureDisplay(getEntryGroup(), getSelection(),
                         getGotoEventSource(), base_plot_group);
    base_display.setShowLabels(false);
    base_display.setScaleFactor(0);

    // base display expander button
    final JButton base_display_button = new JButton("<<");
    final Box base_display_button_box_across = setExpanderButton(base_display_button);
    base_display_button.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)
      {
        if(base_display.isVisible())
          base_display_button.setText(">>");
        else
          base_display_button.setText("<<");
        base_display.setVisible(!base_display.isVisible());
      }
    });
    main_box_panel.add(base_display_button_box_across);
    main_box_panel.add(base_display);

    final boolean show_base_view;

    if(Options.getOptions().getProperty("show_base_view") != null) 
      show_base_view =
        Options.getOptions().getPropertyTruthValue("show_base_view");
    else 
      show_base_view = true;

    if(!show_base_view)
      base_display_button.setText(">>");
    base_display.setVisible(show_base_view);

    final JScrollPane jsp_feature_list;
    feature_list =
      new FeatureList(getEntryGroup(), getSelection(),
                      getGotoEventSource(), base_plot_group);
    jsp_feature_list = new JScrollPane(feature_list);
    feature_list.setFont(default_font);

    // feature list expander button
    final JButton scroll_button = new JButton(">>");
    final Box box_across = setExpanderButton(scroll_button);
    scroll_button.addActionListener(new ActionListener()
    {
      private int hgt = 0;
      public void actionPerformed(ActionEvent event)
      {
        Dimension dim = getSize();
        if(jsp_feature_list.isVisible())
        {
          Dimension dim_box = box_across.getPreferredSize();
          jsp_feature_list.setVisible(false);
          box_across.setPreferredSize(new Dimension(dim.width, dim_box.height)); 
          scroll_button.setText(">>");
          hgt = jsp_feature_list.getSize().height;
        }
        else
        {
          if(hgt == 0)
            hgt = getEntryGroup().getAllFeaturesCount() *
                                    feature_list.getLineHeight();
          
          jsp_feature_list.setPreferredSize(new Dimension(dim.width,hgt));
          jsp_feature_list.setVisible(true);
          scroll_button.setText("<<");
        }
        
        pack();
      }
    });
    main_box_panel.add(box_across);

    if(Options.getOptions().getPropertyTruthValue("show_list"))  
    {
      scroll_button.setText("<<");
      feature_list.setVisible(true);
    }
    else
    {
      scroll_button.setText(">>"); 
      feature_list.setVisible(false);
    }

    mainPanel.add(jsp_feature_list, "Center");
    getContentPane().add(splitPane, BorderLayout.CENTER);
    
    makeMenus(splitPane, jspLookSeq, lookseqPanel);
    pack();

    IconManager.setDockIcon(this, IconManager.ARTEMIS_NAME);

    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
    int screen_height = screen.height;
    int screen_width  = screen.width;

    if(screen_width <= 900 || screen_height <= 800)
      setSize(screen_width * 9 / 10, screen_height * 9 / 10);
    else
      setSize(900, 800);

    final int hgt = getEntryGroup().getAllFeaturesCount() * 
                    feature_list.getLineHeight();
    feature_list.setPreferredSize(new Dimension(getSize().width*4,hgt));
    jsp_feature_list.getVerticalScrollBar().setUnitIncrement(feature_list.getLineHeight());

    Utilities.centreFrame(this);
    
    if(System.getProperty("bam") != null || System.getProperty("bam1") != null)
    {
      SwingWorker worker = new SwingWorker()
      {
        public Object construct()
        {
          EntryEdit.this.setCursor(new Cursor(Cursor.WAIT_CURSOR));
          final String ngs[];
          final int idx;
          
          if(System.getProperty("bam") != null)
          {
            idx = 1;
            ngs = System.getProperty("bam").split("[\\s,]");
          }
          else
          {
            idx = 2;
            ngs = System.getProperty("bam1").split("[\\s,]");
            System.setProperty("bam1", "");
          }
          FileSelectionDialog fileChooser = new FileSelectionDialog(ngs);
          List<String> listBams = fileChooser.getFiles(".*\\.(bam|cram)$");
          final List<String> vcfFiles = fileChooser.getFiles(VCFview.VCFFILE_SUFFIX);
          loadBamAndVcf(listBams, vcfFiles);
          
          for(int i=idx; i<20; i++)
          {
            if(System.getProperty("bam"+i) != null)
            {
              fileChooser = new FileSelectionDialog(
                  System.getProperty("bam"+i).split("[\\s,]"));
              
              List<String> lBams = fileChooser.getFiles(".*\\.(bam|cram)$");
              if(lBams.size() > 0)
                bamView.openBamView(fileChooser.getFiles(".*\\.(bam|cram)$"));
              System.setProperty("bam"+i, "");
            }
          }

          if(System.getProperty("bamClone") != null)
          {
            int nclone = 2;
            try
            {
              nclone = Integer.parseInt(System.getProperty("bamClone"));
            }
            catch(NumberFormatException ne){}
            if(nclone > 10)
              nclone = 10;
            logger4j.debug("No. BamView clones = "+nclone+" bamClone = "+
                           System.getProperty("bamClone"));
            
            for(int i=1;i<nclone;i++)
              bamView.openBamView(listBams);
          }
          System.setProperty("bam", "");
          
          EntryEdit.this.setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
          return null;
        }
      };
      worker.start();
    }
  }

  private Box setExpanderButton(final JButton butt)
  {
    butt.setMargin(new Insets(0,1,0,1));
    butt.setHorizontalAlignment(SwingConstants.LEFT);
    butt.setPreferredSize(new Dimension( butt.getPreferredSize().width,9));
    butt.setBorderPainted(false);
    butt.setFont(new Font("SansSerif", Font.BOLD, 9));
    
    final Box box_across = Box.createHorizontalBox();
    box_across.add(butt);
    box_across.add(Box.createHorizontalGlue());
    return box_across;
  }
  
  protected void resetScrolls()
  {
    base_display.setFirstBase(1);
    base_display.fixScrollbar();
    feature_display.setFirstBase(1);
    feature_display.fixScrollbar();
    one_line_per_entry_display.setFirstBase(1);
    one_line_per_entry_display.fixScrollbar();
  }

  /**
  * Retrieve the feature list object.
  */
  protected FeatureList getFeatureList()
  {
    return feature_list;
  }


  /**
  * Retrieve the base display object.
  */
  protected FeatureDisplay getBaseDisplay()
  {
    return base_display;
  }

  /**
  * Retrieve the one line per entry object.
  */
  public FeatureDisplay getOneLinePerEntryDisplay()
  {
    return one_line_per_entry_display;
  }


  /**
  * Retrieve the base plot display object.
  */
  protected BasePlotGroup getBasePlotGroup()
  {
    return base_plot_group;
  }
  
  protected JPanel getBamPanel()
  {
    return bamPanel;
  }
  
  public BamView getJamView()
  {
    return bamView;
  }
  
  protected JPanel getVcfPanel()
  {
    return vcfPanel;
  }
  
  protected JPanel getVcfView()
  {
    return vcfView;
  }


  /**
  * Retrieve the entry group display object.
  */
  public EntryGroupDisplay getEntryGroupDisplay()
  {
    return group_display;
  }


  /**
  * Retrieve the selection info display object.
  */
  protected SelectionInfoDisplay getSelectionInfoDisplay()
  {
    return selection_info;
  }


  /**
  * Retrieve the feature display object.
  */
  protected FeatureDisplay getFeatureDisplay()
  {
    return feature_display;
  }

  /**
   *  If there are no unsaved changes, close this EntryEdit.  Otherwise ask
   *  the user first.
   **/
  private void closeEntryEdit() 
  {
    if(getEntryGroup().hasUnsavedChanges() &&
       getEntryGroup().refCount() == 1) 
    {
      final YesNoDialog yes_no_dialog =
        new YesNoDialog(EntryEdit.this,
                        "There are unsaved changes - really close?");

      if(!yes_no_dialog.getResult()) 
        return;
    }

    entryEditFinished();
  }

  /**
   *  Return an object that implements the GotoEventSource interface, and is
   *  the controlling object for Goto events associated with this object.
   **/
  public GotoEventSource getGotoEventSource()
  {
    return goto_event_source;
  }

  /**
   *  Return the EntryGroup object that was passed to the constructor.
   **/
  public EntryGroup getEntryGroup() 
  {
    return entry_group;
  }
  
  /**
   * Return the database document entry if present or null.
   * @return
   */
  private DatabaseDocumentEntry getDatabaseDocumentEntry() 
  {
    for(int i=0; i<entry_group.size(); i++)
      if(entry_group.elementAt(i).getEMBLEntry() instanceof DatabaseDocumentEntry)
        return (DatabaseDocumentEntry) entry_group.elementAt(i).getEMBLEntry();
    return null;
  }

  /**
   *  Returns a Selection object containing the selected features/exons.
   **/
  public Selection getSelection() 
  {
    return selection;
  }

  /**
   *  Implementation of the EntryGroupChangeListener interface.  We listen to
   *  EntryGroupChange events so that we can update the File menu when the
   *  EntryGroup changes.
   **/
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

  /**
   *  Implementation of the EntryChangeListener interface.  We listen to
   *  EntryChange events so that we can update the File menu when an Entry
   *  changes.
   **/
  public void entryChanged(final EntryChangeEvent event) 
  {
    if(event.getType() == EntryChangeEvent.NAME_CHANGED)
      makeFileMenu();
  }

  /* private members: */

  /**
   *  This method arranges for the EntryEdit JFrame to go away.  This EntryEdit
   *  object was created by the main program, so the main program must be the
   *  one to delete us.
   **/
  private void entryEditFinished() 
  {
    setVisible(false);
   
    // chado transaction manager
    if(getDatabaseDocumentEntry() != null)
    {
      ChadoTransactionManager ctm = getDatabaseDocumentEntry().getChadoTransactionManager();
      getEntryGroup().removeFeatureChangeListener(ctm);
      getEntryGroup().removeEntryChangeListener(ctm);
    }
    
    if(commitButton != null)
    {
      getEntryGroup().removeFeatureChangeListener(commitButton);
      getEntryGroup().removeEntryChangeListener(commitButton);
      commitButton.close();
    }
    
    getEntryGroup().removeFeatureChangeListener(selection);
    getEntryGroup().removeEntryChangeListener(selection);

    getEntryGroup().removeEntryGroupChangeListener(this);
    getEntryGroup().removeEntryChangeListener(this);

    getEntryGroup().unref();
    
    dispose();
    getEntryGroup().getBases().clearCodonCache();
    //getEntryGroup().getBases().getSequence().clear();
  }

  /**
   *  Write the default Entry in the EntryGroup to the file it came from.
   **/
  private void saveDefaultEntry() 
  {
    if(getEntryGroup().getDefaultEntry() == null) 
      new MessageDialog(EntryEdit.this, "There is no default entry");
    else 
      saveEntry(entry_group.getDefaultEntry(), true, false, true,
                 DocumentEntryFactory.ANY_FORMAT);
  }

  /**
   *  Save the given entry, prompting for a file name if necessary.
   *  @param include_diana_extensions If true then any diana additions to
   *    the embl file format will be included in the output, otherwise they
   *    will be removed.  Also possible problems that would cause an entry to
   *    bounce from the EMBL submission process will be flagged if this is
   *    true.
   *  @param ask_for_name If true then always prompt for a new filename,
   *    otherwise prompt only when the entry name is not set.
   *  @param keep_new_name If ask_for_name is true a file will be written with
   *    the new name the user selects - if keep_new_name is true as well, then
   *    the entry will have it's name set to the new name, otherwise it will
   *    be used for this save and then discarded.
   *  @param destination_type Should be one of EMBL_FORMAT, GENBANK_FORMAT or
   *    ANY_FORMAT.  If ANY_FORMAT then the Entry will be saved in the
   *    same format it was created, otherwise it will be saved in the given
   *    format.
   **/
  void saveEntry(final Entry entry,
                 final boolean include_diana_extensions,
                 final boolean ask_for_name, final boolean keep_new_name,
                 final int destination_type) 
  { 
    saveEntry(entry, include_diana_extensions,
         ask_for_name, keep_new_name,
         destination_type, true); 
  }
  
  private void saveEntry(final Entry entry,
      final boolean include_diana_extensions,
      final boolean ask_for_name, final boolean keep_new_name,
      final int destination_type, final boolean useSwingWorker) 
  {
    if(!include_diana_extensions) 
    {
      if(displaySaveWarnings(entry)) 
        return;
    }

    if(destination_type != DocumentEntryFactory.ANY_FORMAT &&
       entry.getHeaderText() != null) 
    {
      final YesNoDialog yes_no_dialog =
        new YesNoDialog(this, "header section will be lost. continue?");

      if(!yes_no_dialog.getResult()) 
        return;
    }

//  if(!System.getProperty("os.arch").equals("alpha"))
//  {
    if(entry.getEMBLEntry().getSequence() instanceof IndexFastaStream)
    {
      JOptionPane.showMessageDialog(null, 
          entry.getName()+" is an indexed sequence.\n"+
          "This cannot be written out.",
          "Write Option Not Available", 
          JOptionPane.WARNING_MESSAGE);
      return;
    }
    
/*    if(useSwingWorker)
    {
      SwingWorker worker = new SwingWorker()
      {
        public Object construct()
        {
          final EntryFileDialog file_dialog = new EntryFileDialog(
                                                      EntryEdit.this, false);
          file_dialog.saveEntry(entry, include_diana_extensions, ask_for_name,
                                keep_new_name, destination_type);

          return null;
        }
      };
      worker.start();
    }
    else*/
    {
      final EntryFileDialog file_dialog = new EntryFileDialog(this,
                                                              false);
      
      file_dialog.saveEntry(entry, include_diana_extensions, ask_for_name,
                          keep_new_name, destination_type);
    }
  }


  /**
   *  Save the changes to all the Entry objects in the entry_group back to
   *  where the they came from.
   **/
  private void saveAllEntries() 
  {
    SwingWorker worker = new SwingWorker()
    {

      public Object construct()
      {
        final int entry_group_size = entry_group.size();
        for(int entry_index = 0; entry_index < entry_group_size;
            ++entry_index) 
          saveEntry(entry_group.elementAt(entry_index), true, false, true,
                    DocumentEntryFactory.ANY_FORMAT, false);
        return null;
      }
    };
    worker.start();
  }

  /**
   *  Check the given Entry for invalid EMBL features(such as CDS features
   *  without a stop codon) then display a FeatureListFrame list the problem
   *  features.
   *  @return true if and only if the save should be aborted.
   **/
  private boolean displaySaveWarnings(final Entry entry) 
  {
    final FeatureVector invalid_starts = entry.checkFeatureStartCodons();
    final FeatureVector invalid_stops = entry.checkFeatureStopCodons();
    final FeatureVector invalid_keys = entry.checkForNonEMBLKeys();
    final FeatureVector duplicate_features = entry.checkForEMBLDuplicates();
    final FeatureVector overlapping_cds_features =
      entry.checkForOverlappingCDSs();
    final FeatureVector missing_qualifier_features =
      entry.checkForMissingQualifiers();

    // this predicate will filter out those features that aren't in the
    // entry we are trying to save
    final FeaturePredicate predicate = new FeaturePredicate()
    {
      public boolean testPredicate(final Feature feature)
      {
        if(feature.getEntry() == entry) 
          return true;
        else 
          return false;
      }
    };

    String entry_name = entry.getName();

    if(entry_name == null) 
      entry_name = "no name";

    final FilteredEntryGroup filtered_entry_group =
      new FilteredEntryGroup(getEntryGroup(), predicate,
                             "features from " + entry_name);

    if(invalid_starts.size() + invalid_stops.size() +
       invalid_keys.size() + duplicate_features.size() +
       overlapping_cds_features.size() > 0) 
    {

      final YesNoDialog yes_no_dialog =
        new YesNoDialog(this,
                         "warning: some features may have problems. " +
                         "continue with save?");

      if(!yes_no_dialog.getResult()) 
      {
        getSelection().clear();

        if(invalid_starts.size() > 0) 
        {
          getSelection().add(invalid_starts);

          ViewMenu.showBadStartCodons(this,
                                      getSelection(),
                                      filtered_entry_group,
                                      getGotoEventSource(),
                                      base_plot_group);
        }

        if(invalid_stops.size() > 0) 
        {
          getSelection().add(invalid_stops);

          ViewMenu.showBadStopCodons(this, getSelection(),
                                     filtered_entry_group,
                                     getGotoEventSource(),
                                     base_plot_group);
        }

        if(invalid_keys.size() > 0) 
        {
          getSelection().add(invalid_keys);

          ViewMenu.showNonEMBLKeys(this, getSelection(),
                                   filtered_entry_group,
                                   getGotoEventSource(),
                                   base_plot_group);
        }

        if(duplicate_features.size() > 0) 
        {
          getSelection().add(duplicate_features);

          ViewMenu.showDuplicatedFeatures(this, getSelection(),
                                          filtered_entry_group,
                                          getGotoEventSource(),
                                          base_plot_group);
        }

        if(overlapping_cds_features.size() > 0)
        {
          getSelection().add(overlapping_cds_features);

          ViewMenu.showOverlappingCDSs(this, getSelection(),
                                       filtered_entry_group,
                                       getGotoEventSource(),
                                       base_plot_group);
        }

        if(missing_qualifier_features.size() > 0) 
        {
          getSelection().add(missing_qualifier_features);

          ViewMenu.showMissingQualifierFeatures(this, getSelection(),
                                                filtered_entry_group,
                                                getGotoEventSource(),
                                                base_plot_group);
        }
        return true;
      }
    }
    return false;
  }

  /**
   *  Make and add the menus for this component.
   * @param lookseqPanel 
   **/
  private void makeMenus(final JSplitPane splitPane,
                         final JScrollPane jspLookSeq, 
                         final LookSeqPanel lookseqPanel) 
  {
    setJMenuBar(menu_bar);
    makeFileMenu();
    menu_bar.add(file_menu);

    // don't add the menu if this is an applet and we have just one entry
    if(Options.readWritePossible() || getEntryGroup().size() > 1) 
    {
      JMenu entry_group_menu = new EntryGroupMenu(this, getEntryGroup());
      entry_group_menu.setMnemonic(KeyEvent.VK_N);
      menu_bar.add(entry_group_menu);
    }

    SelectMenu select_menu = new SelectMenu(this, getSelection(),
                     getGotoEventSource(), getEntryGroup(),
                     base_plot_group);
    select_menu.setMnemonic(KeyEvent.VK_S);
    menu_bar.add(select_menu);

    ViewMenu view_menu = new ViewMenu(this, getSelection(),
                             getGotoEventSource(), getEntryGroup(),
                             base_plot_group);
    view_menu.setMnemonic(KeyEvent.VK_V);
    menu_bar.add(view_menu);

    JMenu goto_menu = new GotoMenu(this, getSelection(),
                             getGotoEventSource(), getEntryGroup());
    goto_menu.setMnemonic(KeyEvent.VK_O);
    menu_bar.add(goto_menu);

    if(Options.readWritePossible()) 
    {
      EditMenu edit_menu = new EditMenu(this, getSelection(),
                               getGotoEventSource(), getEntryGroup(),
                               base_plot_group, feature_display);
      edit_menu.setMnemonic(KeyEvent.VK_E);
      menu_bar.add(edit_menu);

      AddMenu add_menu = new AddMenu(this, getSelection(), getEntryGroup(),
                             getGotoEventSource(), base_plot_group);
      add_menu.setMnemonic(KeyEvent.VK_C);
      menu_bar.add(add_menu);

      /*JMenu write_menu = new WriteMenu(this, getSelection(), getEntryGroup());
      write_menu.setMnemonic(KeyEvent.VK_W);
      menu_bar.add(write_menu);*/

      JMenu run_menu = new RunMenu(this, getSelection());
      run_menu.setMnemonic(KeyEvent.VK_R);
      menu_bar.add(run_menu);
    }

    JMenu graph_menu = new GraphMenu(this, getEntryGroup(),
                                     base_plot_group,
                                     feature_display, splitPane);
    graph_menu.setMnemonic(KeyEvent.VK_G);
    menu_bar.add(graph_menu);

    final JMenu display_menu = new JMenu("Display");
    display_menu.setMnemonic(KeyEvent.VK_D);
    final JCheckBoxMenuItem show_entry_buttons_item =
                           new JCheckBoxMenuItem("Show Entry Buttons");
    show_entry_buttons_item.setState(true);
    show_entry_buttons_item.addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent event) 
      {
        group_display.setVisible(show_entry_buttons_item.getState());
        // XXX change to revalidate().
        validate();
      }
    });
    display_menu.add(show_entry_buttons_item);

    final JCheckBoxMenuItem show_one_line =
      new JCheckBoxMenuItem("Show One Line Per Entry View", false);
    final JCheckBoxMenuItem show_feature_stack =
        new JCheckBoxMenuItem("Show Feature Stack View", false);
    
    show_one_line.setState(one_line_per_entry_display.isVisible());
    show_one_line.addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent event) 
      {
        if(show_one_line.getState())
          show_feature_stack.setState(false);
        one_line_per_entry_display.setFeatureStackViewFlag(false);
        one_line_per_entry_display.setVisible(show_one_line.getState());
        validate();
      }
    });
    display_menu.add(show_one_line);
    
    show_feature_stack.setState(one_line_per_entry_display.isVisible());
    show_feature_stack.addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent event) 
      {
        if(show_feature_stack.getState())
          show_one_line.setState(false);
        one_line_per_entry_display.setFeatureStackViewFlag(true);
        one_line_per_entry_display.setVisible(show_feature_stack.getState());
        validate();
        one_line_per_entry_display.updateOneLinePerFeatureFlag();
      }
    });
    display_menu.add(show_feature_stack);

    final JCheckBoxMenuItem show_base_display_item =
                                new JCheckBoxMenuItem("Show Base View");

    show_base_display_item.setState(base_display.isVisible());
    show_base_display_item.addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent event) 
      {
        base_display.setVisible(show_base_display_item.getState());
        // XXX change to revalidate().
        validate();
      }
    });
    display_menu.add(show_base_display_item);

    final JCheckBoxMenuItem show_feature_list_item =
                        new JCheckBoxMenuItem("Show Feature List");
    show_feature_list_item.setState(feature_list.isVisible());
    show_feature_list_item.addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent event) 
      {
        feature_list.setVisible(show_feature_list_item.getState());
        // XXX change to revalidate().
        validate();
      }
    });
    display_menu.add(show_feature_list_item);
    
    if (lookseqPanel != null)
    {
      final JCheckBoxMenuItem show_lookseq_item = new JCheckBoxMenuItem(
          "Show lookseq");
      show_lookseq_item.setState(jspLookSeq.isVisible());
      show_lookseq_item.addItemListener(new ItemListener()
      {
        public void itemStateChanged(ItemEvent event)
        {
          if (show_lookseq_item.getState())
          {
            if (lookseqPanel.getQueryStr() == null)
            {
              String urlStr = Options.getOptions().getProperty("lookseq");
              // String urlStr =
              // "http://www.sanger.ac.uk/cgi-bin/teams/team112/lookseq/get_data.pl?";
              String queryStr = "from="
                  + feature_display.getForwardBaseAtLeftEdge()
                  + "&to=8000"
                  + "&chr=MAL1&output=image&width=1024&lane=sample_2a&view=indel&display=|perfect|snps|inversions|pairlinks|potsnps|uniqueness|&debug=0";

              if(Options.getOptions().getProperty("lookseq_chr") != null)
                queryStr = queryStr.replaceFirst(
                    "chr=[^&]+", "chr="+Options.getOptions().getProperty("lookseq_chr").trim());
              
              if(Options.getOptions().getProperty("lookseq_lane") != null)
                queryStr = queryStr.replaceFirst(
                    "lane=[^&]+", "lane="+Options.getOptions().getProperty("lookseq_lane").trim());
              
              lookseqPanel.setUrl(urlStr, queryStr);
            }

            lookseqPanel.setFeatureDisplay(feature_display);
            lookseqPanel.showOptions();
            feature_display.addDisplayAdjustmentListener(lookseqPanel);
          }
          else
            feature_display.removeDisplayAdjustmentListener(lookseqPanel);

          jspLookSeq.setVisible(show_lookseq_item.getState());

          if (show_lookseq_item.getState())
            lookseqPanel.setDisplay(feature_display.getForwardBaseAtLeftEdge(),
                feature_display.getLastVisibleForwardBase(), null);

          validate();
        }
      });
      display_menu.add(show_lookseq_item);
    }

    display_menu.addSeparator();
    final JMenuItem show_Jam_item = new JMenuItem("BAM Alignment");
    show_Jam_item.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)
      {
        if(bamView == null)
          return;
        
        if (!bamPanel.isVisible())
        {
          bamPanel.setVisible(true);
          bamView.setDisplay(feature_display.getFirstVisibleForwardBase(),
                             feature_display.getLastVisibleForwardBase(), null);
          bamView.revalidate();
          one_line_per_entry_display.addDisplayAdjustmentListener(bamView);
          feature_display.addDisplayAdjustmentListener(bamView);
          feature_display.getSelection().addSelectionChangeListener(bamView);
        }
        else
        {
          one_line_per_entry_display.removeDisplayAdjustmentListener(bamView);
          feature_display.removeDisplayAdjustmentListener(bamView);
          feature_display.getSelection().removeSelectionChangeListener(bamView);
          bamPanel.setVisible(false);
        }
        setNGDivider();
      }
    });
    display_menu.add(show_Jam_item);
    
    
    final JMenuItem show_Vcf_item = new JMenuItem("VCF");
    show_Vcf_item.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)
      {
        if(vcfView == null)
          return;
        
        if (!vcfView.isVisible())
        {
          vcfView.setVisible(true);
          vcfView.revalidate();
          feature_display.addDisplayAdjustmentListener(vcfView);
          feature_display.getSelection().addSelectionChangeListener(vcfView);
        }
        else
        {
          feature_display.removeDisplayAdjustmentListener(vcfView);
          feature_display.getSelection().removeSelectionChangeListener(vcfView);
          vcfView.setVisible(false);
        }
        setNGDivider();
      }
    });
    display_menu.add(show_Vcf_item);

    menu_bar.add(display_menu);
  }


  /**
   *  Make a new File menureplacing the current one (if any).
   **/
  private void makeFileMenu() 
  {
    file_menu.removeAll();
    file_menu.setMnemonic(KeyEvent.VK_F);

    if(Options.readWritePossible()) 
    {
      boolean db = false;
      final EntryVector entries = getEntryGroup().getActiveEntries();
      for(int i=0; i<entries.size(); i++)
      {
        Entry entry = entries.elementAt(i);
        if(entry.getEMBLEntry() instanceof DatabaseDocumentEntry)
          db = true;
      }

      if(db && !getEntryGroup().getDefaultEntry().getEMBLEntry().isReadOnly())
      {
        final JMenuItem commit = new JMenuItem("Commit to Database");
        commit.addActionListener(new ActionListener()
        {
          public void actionPerformed(ActionEvent event)
          {
            commitToDatabase(entry_group,
                getDatabaseDocumentEntry().getChadoTransactionManager(), EntryEdit.this, 
                selection, goto_event_source, base_plot_group);
          }
        });
        file_menu.add(commit);
        file_menu.addSeparator();
      }

      JMenuItem popFileManager = new JMenuItem("Show File Manager ...");
      popFileManager.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent event)
        {
          if(ArtemisMain.filemanager == null)
            ArtemisMain.filemanager = new FileManager(EntryEdit.this);
          else
            ArtemisMain.filemanager.setVisible(true);
        }
      });
      file_menu.add(popFileManager);

      // only the standalone version can save or read
      EntrySource filesystem_entry_source = null;
      final int entry_sources_size = entry_sources.size();
     
      for(int source_index = 0; source_index < entry_sources_size;
          ++source_index) 
      {
        final EntrySource this_source =
          entry_sources.elementAt(source_index);

        if(this_source.isFullEntrySource()) 
          continue;

        if(this_source.getSourceName().equals("Filesystem")) 
          filesystem_entry_source = this_source;

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
            readAnEntry(this_source);
          }
        });

        file_menu.add(read_entry);
      }

      JMenu read_features_menu = null;

      if(filesystem_entry_source != null &&
          entry_group != null && entry_group.size() > 0) 
      {
        read_features_menu = new JMenu("Read Entry Into");
        file_menu.add(read_features_menu);
      }

      file_menu.addSeparator();
      
      JMenuItem read_bam_file = new JMenuItem("Read BAM / CRAM / VCF ...");
      read_bam_file.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          FileSelectionDialog fileChooser = new FileSelectionDialog(
              null, false, "BAM / CRAM / VCF View", "BAM / CRAM / VCF");
          List<String> listBams = fileChooser.getFiles(BamView.BAM_SUFFIX);
          List<String> vcfFiles = fileChooser.getFiles(VCFview.VCFFILE_SUFFIX);
          loadBamAndVcf(listBams, vcfFiles);
        }
      });
      file_menu.add(read_bam_file);    
     
      file_menu.addSeparator();

      final JMenuItem save_default =
        new JMenuItem("Save Default Entry");
      save_default.setAccelerator(SAVE_DEFAULT_KEY);
      save_default.addActionListener(new ActionListener() 
      {
        public void actionPerformed(ActionEvent event) 
        {
          saveDefaultEntry();
        }
      });

      file_menu.add(save_default);

      if(entry_group == null || entry_group.size() == 0) 
      {
        // do nothing
      } 
      else 
      {
        final JMenu save_entry_menu = new JMenu("Save An Entry");
        final JMenu save_as_menu    = new JMenu("Save An Entry As");
        final JMenu save_as         = new JMenu("New File");
        final JMenu save_as_embl    = new JMenu("EMBL Format");
        final JMenu save_as_genbank = new JMenu("GENBANK Format");
        final JMenu save_as_genbank_only = new JMenu("Sequin Table Format");
        final JMenu save_as_gff     = new JMenu("GFF Format");
        final JMenu save_embl_only  = new JMenu("EMBL Submission Format");
        final int entry_group_size  = getEntryGroup().size();

        for(int i = 0; i < entry_group_size; ++i) 
        {
          final Entry this_entry = getEntryGroup().elementAt(i);
          String entry_name = this_entry.getName();

          if(entry_name == null) 
            entry_name = "no name";

          final ActionListener save_entry_listener =
            new SaveEntryActionListener(this, this_entry);

          final JMenuItem save_entry_item = new JMenuItem(entry_name);

          save_entry_item.addActionListener(save_entry_listener);

          final ActionListener save_as_listener =
            new SaveEntryAsActionListener(this, this_entry);

          final JMenuItem save_as_item = new JMenuItem(entry_name);

          save_as_item.addActionListener(save_as_listener);

          final ActionListener save_as_embl_listener =
            new SaveEntryAsEMBLActionListener(this, this_entry);

          final JMenuItem save_as_embl_item = new JMenuItem(entry_name);

          save_as_embl_item.addActionListener(save_as_embl_listener);

          final ActionListener save_as_genbank_listener =
            new SaveEntryAsGenbankActionListener(this, this_entry);

          final JMenuItem save_as_genbank_item = new JMenuItem(entry_name);
          save_as_genbank_item.addActionListener(save_as_genbank_listener);

          final JMenuItem save_as_genbank_only_item = new JMenuItem(entry_name);
          save_as_genbank_only_item.addActionListener(new ActionListener()
          {
            public void actionPerformed(ActionEvent e)
            {
              GenbankTblOutputStream.writeEntryAsTbl(this_entry, EntryEdit.this); 
            }
          });
          
          final ActionListener save_as_gff_listener =
            new SaveEntryAsGFFActionListener(this, this_entry);

          final JMenuItem save_as_gff_item = new JMenuItem(entry_name);

          save_as_gff_item.addActionListener(save_as_gff_listener);

          final ActionListener save_embl_only_listener =
            new SaveEntryAsSubmissionActionListener(this, this_entry);

          final JMenuItem save_embl_only_item = new JMenuItem(entry_name);

          save_embl_only_item.addActionListener(save_embl_only_listener);

          if(read_features_menu != null)  
          {
            final ActionListener read_into_listener =
              new ReadFeaturesActionListener(this, filesystem_entry_source,
                                              this_entry);

            final JMenuItem read_into_item = new JMenuItem(entry_name);
            read_into_item.addActionListener(read_into_listener);
            read_features_menu.add(read_into_item);
          }

          save_entry_menu.add(save_entry_item);

          save_as.add(save_as_item);
          save_as_embl.add(save_as_embl_item);
          save_as_genbank.add(save_as_genbank_item);
          save_as_genbank_only.add(save_as_genbank_only_item);
          save_as_gff.add(save_as_gff_item);
          save_embl_only.add(save_embl_only_item);
        }

        save_as_menu.add(save_as);
        save_as_menu.add(save_as_embl);
        save_as_menu.add(save_as_genbank);
        save_as_menu.add(save_as_genbank_only);
        save_as_menu.add(save_as_gff);
        save_as_menu.addSeparator();
        save_as_menu.add(save_embl_only);

        file_menu.add(save_entry_menu);
        file_menu.add(save_as_menu);
      }

      final JMenuItem save_all = new JMenuItem("Save All Entries");
      save_all.addActionListener(new ActionListener() 
      {
        public void actionPerformed(ActionEvent event) 
        {
          saveAllEntries();
        }
      });

      file_menu.add(save_all);
      file_menu.addSeparator();
      
      final WriteMenu write_menu = new WriteMenu(this, getSelection(), getEntryGroup());
      write_menu.setMnemonic(KeyEvent.VK_W);
      file_menu.add(write_menu);
      file_menu.addSeparator();
    }

    final JMenuItem clone_entry_edit = new JMenuItem("Clone This Window");
    clone_entry_edit.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        EntryEdit clone = new EntryEdit(getEntryGroup());
        clone.setVisible(true);
      }
    });

    file_menu.add(clone_entry_edit);
    file_menu.addSeparator();
    
    printMenu();

    file_menu.addSeparator();
    final JMenuItem open_dna_viewer = new JMenuItem("Open in DNAPlotter");
    open_dna_viewer.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        Entry entry = getEntryGroup().getSequenceEntry();

        final DNADraw dna = 
          uk.ac.sanger.artemis.circular.Wizard.getDNADrawFromArtemisEntry(null, 
            getEntryGroup(), entry);
        final String version = dna.getVersion();
        final JFrame f = new JFrame();
        if(version == null)
          f.setTitle("DNAPlotter");
        else
          f.setTitle("DNAPlotter :: "+version);
        dna.setCloseAndDispose(true, f);

        JScrollPane jsp = new JScrollPane(dna);
        jsp.getViewport().setBackground(Color.white);
        f.getContentPane().add(jsp);
        f.setJMenuBar(dna.createMenuBar());
        f.pack();
        Utilities.centreFrame(f);
        f.setVisible(true);
      }
    });
    file_menu.add(open_dna_viewer);
    
    file_menu.addSeparator();
    final JMenuItem prefs = new JMenuItem("Preferences");
    prefs.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)
      {
        final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();

        int nmenus = EntryEdit.super.getJMenuBar().getMenuCount();
        final JTabbedPane shortcut_pane = new JTabbedPane();
        
        Vector<String> sc = new Vector<String>();
        for(int i=0; i<nmenus; i++)
        {
          JMenu menu = EntryEdit.super.getJMenuBar().getMenu(i);
          if(menu instanceof SelectionMenu &&
            ((SelectionMenu)menu).isEditableShortCutMenu() )
            sc.addAll( ((SelectionMenu)menu).getUsedShortCutKeys() );
        }
        
        for(int i=0; i<nmenus; i++)
        {
          JMenu menu = EntryEdit.super.getJMenuBar().getMenu(i);
          if(menu instanceof SelectionMenu &&
            ((SelectionMenu)menu).isEditableShortCutMenu() )
            shortcut_pane.add(menu.getText()+" Menu", ((SelectionMenu)menu).getShortCuts(sc));
        }
        
        shortcut_pane.setPreferredSize(new Dimension(
                           shortcut_pane.getPreferredSize().width+50,
                           screen.height/2));

        final JTabbedPane tabPane = new JTabbedPane();
        
        //
        // Display names
        //
        final Object display_names[] =
          Options.getOptions().getDisplayQualifierNames().toArray();

        String[] description = 
        { "Names of qualifiers to search when attempting to find",
    	  "the display name. These qualifier names are searched in order." };
        ListSelectionPanel displayListSelectionPanel =
              new ListSelectionPanel(entry_group, display_names,
            		  description, true);
     
        //
        //
        //
        tabPane.add("Feature Display Labels", displayListSelectionPanel);

        //
        // Display names
        //
        final Object systematic_names[] =
          Options.getOptions().getSystematicQualifierNames().toArray();

        String[] description2 = 
            { "Names of qualifiers to search when attempting to find",
        	  "the systematic name of a gene."	} ;
        ListSelectionPanel systematicListSelectionPanel =
              new ListSelectionPanel(entry_group, systematic_names,
            		 description2, true);
     
        //
        //
        //
        tabPane.add("Systematic Name Labels", systematicListSelectionPanel);

        //
        // Object Editor - accessing entries via SRS or mfetch
        //
        JTextField cacheSize = null;
        if(System.getProperty("j2ssh") != null)
        {
          boolean remoteMfetch = false;
          if(FileList.isConnected() && !FastaTextPane.isForceUrl())
          {
            FileList fileList = new FileList();
            FileAttributes attr = fileList.stat("/nfs/disk100/pubseq/bin/mfetch");
            if(attr != null)
              remoteMfetch = attr.isFile();
          }
          
          Box yBox = Box.createVerticalBox();
          final String srsUrl = (String)Options.getOptions().getOptionValues("srs_url").elementAt(0);
          final JCheckBox useMfetch = 
            new JCheckBox("Use Sanger Server (Sanger Users Only)", remoteMfetch);
          
          useMfetch.addItemListener(new ItemListener()
          {
            public void itemStateChanged(ItemEvent e)
            {
              if(useMfetch.isSelected())
              {
                FileList fileList = new FileList();
                FileAttributes attr = fileList.stat("/nfs/disk100/pubseq/bin/mfetch");
                if(attr != null)
                {
                  FastaTextPane.setForceUrl(false);
                  FastaTextPane.setRemoteMfetch(useMfetch.isSelected());
                }
              }
              else
              {
                FastaTextPane.setForceUrl(true);
                FastaTextPane.setRemoteMfetch(useMfetch.isSelected());
              }
            }
          });
          final JCheckBox useSrs = new JCheckBox("Use SRS "+srsUrl, !remoteMfetch);
          
          ButtonGroup group = new ButtonGroup();
          group.add(useMfetch);
          group.add(useSrs);
          yBox.add(new JLabel("Select the method for retrieving entry information "+
                              "in the object editor:"));
          yBox.add(useMfetch);
          yBox.add(useSrs);
          yBox.add(Box.createVerticalStrut(5));
          
          // cache
          cacheSize = new JTextField(Integer.toString(BigPane.CACHE_SIZE), 6);
          cacheSize.setMaximumSize(new Dimension(cacheSize.getMaximumSize().width,
                                                 cacheSize.getPreferredSize().height));
          yBox.add(Box.createVerticalStrut(5));
          yBox.add(new JLabel("Size of cache for retrieved entried (max "+
              BigPane.MAX_CACHE_SIZE+"):"));
          yBox.add(cacheSize);
          
          JButton buttClear = new JButton("Empty Cache");
          buttClear.addActionListener(new ActionListener()
          {
            public void actionPerformed(ActionEvent e)
            {
              Arrays.fill(FastaTextPane.cacheHits, null);
            }
          });
          yBox.add(buttClear);
          
          yBox.add(Box.createVerticalGlue());
          tabPane.add("Object Editor", yBox);
        }
        
        tabPane.add("Short Cuts", shortcut_pane);
        
        //
        // Contig ordering options
        //
        final Vector<String> contigKeys = FeatureDisplay.getContigKeys();
        final Vector<String> allPossibleContigKeys = FeatureDisplay.getAllPossibleContigKeys();
        final Vector<String> nonContigKeys = new Vector<String>();
        
        for(int i=0; i<allPossibleContigKeys.size(); i++)
        {
          String keyStr = allPossibleContigKeys.get(i);
          if( !contigKeys.contains(keyStr) )
            nonContigKeys.add(keyStr);
        }

        final DefaultListModel showListModel = new DefaultListModel();
        for(int i=0; i<contigKeys.size(); i++)
          showListModel.addElement(contigKeys.get(i));
        final JList displayList = new JList(showListModel);
        
        
        final DefaultListModel hideListModel = new DefaultListModel();
        for(int i=0; i<nonContigKeys.size(); i++)
          hideListModel.addElement(nonContigKeys.get(i));
        final JList hideList = new JList(hideListModel);
        
        
        final JButton hide_butt = new JButton(">>");
        hide_butt.addActionListener(new ActionListener()
        {
          public void actionPerformed(ActionEvent e)
          {
            while(!displayList.isSelectionEmpty())
            {
              final String hideKey = (String)displayList.getSelectedValue();
              
              nonContigKeys.add(hideKey);
              if(contigKeys.contains(hideKey))
                contigKeys.remove(hideKey);

              hideListModel.add(nonContigKeys.indexOf(hideKey), hideKey);
              showListModel.removeElement(hideKey);
            }
          }
        });

        
        final JPanel panel = new JPanel(new BorderLayout());
        final JPanel contigPanel = new JPanel(new GridBagLayout());
        GridBagConstraints c = new GridBagConstraints();
        
        c.anchor = GridBagConstraints.NORTHWEST;
      
        c.ipadx = 50;
        c.gridx = 0;
        c.gridy = 0;
        contigPanel.add(new JLabel("Contig Ordering Features:"),c);
        c.gridy = 1;
        JScrollPane jspContig = new JScrollPane(displayList);
        Dimension d = new Dimension(300, jspContig.getPreferredSize().height);
        jspContig.setPreferredSize(d);
        contigPanel.add(jspContig,c);
        c.gridy = 2;
        contigPanel.add(hide_butt,c);
        
        final JButton show_butt = new JButton("<<");
        show_butt.addActionListener(new ActionListener()
        {
          public void actionPerformed(ActionEvent e)
          {
            while(!hideList.isSelectionEmpty())
            {
              final String showKey = (String)hideList.getSelectedValue();
              
              if(nonContigKeys.contains(showKey))
                nonContigKeys.remove(showKey);
              contigKeys.add(showKey);
              
              showListModel.add(contigKeys.indexOf(showKey), showKey);
              hideListModel.removeElement(showKey);
            }
          }
        });

        c.gridx = 1;
        c.gridy = 0;
        contigPanel.add(new JLabel("Non-Contig Ordering Features:"),c);
        c.gridy = 1;
        JScrollPane jspNonContig = new JScrollPane(hideList);
        jspNonContig.setPreferredSize(d);
        contigPanel.add(jspNonContig,c);
        c.gridy = 2;
        contigPanel.add(show_butt,c);
 
        panel.add(contigPanel, BorderLayout.CENTER);
        tabPane.add("Contig Tool Options", panel);

        
        /*String urlString = (String)Options.getOptions().getOptionValues("srs_url").elementAt(0);
        Box yBox = Box.createVerticalBox();
        JTextField srsField = new JTextField(urlString);
        yBox.add(srsField);
        yBox.add(Box.createVerticalGlue());
        
        srsField.setSelectionStart(0);
        srsField.setSelectionEnd(urlString.length());
        
        tabPane.add("SRS Site", yBox);*/
        
        JOptionPane.showMessageDialog(null,
                               tabPane,
                               "Preferences",
                               JOptionPane.PLAIN_MESSAGE);
        tabPane.removeAll();
        
        Options.getOptions().setDisplayNameQualifiers(
        		displayListSelectionPanel.getResultString());
        
        if(displayListSelectionPanel.isSaveOption())
          Splash.save_display_name = true;
        
        Options.getOptions().setSystematicQualifierNames(
        		systematicListSelectionPanel.getResultString());
        
        if(systematicListSelectionPanel.isSaveOption())
          Splash.save_systematic_names = true;
        
        //if(!srsField.getText().equals(urlString))
        //  Options.getOptions().setProperty("srs_url", srsField.getText().trim());
        
        if(cacheSize != null && 
           Integer.parseInt(cacheSize.getText()) != BigPane.CACHE_SIZE)
        {
          int cacheSizeValue = Integer.parseInt(cacheSize.getText());
          
          if(cacheSizeValue <= 0)
            return;
          if(cacheSizeValue > BigPane.MAX_CACHE_SIZE)
            BigPane.CACHE_SIZE = BigPane.MAX_CACHE_SIZE;
          else
            BigPane.CACHE_SIZE = Integer.parseInt(cacheSize.getText());
          HitInfo[] cacheHits = FastaTextPane.cacheHits;
        
          FastaTextPane.cacheHits = new HitInfo[BigPane.CACHE_SIZE];
          for(int i = 0; i < cacheHits.length; i++)
          {
            if(i >= FastaTextPane.cacheHits.length)
              break;
            FastaTextPane.cacheHits[i] = cacheHits[i];
          }
        }
     }
    });
    file_menu.add(prefs);

    file_menu.addSeparator();

    final JMenuItem close = new JMenuItem("Close");
    close.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent event) 
      {
        closeEntryEdit();
      }
    });

    file_menu.add(close);
  }
  
  private void loadBamAndVcf(List<String> listBams, List<String> vcfFiles) 
  {
    logger4j.debug("No. BAM FILES="+listBams.size());
    
    if(listBams.size() > 0)
    {
      bamPanel.removeAll();
      try
      {
        bamView = new BamView(listBams, null, feature_display.getMaxVisibleBases(), this,
            feature_display, getEntryGroup().getBases(), bamPanel, null);
        
        
        if(entry_group.getSequenceEntry().getEMBLEntry().getSequence() instanceof IndexFastaStream)
        {
          // add reference sequence selection listeners
          getEntryGroupDisplay().getIndexFastaCombo().addIndexReferenceListener(bamView.getCombo());
          bamView.getCombo().addIndexReferenceListener(getEntryGroupDisplay().getIndexFastaCombo());
        }
        
      }
      catch (Exception ex)
      {
        logger4j.warn("EntryEdit.loadBamAndVcf() "+ex.getMessage());
        
        if(ex.getMessage() != null)
          JOptionPane.showMessageDialog(null, ex.getMessage(), "Error",
            JOptionPane.ERROR_MESSAGE);
        else
          ex.printStackTrace();
        return;
      }

      bamView.getJspView().getVerticalScrollBar().setValue(
          bamView.getJspView().getVerticalScrollBar().getMaximum());
      bamView.setDisplay(feature_display.getFirstVisibleForwardBase(),
          feature_display.getLastVisibleForwardBase(), null);
      
      one_line_per_entry_display.addDisplayAdjustmentListener(bamView);
      feature_display.addDisplayAdjustmentListener(bamView);
      feature_display.getSelection().addSelectionChangeListener(bamView);

      setNGDivider();
    }
    
    logger4j.debug("No. VCF FILES="+vcfFiles.size());
    if (vcfFiles.size() > 0)
    {
      vcfPanel.removeAll();
      vcfView = new VCFview(null, vcfPanel, vcfFiles,
          feature_display.getMaxVisibleBases(), 1, null, null,
          this, feature_display);

      feature_display.addDisplayAdjustmentListener(vcfView);
      feature_display.getSelection().addSelectionChangeListener(vcfView);

      if(entry_group.getSequenceEntry().getEMBLEntry().getSequence() instanceof IndexFastaStream)
      {
        // add reference sequence selection listeners
        getEntryGroupDisplay().getIndexFastaCombo().addIndexReferenceListener(vcfView.getCombo());
        vcfView.getCombo().addIndexReferenceListener(getEntryGroupDisplay().getIndexFastaCombo());
      }
  
      setNGDivider();
    }
  }

  /**
   * Handle the split panes divider positions for BamView and VcfView.
   */
  public void setNGDivider() 
  {
    if( (bamPanel.getComponents().length > 0 && bamPanel.isVisible()) &&
        (vcfPanel.getComponents().length > 0 && vcfView.isVisible()))
    {
      ngSplitPane.setVisible(true);
      lowerSplitPane.setDividerSize(3);
      lowerSplitPane.setDividerLocation(0.4d);
      
      ngSplitPane.setResizeWeight(0.5);
      ngSplitPane.setDividerSize(3);
      ngSplitPane.setDividerLocation(0.5);
      
      logger4j.debug("BAM & VCF visible");
    }
    else if(vcfPanel.getComponents().length > 0 && vcfView.isVisible())
    {
      ngSplitPane.setVisible(true);
      lowerSplitPane.setDividerSize(3);
      lowerSplitPane.setDividerLocation(0.25d);
      
      ngSplitPane.setResizeWeight(0);
      ngSplitPane.setDividerSize(0);
      ngSplitPane.setDividerLocation(0.);
      logger4j.debug("VCF visible");
    }
    else if(bamPanel.getComponents().length > 0 && bamPanel.isVisible()) 
    {
      ngSplitPane.setVisible(true);
      lowerSplitPane.setDividerSize(3);
      lowerSplitPane.setDividerLocation(0.3d);
      
      ngSplitPane.setResizeWeight(1);
      ngSplitPane.setDividerSize(1);
      ngSplitPane.setDividerLocation(1.d);
      logger4j.debug("BAM visible");
    }
    else
    {
      lowerSplitPane.setResizeWeight(0);
      lowerSplitPane.setDividerSize(0);
      lowerSplitPane.setDividerLocation(0);
      logger4j.debug("BAM & VCF not visible");
    }
  }
  
  private void printMenu()
  {
    JMenuItem printImage = new JMenuItem("Save As Image Files (png/svg)...");
    printImage.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        PrintArtemis part = new PrintArtemis(EntryEdit.this);
        part.print();
      }
    });
    file_menu.add(printImage);
    
    JMenuItem printPS = new JMenuItem("Print...");
    printPS.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        PrintArtemis part = new PrintArtemis(EntryEdit.this);
        part.validate();
        part.doPrintActions();
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
        PrintArtemis part = new PrintArtemis(EntryEdit.this);
        part.printPreview();
      }
    });

  }
  
  private static boolean canCommitChecks(final EntryGroup entry_group,
                                         final ChadoTransactionManager ctm,
                                         final JFrame frame,
                                         final Selection selection,
                                         final GotoEventSource getGotoEventSource,
                                         final BasePlotGroup  base_plot_group)
  {
    if(entry_group.getDefaultEntry() == null)
    {
      JOptionPane.showMessageDialog(frame, 
          "No default entry selected.\n"+
          "Select the default from the entry bar at the top.");
      return false; 
    }
    else if(!ctm.hasTransactions())
    {
      JOptionPane.showMessageDialog(frame, 
          "No changes to commit to the database.");
      return false;
    }
    
    if(System.getProperty("nocommit") == null ||
       System.getProperty("nocommit").equals("false"))
    {
      int select = JOptionPane.showConfirmDialog(frame, 
          "Commit "+ctm.getTransactionCount()+" change(s) to the database?", 
          "Commit", JOptionPane.OK_CANCEL_OPTION);
      if(select == JOptionPane.CANCEL_OPTION)
        return false;
    }
    
    if(!isUniqueID(entry_group, ctm, selection, 
          getGotoEventSource, base_plot_group))
      return false;
    return true;
  }

  public static boolean commitToDatabase(final EntryGroup entry_group,
      final ChadoTransactionManager ctm,
      final JFrame frame,
      final Selection selection,
      final GotoEventSource getGotoEventSource,
      final BasePlotGroup  base_plot_group)
  {
    return commitToDatabase(entry_group, ctm, frame, selection,
                   getGotoEventSource, base_plot_group, false);
  }
  
  /**
   * Commit back to the database.
   * @param entry_group
   * @param ctm
   * @param frame
   * @param selection
   * @param getGotoEventSource
   * @param base_plot_group
   * @param force force the commit - i.e. commit everything that doesn't
   * throw an exception
   * @return
   */
  private static boolean commitToDatabase(final EntryGroup entry_group,
                                          final ChadoTransactionManager ctm,
                                          final JFrame frame,
                                          final Selection selection,
                                          final GotoEventSource getGotoEventSource,
                                          final BasePlotGroup  base_plot_group,
                                          final boolean force)
  {
    try
    {
      if(!canCommitChecks(entry_group, ctm, frame, selection,
                        getGotoEventSource, base_plot_group))
        return false;
      
      final Document dbDoc =
          ((DocumentEntry)entry_group.getDefaultEntry().getEMBLEntry()).getDocument();

      if(dbDoc instanceof DatabaseDocument)
      {
        frame.setCursor(new Cursor(Cursor.WAIT_CURSOR));
        ChadoTransactionManager.commit((DatabaseDocument)dbDoc, force, ctm);
        
        final String nocommit = System.getProperty("nocommit");
        if(!force && ctm.hasTransactions() &&
           (nocommit == null || nocommit.equals("false")))
        {
          int forceCommit = JOptionPane.showConfirmDialog(frame, 
              "Force commit (saving any changes that can be\n"+ 
              "committed to the database and ignoring those that fail)?", 
              "Force Database Commit", JOptionPane.OK_CANCEL_OPTION);
          
          if(forceCommit == JOptionPane.OK_OPTION)
            ChadoTransactionManager.commit((DatabaseDocument)dbDoc, true, ctm);
        }
        frame.setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
      }
      else
        new MessageDialog(frame,
                 "No database associated with the default entry");
    }
    catch(NullPointerException npe)
    {
      frame.setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
      new MessageDialog(frame,
                 "Problem writing to the database - check log window");
      npe.printStackTrace();
      logger4j.error(npe.getMessage()); 
    }
    
    if(ctm.hasTransactions())
      return false;
    return true;
  }
  

  /**
   * Test to ensure ID (chado uniquename) is unique.
   * @param entry_group
   * @return
   */
  private static boolean isUniqueID(final EntryGroup entry_group,
                                    final ChadoTransactionManager ctm, 
                                    final Selection selection,
                                    final GotoEventSource getGotoEventSource,
                                    final BasePlotGroup  base_plot_group)
  {
    // only need to check if the Feature table has been changed
    List<String> changed_features = ctm.getFeatureInsertUpdate();
    if(changed_features == null)
      return true;
    
    // filter out non-db entries
    FilteredEntryGroup db_entry_group = new FilteredEntryGroup(entry_group, 
                                            new FeatureDatabasePredicate(), 
                                           "Database Entries");
    FeatureVector features = db_entry_group.getAllFeatures();
    FeatureVector duplicateIDs = new FeatureVector();
    for(int i=0; i<features.size()-1; i++)
    {
      Feature ifeature = features.elementAt(i);
      String id;
      try
      {
        id = (String)ifeature.getQualifierByName("ID").getValues().get(0);
        if(!changed_features.contains(id))
          continue;
      }
      catch(InvalidRelationException e)
      {
        continue;
      }
      catch(NullPointerException npe)
      {
        npe.printStackTrace();
        logger4j.error(ifeature.getLocation().toStringShort());
        continue;
      }
      
      FeaturePredicate predicate =
        new FeatureKeyQualifierPredicate(null, "ID", id, 
                                         false, true);
      
      for(int j=i+1; j<features.size(); j++)
      {
        Feature jfeature = features.elementAt(j);
        if(predicate.testPredicate(jfeature))
        {
          duplicateIDs.add(ifeature);
          duplicateIDs.add(jfeature);
        }
      }

    }
    
    if(duplicateIDs.size() > 0)
    {
      String filtern_name = "Features with non-unique feature ID's";
      selection.set(duplicateIDs);
      final FilteredEntryGroup filtered_entry_group =
        new FilteredEntryGroup(entry_group, duplicateIDs, 
                               filtern_name);
      
      final FeatureListFrame feature_list_frame =
        new FeatureListFrame(filtern_name, selection,
            getGotoEventSource, filtered_entry_group,
                             base_plot_group);

      feature_list_frame.setVisible (true);
      return false;
    }
    return true;
  }

  /**
   *  Read an entry
   **/
  private void readAnEntry(final EntrySource this_source)
  {
    try
    {
      final Entry new_entry = this_source.getEntry(entry_group.getBases(),
                                                   true);
      if(new_entry != null)
        getEntryGroup().add(new_entry);
    }
    catch(final OutOfRangeException e)
    {
      new MessageDialog(EntryEdit.this,
                     "read failed: one of the features " +
                     "in the entry has an out of " +
                     "range location: " +
                     e.getMessage());
    }
    catch(final IOException e)
    {
      new MessageDialog(EntryEdit.this,
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
  
  class CommitButton extends JButton
        implements FeatureChangeListener, SequenceChangeListener, EntryChangeListener
  {
    private static final long serialVersionUID = 1L;
    private Color DEFAULT_FOREGROUND;
    private CommitFrame commitFrame;
    private final ChadoTransactionManager ctm;
    
    public CommitButton()
    {
      super("Commit");
      
      ctm = getDatabaseDocumentEntry().getChadoTransactionManager();
      setBackground(EntryGroupDisplay.background_colour);
      addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          if(commitToDatabase(entry_group, ctm, EntryEdit.this, 
                selection, goto_event_source, base_plot_group, false))
          {
            setForeground(DEFAULT_FOREGROUND);
            setFont(getFont().deriveFont(Font.PLAIN));
            
            if(commitFrame != null)
              commitFrame.setList(ctm);
          }
        }    
      });
      
      addMouseListener(new MouseAdapter()
      {
        public void mouseReleased(MouseEvent e)
        {
          if(e.getButton() == MouseEvent.BUTTON3)
          {
            commitFrame =
              new CommitFrame(ctm, entry_group, EntryEdit.this, 
                selection, goto_event_source, base_plot_group);
          }
        }
      });
      
      setMargin(new Insets(0,2,0,2));
      setToolTipText("commit to database");
      
      DEFAULT_FOREGROUND = super.getForeground();
    }
    
    public String getToolTipText()
    {
      if(ctm.hasTransactions())
        return Integer.toString(ctm.getTransactionCount());
      else
        return "no transaction to commit";
    }

    public void featureChanged(FeatureChangeEvent event)
    {
      if(commitFrame != null)
        commitFrame.setList(ctm);
      if(ctm.hasTransactions())
      {
        setFont(getFont().deriveFont(Font.BOLD));
        setForeground(Color.red);
      }
    }

    public void sequenceChanged(SequenceChangeEvent event)
    {
      if(commitFrame != null)
        commitFrame.setList(ctm);
      if(ctm.hasTransactions())
      {
        setForeground(Color.red);
        setFont(getFont().deriveFont(Font.BOLD)); 
      }
    }

    public void entryChanged(EntryChangeEvent event)
    {
      if(commitFrame != null)
        commitFrame.setList(ctm);
      if(ctm.hasTransactions())
      {
        setForeground(Color.red);
        setFont(getFont().deriveFont(Font.BOLD)); 
      }
    }
    
    protected void close()
    {
      if(commitFrame != null)
        commitFrame.dispose();
    }
  }

}

/**
 *  This is an EntryActionListener that will get an entry from entry_source
 *  and then copy them to destination_entry when actionPerformed() is called.
 **/
class ReadFeaturesActionListener extends EntryActionListener 
{
  final EntrySource entry_source;

  ReadFeaturesActionListener(final EntryEdit entry_edit,
                             final EntrySource entry_source,
                             final Entry destination_entry) 
  {
    super(entry_edit, destination_entry);
    this.entry_source = entry_source;
  }

  public void actionPerformed(final ActionEvent event) 
  {
    try 
    {
      if(getEntry().isReadOnly()) 
      {
        final String message =
          "the default entry is read only - cannot continue";
        new MessageDialog(getEntryEdit(), message);
      }

      final Entry source_entry =
        entry_source.getEntry(getEntryEdit().getEntryGroup().getBases(),
                              true);

      for(int i = 0; i < source_entry.getFeatureCount(); ++i) 
      {
        final Feature this_feature = source_entry.getFeature(i);
        try 
        {
          this_feature.copyTo(getEntry());
        }
        catch(OutOfRangeException e) 
        {
          throw new Error("internal error - unexpected exception: " + e);
        }
        catch(EntryInformationException e)
        {
          final String message =
            "couldn't move one of the features(" +
            this_feature.getIDString() + "): " + e.getMessage();
          new MessageDialog(getEntryEdit(), message);
        }
        catch(ReadOnlyException e) 
        {
          final String message =
            "the default entry is read only - cannot continue";
          new MessageDialog(getEntryEdit(), message);
          return;
        }
      }
    } 
    catch(OutOfRangeException e) 
    {
      new MessageDialog(getEntryEdit(),
                         "read failed: one of the features in " +
                         "the source entry has an out of range location");
    } 
    catch(IOException e) 
    {
      new MessageDialog(getEntryEdit(),
                         "read failed due to an IO error: " +
                         e.getMessage());
    }
    catch(NullPointerException e)
    {
      new MessageDialog(getEntryEdit(),
                         "read failed due to a null pointer error: " +
                         e.getMessage());
    }

  }

}

/**
 *  This is an EntryActionListener that will call saveEntry() when
 *  actionPerformed() is called.
 **/
class SaveEntryActionListener extends EntryActionListener 
{
  SaveEntryActionListener(final EntryEdit entry_edit,
                           final Entry entry) 
  {
    super(entry_edit, entry);
  }

  public void actionPerformed(final ActionEvent event) 
  {
    getEntryEdit().saveEntry(getEntry(), true, false, true,
                               DocumentEntryFactory.ANY_FORMAT);
  }
}

/**
 *  This is an EntryActionListener that will call saveEntry() when
 *  actionPerformed() is called.
 **/
class SaveEntryAsActionListener extends EntryActionListener 
{
  SaveEntryAsActionListener(final EntryEdit entry_edit,
                             final Entry entry) 
  {
    super(entry_edit, entry);
  }

  public void actionPerformed(final ActionEvent event) 
  {
    getEntryEdit().saveEntry(getEntry(), true, true, true,
                               DocumentEntryFactory.ANY_FORMAT);
  }
}

/**
 *  This is an EntryActionListener that will call saveEntry() when
 *  actionPerformed() is called.  The output file type will be EMBL.
 **/
class SaveEntryAsEMBLActionListener extends EntryActionListener 
{
  SaveEntryAsEMBLActionListener(final EntryEdit entry_edit,
                                 final Entry entry) 
  {
    super(entry_edit, entry);
  }

  public void actionPerformed(final ActionEvent event) 
  {
    getEntryEdit().saveEntry(getEntry(), true, true, false,
                               DocumentEntryFactory.EMBL_FORMAT);
  }
}

/**
 *  This is an EntryActionListener that will call saveEntry() when
 *  actionPerformed() is called.  The output file type will be GENBANK.
 **/
class SaveEntryAsGenbankActionListener extends EntryActionListener 
{
  SaveEntryAsGenbankActionListener(final EntryEdit entry_edit,
                                   final Entry entry) 
  {
    super(entry_edit, entry);
  }

  public void actionPerformed(final ActionEvent event) 
  {
    getEntryEdit().saveEntry(getEntry(), true, true, false,
                             DocumentEntryFactory.GENBANK_FORMAT);
  }
}

/**
 *  This is an EntryActionListener that will call saveEntry() when
 *  actionPerformed() is called.  The output file type will be GFF, with the
 *  sequence(if any) in FASTA format.
 **/
class SaveEntryAsGFFActionListener extends EntryActionListener 
{
  SaveEntryAsGFFActionListener(final EntryEdit entry_edit,
                                    final Entry entry) 
  {
    super(entry_edit, entry);
  }

  public void actionPerformed(final ActionEvent event) 
  {
    getEntryEdit().saveEntry(getEntry(), true, true, false,
                             DocumentEntryFactory.GFF_FORMAT);
  }
}

/**
 *  This is an EntryActionListener that will call saveEntry() when
 *  actionPerformed() is called.
 **/
class SaveEntryAsSubmissionActionListener extends EntryActionListener 
{
  SaveEntryAsSubmissionActionListener(final EntryEdit entry_edit,
                                       final Entry entry) 
  {
    super(entry_edit, entry);
  }

  public void actionPerformed(final ActionEvent event) 
  {
    getEntryEdit().saveEntry(getEntry(), false, true, false,
                             DocumentEntryFactory.EMBL_FORMAT);
  }
}

class FeatureDatabasePredicate implements FeaturePredicate
{

  public FeatureDatabasePredicate()
  {
  }

  public boolean testPredicate(final Feature feature)
  {
    if(feature.getEntry().getEMBLEntry() instanceof DatabaseDocumentEntry)
      return true;
    else
      return false;
  }

}


