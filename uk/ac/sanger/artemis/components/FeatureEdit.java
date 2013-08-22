/* FeatureEdit.java
 *
 * created: Tue Dec  1 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/FeatureEdit.java,v 1.70 2009-09-24 15:01:27 tjc Exp $
 **/

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.util.*;
import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.sequence.MarkerRange;

import uk.ac.sanger.artemis.io.ChadoCanonicalGene;
import uk.ac.sanger.artemis.io.DatabaseDocumentEntry;
import uk.ac.sanger.artemis.io.DocumentEntry;
import uk.ac.sanger.artemis.io.GFFDocumentEntry;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.OutOfDateException;
import uk.ac.sanger.artemis.io.LocationParseException;
import uk.ac.sanger.artemis.io.PartialSequence;
import uk.ac.sanger.artemis.io.QualifierLazyLoading;
import uk.ac.sanger.artemis.io.QualifierParseException;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.io.RangeVector;
import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.io.KeyVector;
import uk.ac.sanger.artemis.io.Location;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.EntryInformationException;
import uk.ac.sanger.artemis.io.StreamQualifier;
import uk.ac.sanger.artemis.io.QualifierInfo;
import uk.ac.sanger.artemis.io.ValidateFeature;

import uk.ac.sanger.artemis.components.ProgressThread;
import uk.ac.sanger.artemis.components.genebuilder.BasicGeneBuilderFrame;
import uk.ac.sanger.artemis.components.genebuilder.GeneBuilderFrame;
import uk.ac.sanger.artemis.components.genebuilder.GeneEditorPanel;
import uk.ac.sanger.artemis.components.genebuilder.GeneUtils;
import uk.ac.sanger.artemis.components.genebuilder.ProteinMapPanel;
import uk.ac.sanger.artemis.components.genebuilder.ReferencesPanel;
import uk.ac.sanger.artemis.components.genebuilder.cv.CVPanel;
import uk.ac.sanger.artemis.components.genebuilder.gff.PropertiesPanel;
import uk.ac.sanger.artemis.components.genebuilder.ortholog.MatchPanel;

import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.util.Date;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Set;
import java.util.Vector;
import java.util.Collections;
import java.util.Comparator;

import javax.swing.*;


/**
 *  FeatureEdit class
 *
 *  @author Kim Rutherford
 *  @version $Id: FeatureEdit.java,v 1.70 2009-09-24 15:01:27 tjc Exp $
 **/
public class FeatureEdit extends JPanel
                         implements EntryChangeListener, FeatureChangeListener 
{
  private static final long serialVersionUID = 1L;


  /** The choice of feature keys - created in createComponents(). */
  private KeyChoice key_choice;

  /** The choice of qualifiers - created in createComponents(). */
  private QualifierChoice qualifier_choice = null;

  private final static int LOCATION_TEXT_WIDTH = 80;

  /** The location text - set by updateLocation(). */
  private JTextField location_text = new JTextField(LOCATION_TEXT_WIDTH);

  /** When pressed - discard changes and dispose of the component. */
  private JButton cancel_button = new JButton("Cancel");

  /** When pressed - apply changes and keep the component open. */
  private JButton apply_button = new JButton("Apply");

  /** Edit area for qualifiers - created by createComponents(). */
  private QualifierTextArea qualifier_text_area;

  /** The Feature this object is displaying. */
  private Feature edit_feature;

  /**
   *  The GotoEventSource that was passed to the constructor - used for the
   *  "Goto Feature" button.
   **/
  private GotoEventSource goto_event_source;

  /** Entry containing the Feature this object is displaying. */
  private Entry edit_entry;

  /** EntryGroup that contains this Feature (passed to the constructor). */
  private EntryGroup entry_group;

  /**
   *  The datestamp of the RWCorbaFeature when updateFromFeature() was last
   *  called or null if this is not a RWCorbaFeature.
   **/
  private Date datestamp = null;

  /** The Selection that was passed to the constructor. */
  private Selection selection;

  /**
   *  The contents of the QualifierTextArea before the user edits anything.
   *  This is used to work out if anything has changed since the creation of
   *  the FeatureEdit.
   **/
  private final String orig_qualifier_text;

  private JFrame frame;
  
  /** controlled vocabulary tab */
  private CVPanel cvForm;
  
  /** GFF tab */
  private PropertiesPanel propertiesPanel;
  
  /** similarity/ortholog/paralog tab */
  private MatchPanel matchForm;
  
  private ReferencesPanel refPanel;
  
  private EntryInformation entry_information;
  
  private static boolean isTabbedView = false;
  
  private GeneEditorPanel editorPanel;
  
  private static UserDefinedQualifiers userDefinedQualifierFrame;
  
  /**
   *  Create a new FeatureEdit object from the given Feature.
   *  @param entry_group The EntryGroup that contains this Feature.
   *  @param selection The Selection operate on.  The operations are "Remove
   *    Range" and "Grab Range"
   *  @param goto_event_source The object the we will call gotoBase() on.
   **/
  public FeatureEdit(final Feature edit_feature,
                     final EntryGroup entry_group,
                     final Selection selection,
                     final GotoEventSource goto_event_source,
                     final JFrame frame) 
  {
    this.frame = frame;
    this.edit_feature = edit_feature;
    this.edit_entry   = edit_feature.getEntry();
    this.entry_group  = entry_group;
    this.selection    = selection;
    this.goto_event_source = goto_event_source;
    this.entry_information = edit_feature.getEntry().getEntryInformation();

    setLayout(new BorderLayout());
    createComponents();
    updateFromFeature();

    orig_qualifier_text = qualifier_text_area.getText();
  
    if(edit_feature.getEntry() != null)
    {
      edit_feature.getEntry().addEntryChangeListener(this);
      edit_feature.addFeatureChangeListener(this);

      frame.addWindowListener(new WindowAdapter()
      {
        public void windowClosing(WindowEvent event)
        {
          stopListening();
          frame.dispose();
          if(userDefinedQualifierFrame != null)
            userDefinedQualifierFrame.setVisible(false);
        }
      });
    }
    
    qualifier_text_area.requestFocus();
  }
  
  private boolean isPartialSequence()
  {
    return entry_group.getSequenceEntry().getBases().getSequence() instanceof PartialSequence;
  }
  
  /**
   * Set the feature to edit
   * @param edit_feature
   * @param isSet
   */
  public void setActiveFeature(final Feature edit_feature,
                               final boolean isSet)
  {
    if(!isPartialSequence() && isSet)
      setFeature();
    this.edit_feature = edit_feature;
    this.edit_entry   = edit_feature.getEntry();
    
    if(edit_feature.getEntry() != null)
    {
      edit_feature.getEntry().addEntryChangeListener(this);
      edit_feature.addFeatureChangeListener(this);
    }
    updateFromFeature();
  }

  /**
   *  Remove this object as a feature and entry change listener.
   **/
  public void stopListening() 
  {
    getEntry().removeEntryChangeListener(this);
    getFeature().removeFeatureChangeListener(this);
    if(cvForm != null)
      getFeature().removeFeatureChangeListener(cvForm);
    if(propertiesPanel != null)
      getFeature().removeFeatureChangeListener(propertiesPanel);
    if(matchForm != null)
      getFeature().removeFeatureChangeListener(matchForm);
  }

  /**
   *  Implementation of the EntryChangeListener interface.  We listen to
   *  EntryChange events so we can notify the user if of this component if the
   *  feature gets deleted.
   **/
  public void entryChanged(EntryChangeEvent event) 
  {
    switch(event.getType())
    {
      case EntryChangeEvent.FEATURE_DELETED:
        if(event.getFeature() == edit_feature) 
        {
          stopListening();
          frame.dispose();
        }
        break;
      default:
        // do nothing;
        break;
    }
  }

  /**
   *  Add an ActionListener to the Cancel JButton of this FeatureEdit.
   **/
  protected void addCancelActionListener(final ActionListener l) 
  {
    cancel_button.addActionListener(l);
  }

  /**
   *  Remove an ActionListener from the Cancel JButton of this FeatureEdit.
   **/
  protected void removeCancelActionListener(final ActionListener l) 
  {
    cancel_button.removeActionListener(l);
  }

  /**
   *  Add an ActionListener to the Apply JButton of this FeatureEdit.
   **/
  protected void addApplyActionListener(final ActionListener l) 
  {
    apply_button.addActionListener(l);
  }

  /**
   *  Implementation of the FeatureChangeListener interface.  We need to
   *  listen to feature change events from the Features in this object so that
   *  we can update the display.
   *  @param event The change event.
   **/
  public void featureChanged(FeatureChangeEvent event) 
  {
    getFeature().resetColour();
    // re-read the information from the feature
    switch(event.getType()) 
    {
      case FeatureChangeEvent.LOCATION_CHANGED:
        updateLocation();
        break;
      case FeatureChangeEvent.SEGMENT_CHANGED:
        updateLocation();
        break;
      case FeatureChangeEvent.KEY_CHANGED:
        updateKey();
        break;
      case FeatureChangeEvent.QUALIFIER_CHANGED:
        if(qualifier_text_area.getText().equals(orig_qualifier_text)) 
          updateFromFeature();
        else
        {
          if(!event.getFeature().getIDString().equals(getFeature().getIDString()))
            break;
          
          final String message =
            "warning: the qualifiers have changed outside the editor - " +
            "view now?";

          final YesNoDialog yes_no_dialog =
            new YesNoDialog(frame, message);

          if(yes_no_dialog.getResult()) 
            new FeatureViewer(getFeature());
        }
        break;
      default:
        updateFromFeature();
        break;
    }
  }


  /**
   *  Create all the components for this FeatureEdit component.
   **/
  private void createComponents()
  {
    qualifier_text_area = new QualifierTextArea();
    //qualifier_text_area.setWrapStyleWord(true);

    FlowLayout flowLayoutZeroHVgap = new FlowLayout(FlowLayout.LEADING, 0, 0);

    key_choice =
      new KeyChoice(getEntryInformation(),getFeature().getKey());
    key_choice.setLayout(flowLayoutZeroHVgap);
    final JPanel key_and_qualifier_panel = new JPanel();
    location_text.setBackground(Color.white);

    final JPanel key_panel = new JPanel(flowLayoutZeroHVgap);
    final JLabel locLabel = new JLabel("Location: ");
    final JLabel keyLabel = new JLabel("Key: ");
    keyLabel.setHorizontalAlignment(SwingConstants.RIGHT);
    keyLabel.setPreferredSize(locLabel.getPreferredSize());
    
    key_panel.add(keyLabel);
    key_panel.add(key_choice);

    key_and_qualifier_panel.setLayout(new BorderLayout());
    key_and_qualifier_panel.add(key_panel, "West");

    boolean isGFF = false;
    if(edit_feature.getEmblFeature().getEntry() instanceof GFFDocumentEntry)
      isGFF = true;
    
    qualifier_choice = new QualifierChoice(getEntryInformation(),
                                  key_choice.getSelectedItem(),null,
                                  isGFF);
    
    final JPanel qualifier_panel = new JPanel(new FlowLayout(FlowLayout.TRAILING,0,0));
    final JButton qualifier_add_button = new JButton("Add Qualifier:");

    qualifier_panel.add(qualifier_add_button);
    qualifier_panel.add(qualifier_choice);

    key_and_qualifier_panel.add(qualifier_panel, "East");

    key_choice.addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent _) 
      {
        qualifier_choice.setKey(key_choice.getSelectedItem());
      }
    });

    qualifier_add_button.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent e)
      {
        final String qualifier_name =
          (String)qualifier_choice.getSelectedItem();

        QualifierInfo qualifier_info =
          getEntryInformation().getQualifierInfo(qualifier_name);

        if(qualifier_info == null) 
          qualifier_info = new QualifierInfo(qualifier_name,
                              QualifierInfo.OPTIONAL_QUOTED_TEXT,
                              null, null, false);

        qualifier_text_area.append("/" + qualifier_name);

        switch(qualifier_info.getType()) 
        {
          case QualifierInfo.QUOTED_TEXT:
            if(qualifier_name.equals("GO")) 
            {
              // special case for /GO
              final java.util.Calendar calendar =
                      java.util.Calendar.getInstance();
                
              final Date current_time = calendar.getTime();
            
              final java.text.SimpleDateFormat date_format =
                  new java.text.SimpleDateFormat("yyyyMMdd");
            
              final StringBuffer result_buffer = new StringBuffer();
            
              date_format.format(current_time, result_buffer,
                                 new java.text.FieldPosition(java.text.DateFormat.DATE_FIELD));
            
              final String go_string = "aspect=; term=; GOid=GO:; "+
                                       "evidence=ISS; db_xref=GO_REF:0000001; " +
                                       "with=UniProtKB:; date=" + result_buffer;
              qualifier_text_area.append("=\"" + go_string + "\"");
            } 
            else if(qualifier_name.equals("controlled_curation"))
            {
              final java.util.Calendar calendar =
                      java.util.Calendar.getInstance();

              final Date current_time = calendar.getTime();

              final java.text.SimpleDateFormat date_format =
                  new java.text.SimpleDateFormat("yyyyMMdd");

              final StringBuffer result_buffer = new StringBuffer();

              date_format.format(current_time, result_buffer,
                                 new java.text.FieldPosition(java.text.DateFormat.DATE_FIELD));

              final String cc_string = "term=; db_xref=; date="+ result_buffer;
              qualifier_text_area.append("=\"" + cc_string + "\"");
            }
            else 
              qualifier_text_area.append ("=\"\"");
            break;

          case QualifierInfo.NO_VALUE:
          case QualifierInfo.OPTIONAL_QUOTED_TEXT:
            break;

          default:
            qualifier_text_area.append("=");
        }

        qualifier_text_area.append("\n");
      }
    });

    final JPanel middle_panel = new JPanel();
    middle_panel.setLayout(new BorderLayout(0,0));

    final JPanel lower_panel = new JPanel();
    lower_panel.setLayout(new BorderLayout(0,0));

    //final JPanel outer_location_button_panel = new JPanel(flowLayoutZeroHgap);
    //outer_location_button_panel.setLayout(new BorderLayout(0,0));

    final JToolBar location_button_panel = new JToolBar();
    location_button_panel.setRollover(true);
    //outer_location_button_panel.add(location_button_panel, "West");
    lower_panel.add(location_button_panel, "North");
    
    final JPanel location_panel = new JPanel(new BorderLayout(0,0));
    location_panel.add(locLabel, "West");
    location_panel.add(location_text, "Center");

    final JButton complement_button = new JButton("Complement");
    complement_button.setToolTipText("Complement position");
    location_button_panel.add(complement_button);
    complement_button.addActionListener(new ActionListener () 
    {
      public void actionPerformed(ActionEvent e) 
      {
        complementLocation();
      }
    });
    
    if(GeneUtils.isDatabaseEntry(getFeature().getEmblFeature()))
    {
      final JButton refresh_button = new JButton("Refresh");
      location_button_panel.add(refresh_button);
      refresh_button.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          refresh();
        }
      });
    }

    final JButton grab_button = new JButton("Grab Range");
    grab_button.setToolTipText("Add selected base range from feature coordinates");
    location_button_panel.add(grab_button);
    grab_button.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        grabSelectedRange();
      }
    });

    final JButton remove_button = new JButton("Remove Range");
    remove_button.setToolTipText("Remove selected base range from feature coordinates");
    location_button_panel.add(remove_button);
    remove_button.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent e)
      {
        removeSelectedRange();
      }
    });

    final JButton goto_button = new JButton("Goto Feature");
    goto_button.setToolTipText("Goto and select this feature");
    location_button_panel.add(goto_button);
    goto_button.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent e)
      {
        goto_event_source.gotoBase(getFeature().getFirstBaseMarker());
        getSelection().set(getFeature());
      }
    });

/*    final JButton select_button = new JButton("Select Feature");
    location_button_panel.add(select_button);
    select_button.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e) 
      {
        getSelection().set(getFeature());
      }
    });*/

    if(Options.getOptions().getPropertyTruthValue("sanger_options"))
    {
      // a PSU only hack 
      final JButton tidy_button = new JButton("Tidy");
      location_button_panel.add(tidy_button);
      tidy_button.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e) 
        {
          try 
          {
            tidy();
            tidyGO();
          } 
          catch(QualifierParseException exception) 
          {
            new MessageDialog(frame,
                 "Cannot tidy - qualifier error: " +
                  exception.getMessage());
          }
        }
      });
    }
    
    final JButton transferAnnotationBbutton = new JButton("TAT");
    location_button_panel.add(transferAnnotationBbutton);
    transferAnnotationBbutton.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e) 
      {
      	new TransferAnnotationTool(getFeature(), entry_group, matchForm);
      }
    });


    if(Options.isUnixHost())
    {
      JButton oo_edit_button = new JButton("ObjectEdit");
      location_button_panel.add(oo_edit_button);
      final uk.ac.sanger.artemis.editor.BigPane bp =
                   new uk.ac.sanger.artemis.editor.BigPane(getEntryInformation());

      oo_edit_button.addActionListener(new ActionListener ()
      {
        public void actionPerformed(ActionEvent e)
        {
          String qualifier_txt = qualifier_text_area.getText();
        
          File base_dir = getBaseDirectoryFromEntry(edit_entry);
          String baseDirStr = "";

          if(base_dir != null)
            baseDirStr = base_dir.getAbsolutePath() + 
                              System.getProperty("file.separator");

          StringReader strRead = new StringReader(qualifier_txt);
          BufferedReader buff = new BufferedReader(strRead);
          String line;
          final Hashtable<String, Vector<String>> dataFile = new Hashtable<String, Vector<String>>();
          try
          {
            int ind;
            while((line = buff.readLine()) != null)
            {
              if(line.startsWith("/fasta_file="))
              {
                if((ind = line.indexOf(':'))>-1)
                  line = baseDirStr+line.substring(ind+1);
                else
                  line = baseDirStr+line.substring(13);
                
                ind = line.lastIndexOf("\"");
                if(ind > -1)
                  line = line.substring(0, ind);
                
                Vector<String> v;
                if(dataFile.containsKey("fasta"))
                  v = dataFile.get("fasta");
                else
                  v = new Vector<String>();
                v.add(line);
                dataFile.put("fasta",v);
              }
              else if(line.startsWith("/blastp_file="))
              {
                if((ind = line.indexOf(':'))>-1)
                  line = baseDirStr+line.substring(ind+1);
                else
                  line = baseDirStr+line.substring(14);
                
                ind = line.lastIndexOf("\"");
                if(ind > -1)
                  line = line.substring(0, ind);
                
                Vector<String> v;
                if(dataFile.containsKey("blastp"))
                  v = dataFile.get("blastp");
                else
                  v = new Vector<String>();
                v.add(line);
                dataFile.put("blastp",v);
              }
              else if(line.startsWith("/blastp+go_file="))
              {
                if((ind = line.indexOf(':'))>-1)
                  line = line.substring(ind+1);
                else
                  line = baseDirStr+line.substring(17);
                ind = line.lastIndexOf("\"");
                if(ind > -1)
                  line = line.substring(0, ind);
                
                Vector<String> v;
                if(dataFile.containsKey("blastp+go"))
                  v = dataFile.get("blastp+go");
                else
                  v = new Vector<String>();
                v.add(line);
                dataFile.put("blastp+go",v);
              }   
            }
          }
          catch(IOException ioe){}
 
          FeatureEdit.this.setCursor(new Cursor(Cursor.WAIT_CURSOR)); 
          final ProgressThread progress = new ProgressThread(null,
                                        "Loading Data....");
          progress.start();
          SwingWorker ooEd = new SwingWorker()
          {
            public Object construct()
            {
              // find overlaping features
              final Location this_loc = edit_feature.getLocation();
              final int this_start = this_loc.getFirstBase();
              final int this_end   = this_loc.getLastBase();

              FeaturePredicate predicate = new FeaturePredicate()
              {
                public boolean testPredicate (final Feature feature) 
                {
                  final Location loc = feature.getLocation();

                  final int start = loc.getFirstBase();
                  final int end   = loc.getLastBase();

                  if((start > this_start &&
                      start < this_end) ||
                     (end > this_start &&
                      end < this_end))
                  {
                    final String note = feature.getNote();
                    if(note != null &&
                       note.indexOf("Pfam")>-1)
                      return true;
                  }
                  
                  return false;
                }
              };

              final FeatureVector overlapFeatures = new FeatureVector();
              final FeatureEnumeration featureEnum = entry_group.features();
              while(featureEnum.hasMoreFeatures()) 
              {
                Feature this_feature = featureEnum.nextFeature();
                if(predicate.testPredicate(this_feature))
                  overlapFeatures.add(this_feature);
              }

              // show object editor
              try
              {
                if(dataFile.size() > 0)
                  bp.set(dataFile, qualifier_text_area, overlapFeatures,
                         edit_feature, matchForm, cvForm);
                else
                  JOptionPane.showMessageDialog(null,"No results files.",
                      "Warning",
                      JOptionPane.WARNING_MESSAGE); 
              }
              catch(ArrayIndexOutOfBoundsException e)
              {
                JOptionPane.showMessageDialog(null,"No results files.",
                        "Warning",
                        JOptionPane.WARNING_MESSAGE);
              }

              progress.finished();
              FeatureEdit.this.setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
              return null;
            }
          };
          ooEd.start();   
        }
      });
    }

    final JButton userQualifiers = new JButton("User Qualifiers");
    userQualifiers.setToolTipText("User defined qualifier selection");
    location_button_panel.add(userQualifiers);
    userQualifiers.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e) 
      {
        if(userDefinedQualifierFrame == null) 
        {
          userDefinedQualifierFrame = new UserDefinedQualifiers();
          userDefinedQualifierFrame.setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
        }

        userDefinedQualifierFrame.pack();
        final JFrame topFrame = 
              (JFrame) SwingUtilities.getWindowAncestor(FeatureEdit.this);
        Point p = topFrame.getLocationOnScreen();
        p.x -= userDefinedQualifierFrame.getWidth();
        if(p.x < 10)
          p.x = 10;
        userDefinedQualifierFrame.setLocation(p);

        userDefinedQualifierFrame.setQualifierTextArea(qualifier_text_area);
        userDefinedQualifierFrame.setSelection(selection);
        userDefinedQualifierFrame.setVisible(true);
      }
    });
    
    middle_panel.add(location_panel, "North");
    add(key_and_qualifier_panel, "North");

    cancel_button.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e) 
      {
        if(edit_feature.getEntry() != null)
          stopListening();
        frame.dispose();
        if(userDefinedQualifierFrame != null)
          userDefinedQualifierFrame.setVisible(false);
      }
    });

    final JButton ok_button = new JButton("OK");
    if(!getFeature().isReadOnly())
    {
      ok_button.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          if(setFeature()) 
          {
            stopListening();
            
            if(propertiesPanel != null)
              propertiesPanel.updateSettings();
            
            frame.dispose();
            if(userDefinedQualifierFrame != null)
              userDefinedQualifierFrame.setVisible(false);
          }
        }
      });

      apply_button.addActionListener(new ActionListener() 
      {
        public void actionPerformed(ActionEvent e) 
        {
          setFeature();
        }
      });
    }

    final FlowLayout flow_layout =
                 new FlowLayout(FlowLayout.CENTER, 18, 1);

    final JPanel ok_cancel_update_panel = new JPanel(flow_layout);
    Box fillerBox = Box.createHorizontalBox();
    if(GeneUtils.isDatabaseEntry(getFeature().getEmblFeature()))
    {
      cvForm = new CVPanel(getFeature());
      cvForm.setBackground(Color.WHITE);

      matchForm = new MatchPanel(getFeature(), 
          (DocumentEntry)getFeature().getEmblFeature().getEntry());
      matchForm.setBackground(Color.WHITE);
      
      propertiesPanel = new PropertiesPanel(getFeature());
      propertiesPanel.setBackground(Color.WHITE);

      refPanel = new ReferencesPanel(getFeature());
      refPanel.setBackground(Color.WHITE);
      
      addGffAnnotationView(lower_panel);
      
      final JCheckBox tabbedView = new JCheckBox("Tab View", isTabbedView);
      tabbedView.addItemListener(new ItemListener()
      {
        public void itemStateChanged(ItemEvent e)
        {
          isTabbedView = tabbedView.isSelected();
          addGffAnnotationView(lower_panel);
          lower_panel.revalidate();
          lower_panel.repaint();
        }
      });


      final JCheckBox oneView = new JCheckBox("Overview", false);
      oneView.addItemListener(new ItemListener()
      {
        public void itemStateChanged(ItemEvent e)
        {
          if(setFeature()) 
          {
            stopListening();
            
            if(propertiesPanel != null)
              propertiesPanel.updateSettings();
            
            frame.dispose();
            System.setProperty("basic", "true");
            new BasicGeneBuilderFrame(getFeature(), entry_group,
                selection, null);
          }
        }
      });

      if(((GFFStreamFeature)getFeature().getEmblFeature()).getChadoGene() != null)
        ok_cancel_update_panel.add(oneView);
      ok_cancel_update_panel.add(tabbedView);
      fillerBox.add(Box.createHorizontalStrut( 
          tabbedView.getPreferredSize().width ));
    }
    else
      lower_panel.add(new JScrollPane(qualifier_text_area), "Center");
    
    if(!getFeature().isReadOnly()) 
      ok_cancel_update_panel.add(ok_button);

    ok_cancel_update_panel.add(cancel_button);

    if(!getFeature().isReadOnly()) 
      ok_cancel_update_panel.add(apply_button);
    ok_cancel_update_panel.add(fillerBox);
    
    add(ok_cancel_update_panel, "South");
    
    middle_panel.add(lower_panel, "Center");

    add(middle_panel, "Center");
  }

  /**
   * Refresh the annotation for the active feature
   */
  private void refresh()
  {
    // refresh from the database
    final DatabaseDocument originalDocument =
      (DatabaseDocument)((DocumentEntry)edit_feature.getEmblFeature().getEntry()).getDocument();

    final Set<String> uniquenames = ((GFFStreamFeature)edit_feature.getEmblFeature()).getSegmentRangeStore().keySet();
    final Iterator<String> it = uniquenames.iterator();
    final String uniquename = it.next();
    final DatabaseDocument newDocument = new DatabaseDocument(originalDocument,
        uniquename, null, true, null);
    newDocument.setLazyFeatureLoad(false);
    newDocument.setReadChildren(false);
    
    try
    {
      DatabaseDocumentEntry dbentry = new DatabaseDocumentEntry(newDocument, null);
      uk.ac.sanger.artemis.io.Feature databaseFeature = dbentry.getAllFeatures().featureAt(0);
      
      // compare timelastmodified
      Qualifier qualifier = edit_feature.getQualifierByName("timelastmodified");
      String active_timelastmodified = (String)qualifier.getValues().get(0);
      qualifier = databaseFeature.getQualifierByName("timelastmodified");
      String database_timelastmodified = (String)qualifier.getValues().get(0);
      
      if(active_timelastmodified.equals(database_timelastmodified))
      {
        JOptionPane.showMessageDialog(this, 
            "No new changes found for the feature\n"+
            uniquename+"\n"+
            "in the database since:\n"+database_timelastmodified, 
            "No Updates", JOptionPane.INFORMATION_MESSAGE);
        return;
      }
      else
      {
        JOptionPane.showMessageDialog(this, 
            "Changes found for the feature\n"+
            uniquename+"\n"+
            "in the database at:\n"+database_timelastmodified, 
            "Changes Found", JOptionPane.INFORMATION_MESSAGE);
      }
      
      final QualifierVector db_qv = databaseFeature.getQualifiers();
      final QualifierVector qv = edit_feature.getQualifiers();
      final QualifierVector new_qv = new QualifierVector();
      
      for(int i=0; i<qv.size(); i++)
      {
        Qualifier q = (Qualifier)qv.get(i);
        if(q.getName().equals("Parent") ||
           q.getName().equals("Derives_from") ||
           q.getName().equals("ID")  )
          new_qv.addQualifierValues(q);
      }
      
      
      for(int i=0; i<db_qv.size(); i++)
      {
        Qualifier q = (Qualifier)db_qv.get(i);
        if(q.getName().equals("Parent") ||
           q.getName().equals("Derives_from") ||
           q.getName().equals("ID")  )
          continue;
        new_qv.add(q);  
      }
      
      edit_feature.getQualifiers().removeAllElements();
      edit_feature.getQualifiers().addAll(new_qv);
      updateQualifiers();
    }
    catch(EntryInformationException e) {}
    catch(IOException e) {}    
  }
  
  /**
   * Add the annotation view as tabbed or in a single pane
   * @param lower_panel
   */
  private void addGffAnnotationView(final JPanel lower_panel)
  {
    Component c[] = lower_panel.getComponents();
    
    if(isTabbedView)
    {
      for(int i=0; i<c.length; i++)
      {
        if(c[i] instanceof JScrollPane)
          lower_panel.remove(c[i]);
      }
      // tabbed pane of core and cv annotaion
      JTabbedPane tabbedPane = new JTabbedPane();
      
      JScrollPane jspGff = new JScrollPane(propertiesPanel);
      propertiesPanel.setVisible(true);    // ensure visible
      //jspGff.setPreferredSize(jspCore.getPreferredSize());
      tabbedPane.add("Properties", jspGff);
      
      JScrollPane jspCore = new JScrollPane(qualifier_text_area);
      tabbedPane.add("Core", jspCore);
      JScrollPane jspCV = new JScrollPane(cvForm);
      cvForm.setVisible(true);      // ensure visible
      jspCV.setPreferredSize(jspCore.getPreferredSize());
      tabbedPane.add("CV", jspCV);
      
      
      JScrollPane jspRef = new JScrollPane(refPanel);
      refPanel.setVisible(true);      // ensure visible
      jspRef.setPreferredSize(jspCore.getPreferredSize());
      tabbedPane.add("References", jspRef);
      
      JScrollPane jspOrtholog = new JScrollPane(matchForm);
      matchForm.setVisible(true);   // ensure visible
      jspOrtholog.setPreferredSize(getPreferredSize());
      tabbedPane.add("Match", jspOrtholog);
      
      lower_panel.add(tabbedPane, "Center");
    }
    else
    {
      for(int i=0; i<c.length; i++)
      {
        if(c[i] instanceof JTabbedPane)
          lower_panel.remove(c[i]);
      }
      
      editorPanel = new GeneEditorPanel(qualifier_text_area, cvForm,
          refPanel, matchForm, propertiesPanel);
      JScrollPane jsp = new JScrollPane(editorPanel);
          
      jsp.setPreferredSize(
          new Dimension(qualifier_text_area.getPreferredSize().width,
                        (int)(qualifier_text_area.getPreferredSize().height*1.5)));
      // single panel containing annotation forms
      lower_panel.add(jsp, "Center");
    }
  }
  
  /**
   *  Return the dirtectory that the given entry was read from.
   **/
  private File getBaseDirectoryFromEntry(final Entry entry)
  {
    final uk.ac.sanger.artemis.io.Entry embl_entry = entry.getEMBLEntry();

    if(embl_entry instanceof DocumentEntry)
    {
      final DocumentEntry document_entry =(DocumentEntry) embl_entry;

      if(document_entry.getDocument() instanceof FileDocument)
      {
        final FileDocument file_document =
         (FileDocument) document_entry.getDocument();

        if(file_document.getFile().getParent() != null)
          return new File(file_document.getFile().getParent());
      }
    }

    return null;
  }

  /**
   *  Read the key, location and qualifier information from the feature and
   *  update the components.
   **/
  private void updateFromFeature() 
  {
    datestamp = getFeature().getDatestamp();

    updateKey();
    updateLocation();
    updateQualifiers();
  }

  /**
   *  Read the location from the feature and update the location field.
   **/
  private void updateLocation() 
  {
    location_text.setText(getFeature().getLocation().toStringShort());
  }

  /**
   *   Complement the current location_text.
   **/
  private void complementLocation() 
  {
    if(GeneUtils.isDatabaseEntry(getFeature().getEmblFeature()))
    {
      final ChadoCanonicalGene chadoGene = 
        ((GFFStreamFeature)getFeature().getEmblFeature()).getChadoGene();
      GeneUtils.complementGeneModel(chadoGene);
      ((GeneBuilderFrame)frame).dispose(true);
      return;  
    }
    
    if(rationalizeLocation()) 
    {
      if(location_text.getText().startsWith("complement(")) 
      {
        final String new_text = location_text.getText().substring(11);
        if (new_text.endsWith(")")) 
        {
          final String new_text2 =
            new_text.substring(0, new_text.length () - 1);
          location_text.setText(new_text2);
        } 
        else 
          location_text.setText(new_text);
      }
      else 
      {
        final String new_text = location_text.getText ();
        location_text.setText("complement(" + new_text + ")");
      }
    } 
    else 
      new MessageDialog(frame, "complement failed - " +
                        "current location cannot be parsed");
  }

  /**
   *  Tidy the qualifiers by removing any indication that the annotation has
   *  been transferred from another gene.  Remove the "transferred_" part of
   *  all qualifier names and the "[[FROM sc_001234 abc1]]" bit of each
   *  corresponding qualifier value.
   **/
  private void tidy() throws QualifierParseException
  {
    final StringBuffer buffer = new StringBuffer();
    final QualifierVector qualifiers =
      qualifier_text_area.getParsedQualifiers(getEntryInformation());

    for(Qualifier this_qualifier: qualifiers) 
    {
      final QualifierInfo qualifier_info =
            getEntryInformation().getQualifierInfo(this_qualifier.getName());
      final StringVector qualifier_strings =
          StreamQualifier.toStringVector(qualifier_info, this_qualifier);
      
      for(String qualifier_string: qualifier_strings)
        buffer.append(tidyHelper(qualifier_string) + "\n");
    }
    qualifier_text_area.setText(buffer.toString());
  }

  private void tidyGO() throws QualifierParseException
  {
    String qualifier_txt = qualifier_text_area.getText();
    BufferedReader buff = new BufferedReader(new StringReader(qualifier_txt));
    final Vector<String> qual_str = new Vector<String>();
    try
    {
      String line;
      while((line = buff.readLine()) != null)
        qual_str.add(line);
    }
    catch(IOException ioe){}

    Comparator<String> comparator = new Comparator<String>()
    {
      public int compare (String fst, String snd)
      {
        if( !fst.startsWith("/GO") ||
            !snd.startsWith("/GO") )
          return 0;
        return fst.compareTo(snd);
      }
    };
    Collections.sort(qual_str, comparator);
    
    StringBuffer buffer = new StringBuffer();
    for(String qualifier_string: qual_str)
      buffer.append(tidyHelper(qualifier_string) + "\n");

    qualifier_text_area.setText(buffer.toString());
  }


  /**
   *  Perform the action of tidy() on the String version of one qualifier.
   **/
  private String tidyHelper(final String qualifier_string)
  {
    final String temp_string;

    if (qualifier_string.startsWith("/transferred_")) 
      temp_string = "/" + qualifier_string.substring(13);
    else 
      temp_string = qualifier_string;

    final int start_index = temp_string.indexOf ("<<FROM ");
    final int end_index = temp_string.indexOf (">>");

    if(start_index != -1 && end_index != -1) 
    {
      if(temp_string.length() > end_index + 2 &&
         temp_string.charAt(end_index + 2) == ' ') 
        return temp_string.substring(0, start_index) +
               temp_string.substring(end_index + 3);
      else 
        return temp_string.substring(0, start_index) +
               temp_string.substring(end_index + 2);
    } 
    else 
      return temp_string;
  }

  /**
   *  Add the currently selected range to location_text.
   **/
  private void grabSelectedRange() 
  {
    if(!rationalizeLocation ()) 
    {
      new MessageDialog(frame,
                        "grab failed - current location cannot be parsed");
      return;
    }

    final Range selected_range = getSelection().getSelectionRange();

    if(selected_range == null) 
    {
      new MessageDialog(frame, "grab failed - nothing is selected");
      return;
    }

    // save it in case it gets mangled
    final String old_location_text = location_text.getText();

    if (old_location_text.endsWith ("))")) 
    {
      final String new_text =
        old_location_text.substring(0, old_location_text.length () - 2);

      location_text.setText(new_text + "," + selected_range.getStart () +
                            ".." + selected_range.getEnd () + "))");
    } 
    else
    {
      if(old_location_text.endsWith (")")) 
      {
        final String new_text =
          old_location_text.substring(0, old_location_text.length () - 1);

        location_text.setText(new_text + "," + selected_range.getStart () +
                              ".." + selected_range.getEnd () + ")");
      } 
      else
        location_text.setText(old_location_text + "," +
                              selected_range.getStart () +
                              ".." + selected_range.getEnd ());
    }

    if(!rationalizeLocation())
    {
      location_text.setText (old_location_text);
      new MessageDialog(frame, "grab failed - location cannot be parsed after " +
                        "grabbing");
    }
  }

  /**
   *  Remove the currently selected range of bases from location_text.
   **/
  private void removeSelectedRange()
  {
    if(!rationalizeLocation())
    {
      new MessageDialog(frame,
                        "grab failed - current location cannot be parsed");
      return;
    }

    final MarkerRange selected_marker_range =
                                           getSelection().getMarkerRange();

    if(selected_marker_range == null) 
    {
      new MessageDialog(frame, "remove range failed - no bases are selected");
      return;
    }

    final Range selected_range = selected_marker_range.getRawRange();

    if (selected_marker_range.getStrand() != getFeature().getStrand()) 
    {
      new MessageDialog(frame, "remove range failed - you need to select " +
                         "some bases on the other strand");
      return;
    }

    final Location location;

    try 
    {
      location = new Location(location_text.getText ());
    }
    catch (LocationParseException e) 
    {
      // this shouldn't happen because we called rationalizeLocation ()
      throw new Error("internal error - unexpected exception: " + e);
    }

    final Range location_total_range = location.getTotalRange();

    if(!selected_range.overlaps(location_total_range))
    {
      new MessageDialog(frame, "remove range failed - the range you " +
                        "selected does not overlap the feature");
      return;
    }

    if(selected_range.contains(location_total_range)) 
    {
      new MessageDialog(frame, "remove range failed - the range you " +
                        "selected overlaps the whole feature");
      return;
    }

    final RangeVector location_ranges = location.getRanges();
    final boolean location_is_complemented = location.isComplement();

    final RangeVector new_ranges = new RangeVector();

    // if the selected_range completely covers a range remove the
    // range. otherwise if the selected_range is completely within one of the
    // ranges two new ranges are created.  if the selected_range is not
    // completely contained then the appropriate end of the range is truncated
    for(int i = 0; i < location_ranges.size(); ++i) 
    {
      final Range this_range = (Range)location_ranges.elementAt(i);

      if(selected_range.overlaps(this_range)) 
      {
        try 
        {
          if(this_range.contains(selected_range) &&
             this_range.getStart() != selected_range.getStart() &&
             this_range.getEnd() != selected_range.getEnd()) 
          {
            // chop a piece out of the middle and make two new ranges
            final Range new_start_range =
                    this_range.change(selected_range.getEnd() + 1,
                                      this_range.getEnd());
            new_ranges.add(new_start_range);
            final Range new_end_range =
              this_range.change(this_range.getStart(),
                                selected_range.getStart() - 1);
            new_ranges.add(new_end_range);
          } 
          else
          {
            if(selected_range.contains(this_range)) {
              // delete (ie. don't copy) the range
            } 
            else
            {
              if(this_range.getStart() < selected_range.getStart()) 
              {
                // truncate the end of the range
                final Range new_start_range =
                  this_range.change(this_range.getStart(),
                                    selected_range.getStart() - 1);
                new_ranges.add(new_start_range);
              } 
              else
              {
                if(this_range.getEnd() > selected_range.getEnd())
                {
                  // truncate the start of the range
                  final Range new_end_range =
                    this_range.change(selected_range.getEnd() + 1,
                                      this_range.getEnd());
                  new_ranges.add(new_end_range);
                } 
                else
                  throw new Error ("internal error - can't remove range");
              }
            }
          }
        }
        catch(OutOfRangeException e)
        {
          throw new Error ("internal error - unexpected exception: " + e);
        }
      }
      else  
        new_ranges.add(this_range); // copy it unchanged
    }

    final Location new_location =
      new Location(new_ranges, location_is_complemented);

    location_text.setText(new_location.toStringShort());
  }

  /**
   *  Attempt to parse the current location_text as a Location.  If it can be
   *  parsed it will be canonicalized (ie. the complement, if any, will be
   *  outermost).  Returns true if and only if the location_text could be
   *  parsed.
   **/
  private boolean rationalizeLocation()
  {
    try 
    {
      final Location location = new Location(location_text.getText());
      location_text.setText(location.toStringShort());
      return true;
    }
    catch(LocationParseException e) 
    {
      return false;
    }
  }


  /**
   *  Edit the qualifiers of this Feature in an external editor.  The
   *  qualifiers will be set when the editor finishes.  This method works by
   *  writing the qualifiers to a temporary file and the sequence of the
   *  feature to a different file.
   *  @param editor_extra_args Extra arguments to pass to the editor.  null
   *    means there are no extra args.
   **/
  /*private void externalEdit(final String[] editor_extra_args) 
  {
    try 
    {
      final String pre_edit_text = qualifier_text_area.getText();

      // write to a temporary file
      final Date current_time = calendar.getTime();

      final String temp_file_name =
               "/tmp/artemis_temp." + current_time.getTime();

      final File temp_file = new File(temp_file_name);

      final FileWriter out_writer    = new FileWriter(temp_file);
      final PrintWriter print_writer = new PrintWriter(out_writer);

      print_writer.write(qualifier_text_area.getText());
      print_writer.close();
      out_writer.close();

      final File sequence_temp_file = new File(temp_file_name + ".seq");
      final FileWriter sequence_out_writer =
                                     new FileWriter(sequence_temp_file);
      final PrintWriter sequence_print_writer =
                                   new PrintWriter(sequence_out_writer);

      getFeature().writeBasesOfFeature(sequence_print_writer);
      sequence_print_writer.close();
      sequence_out_writer.close();

      final String editor_path =
        Options.getOptions().getProperty("external_editor");

      final String[] process_args;

      if(editor_extra_args == null) 
      {
        process_args = new String[1];
        process_args[0] = temp_file.getCanonicalPath();
      } 
      else
      {
        process_args = new String[editor_extra_args.length + 1];
        System.arraycopy(editor_extra_args, 0, process_args, 0,
                         editor_extra_args.length);
        process_args[process_args.length - 1] = temp_file.getCanonicalPath();
      }


      System.out.println(editor_path);
      for(int i=0;i<process_args.length;i++)
        System.out.println(process_args[i]);

      final Process process =
        ExternalProgram.startProgram(editor_path, process_args);

      final ProcessWatcher process_watcher =
                                new ProcessWatcher(process, "editor", false);

      final Thread watcher_thread = new Thread(process_watcher);
      watcher_thread.start();

      final ProcessWatcherListener listener = new ProcessWatcherListener()
      {
        public void processFinished(final ProcessWatcherEvent event)
        {
          try 
          {
            final FileReader file_reader = new FileReader(temp_file);
            final BufferedReader buffered_reader = 
                                     new BufferedReader(file_reader);

            final StringBuffer buffer = new StringBuffer();
            String line;

            while((line = buffered_reader.readLine()) != null) 
              buffer.append(line + "\n");

            //ensure current qualifier text has not changed
            if(!qualifier_text_area.getText().equals(pre_edit_text))
            {
              final String message =
                  "the qualifiers have changed - apply changes from the " +
                  "external editor?";

              final YesNoDialog yes_no_dialog =
                  new YesNoDialog(frame, message);

              if(!yes_no_dialog.getResult())
                return;
            }

            qualifier_text_area.setText(buffer.toString());
            temp_file.delete();
            sequence_temp_file.delete();

            return;
          }
          catch(IOException e) 
          {
            new MessageDialog(frame, "an error occured while " +
                              "reading from the editor: " + e);
          }
        }
      };

      process_watcher.addProcessWatcherListener(listener);
    }
    catch(IOException e) 
    {
      new MessageDialog(frame, "error while starting editor: " + e);
    } 
    catch(ExternalProgramException e) 
    {
      new MessageDialog(frame, "error while starting editor: " + e);
    }
  }*/

  /**
   *  Read the qualifiers from the feature and update the qualifier JTextArea.
   **/
  private void updateQualifiers() 
  {
    if(GeneUtils.isDatabaseEntry(getFeature().getEmblFeature()))
      GeneUtils.addLazyQualifiers((GFFStreamFeature)getFeature().getEmblFeature());
    
    qualifier_text_area.setText(getQualifierString());
    
    if(GeneUtils.isDatabaseEntry(getFeature().getEmblFeature()))
    {  
      // load synonym
      if(cvForm != null)
        cvForm.updateFromFeature(getFeature());
      
      if(refPanel != null)
        refPanel.updateFromFeature(getFeature());

      if(propertiesPanel != null)
        propertiesPanel.updateFromFeature(getFeature());

      if(matchForm != null)
        matchForm.updateFromFeature(getFeature());

      if(!isTabbedView && editorPanel != null)
        editorPanel.updatePanelState();
    }
  }

  /**
   *  Return a string containing one qualifier per line.  These are the
   *  original qualifiers, not the qualifiers from the qualifier_text_area.
   **/
  private String getQualifierString() 
  {
    final StringBuffer buffer = new StringBuffer();
    final QualifierVector qualifiers = getFeature().getQualifiers();       
    
    for(Qualifier this_qualifier: qualifiers) 
    {
      //
      // strip out CV qualifiers
      //
      if( (cvForm != null && CVPanel.isCvTag(this_qualifier)) ||
          (refPanel != null && ReferencesPanel.isReferenceTag(this_qualifier)) ||
          (propertiesPanel != null && PropertiesPanel.isPropertiesTag(this_qualifier, getFeature())) ||
          (matchForm != null && MatchPanel.isMatchTag(this_qualifier)) ||
          (propertiesPanel != null && ProteinMapPanel.isProteinMapElement(this_qualifier)) )
        continue;
      
      if(this_qualifier instanceof QualifierLazyLoading)
        ((QualifierLazyLoading)this_qualifier).setForceLoad(true);

      final QualifierInfo qualifier_info =
                       getEntryInformation().getQualifierInfo(this_qualifier.getName());

      final StringVector qualifier_strings =
                       StreamQualifier.toStringVector(qualifier_info, this_qualifier);

      for(String qualStr: qualifier_strings)
        buffer.append(qualStr + "\n");
    }

    return buffer.toString();
  }


  
  /**
   *  Set the key, location and qualifiers of the feature to be the same as
   *  what values currently shown in the components.
   *  @return true if and only if action succeeds.  It may fail because of an
   *    illegal location or qualifier, in which case a message will be
   *    displayed before returning.
   **/
  private boolean setFeature() 
  {
    final Key key = key_choice.getSelectedItem();
    final KeyVector possible_keys = getEntryInformation().getValidKeys();

    if(possible_keys != null && !possible_keys.contains(key)) 
    {
      final YesNoDialog dialog =
        new YesNoDialog(frame, "Add this new key: " + key + "?");

      if(dialog.getResult()) // yes
        getEntryInformation ().addKey (key);
      else 
        return false;
    }

    final Location location;
    try 
    {
      location = new Location(location_text.getText());
    }
    catch(LocationParseException exception) 
    {
      final String error_string = exception.getMessage ();
      System.out.println(error_string);
      new MessageDialog(frame, "Cannot apply changes because of location error: " +
                        error_string);
      return false;
    }


    final QualifierVector qualifiers;

    try 
    {
      qualifiers =
        qualifier_text_area.getParsedQualifiers(getEntryInformation ());
      
      updateGffIds(qualifiers);
      
      // if using controlled vocab form
      if(cvForm != null)
      {
        QualifierVector cvQualifiers = cvForm.getCvQualifiers();
        if(cvQualifiers != null && cvQualifiers.size() > 0)
          qualifiers.addAll(cvQualifiers);
      }
      
      if(refPanel != null)
      {
        QualifierVector refQualifiers = refPanel.getQualifiers();
        if(refQualifiers != null && refQualifiers.size() > 0)
          qualifiers.addAll(refQualifiers);
      }
      
      if(propertiesPanel != null)
      {
        QualifierVector gffQualifiers = propertiesPanel.getGffQualifiers(getFeature());
        if(gffQualifiers != null && gffQualifiers.size() > 0)
          qualifiers.addAll(gffQualifiers);
        
        QualifierVector mapQualifiers = ProteinMapPanel.getProteinMapQualifiers(getFeature());
        if(mapQualifiers != null && mapQualifiers.size() > 0)
          qualifiers.addAll(mapQualifiers);
      }
      
      if(matchForm != null)
      {
        QualifierVector orthologQualifiers = matchForm.getMatchQualifiers();
        if(orthologQualifiers != null && orthologQualifiers.size() > 0)
          qualifiers.addAll(orthologQualifiers);
      }
      
      final String goErrs = ValidateFeature.validateGO(qualifiers, getEntryInformation());
      if(goErrs.length()>0)
      {
        Object[] options = { "CANCEL", "CONTINUE" };
        int opt = JOptionPane.showOptionDialog(null, goErrs, "GO errors",
            JOptionPane.DEFAULT_OPTION, JOptionPane.WARNING_MESSAGE,
            null, options, options[0]);

        if(opt == 0)
          return false;
      }
      //if(similarityTextArea != null)
      //  similarityTextArea.checkForChanges();
    }
    catch(QualifierParseException exception) 
    {
      final String error_string = exception.getMessage();
      System.out.println(error_string);
      new MessageDialog(frame, "Cannot apply changes because of a qualifier " +
                        "error: " + error_string);
      return false;
    }

    try 
    {
      entry_group.getActionController().startAction();

      try 
      {
        getFeature().set(datestamp, key, location, qualifiers);
      }
      catch(OutOfDateException e) 
      {
        final YesNoDialog dialog =
          new YesNoDialog(frame, "the feature has changed since the edit " +
                          "window was opened, continue?");

        if(dialog.getResult())  // yes - ignore the datestamp
          getFeature().set(key, location, qualifiers);
        else 
          return false;
      }
      catch(java.lang.Error err)
      {
        err.printStackTrace();

        if(err.getMessage().indexOf("InvalidRelationException")>-1)
        {
          JScrollPane jsp = new JScrollPane(new JLabel(err.getMessage()));
          jsp.setPreferredSize(new Dimension(200,100));
          JOptionPane.showMessageDialog(null, jsp, 
    				  "Error", JOptionPane.ERROR_MESSAGE);
        }
        return false;
      }
    } 
    catch(EntryInformationException e) 
    {
      final String error_string = e.getMessage();
      new MessageDialog(frame, "Cannot apply changes: " + error_string);

      return false;
    } 
    catch(OutOfRangeException e) 
    {
      new MessageDialog(frame, "Cannot apply changes - the location is out of " +
                        "range for this sequence");
      return false;
    } 
    catch(ReadOnlyException e) 
    {
      new MessageDialog(frame, "Cannot apply changes - the feature is " +
                        "read only");
      return false;
    } 
    finally
    {
      entry_group.getActionController ().endAction ();
    }

    dribble();

    return true;
  }

  /**
   * Propagate any changes to GFF ID qualifiers through the 
   * gene model
   * @param qualifiers
   */
  private void updateGffIds(QualifierVector qualifiers)
  {
    if( !GeneUtils.isDatabaseEntry(getFeature().getEmblFeature()) &&
        getFeature().getEmblFeature() instanceof GFFStreamFeature )
    {
      final GFFStreamFeature gffFeature = (GFFStreamFeature)getFeature().getEmblFeature();
      if(gffFeature.getChadoGene() != null)
      {
        try
        {
          final String newName = ((String) (qualifiers.getQualifierByName("ID").getValues().get(0))).trim();
          final String oldName = ((String) (gffFeature.getQualifierByName("ID").getValues().get(0))).trim();
         
          if(!newName.equals(oldName))
          {
            int val = JOptionPane.showConfirmDialog(null, 
                "Change name of children based on this ["+newName+"]?", 
                "Name Change", JOptionPane.OK_CANCEL_OPTION);
            if(val == JOptionPane.CANCEL_OPTION)
              return;

            final Set<uk.ac.sanger.artemis.io.Feature> children = 
              gffFeature.getChadoGene().getChildren(gffFeature);
            GeneUtils.propagateId(gffFeature, newName, children);
            GeneUtils.fixParentQualifier(oldName, newName, children);
            
            final Iterator<uk.ac.sanger.artemis.io.Feature> it = children.iterator();
            while(it.hasNext())
            {
              final GFFStreamFeature child = (GFFStreamFeature)it.next();
              if( child.getSegmentRangeStore().size() == 1 && 
                 !child.getKey().getKeyString().equals("CDS") &&
                  child.getKey().getKeyString().indexOf("exon") == -1)
                child.setSegmentRangeStore(null);
            }
          }
        }
        catch(Exception e){ }
      }
    }
  }
  
  /**
   *  Return the Feature we are editing as passed to the constructor.
   **/
  public Feature getFeature() 
  {
    return edit_feature;
  }

  /**
   *  On Unix machines this method will append the text of the feature to a
   *  file in a current directory called .dribble + <the entry name>
   **/
  private void dribble()
  {
    if(!Options.isUnixHost()) 
      return;

    final String dribble_file_name;

    if(getEntry().getName() != null) 
      dribble_file_name = ".dribble." + getEntry().getName();
    else 
      dribble_file_name = ".dribble.no_name";

    try  
    {
      final Writer writer = new FileWriter(dribble_file_name, true);
      getFeature().writeNative(writer);
      writer.flush();
      writer.close();
    } 
    catch(IOException e) 
    {
      System.err.println("IO exception while accessing " + dribble_file_name +
                         ": " + e.getMessage());
    }
  }

  /**
   *  Return the Entry that contains the Feature this object is displaying.
   **/
  private Entry getEntry()
  {
    return edit_entry;
  }

  /**
   *  Read the key from the feature and update the key chooser.
   **/
  private void updateKey() 
  {
    final Key feature_key = getFeature().getKey();
    key_choice.setKey(feature_key);
  }

  /**
   *  Return the EntryInformation object of the entry containing the feature.
   **/
  public EntryInformation getEntryInformation() 
  {
    if(entry_information == null)
      entry_information = getEntry().getEntryInformation();
    return entry_information;
  }

  /**
   *  Return the Selection that was passed to the constructor.
   **/
  private Selection getSelection()
  {
    return selection;
  }

  public static boolean isTabbedView()
  {
    return isTabbedView;
  }

  public static void setTabbedView(boolean isTabbedView)
  {
    FeatureEdit.isTabbedView = isTabbedView;
  }
  
  /**
   * Set whether the feature is obsolete (database mode).
   * @param obsoleteChanged
   */
  public void setObsoleteChanged(boolean obsoleteChanged)
  {
    propertiesPanel.setObsoleteChanged(obsoleteChanged);
  }

  public QualifierTextArea getQualifierTextArea()
  {
    return qualifier_text_area;
  }

}
