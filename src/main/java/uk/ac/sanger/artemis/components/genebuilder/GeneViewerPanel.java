/* GeneViewerPanel
 *
 * created: 2006
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2006  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/genebuilder/GeneViewerPanel.java,v 1.85 2009-06-18 14:59:05 tjc Exp $
 */

package uk.ac.sanger.artemis.components.genebuilder;

import javax.swing.*;

import java.awt.*;
import java.awt.event.*;
import java.util.List;
import java.util.Hashtable;
import java.util.Set;
import java.util.Iterator;
import java.util.Vector;

import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.FeatureSegment;
import uk.ac.sanger.artemis.FeatureSegmentVector;
import uk.ac.sanger.artemis.LastSegmentException;
import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.Selection;
import uk.ac.sanger.artemis.FeatureVector;

import uk.ac.sanger.artemis.components.SegmentBorder;
import uk.ac.sanger.artemis.io.DatabaseInferredFeature;
import uk.ac.sanger.artemis.io.EntryInformationException;
import uk.ac.sanger.artemis.io.Feature;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.ChadoCanonicalGene;
import uk.ac.sanger.artemis.io.InvalidRelationException;
import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.io.Location;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.io.RangeVector;
import uk.ac.sanger.artemis.sequence.Marker;
import uk.ac.sanger.artemis.sequence.MarkerRange;
import uk.ac.sanger.artemis.sequence.Strand;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.util.ReadOnlyException;

public class GeneViewerPanel extends MapPanel
{
  
  /** */
  private static final long serialVersionUID = 1L;

  /** pop up menu */
  private JPopupMenu popup;

  private float fraction;
  
  private int start;
  
  private MarkerRange click_range = null;
  
  /**
   *  The last(FeatureSegment) Marker that the user clicked on.  This is used
   *  for dragging the ends of segments.
   **/
  private Marker click_segment_marker = null;

  /**
   *  This is true if click_segment_marker is the Marker at the start of
   *  segment false otherwise.  The value is only useful if
   *  click_segment_marker is set.
   **/
  private boolean click_segment_marker_is_start_marker = false;

  /**
   *  When a FeatureSegment Marker drag starts, this is set to the Marker at
   *  the other end of the segment.  This is used to check that the drag has
   *  not move the Marker too far(past the end of the segment).
   **/
  private Marker other_end_of_segment_marker = null;
  
  /** This is a record of the feature last clicked */
  private uk.ac.sanger.artemis.Feature clicked_feature;
  
  private Point last_cursor_position;
  
  /**
   *  The shortcut for Delete Selected Features.
   **/
  final static KeyStroke DELETE_FEATURES_KEY =
    KeyStroke.getKeyStroke(KeyEvent.VK_DELETE,
                           Toolkit.getDefaultToolkit().getMenuShortcutKeyMask());
  final static public int DELETE_FEATURES_KEY_CODE = KeyEvent.VK_DELETE;
  
  final static KeyStroke CREATE_FEATURES_KEY =
    KeyStroke.getKeyStroke (KeyEvent.VK_C,
                            Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()); 

  final static public int CREATE_FEATURES_KEY_CODE = KeyEvent.VK_C;
  
  private GeneBuilderFrame gene_builder;
  
  public GeneViewerPanel(final ChadoCanonicalGene chado_gene,
                         final Selection selection,
                         final EntryGroup entry_group,
                         final GeneBuilderFrame gene_builder,
                         final JLabel status_line)
  {
    this.chado_gene   = chado_gene;
    this.selection    = selection;
    this.gene_builder = gene_builder;
    
    Dimension dim = new Dimension(400,350);
    setPreferredSize(dim);
    setBackground(Color.white);
    
//  Popup menu
    addMouseListener(new PopupListener());
    popup = new JPopupMenu();
    createMenus(popup, entry_group);
    
    // Listen for mouse motion events so that we can select ranges of bases.
    addMouseMotionListener(new MouseMotionAdapter()
    {
      public void mouseDragged(MouseEvent event)
      {
        if(event.isPopupTrigger() || entry_group == null ||
           event.getButton() == MouseEvent.BUTTON3)
          return;

        int select_start = (int) ((event.getX() - border) / fraction) + start;
        
        final Strand strand = ((uk.ac.sanger.artemis.Feature) (chado_gene.getGene()
            .getUserData())).getStrand();

        if(!strand.isForwardStrand())
          select_start = strand.getBases().getComplementPosition(select_start);

        if(click_segment_marker != null)
        {
          final int new_position;
          MarkerRange drag_range;
          try
          {
            drag_range = new MarkerRange(strand, select_start, select_start + 1);
          }
          catch(OutOfRangeException e1)
          {
            e1.printStackTrace();
            return;
          }

          if(click_segment_marker_is_start_marker)
          {
            // the Marker is at the start of the segment
            new_position = drag_range.getStart().getPosition();

            // don't go past the other end of the segment
            if(new_position > other_end_of_segment_marker.getPosition())
              return;
          }
          else
          {
            new_position = drag_range.getEnd().getPosition();

            // don't go past the other end of the segment
            if(new_position < other_end_of_segment_marker.getPosition())
              return;
          }

          try
          {
            click_segment_marker.setPosition(new_position);
            gene_builder.setActiveFeature(clicked_feature, false);
            return;
          }
          catch(OutOfRangeException e)
          {
            throw new Error("internal error - unexpected OutOfRangeException");
          }
        }

        final MarkerRange selected_range = selection.getMarkerRange();

        try
        {
          MarkerRange drag_range = new MarkerRange(strand, select_start,
              select_start + 1);

          // final MarkerRange new_marker_range;
          if(selected_range == null || click_range == null)
          {
            click_range = drag_range;
            status_line.setText("");
          }
          else
          {
            click_range = selected_range.combineRanges(drag_range, true);
            status_line.setText(selected_range.getRawRange().getStart() + ".."
                + selected_range.getRawRange().getEnd());
          }

          last_cursor_position = event.getPoint();
          selection.setMarkerRange(click_range);

          repaint();
        }
        catch(OutOfRangeException e)
        {
          e.printStackTrace();
        }
      }
    });
    
  }

  /**
   * Create menu for gene editor
   * 
   * @param menu
   * @param entry_group
   */
  protected void createMenus(JComponent menu,
                             final EntryGroup entry_group)
  {
    JMenuItem deleteMenu = new JMenuItem("Delete Selected Features");
    deleteMenu.setAccelerator(DELETE_FEATURES_KEY);
    deleteMenu.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)  
      {
        FeatureVector features = selection.getAllFeatures();
        boolean delete = deleteFeatures(features, chado_gene);
        
        try
        {
          if(delete)
          {
            for(int i = 0; i < features.size(); i++)
              GeneUtils.deleteAllFeature(features.elementAt(i), chado_gene);
            repaint();
          }
          else
            gene_builder.getFeatureEdit().setObsoleteChanged(true);
        }
        catch(ReadOnlyException e)
        {
          JOptionPane.showMessageDialog(null, 
              e.getMessage(), "Read Only", 
              JOptionPane.WARNING_MESSAGE);
        } 
      }
    });
    menu.add(deleteMenu);
    
    
    final JMenuItem deleteSegmentMenu = new JMenuItem("Delete Selected Exon");
    deleteSegmentMenu.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)  
      {
        final FeatureSegmentVector features = selection.getAllSegments();

        if(features == null || features.size() < 1)
        {
          JOptionPane.showMessageDialog(null, 
              "Select an exon and try again.", 
              "Exon to delete not found!",
              JOptionPane.ERROR_MESSAGE);
          return;
        }
        else
        {
          for(int i=0; i<features.size(); i++)
          {
            final uk.ac.sanger.artemis.Feature f = features.elementAt(i).getFeature();
            if(!f.getKey().getKeyString().equals(DatabaseDocument.EXONMODEL) &&
               !f.getKey().getKeyString().equals("pseudogenic_exon") )
            {
              JOptionPane.showMessageDialog(null, 
                  "Other feature types are selected.\nSelect an exon and try again.", 
                  "Select exon to delete!",
                  JOptionPane.ERROR_MESSAGE);
              return; 
            }
          }
        }
        
        uk.ac.sanger.artemis.FeatureSegment segment = null;
        int option = JOptionPane.showConfirmDialog(null, 
            "Delete Selected Exon", 
            "Delete Selected Exon", 
            JOptionPane.OK_CANCEL_OPTION);
        
        if(option == JOptionPane.CANCEL_OPTION)
          return;
        try
        {
          uk.ac.sanger.artemis.Feature feature = selection.getAllFeatures().elementAt(0);
          for(int i = 0; i < features.size(); i++)
          {
            segment = features.elementAt(i);
            segment.removeFromFeature();
            selection.remove(segment);
          }
          
          //selection.add(feature);
          gene_builder.setActiveFeature(feature, false);
          repaint();
        }
        catch(ReadOnlyException e)
        {
          e.printStackTrace();
        }
        catch(LastSegmentException e)
        {
          try
          {
            GeneUtils.deleteAllFeature(segment.getFeature(), chado_gene);
          }
          catch(ReadOnlyException e1)
          {
            // TODO Auto-generated catch block
            e1.printStackTrace();
          }
        }
      }
    });
    menu.add(deleteSegmentMenu);
    
    menu.add(new JSeparator());
    final JMenuItem createTranscript = new JMenuItem("Create transcript");
    createTranscript.setAccelerator(CREATE_FEATURES_KEY);
    createTranscript.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)  
      {
        createTranscript(chado_gene, entry_group);
        repaint();
      }
    });
    menu.add(createTranscript);
    
    
    final JMenuItem duplicateTranscript = new JMenuItem("Duplicate selected transcript");
    duplicateTranscript.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)  
      {
        duplicateTranscript(entry_group);
        repaint();
      }
    });
    menu.add(duplicateTranscript);
    
    
    final JMenu createFeatureMenu = new JMenu("Add to transcript in selected range");
    
    final JMenuItem createExon = new JMenuItem("exon");
    createExon.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)  
      {
        if(last_cursor_position == null)
          return;
        Feature transcript = getTranscriptAt(last_cursor_position);
        
        if(transcript == null)
        {
          JOptionPane.showMessageDialog(null, 
              "Select a single transcript to add an exon to and try again.", 
              "Transcript Not Found",
              JOptionPane.ERROR_MESSAGE);
          return;
        }
        
        String uniquename  = getQualifier(transcript, "ID");
        
        final List exons;
        final Key exonKey;
        if(chado_gene.getGene().getKey().getKeyString().equals("pseudogene"))
        {
          exons = chado_gene.getSpliceSitesOfTranscript(uniquename, "pseudogenic_exon");
          exonKey = new Key("pseudogenic_exon");
        }
        else 
        {
          exons = chado_gene.getSpliceSitesOfTranscript(uniquename, DatabaseDocument.EXONMODEL);
          exonKey = new Key(DatabaseDocument.EXONMODEL);
          /*
          exons = chado_gene.getSpliceSitesOfTranscript(uniquename, "exon");
          exonKey = new Key("exon");
          */
        }
        
        GFFStreamFeature embl_exon = null;
        if(exons != null && exons.size() > 0)
          embl_exon = (GFFStreamFeature)exons.get(0);
        
        Range range_selected = selection.getSelectionRange();
    
        addExonFeature(chado_gene, entry_group, embl_exon, 
                       range_selected, uniquename, selection, exonKey, gene_builder);
      }
    });
    createFeatureMenu.add(createExon);
    
    
    
    JMenuItem createFeature5Utr = new JMenuItem("5'UTR");
    createFeature5Utr.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)  
      {
        addDnaFeature(last_cursor_position, selection, 
                      entry_group, new Key("five_prime_UTR"));
      }
    });
    
    
    JMenuItem createFeature3Utr = new JMenuItem("3'UTR");
    createFeature3Utr.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)  
      {
        addDnaFeature(last_cursor_position, selection, 
                      entry_group, new Key("three_prime_UTR"));
      }
    });
    
    
    JMenuItem createFeatureDna = new JMenuItem("region");
    createFeatureDna.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)  
      {
        addDnaFeature(last_cursor_position, selection, 
                      entry_group, new Key("region"));
      }
    });
    
    
    
    JMenuItem createFeatureProtein = new JMenuItem("polypeptide");
    createFeatureProtein.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)  
      {
        if(selection.getAllFeatures().size() < 1)
        {
          JOptionPane.showMessageDialog(null, 
              "Select a single transcript and try again.", 
              "Transcript Selection",
              JOptionPane.ERROR_MESSAGE);
          return;
        }
        uk.ac.sanger.artemis.Feature transcript = selection.getAllFeatures().elementAt(0);
        final String transcriptName;
        try
        {
          transcriptName = (String)transcript.getQualifierByName("ID").getValues().get(0);
        }
        catch(InvalidRelationException e1)
        {
          e1.printStackTrace();
          return;
        }
        
        if(!chado_gene.isTranscript(transcriptName))
        {
          JOptionPane.showMessageDialog(null, 
              "Select a single transcript and try again.", 
              "Transcript Selection",
              JOptionPane.ERROR_MESSAGE);
          return;
        }
          
        addProteinFeature(chado_gene, entry_group, transcriptName, transcript);
      }
    });
    createFeatureMenu.add(createFeatureProtein);
    
    if(!DatabaseDocument.CHADO_INFER_CDS)
    {
      createFeatureMenu.add(createFeature5Utr);
      createFeatureMenu.add(createFeature3Utr);
    }
    createFeatureMenu.add(createFeatureDna);
    
    menu.add(createFeatureMenu);
    
    menu.add(new JSeparator());
    JMenuItem adjustCoords = new JMenuItem("Fix selected transcript coordinates");
    adjustCoords.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)  
      {
        FeatureVector features = selection.getAllFeatures();
        if(features.size() != 1)
        {
          JOptionPane.showMessageDialog(null, 
              "Select a single transcript and try again.", 
              "Transcript Selection",
              JOptionPane.ERROR_MESSAGE);
          return;
        }
 
        final uk.ac.sanger.artemis.Feature transcript = features.elementAt(0);
        GeneUtils.checkTranscriptBoundary(transcript, chado_gene);
        gene_builder.dispose(true);
      }
    });
    menu.add(adjustCoords);
    
    
    JMenuItem adjustGeneCoords = new JMenuItem("Fix gene model coordinates");
    adjustGeneCoords.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)  
      {
        GeneUtils.checkGeneBoundary(chado_gene);
        gene_builder.dispose(true);
      }
    });
    menu.add(adjustGeneCoords);
    
    menu.add(new JSeparator());
    
    final JMenuItem convertPsuedoMenu = new JMenuItem("Convert gene to/from pseudogene");
    convertPsuedoMenu.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)  
      {
        try
        {
          GeneUtils.convertPseudo(chado_gene);
          gene_builder.dispose(true);
        }
        catch(ReadOnlyException e)
        {
          e.printStackTrace();
        }
        catch(EntryInformationException e)
        {
          e.printStackTrace();
        }
        catch(OutOfRangeException e)
        {
          e.printStackTrace();
        }
      }
    });
    menu.add(convertPsuedoMenu);
  }
  
  /**
   * Delete or make features obsolete
   */
  protected static Boolean deleteFeatures(FeatureVector features,
                                          ChadoCanonicalGene chado_gene)
  {   
    final String feature_count_string;
    if (features.size () < 2)
      feature_count_string = "the selected feature";
    else
      feature_count_string = features.size () + " features";

    Box boption = Box.createVerticalBox();

    JCheckBox delete = new JCheckBox("permanently delete", 
        !Options.getOptions().getPropertyTruthValue("set_obsolete_on_delete"));
    boption.add(new JLabel("Make "+feature_count_string+" obsolete?"));
    boption.add(delete);
    int option = JOptionPane.showConfirmDialog(null,
                  boption, "Make obsolete", 
                  JOptionPane.OK_CANCEL_OPTION,
                  JOptionPane.QUESTION_MESSAGE);
    if(option == JOptionPane.CANCEL_OPTION)
      return null;
    return delete.isSelected();
  }
  
  public static uk.ac.sanger.artemis.Feature 
                     createTranscript(final ChadoCanonicalGene chadoGene,
                                      final EntryGroup entry_group)
  {
    return createTranscript(chadoGene, entry_group, chadoGene.getGene().getLocation());
  }
  
  public static uk.ac.sanger.artemis.Feature 
                     createTranscript(final ChadoCanonicalGene chadoGene,
                                      final EntryGroup entry_group,
                                      final Location location)
  {
    try
    {
      String gene_name = 
        (String)chadoGene.getGene().getQualifierByName("ID").getValues().get(0);
      
      String ID = chadoGene.autoGenerateTanscriptName(DatabaseDocument.TRANSCRIPT);
      
      QualifierVector qualifiers = new QualifierVector();
      qualifiers.add(new Qualifier("Parent", gene_name));
      
      if(ID != null)
        qualifiers.add(new Qualifier("ID", ID));
      
      final Key transcriptKey;
      if(chadoGene.getGene().getKey().getKeyString().equals("pseudogene"))
        transcriptKey = new Key("pseudogenic_transcript");
      else
        transcriptKey = new Key(DatabaseDocument.TRANSCRIPT);
      
      uk.ac.sanger.artemis.Feature feature = createFeature(
                    location, entry_group, transcriptKey,
                    qualifiers);
      
      ((GFFStreamFeature)(feature.getEmblFeature())).setChadoGene(chadoGene);
      chadoGene.addTranscript(feature.getEmblFeature());
      return feature;
    }
    catch(InvalidRelationException e)
    {
      e.printStackTrace();
    }
    return null;
  }
  
  
  /**
   * 
   */
  public void paintComponent(Graphics g)
  {
    super.paintComponent(g);
    Graphics2D g2d = (Graphics2D)g;
    Feature embl_gene = (Feature)chado_gene.getGene();
    uk.ac.sanger.artemis.Feature gene =
      (uk.ac.sanger.artemis.Feature)embl_gene.getUserData();
    
    setFont(uk.ac.sanger.artemis.Options.getOptions().getFont());

    start   = embl_gene.getFirstBase();
    int end = embl_gene.getLastBase();
    g2d.setColor( gene.getColour() );
    int ypos = border;
    
    // draw gene
    g2d.drawString(gene.getIDString()+
        ( GeneUtils.isObsolete((GFFStreamFeature)embl_gene) ? " (obsolete)" : "") , border, ypos);
    drawFeature(g2d, border, 
                getSize().width - border, 
                ypos, gene.getColour(), 1, 
                selection.contains(gene), 2.f,
                getFontHeight());
    
    List transcripts = chado_gene.getTranscripts();   
    fraction = (float)(getSize().width - (2*border))/
               (float)(end-start);
    
    ypos += border*2;
    
    for(int i=0; i<transcripts.size(); i++)
    {
      drawTranscriptOnLine(g2d, (Feature)transcripts.get(i), 
                           start, end, ypos, 
                           fraction, selection, chado_gene,
                           getFontHeight());
      
      if(i != transcripts.size()-1)
        ypos += getTranscriptSize();
    }
    
    // draw mouse drag selection
    if(selection.getMarkerRange() != null &&
       last_cursor_position != null)
    {
      Range range = selection.getSelectionRange();
      
      int ntranscript = (last_cursor_position.y - (border*3))/getTranscriptSize();
      if(ntranscript < transcripts.size())
      {
        int select_start = border+(int)((range.getStart()-start)*fraction);
        int select_end   = border+(int)((range.getEnd()-start)*fraction);
        ypos = (border*5)+(getTranscriptSize()*ntranscript);
        drawFeature(g2d, select_start, select_end, 
                    ypos, Color.YELLOW, 1.5f,
                    false, 2.f, getFontHeight());
      }
    }
    setPreferredSize(new Dimension(getSize().width, ypos+border));
  }
  
  /**
   * Duplicate the selected transcript and children
   * @param entry_group
   */
  private void duplicateTranscript(final EntryGroup entry_group)
  {
    FeatureVector features = selection.getAllFeatures();
    Feature transcript = null;
    
    if(features.size() == 1)
      transcript = features.elementAt(0).getEmblFeature();
    
    if(transcript == null || 
       transcript.getKey().getKeyString().indexOf(DatabaseDocument.TRANSCRIPT) == -1)
    {
      JOptionPane.showMessageDialog(null, 
          "Select a single transcript and try again.", 
          "Transcript Selection",
          JOptionPane.ERROR_MESSAGE);
      return;
    }

    final uk.ac.sanger.artemis.Feature newTranscript =
      createTranscript(chado_gene, entry_group, transcript.getLocation());
    String newTranscriptName = getQualifier(newTranscript.getEmblFeature(), "ID");
    last_cursor_position = getPointFromTranscriptName(newTranscriptName);
    
    Set childFeatures = chado_gene.getChildren(transcript);
    Iterator it = childFeatures.iterator();
    while(it.hasNext())
    {
      Feature f = (Feature)it.next();
      if(f.getKey().getKeyString().equals(DatabaseDocument.EXONMODEL) ||
         f.getKey().getKeyString().equals("pseudogenic_exon") )
      {
        GFFStreamFeature gff_exon = null;
        RangeVector ranges = f.getLocation().getRanges();
        
        for(int i=0; i<ranges.size(); i++)
          gff_exon = addExonFeature(chado_gene, entry_group, gff_exon, 
              (Range)ranges.get(i), newTranscriptName, selection, f.getKey(), gene_builder);
      }
      else if(!f.getKey().equals("polypeptide"))
      {
        selection.clear();
        selection.add((uk.ac.sanger.artemis.Feature) f.getUserData());
        addDnaFeature(last_cursor_position, selection, entry_group, f.getKey());
      }
    }

    addProteinFeature(chado_gene, entry_group, newTranscriptName, 
                      newTranscript); 
  }
  
  /**
   * Create and return a new Artemis feature.
   * @param new_location  location of the feature
   * @param entry_group   group to create feature in
   * @param key           the new features key
   * @param qualifiers    the new features qualifiers
   * @return
   */
  private static uk.ac.sanger.artemis.Feature createFeature(
                             final Location new_location,
                             final EntryGroup entry_group,
                             final Key key,
                             final QualifierVector qualifiers)
  {
    final Entry default_entry = entry_group.getDefaultEntry();
    uk.ac.sanger.artemis.Feature new_feature = null;

    try 
    {          
      new_feature = default_entry.createFeature(key, 
                                      new_location, qualifiers);
      return new_feature;
    } 
    catch (EntryInformationException e) 
    {
      e.printStackTrace();
    }
    catch(ReadOnlyException e)
    {
      e.printStackTrace();
    }
    catch(OutOfRangeException e)
    {
      e.printStackTrace();
    }
    return new_feature;
  }
 
  
  /**
   * Return the closest transcript feature from a given point on the
   * panel.
   * @param p 
   * @return
   */
  private Feature getTranscriptAt(Point p)
  {
    List transcripts = chado_gene.getTranscripts();
    
    int ntranscript = (p.y - (border*3))/getTranscriptSize();
    if(ntranscript < transcripts.size())
      return (Feature)transcripts.get(ntranscript);
    
    return (Feature)transcripts.get(transcripts.size()-1);
  }
  
  /**
   * Given a point on the panel find the feature drawn at that
   * location.
   * @param p
   * @return
   */
  private Object getFeatureAt(Point p)
  {
    if(p.y <= border+getFontHeight())
      return chado_gene.getGene();
    
    List transcripts = chado_gene.getTranscripts();
    
    for(int i=0; i<transcripts.size(); i++)
    {
      if(p.y >= (border*3)+(getTranscriptSize()*i) &&
         p.y <= (border*3)+(getTranscriptSize()*i)+getFontHeight())
      {
        return (Feature)transcripts.get(i);
      }
      else if(p.y >= (border*3)+(getTranscriptSize()*i)+getFontHeight() &&
              p.y <= (border*3)+(getTranscriptSize()*(i+1)))
      {
        Feature feature = (Feature)transcripts.get(i);
        String transcript_name = getQualifier(feature, "ID");
        List splicedFeatures = chado_gene.getSplicedFeaturesOfTranscript(transcript_name);
        
        if(splicedFeatures != null)
        {
          for(int j = 0; j < splicedFeatures.size(); j++)
          {
            FeatureSegmentVector segments = ((uk.ac.sanger.artemis.Feature) ((Feature) 
                splicedFeatures.get(j)).getUserData()).getSegments();

            if(segments.size() == 0)
              return (Feature)splicedFeatures.get(j);

            for(int k = 0; k < segments.size(); k++)
            {
              FeatureSegment segment = segments.elementAt(k);
              Range segment_range = segment.getRawRange();

              int segment_start = border
                  + (int) ((segment_range.getStart() - start) * fraction);
              int segment_end = border
                  + (int) ((segment_range.getEnd() - start) * fraction);
              if(p.x >= segment_start && p.x <= segment_end)
                return segment;
            }
          }
        }
        
        final List utr_3 = chado_gene.get3UtrOfTranscript(transcript_name);
        final List utr_5 = chado_gene.get5UtrOfTranscript(transcript_name);
        final List utrs = new Vector();
        if(utr_3 != null)
          utrs.addAll(utr_3);
        if(utr_5 != null)
          utrs.addAll(utr_5);
        
        for(int j = 0; j < utrs.size(); j++)
        {
          Feature utr = (Feature) utrs.get(j);
          Range range = utr.getLocation().getTotalRange();
          int utr_start = border
              + (int) ((range.getStart() - start) * fraction);
          int utr_end = border + (int) ((range.getEnd() - start) * fraction);
          if(p.x >= utr_start && p.x <= utr_end)
            return utr;
        }
        
        
        // anything else
        List others = chado_gene.getOtherFeaturesOfTranscript(transcript_name);
        if(others != null)
        {
          for(int j=0; j<others.size(); j++)
          {
            Feature other = (Feature)others.get(j);
            Range range = other.getLocation().getTotalRange();
            int r_start = border
                + (int) ((range.getStart() - start) * fraction);
            int r_end = border
                + (int) ((range.getEnd() - start) * fraction);
            if(p.x >= r_start && p.x <= r_end)
              return other;
          }
        }
      }
    }
         
    return null;
  }
  
  /**
   * Get the Point for a given transcript name
   * @param transcriptName
   * @return
   */
  private Point getPointFromTranscriptName(final String transcriptName) 
  {
    final List transcripts = chado_gene.getTranscripts();
    
    for(int i=0; i<transcripts.size(); i++)
    {
      String name = GeneUtils.getUniqueName((Feature)transcripts.get(i));
      if( name.equals(transcriptName) )
        return new Point(1, (border*3)+(getTranscriptSize()*i));
    }
    
    return null;
  }

  
  /**
   * Draw the transcript and child features.
   * @param g2d
   * @param embl_transcript
   * @param start
   * @param end
   * @param ypos
   * @param fraction
   */
  protected static void drawTranscriptOnLine(Graphics2D g2d, Feature embl_transcript, 
                                    final int start, final int end, int ypos, 
                                    float fraction, Selection selection, 
                                    ChadoCanonicalGene chado_gene,
                                    int fontHeight)
  {
    BasicStroke stroke = new BasicStroke(48.f);
    g2d.setStroke(stroke);
    
    uk.ac.sanger.artemis.Feature transcript = 
       (uk.ac.sanger.artemis.Feature)embl_transcript.getUserData();

    int t_start = border+(int)((embl_transcript.getFirstBase()-start)*fraction);
    int t_end   = border+(int)((embl_transcript.getLastBase()-start)*fraction);

    g2d.setColor( transcript.getColour() );

    g2d.drawString(transcript.getIDString() +
        ( GeneUtils.isObsolete((GFFStreamFeature)embl_transcript) ? " (obsolete)" : ""), border, ypos);
    drawFeature(g2d, t_start, t_end, 
                ypos, transcript.getColour(), 1, 
                selection.contains(transcript), 2.f, fontHeight);

    //List exons = chado_gene.getSpliceSitesOfTranscript(
    //    getQualifier( embl_transcript, "ID" ), "exon");

    Set spliceSiteTypes = chado_gene.getSpliceTypes(
        getQualifier( embl_transcript, "ID" ));
    ypos += border*2;
    
    if(spliceSiteTypes != null)
    {
      Iterator it = spliceSiteTypes.iterator();
      
      while(it.hasNext())
      {
        final String type = (String)it.next();
        List splicedFeatures = chado_gene.getSpliceSitesOfTranscript(
            getQualifier( embl_transcript, "ID" ), type);
        
        boolean last_segment = false;

        // build from artemis objects
        for(int i = 0; i < splicedFeatures.size(); i++)
        {
          int last_ex_start = 0;
          int last_ex_end = 0;
          int last_ypos = 0;

          Feature embl_exon = (Feature) splicedFeatures.get(i);
          if(embl_exon instanceof DatabaseInferredFeature)
            continue;
          
          uk.ac.sanger.artemis.Feature exon = (uk.ac.sanger.artemis.Feature) embl_exon
              .getUserData();

          RangeVector ranges = exon.getLocation().getRanges();
          FeatureSegmentVector segments = null;

          try
          {
            segments = exon.getSegments();
          }
          catch(NullPointerException npe)
          {
          }

          float selected_size;
          for(int j = 0; j < ranges.size(); j++)
          {
            Range range = (Range) ranges.get(j);

            int ex_start = border
                + (int) ((range.getStart() - start) * fraction);
            int ex_end = border + (int) ((range.getEnd() - start) * fraction);

            if(exon.getColour() != null)
              g2d.setColor(exon.getColour());

            if(j == ranges.size() - 1)
              last_segment = true;

            selected_size = 2.f;
            Color borderColor = Color.BLACK;
            if(segments != null)
            {
              for(int k = 0; k < segments.size(); k++)
              {
                FeatureSegment segment = segments.elementAt(k);
                if(range.equals(segment.getRawRange())
                    && selection.contains(segment))
                {
                  selected_size = 4.f;
                  borderColor = SegmentBorder.HIGHLIGHT_BORDER_COLOUR;
                }
              }
            }

            drawExons(g2d, ex_start, ex_end, last_ex_start, last_ex_end,
                last_ypos, 0, ypos, exon.getColour(), borderColor, 1.5f, exon
                    .isForwardFeature(), last_segment,
                selection.contains(exon), selected_size, fontHeight);

            last_ex_end = ex_end;
            last_ex_start = ex_start;
            last_ypos = ypos;
          }
        }
      }
    }
    
    // draw utr's
    String transcript_id = getQualifier( embl_transcript, "ID" );
    List embl_utr = chado_gene.get3UtrOfTranscript(
        transcript_id);
    
    if(embl_utr != null)
      drawFeatureList(g2d, embl_utr, ypos, selection, fontHeight, start, fraction);
    
    embl_utr = chado_gene.get5UtrOfTranscript(
        transcript_id);
    
    if(embl_utr != null)
      drawFeatureList(g2d, embl_utr, ypos, selection, fontHeight, start, fraction);
    
    // draw other transcript child features
    List embl_other = chado_gene.getOtherFeaturesOfTranscript(transcript_id);
    if(embl_other != null)
      drawFeatureList(g2d, embl_other, ypos, selection, fontHeight, start, fraction);
  }
  
  /**
   * Method to draw UTR features
   * @param g2d
   * @param embl_utr
   * @param ypos
   */
  private static void drawFeatureList(final Graphics2D g2d,
                       final List feature_list,
                       final int ypos, 
                       Selection selection, 
                       int fontHeight,
                       int start,
                       float fraction)
  {
    for(int i=0; i<feature_list.size(); i++)
    {
      Feature embl_feature = (Feature)feature_list.get(i);
      uk.ac.sanger.artemis.Feature feature = 
        (uk.ac.sanger.artemis.Feature)embl_feature.getUserData();
      RangeVector ranges = embl_feature.getLocation().getRanges();
  
      if(feature == null)
        feature = new uk.ac.sanger.artemis.Feature(embl_feature);
      
      
      for(int j = 0; j < ranges.size(); j++)
      {
        Range range = (Range) ranges.get(j);

        int r_start = border + (int) ((range.getStart() - start) * fraction);
        int r_end   = border + (int) ((range.getEnd() - start) * fraction);

        drawFeature(g2d, r_start, r_end, ypos, feature.getColour(), 1.5f,
                    selection.contains(feature), 2.f, fontHeight);
      }
    }
  }
  
  /**
   * Draw frame lines
   * @param g2d
   * @param ypos
   * @param start
   * @param end
   * @param fraction
   */
  /*private void drawFrameLines(Graphics2D g2d, int ypos,
                              int start, int end, float fraction)
  {
    int offset;
    g2d.setStroke(new BasicStroke(getFontHeight()));
    for(int k=0; k<8; k++)
    {
      offset = (k * getFontHeight() * 2) + (getFontHeight()/2);
      if(k == 3 || k == 4)
        g2d.setColor( Color.LIGHT_GRAY );
      else
        g2d.setColor(light_grey);
      g2d.drawLine(border, ypos+offset, 
                   (int)((end-start)*fraction)+border, ypos+offset);
    }
  }*/
  
  /**
   * Draw exon features
   * @param g2d
   * @param ex_start
   * @param ex_end
   * @param last_ex_start
   * @param last_ex_end
   * @param last_ypos
   * @param offset
   * @param ypos
   * @param exon_colour
   * @param size
   * @param isForward
   * @param last_segment
   * @param selected
   * @param selected_size
   */
  private static void drawExons(Graphics2D g2d, int ex_start, int ex_end, 
                         int last_ex_start, int last_ex_end, int last_ypos,
                         int offset, int ypos, Color exon_colour,
                         Color borderColor, float size,
                         boolean isForward, boolean last_segment,
                         boolean selected, float selected_size, int fontHeight)
  {   
    drawFeature(g2d, ex_start, ex_end, ypos+offset, exon_colour, borderColor,
        size, selected, selected_size, fontHeight);
    
    // draw connections
    if(last_ex_end != 0 ||
       last_ex_start != 0)
    {
      BasicStroke stroke = new BasicStroke(1.f);
      g2d.setStroke(stroke);
      
      int ymid;
      if(last_ypos < ypos+offset)
        ymid = last_ypos;
      else
        ymid = ypos+offset; 

      if(isForward)
      {
        g2d.drawLine(last_ex_end, last_ypos, 
                     last_ex_end+((ex_start-last_ex_end)/2), ymid-fontHeight/2);
        g2d.drawLine(last_ex_end+((ex_start-last_ex_end)/2), ymid-fontHeight/2, 
                     ex_start, ypos+offset); 
      }
      else
      { 
        g2d.drawLine(last_ex_start, last_ypos, 
                     last_ex_start+((ex_end-last_ex_start)/2), ymid-fontHeight/2);
        g2d.drawLine(last_ex_start+((ex_end-last_ex_start)/2), ymid-fontHeight/2, 
                     ex_end, ypos+offset); 
      }  
    }
  
    // draw arrow
    if(last_segment)
    {
      BasicStroke stroke = new BasicStroke(1.f);
      g2d.setStroke(stroke);
      
      if(isForward)
      {      

        g2d.drawLine(ex_end, ypos, 
                     ex_end+fontHeight/2, ypos+(int)((fontHeight*size)/2));
        g2d.drawLine(ex_end+fontHeight/2, ypos+(int)((fontHeight*size)/2),
                     ex_end, ypos+(int)((fontHeight*size)));
      }
      else
      {
        g2d.drawLine(ex_start, ypos, 
                     ex_start-fontHeight/2, ypos+(int)((fontHeight*size)/2));
        g2d.drawLine(ex_start-fontHeight/2, ypos+(int)((fontHeight*size)/2),
                     ex_start, ypos+(int)((fontHeight*size)));
      }
    }
  }
  

  
  /**
   * Get the <code>Color</code> for a feature from its colour attribute.
   * @param feature
   * @return
   */
  /*private Color getColorFromAttributes(org.gmod.schema.sequence.Feature feature)
  {
    List properties = new Vector(feature.getFeatureProps());
    for(int i=0; i<properties.size(); i++)
    {
      FeatureProp property = (FeatureProp)properties.get(i);
      
      if(property.getCvTerm().getName().equals("colour") ||
         property.getCvTerm().getName().equals("color") )
        return Options.getOptions().getColorFromColourNumber(Integer.parseInt(property.getValue()));
    }  
    return Color.CYAN;
  }*/
  
  /**
   * Get the frame id for a feature segment
   * @param chado_gene  the chado representation of the gene model
   * @param loc         feature location of the feature
   * @param nexon       number of the exon
   * @param exons       List of exons
   * @return frame id
   */
  /*private int getFrameID(ChadoCanonicalGene chado_gene, 
                         FeatureLoc loc, 
                         int nexon, List exons)
  {
    final int position_on_strand;
    
    if(loc.getStrand().shortValue() == -1)
      position_on_strand = chado_gene.getSeqlen()-loc.getFmax().intValue();
    else
      position_on_strand = loc.getFmin().intValue();
    
    // this will be 0, 1 or 2 depending on which frame the segment is in
    final int start_base_modulo =
      (position_on_strand + getFrameShift(nexon, exons, chado_gene, loc)) % 3;

    if(loc.getStrand().shortValue() == 1)
    {
      switch (start_base_modulo)
      {
      case 0:
        return FeatureSegment.FORWARD_FRAME_1;
      case 1:
        return FeatureSegment.FORWARD_FRAME_2;
      case 2:
        return FeatureSegment.FORWARD_FRAME_3;
      }
    } 
    else
    {
      switch (start_base_modulo)
      {
      case 0:
        return FeatureSegment.REVERSE_FRAME_1;
      case 1:
        return FeatureSegment.REVERSE_FRAME_2;
      case 2:
        return FeatureSegment.REVERSE_FRAME_3;
      }
    }

    return FeatureSegment.NO_FRAME;
  }*/
  
  
  /**
   *  Returns 0, 1 or 2 depending on which translation frame this segment is
   *  in.  A frame shift of zero means that the bases should be translated
   *  starting at the start position of this segment, 1 means start
   *  translating one base ahead of the start position and 2 means start
   *  translating two bases ahead of the start position.
   **/
  /*private int getFrameShift(int nexon, List exons, 
                            ChadoCanonicalGene chado_gene,
                            FeatureLoc loc) 
  {
    // find the number of bases in the segments before this one
    int base_count = 0;
    int direction  = 0;
    
    for(int i = 0; i < exons.size(); ++i) 
    {
      org.gmod.schema.sequence.Feature this_feature = 
        (org.gmod.schema.sequence.Feature)exons.get(i);
      FeatureLoc featureLoc = uk.ac.sanger.artemis.util.DatabaseDocument.getFeatureLoc(
          new Vector(this_feature.getFeatureLocsForFeatureId()), chado_gene.getSrcfeature_id());
      
      int this_direction;
      if(featureLoc.getStrand().shortValue() == 1)
        this_direction = 1;
      else
        this_direction = -1;

      if(i == nexon) 
      {
        if(i != 0 && this_direction != direction)
          base_count = 0;

        break;
      }
      else 
      {
        if(i == 0)
          direction = this_direction;
        else if(this_direction != direction)
          base_count = 0;

        base_count += featureLoc.getFmax().intValue()-featureLoc.getFmin().intValue();
      }
    }
    
    int codon_start = loc.getPhase().shortValue();
    int mod_value   = (base_count + 3 - codon_start) % 3;

    if(mod_value == 1) 
      return 2;
    else if(mod_value == 2)
      return 1;
    else 
      return 0;
  }*/

 
  protected static uk.ac.sanger.artemis.Feature addFeature(
                          final Range range,
                          final String transcriptName,
                          String featureName,
                          final boolean isComplement,
                          final boolean isParent,
                          final ChadoCanonicalGene chado_gene,
                          final EntryGroup entry_group,
                          final Key key, final String keyName) 
          throws OutOfRangeException
  {
    if(range == null)
    {
      JOptionPane.showMessageDialog(null, 
          "Select a range and try again.", 
          "Range Selection",
          JOptionPane.ERROR_MESSAGE);
      return null;
    }
    
    if(featureName == null)
    {
      final String autoId = 
        chado_gene.autoGenerateFeatureName(transcriptName, keyName);
      featureName = JOptionPane.showInputDialog(null,
          "Provide a feature unique name : ", autoId);
    
      if(featureName == null || featureName.equals(""))
        return null;
    }
    
    QualifierVector qualifiers = new QualifierVector();
    
    qualifiers.add(new Qualifier("ID", featureName.trim()));

    if(isParent)
    {
      //key = new Key("region");
      qualifiers.add(new Qualifier("Parent", transcriptName));
    }
    else
    {
      //key = new Key("polypeptide");
      qualifiers.add(new Qualifier("Derives_from", transcriptName));
    }
    
    //final Entry default_entry = entry_group.getDefaultEntry();
    //final Key default_key =
    //  default_entry.getEntryInformation().getDefaultKey();
    
    uk.ac.sanger.artemis.Feature newFeature = createFeature(
        new Location(new RangeVector(range), isComplement),
        entry_group, key,
        qualifiers);
    
    if(GeneUtils.isHiddenFeature(key.getKeyString()))
      ((GFFStreamFeature)(newFeature.getEmblFeature())).setVisible(false);
      
    ((GFFStreamFeature)(newFeature.getEmblFeature())).setChadoGene(chado_gene);
    
    if(isParent)
      chado_gene.addOtherFeatures(transcriptName, 
          newFeature.getEmblFeature());
    else
      chado_gene.addProtein(transcriptName, 
          newFeature.getEmblFeature());
    return newFeature;
  }
  
  
  private void addDnaFeature(final Point last_cursor_position,
                             final Selection selection,
                             final EntryGroup entry_group,
                             final Key key)
  {
    if(last_cursor_position == null)
      return;
    
    final String name;
    if(key.getKeyString().equals("five_prime_UTR"))
      name = "5UTR";
    else if(key.getKeyString().equals("three_prime_UTR"))
      name = "3UTR";
    else
      name = key.getKeyString();
    
    Feature transcript = getTranscriptAt(last_cursor_position);
    String transcriptName  = getQualifier(transcript, "ID");
    Range range_selected = selection.getSelectionRange();
    
    
    try
    {
      uk.ac.sanger.artemis.Feature newFeature = 
          addFeature(range_selected, transcriptName, null,
                     transcript.getLocation().isComplement(), 
                     true, chado_gene, entry_group, key, name);
      if(newFeature != null)
      {
        selection.clear();
        selection.add(newFeature);
        gene_builder.setActiveFeature(newFeature, false);
      }
    }
    catch(OutOfRangeException e)
    {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }  
  }
  
  public static uk.ac.sanger.artemis.Feature 
                     addProteinFeature(final ChadoCanonicalGene chadoGene,
                                       final EntryGroup entry_group,
                                       final String transcriptName,
                                       final uk.ac.sanger.artemis.Feature transcript)
  {
    final List exons;
    if(chadoGene.getGene().getKey().getKeyString().equals("pseudogene"))
      exons = chadoGene.getSpliceSitesOfTranscript(transcriptName, "pseudogenic_exon");
    else 
      exons = chadoGene.getSpliceSitesOfTranscript(transcriptName, DatabaseDocument.EXONMODEL);
   
    final Range range_selected;
    if(exons != null && exons.size() > 0)
      range_selected = ((GFFStreamFeature)exons.get(0)).getLocation().getTotalRange();
    else
      range_selected = transcript.getLocation().getTotalRange();
    
    final String pepName = chadoGene.autoGeneratePepName(transcriptName);
    
    try
    {
      return addFeature(range_selected, transcriptName, pepName,
          transcript.getLocation().isComplement(), false, chadoGene, entry_group,
          new Key("polypeptide"), "pep");
    }
    catch(OutOfRangeException e)
    {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }  
    return null;
  }
  
  /**
   * Create an exons feature
   * @param feature
   * @param range
   * @param transcript_name
   */
  public static GFFStreamFeature addExonFeature(
                              final ChadoCanonicalGene chadoGene,
                              final EntryGroup entry_group,
                              final GFFStreamFeature feature, Range range,
                              final String transcript_name,
                              final Selection selection,
                              final Key exonKey, 
                              final GeneBuilderFrame gene_builder)
  {
    try
    {
      if(feature == null)
      {
        QualifierVector qualifiers = new QualifierVector();
        qualifiers.add(new Qualifier("Parent", transcript_name));
      
        final String ID = chadoGene.autoGenerateSplicedFeatureName(transcript_name);

        if(ID != null)
          qualifiers.add(new Qualifier("ID", ID));
        
        uk.ac.sanger.artemis.Feature exon = createFeature(
            new Location(new RangeVector(range), 
            chadoGene.getGene().getLocation().isComplement()),
            entry_group, exonKey,
            qualifiers);
      
        GFFStreamFeature gff_exon = (GFFStreamFeature)exon.getEmblFeature();
        chadoGene.addSplicedFeatures(transcript_name, gff_exon);
        ((GFFStreamFeature)(exon.getEmblFeature())).setChadoGene(chadoGene);
        
        if(ID != null)
        {
          final Hashtable id_range_store = new Hashtable();
          id_range_store.put(ID, range);
          gff_exon.setSegmentRangeStore(id_range_store);
        }
        return gff_exon;
      }
      else
      {
        // add new ID
        final RangeVector ranges = new RangeVector();
        ranges.add(range);
        GeneUtils.addSegment(feature, ranges, transcript_name);
        final QualifierVector old_qualifiers = feature.getQualifiers().copy();
        ((uk.ac.sanger.artemis.Feature)feature.getUserData()).addSegment(range, old_qualifiers);
        
        if(gene_builder != null)
          gene_builder.setActiveFeature((uk.ac.sanger.artemis.Feature)feature.getUserData(), false);
      }
    }
    catch(InvalidRelationException e)
    {
      e.printStackTrace();
    }
    catch(ReadOnlyException e)
    {
      e.printStackTrace();
    }
    catch(EntryInformationException e)
    {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    return feature;
  }
  
  

  
  private static String getQualifier(Feature feature, String name) 
  {
    Qualifier qualifier = null;
    try
    {
      qualifier = feature.getQualifierByName(name);
    }
    catch(InvalidRelationException e)
    {
      e.printStackTrace();
    }
    if(qualifier == null)
      return null;

    return (String) (qualifier.getValues().get(0));
  }
  
  /**
   * Popup listener
   */
  class PopupListener extends MouseAdapter
  {
    public void mousePressed(MouseEvent e)
    {
      click_segment_marker = null;
      maybeShowPopup(e);
    }

    public void mouseReleased(MouseEvent e)
    {
      if(e.isPopupTrigger())
        maybeShowPopup(e);
      click_segment_marker = null;
    }

    private void maybeShowPopup(MouseEvent e)
    {
      if(e.isPopupTrigger())
      {
        popup.show(e.getComponent(),
                e.getX(), e.getY());
      }
      else if(e.getButton() != MouseEvent.BUTTON3 &&
              e.getClickCount() == 1)
      {
        if(selection.getMarkerRange() != null &&
           e.isShiftDown())
          return;
        
        click_range = null;
        
        if(!e.isShiftDown())
          selection.clear();
        Object feature = getFeatureAt(e.getPoint());
        if(feature == null)
          return;      
        
        if(e.isShiftDown())
        {
          if(feature instanceof Feature)
            selection.add(
                (uk.ac.sanger.artemis.Feature)((Feature)feature).getUserData());
          else
            selection.add((FeatureSegment)feature);
        }
        else
        {
          if(feature instanceof Feature)
          {
            selection.set(
                (uk.ac.sanger.artemis.Feature)((Feature)feature).getUserData());
            
            boolean isSet = true;
            if(((Feature)feature).isReadOnly())
              isSet = false;
            
            gene_builder.setActiveFeature(
                (uk.ac.sanger.artemis.Feature)((Feature)feature).getUserData(), isSet);
          }
          else
          {
            selection.set((FeatureSegment)feature);
            gene_builder.setActiveFeature(((FeatureSegment)feature).getFeature(), true);
          }
        }

        
        // the following is used by mouseDragged() to edit feature locations
        if(feature instanceof Feature)
          clicked_feature = (uk.ac.sanger.artemis.Feature)((Feature)feature).getUserData();
        else
          clicked_feature = ((FeatureSegment)feature).getFeature();
        click_segment_marker = null;
        
        if(Options.getOptions().canDirectEdit() &&
           !clicked_feature.isReadOnly() &&
           clicked_feature.getEntry().getBases() != null)
         {
           // search the feature to find if the start Marker of any of the
           // segments matches the start Marker of new_click_range or if an end
           // Marker matches the end Marker of new_click_range

          final FeatureSegmentVector segments = clicked_feature.getSegments();

          for(int k = 0; k < segments.size(); k++)
          {
            FeatureSegment this_segment = segments.elementAt(k);
            Range segment_range = this_segment.getRawRange();

            int segment_start = border
                + (int) ((segment_range.getStart() - start) * fraction);
            int segment_end = border
                + (int) ((segment_range.getEnd() - start) * fraction);
            
            // allow very ends of features to be dragable
            if( (e.getPoint().x >= segment_start && e.getPoint().x < segment_start+3 && e.getPoint().x <= segment_end) ||
                (e.getPoint().x <= segment_end   && e.getPoint().x > segment_end-3   && e.getPoint().x >= segment_start) )
            {
              int segment_mid = (segment_start+segment_end)/2;
              
              if((e.getPoint().x < segment_mid &&  this_segment.isForwardSegment()) ||
                 (e.getPoint().x > segment_mid && !this_segment.isForwardSegment()))
              {
                click_segment_marker = this_segment.getStart(); 
                click_segment_marker_is_start_marker = true;
                other_end_of_segment_marker = this_segment.getEnd();
              }
              else
              {
                click_segment_marker = this_segment.getEnd();
                click_segment_marker_is_start_marker = false;
                other_end_of_segment_marker = this_segment.getStart();
              }
            }
          }
          
        }
        repaint();
      }
    }
  }
 
}