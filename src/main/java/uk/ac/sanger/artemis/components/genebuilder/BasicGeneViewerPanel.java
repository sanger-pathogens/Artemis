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
import java.util.Vector;

import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.FeatureSegment;
import uk.ac.sanger.artemis.FeatureSegmentVector;
import uk.ac.sanger.artemis.LastSegmentException;
import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.Selection;
import uk.ac.sanger.artemis.FeatureVector;

import uk.ac.sanger.artemis.io.EntryInformationException;
import uk.ac.sanger.artemis.io.Feature;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.ChadoCanonicalGene;
import uk.ac.sanger.artemis.io.InvalidRelationException;
import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.sequence.Marker;
import uk.ac.sanger.artemis.sequence.MarkerRange;
import uk.ac.sanger.artemis.sequence.Strand;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.util.ReadOnlyException;

public class BasicGeneViewerPanel extends MapPanel
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
  
  private BasicGeneBuilderFrame gbFrame;
  
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
  
  final private BasicGeneBuilderFrame gene_builder;
  
  public BasicGeneViewerPanel(final BasicGeneBuilderFrame gene_builder,
                              final ChadoCanonicalGene chado_gene,
                              final Selection selection,
                              final EntryGroup entry_group,
                              final BasicGeneBuilderFrame gbFrame,
                              final JLabel status_line)
  {
    this.chado_gene   = chado_gene;
    this.selection    = selection;
    this.gbFrame      = gbFrame;
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
        
        final Strand strand = ((uk.ac.sanger.artemis.Feature) 
            (chado_gene.getGene().getUserData())).getStrand();

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
            gene_builder.setActiveFeature(clicked_feature);
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
        boolean delete = GeneViewerPanel.deleteFeatures(features, chado_gene);
        
        try
        {
          if(delete)
          {
            for(int i = 0; i < features.size(); i++)
              GeneUtils.deleteAllFeature(features.elementAt(i), chado_gene);
          }
          else
          {
            gbFrame.setObsoleteChanged(true, features);
          }
          gbFrame.dispose(true);
        }
        catch(NullPointerException npe)
        {
          // can't reopen
          gbFrame.dispose(false);
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
          for(int i = 0; i < features.size(); i++)
          {
            segment = features.elementAt(i);
            segment.removeFromFeature();
            selection.remove(segment);
          }
          
          //selection.add(feature);
          //gene_builder.setActiveFeature(feature);
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
      
    final JMenu createFeatureMenu = new JMenu("Add to transcript in selected range");
    
    final JMenuItem createExon = new JMenuItem("exon");
    createExon.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)  
      {
        if(last_cursor_position == null)
          return;
        Feature transcript = gbFrame.getSelectedTranscriptFeature().getEmblFeature();
        
        if(transcript == null)
        {
          JOptionPane.showMessageDialog(null, 
              "Select a single transcript to add an exon to and try again.", 
              "Transcript Not Found",
              JOptionPane.ERROR_MESSAGE);
          return;
        }
        
        String uniquename  = GeneUtils.getUniqueName(transcript);
        
        final List<Feature> exons;
        final Key exonKey;
        if(chado_gene.getGene().getKey().getKeyString().equals("pseudogene"))
        {
          exons = chado_gene.getSpliceSitesOfTranscript(uniquename, "pseudogenic_exon");
          exonKey = new Key("pseudogenic_exon");
        }
        else 
        {
          exons = chado_gene.getSpliceSitesOfTranscript(
              uniquename, DatabaseDocument.EXONMODEL);
          exonKey = new Key(DatabaseDocument.EXONMODEL);
        }
        
        GFFStreamFeature embl_exon = null;
        if(exons != null && exons.size() > 0)
          embl_exon = (GFFStreamFeature)exons.get(0);
        
        Range range_selected = selection.getSelectionRange();
        GeneViewerPanel.addExonFeature(chado_gene, entry_group, embl_exon, 
                       range_selected, uniquename, selection, exonKey, null);
        repaint();
      }
    });
    createFeatureMenu.add(createExon);
    
    
    
    JMenuItem createFeature5Utr = new JMenuItem("5'UTR");
    createFeature5Utr.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)  
      {
        // TODO
        //addDnaFeature(last_cursor_position, selection, 
        //              entry_group, new Key("five_prime_UTR"));
      }
    });
    
    
    JMenuItem createFeature3Utr = new JMenuItem("3'UTR");
    createFeature3Utr.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)  
      {
        // TODO
        // addDnaFeature(last_cursor_position, selection, 
        //              entry_group, new Key("three_prime_UTR"));
      }
    });
    
    
    JMenuItem createFeatureDna = new JMenuItem("region");
    createFeatureDna.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)  
      {
        // TODO
        //addDnaFeature(last_cursor_position, selection, 
        //              entry_group, new Key("region"));
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
          
     // TODO
        //addProteinFeature(chado_gene, entry_group, transcriptName, transcript);
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
     // TODO  gene_builder.dispose(true);
      }
    });
    menu.add(adjustCoords);
    
    
    JMenuItem adjustGeneCoords = new JMenuItem("Fix gene model coordinates");
    adjustGeneCoords.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)  
      {
        GeneUtils.checkGeneBoundary(chado_gene);
        
     // TODO
     //   gene_builder.dispose(true);
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
    GeneViewerPanel.drawFeature(g2d, border, 
                getSize().width - border, 
                ypos, gene.getColour(), 1, 
                selection.contains(gene), 2.f, getFontHeight());
    
    List<Feature> transcripts = chado_gene.getTranscripts();   
    fraction = (float)(getSize().width - (2*border))/
               (float)(end-start);
    
    ypos += border*2;
    
    try
    {
      Feature transcript = gbFrame.getSelectedTranscriptFeature().getEmblFeature();
      GeneViewerPanel.drawTranscriptOnLine(g2d, transcript, 
          start, end, ypos, 
          fraction, selection, chado_gene,
          getFontHeight());
    }
    catch(NullPointerException npe)
    {
      return;
    }

    ypos += getTranscriptSize();

    setPreferredSize(new Dimension(getSize().width, ypos+border));
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
  }

  
  /**
   * Given a point on the panel find the feature drawn at that
   * location.
   * @param p
   * @return
   */
  private Object getFeatureAt(Point p)
  {
    if (p.y <= border + getFontHeight())
      return chado_gene.getGene();

    Feature transcript = gbFrame.getSelectedTranscriptFeature().getEmblFeature();

    if (p.y >= (border * 3) && p.y <= (border * 3) + getFontHeight())
      return transcript;
    else if (p.y >= (border * 3) + getFontHeight()
        && p.y <= (border * 3) + (getTranscriptSize()))
    {
      String transcript_name = GeneUtils.getUniqueName(transcript);
      List<Feature> splicedFeatures = chado_gene
          .getSplicedFeaturesOfTranscript(transcript_name);

      if (splicedFeatures != null)
      {
        for (int j = 0; j < splicedFeatures.size(); j++)
        {
          FeatureSegmentVector segments = ((uk.ac.sanger.artemis.Feature) (splicedFeatures
              .get(j)).getUserData()).getSegments();

          if (segments.size() == 0)
            return (Feature) splicedFeatures.get(j);

          for (int k = 0; k < segments.size(); k++)
          {
            FeatureSegment segment = segments.elementAt(k);
            Range segment_range = segment.getRawRange();

            int segment_start = border
                + (int) ((segment_range.getStart() - start) * fraction);
            int segment_end = border
                + (int) ((segment_range.getEnd() - start) * fraction);
            if (p.x >= segment_start && p.x <= segment_end)
              return segment;
          }
        }
      }

      final List<Feature> utr_3 = chado_gene
          .get3UtrOfTranscript(transcript_name);
      final List<Feature> utr_5 = chado_gene
          .get5UtrOfTranscript(transcript_name);
      final List<Feature> utrs = new Vector<Feature>();
      if (utr_3 != null)
        utrs.addAll(utr_3);
      if (utr_5 != null)
        utrs.addAll(utr_5);

      for (int j = 0; j < utrs.size(); j++)
      {
        Feature utr = (Feature) utrs.get(j);
        Range range = utr.getLocation().getTotalRange();
        int utr_start = border + (int) ((range.getStart() - start) * fraction);
        int utr_end = border + (int) ((range.getEnd() - start) * fraction);
        if (p.x >= utr_start && p.x <= utr_end)
          return utr;
      }

      // anything else
      List<Feature> others = chado_gene
          .getOtherFeaturesOfTranscript(transcript_name);
      if (others != null)
      {
        for (int j = 0; j < others.size(); j++)
        {
          Feature other = others.get(j);
          Range range = other.getLocation().getTotalRange();
          int r_start = border + (int) ((range.getStart() - start) * fraction);
          int r_end = border + (int) ((range.getEnd() - start) * fraction);
          if (p.x >= r_start && p.x <= r_end)
            return other;
        }
      }
    }
 
    return null;
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
            selection.set(
                (uk.ac.sanger.artemis.Feature)((Feature)feature).getUserData());
          else
            selection.set((FeatureSegment)feature);
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
        
        gene_builder.setActiveFeature(clicked_feature);
        repaint();
      }
    }
  }
 
}