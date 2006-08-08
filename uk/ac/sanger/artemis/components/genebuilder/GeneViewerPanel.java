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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/genebuilder/GeneViewerPanel.java,v 1.20 2006-08-08 10:20:23 tjc Exp $
 */

package uk.ac.sanger.artemis.components.genebuilder;

import javax.swing.*;

import java.awt.*;
import java.awt.event.*;
import java.util.List;
import java.util.Collections;
import java.util.Hashtable;
import java.util.Vector;
import java.awt.geom.RoundRectangle2D;

import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.FeatureSegment;
import uk.ac.sanger.artemis.FeatureSegmentVector;
import uk.ac.sanger.artemis.LastSegmentException;
import uk.ac.sanger.artemis.Selection;
import uk.ac.sanger.artemis.FeatureVector;

import uk.ac.sanger.artemis.chado.ChadoFeature;
import uk.ac.sanger.artemis.chado.ChadoFeatureLoc;
import uk.ac.sanger.artemis.chado.ChadoFeatureProp;
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
import uk.ac.sanger.artemis.sequence.MarkerRange;
import uk.ac.sanger.artemis.sequence.Strand;
import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.util.ReadOnlyException;
import uk.ac.sanger.artemis.Options;

public class GeneViewerPanel extends JPanel
{
  
  private ChadoCanonicalGene chado_gene;
  private int border = 15;
  /** Used to colour the frames. */
  private Color light_grey = new Color(240, 240, 240);
  /** pop up menu */
  private JPopupMenu popup;
  /** overlay transcript features */
  private boolean overlay_transcripts = false;
  
  private Selection selection;

  private float fraction;
  
  private int start;
  //private int end;
  
  private MarkerRange click_range = null;
  
  private Point last_cursor_position;
  
  private EntryGroup entry_group;
  
  /**
   *  The shortcut for Delete Selected Features.
   **/
  final static KeyStroke DELETE_FEATURES_KEY =
    KeyStroke.getKeyStroke(KeyEvent.VK_DELETE,
                           Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()); // InputEvent.CTRL_MASK);
  final static public int DELETE_FEATURES_KEY_CODE = KeyEvent.VK_DELETE;
  
  final static KeyStroke CREATE_FEATURES_KEY =
    KeyStroke.getKeyStroke (KeyEvent.VK_C,
                            Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()); //InputEvent.CTRL_MASK);

  final static public int CREATE_FEATURES_KEY_CODE = KeyEvent.VK_C;
  
  
  public GeneViewerPanel(final ChadoCanonicalGene chado_gene,
                         final Selection selection,
                         final EntryGroup entry_group)
  {
    this.chado_gene  = chado_gene;
    this.selection   = selection;
    this.entry_group = entry_group;
    
    Dimension dim = new Dimension(400,400);
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
        if(event.isPopupTrigger())
          return;
        
        final MarkerRange selected_range = selection.getMarkerRange();
        
        int select_start = (int)((event.getX() - border)/fraction)+start;
        Strand strand = 
          ((uk.ac.sanger.artemis.Feature)(chado_gene.getGene().getUserData())).getStrand();
          
        if(!strand.isForwardStrand())
          select_start = strand.getBases().getComplementPosition(select_start);
        
        try
        {
          MarkerRange drag_range = 
            new MarkerRange(strand, select_start, select_start+1);
          
          //final MarkerRange new_marker_range;
          if(selected_range == null || click_range == null) 
            click_range = drag_range;
          else 
            click_range = selected_range.combineRanges(drag_range, true);
          
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
        
        int option = JOptionPane.showConfirmDialog(null, 
            "Delete selected features", 
            "Delete selected features", 
            JOptionPane.OK_CANCEL_OPTION);
        
        if(option == JOptionPane.CANCEL_OPTION)
          return;
        try
        {
          for(int i = 0; i < features.size(); i++)
          {
            uk.ac.sanger.artemis.Feature feature = features.elementAt(i);
            Vector children = chado_gene.getChildren(feature.getEmblFeature());
            deleteFeature(feature);
            chado_gene.deleteFeature(feature.getEmblFeature());    
            
            Feature embl_feature;
            for(int j=0; j<children.size(); j++)
            {  
              embl_feature = (Feature)children.get(j);
              deleteFeature((uk.ac.sanger.artemis.Feature)embl_feature.getUserData());
              chado_gene.deleteFeature(embl_feature);
            }       
          }
          repaint();
        }
        catch(ReadOnlyException e)
        {
          e.printStackTrace();
        }
        
      }
    });
    menu.add(deleteMenu);
    
    
    JMenuItem deleteSegmentMenu = new JMenuItem("Delete Selected Exon");
    deleteSegmentMenu.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)  
      {
        FeatureSegmentVector features = selection.getAllSegments();
        
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
            uk.ac.sanger.artemis.FeatureSegment segment = 
              features.elementAt(i);
            segment.removeFromFeature();
          }
          repaint();
        }
        catch(ReadOnlyException e)
        {
          e.printStackTrace();
        }
        catch(LastSegmentException e)
        {
          e.printStackTrace();
        }
      }
    });
    menu.add(deleteSegmentMenu);
    
    menu.add(new JSeparator());
    JMenuItem createTranscript = new JMenuItem("Create transcript");
    createTranscript.setAccelerator(CREATE_FEATURES_KEY);
    createTranscript.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)  
      {
        try
        {
          String gene_name = 
            (String)chado_gene.getGene().getQualifierByName("ID").getValues().get(0);
          
          String ID = chado_gene.autoGenerateTanscriptName("mRNA");
          
          QualifierVector qualifiers = new QualifierVector();
          qualifiers.add(new Qualifier("Parent", gene_name));
          
          if(ID != null)
            qualifiers.add(new Qualifier("ID", ID));
          
          uk.ac.sanger.artemis.Feature feature = createFeature(
                        chado_gene.getGene().getLocation(),
                        entry_group, new Key("mRNA"),
                        qualifiers);
          
          chado_gene.addTranscript(feature.getEmblFeature());
        }
        catch(InvalidRelationException e)
        {
          // TODO Auto-generated catch block
          e.printStackTrace();
        } 

        repaint();
      }
    });
    menu.add(createTranscript);
    
    JMenuItem createExon = new JMenuItem("Add exon from selected range");
    createExon.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)  
      {
        if(last_cursor_position == null)
          return;
        Feature transcript = getTranscriptAt(last_cursor_position);
        String uniquename  = getFeatureID(transcript);
        
        List exons = chado_gene.getExonsOfTranscript(uniquename);
        GFFStreamFeature embl_exon = null;
        if(exons != null)
          embl_exon = (GFFStreamFeature)exons.get(0);
        Range range_selected = selection.getSelectionRange();
    
        addExonFeature(embl_exon, range_selected, uniquename);
      }
    });
    menu.add(createExon);
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
    g2d.drawString(gene.getIDString(), border, ypos);
    drawFeature(g2d, border, 
                getSize().width - border, 
                ypos, gene.getColour(), 1, 
                selection.contains(gene), 2.f);
    
    List transcripts = chado_gene.getTranscripts();   
    fraction = (float)(getSize().width - (2*border))/
               (float)(end-start);
    
    ypos += border*2;
    
    for(int i=0; i<transcripts.size(); i++)
    {
      /*
      // draw frame lines
   
      if(!overlay_transcripts || i == 0)
        drawFrameLines(g2d, ypos, 
            start, end, fraction);
    
      drawTranscriptOnFrameLine(g2d, (Feature)transcripts.get(i), 
          start, end, ypos, 
          fraction);
   
      if(!overlay_transcripts)
        ypos += 9 * getFontHeight() * 2;
      */
       
      drawTranscriptOnLine(g2d, (Feature)transcripts.get(i), 
                           start, end, ypos, 
                           fraction);
      
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
                    ypos, Color.YELLOW, 2,
                    false, 2.f);
      }
    }
    setPreferredSize(new Dimension(getSize().width, ypos+border));
  }
  
  
  private uk.ac.sanger.artemis.Feature createFeature(
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
   * Macro for getting the size of the transcipt and
   * exon image.
   * @return
   */
  protected int getTranscriptSize()
  {
    return (2 * border) + (getFontHeight() * 4);  
  }
  
  protected int getViewerBorder()
  {
    return border; 
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
    return null;
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
        String transcript_name = getFeatureID(feature);
        List exons = chado_gene.getExonsOfTranscript( transcript_name );
        
        if(exons != null)
        {
          FeatureSegmentVector segments = ((uk.ac.sanger.artemis.Feature) 
              ((Feature) exons.get(0)).getUserData()).getSegments();

          if(segments.size() == 0)
            return (Feature) exons.get(0);

          for(int j = 0; j < segments.size(); j++)
          {
            FeatureSegment segment = segments.elementAt(j);
            Range segment_range = segment.getRawRange();

            int segment_start = border
                + (int) ((segment_range.getStart() - start) * fraction);
            int segment_end = border
                + (int) ((segment_range.getEnd() - start) * fraction);
            if(p.x >= segment_start && p.x <= segment_end)
              return segment;
          }
        }
        
        List utr_3 = chado_gene.get3UTRTranscript(transcript_name);
        List utr_5 = chado_gene.get5UTRTranscript(transcript_name);
        List utrs = null;
        if(utr_3 != null)
          utrs = utr_3;
        if(utr_5 != null)
        {
          if(utrs == null)
            utrs = utr_5;
          else
            utrs.addAll(utr_5);
        }
        
        if(utrs != null)
        {
          for(int j=0; j<utrs.size(); j++)
          {
            Feature utr = (Feature)utrs.get(j);
            Range range = utr.getLocation().getTotalRange();
            int utr_start = border
                + (int) ((range.getStart() - start) * fraction);
            int utr_end = border
                + (int) ((range.getEnd() - start) * fraction);
            if(p.x >= utr_start && p.x <= utr_end)
              return utr;
          }
        }
      }
    }
         
    return null;
  }
  
  /**
   * Draw the features on frame lines.
   * @param g2d
   * @param embl_transcript
   * @param start
   * @param end
   * @param ypos
   * @param fraction
   */
  private void drawTranscriptOnFrameLine(Graphics2D g2d, Feature embl_transcript, 
                                         int start, int end, int ypos, 
                                         float fraction)
  {

    uk.ac.sanger.artemis.Feature transcript = 
       (uk.ac.sanger.artemis.Feature)embl_transcript.getUserData();
    
    int t_start = border+(int)((embl_transcript.getFirstBase()-start)*fraction);
    int t_end   = border+(int)((embl_transcript.getLastBase()-start)*fraction);
    
    g2d.setColor( transcript.getColour() );
    int nframe;
    if(!embl_transcript.getLocation().isComplement())
      nframe = (3 * getFontHeight() * 2) ;
    else
      nframe = (4 * getFontHeight() * 2) ;
    
    g2d.drawString(transcript.getIDString(), border, ypos+nframe);
    drawFeature(g2d, t_start, t_end, 
                ypos+nframe, transcript.getColour(), 1, 
                selection.contains(transcript), 2.f);
    
    List exons = chado_gene.getExonsOfTranscript(
        getFeatureID( embl_transcript ));
        //(String)embl_transcript.getQualifierByName("ID").getValues().get(0));
    
    if(exons == null)
    {
      if(!overlay_transcripts)
        ypos += 9 * getFontHeight() * 2;
      return;
    }
    
    int offset = 0;
    boolean last_segment = false;
    
    if(exons.get(0) instanceof ChadoFeature)
    {
      int last_ex_start = 0;
      int last_ex_end   = 0;
      int last_ypos     = 0;
      
      ChadoFeature start_exon = (ChadoFeature)exons.get(0);
      ChadoFeatureLoc loc = ChadoFeature.getFeatureLoc(
           start_exon.getFeaturelocsForFeatureId(), chado_gene.getSrcfeature_id());
      
      if(loc.getStrand() == -1)
      {
        ChadoFeatureLoc loc_last = ChadoFeature.getFeatureLoc(
            ((ChadoFeature)exons.get(exons.size()-1)).getFeaturelocsForFeatureId(),
            chado_gene.getSrcfeature_id());
            
        if(loc.getFmin() < loc_last.getFmin())
          Collections.reverse(exons);
      }
      
      for(int j=0; j<exons.size(); j++)
      {           
        ChadoFeature exon = (ChadoFeature)exons.get(j);
        loc = ChadoFeature.getFeatureLoc(
            exon.getFeaturelocsForFeatureId(), chado_gene.getSrcfeature_id());
        
        int ex_start = border+(int)((loc.getFmin()+1-start)*fraction);
        int ex_end   = border+(int)((loc.getFmax()-start)*fraction);
           
        Color exon_col = getColorFromAttributes(exon);
     
        offset = getFrameID(chado_gene, loc, j, exons) * getFontHeight() * 2;
        
        boolean isForward = false;
        if(loc.getStrand() == 1)
          isForward = true;
        
        if(j == exons.size()-1)
          last_segment = true;
        
        drawExons(g2d, ex_start, ex_end, 
                  last_ex_start, last_ex_end, last_ypos,
                  offset, ypos, exon_col,
                  1, isForward, last_segment, false, 2.f);
        
        last_ex_end   = ex_end;
        last_ex_start = ex_start;
        last_ypos   = ypos+offset;
      }

      return;
    }
    
    
    // build from artemis objects
    for(int j=0; j<exons.size(); j++)
    {
      int last_ex_start = 0;
      int last_ex_end   = 0;
      int last_ypos     = 0;
      
      Feature embl_exon = (Feature)exons.get(j);
      
      uk.ac.sanger.artemis.Feature exon = 
        (uk.ac.sanger.artemis.Feature)embl_exon.getUserData();
        
      FeatureSegmentVector segments = exon.getSegments();

      for(int k=0; k<segments.size(); k++)
      {
        FeatureSegment segment = segments.elementAt(k);
        
        Range range = segment.getRawRange();
        offset = segment.getFrameID() * getFontHeight() * 2;
        
        int ex_start = border+(int)((range.getStart()-start)*fraction);
        int ex_end   = border+(int)((range.getEnd()-start)*fraction);

        if(exon.getColour() != null)
          g2d.setColor( exon.getColour() );
             

        if(k == segments.size()-1)
          last_segment = true;
        
        drawExons(g2d, ex_start, ex_end, 
                 last_ex_start, last_ex_end, last_ypos,
                 offset, ypos, exon.getColour(),
                 1, segment.isForwardSegment(), last_segment,
                 selection.contains(exon), 2.f);
        
        last_ex_end   = ex_end;
        last_ex_start = ex_start;
        last_ypos   = ypos+offset; 
      }
    }
  
  }
  
  /**
   * 
   * @param g2d
   * @param embl_transcript
   * @param start
   * @param end
   * @param ypos
   * @param fraction
   */
  private void drawTranscriptOnLine(Graphics2D g2d, Feature embl_transcript, 
                                    final int start, final int end, int ypos, 
                                    float fraction)
  {
    BasicStroke stroke = new BasicStroke(48.f);
    g2d.setStroke(stroke);
    
    uk.ac.sanger.artemis.Feature transcript = 
       (uk.ac.sanger.artemis.Feature)embl_transcript.getUserData();

    int t_start = border+(int)((embl_transcript.getFirstBase()-start)*fraction);
    int t_end   = border+(int)((embl_transcript.getLastBase()-start)*fraction);

    g2d.setColor( transcript.getColour() );

    g2d.drawString(transcript.getIDString(), border, ypos);
    drawFeature(g2d, t_start, t_end, 
                ypos, transcript.getColour(), 1, 
                selection.contains(transcript), 2.f);

    List exons = chado_gene.getExonsOfTranscript(
        getFeatureID( embl_transcript ));

    ypos += border*2;
    
    if(exons != null)
    {
      boolean last_segment = false;
    
      // build from artemis objects
      for(int i=0; i<exons.size(); i++)
      {
        int last_ex_start = 0;
        int last_ex_end = 0;
        int last_ypos = 0;

        Feature embl_exon = (Feature) exons.get(i);

        uk.ac.sanger.artemis.Feature exon = 
          (uk.ac.sanger.artemis.Feature)embl_exon.getUserData();

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

          int ex_start = border + (int)((range.getStart() - start) * fraction);
          int ex_end   = border + (int)((range.getEnd() - start) * fraction);

          if(exon.getColour() != null)
            g2d.setColor(exon.getColour());

          if(j == ranges.size() - 1)
            last_segment = true;

          selected_size = 2.f;
          if(segments != null)
          {
            for(int k = 0; k < segments.size(); k++)
            {
              FeatureSegment segment = segments.elementAt(k);
              if(range.equals(segment.getRawRange())
                  && selection.contains(segment))
                selected_size = 4.f;
            }
          }

          drawExons(g2d, ex_start, ex_end, last_ex_start, last_ex_end,
              last_ypos, 0, ypos, exon.getColour(), 2, exon.isForwardFeature(),
              last_segment, selection.contains(exon), selected_size);

          last_ex_end = ex_end;
          last_ex_start = ex_start;
          last_ypos = ypos;
        }
      }
    }
    
    // draw utr's
    List embl_utr = chado_gene.get3UTRTranscript(
        getFeatureID( embl_transcript ));
    
    if(embl_utr != null)
      drawUTR(g2d, embl_utr, ypos);
    
    embl_utr = chado_gene.get5UTRTranscript(
        getFeatureID( embl_transcript ));
    
    if(embl_utr != null)
      drawUTR(g2d, embl_utr, ypos);
    
  }
  
  /**
   * Method to draw UTR features
   * @param g2d
   * @param embl_utr
   * @param ypos
   */
  private void drawUTR(final Graphics2D g2d,
                       final List utrs,
                       final int ypos)
  {
    for(int i=0; i<utrs.size(); i++)
    {
      Feature embl_utr = (Feature)utrs.get(i);
      uk.ac.sanger.artemis.Feature utr = 
        (uk.ac.sanger.artemis.Feature)embl_utr.getUserData();
      RangeVector ranges = utr.getLocation().getRanges();
      for(int j = 0; j < ranges.size(); j++)
      {
        Range range = (Range) ranges.get(j);

        int utr_start = border + (int) ((range.getStart() - start) * fraction);
        int utr_end = border + (int) ((range.getEnd() - start) * fraction);

        drawFeature(g2d, utr_start, utr_end, ypos, utr.getColour(), 2,
            selection.contains(utr), 2.f);
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
  private void drawFrameLines(Graphics2D g2d, int ypos,
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
  }
  
  private void drawExons(Graphics2D g2d, int ex_start, int ex_end, 
                         int last_ex_start, int last_ex_end, int last_ypos,
                         int offset, int ypos, Color exon_colour,
                         int size,
                         boolean isForward, boolean last_segment,
                         boolean selected, float selected_size)
  {   
    drawFeature(g2d, ex_start, ex_end, 
                ypos+offset, exon_colour, size, selected, selected_size);
    
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
                     last_ex_end+((ex_start-last_ex_end)/2), ymid-getFontHeight()/2);
        g2d.drawLine(last_ex_end+((ex_start-last_ex_end)/2), ymid-getFontHeight()/2, 
                     ex_start, ypos+offset); 
      }
      else
      { 
        g2d.drawLine(last_ex_start, last_ypos, 
                     last_ex_start+((ex_end-last_ex_start)/2), ymid-getFontHeight()/2);
        g2d.drawLine(last_ex_start+((ex_end-last_ex_start)/2), ymid-getFontHeight()/2, 
                     ex_end, ypos+offset); 
      }  
    }
  
    // draw arrow
    if(last_segment)
    {
      if(isForward)
      {      

        g2d.drawLine(ex_end, ypos, 
                     ex_end+getFontHeight()/2, ypos+(getFontHeight()*size)/2);
        g2d.drawLine(ex_end+getFontHeight()/2, ypos+(getFontHeight()*size)/2,
                     ex_end, ypos+(getFontHeight()*size));
      }
      else
      {
        g2d.drawLine(ex_start, ypos, 
                     ex_start-getFontHeight()/2, ypos+(getFontHeight()*size)/2);
        g2d.drawLine(ex_start-getFontHeight()/2, ypos+(getFontHeight()*size)/2,
                     ex_start, ypos+(getFontHeight()*size));
      }
    }
  }
  
  /**
   * Draw rectangular box for a feature.
   * @param g2d
   * @param start   start of feature
   * @param end     end of feature
   * @param ypos    y position
   * @param colour  feature colour
   * @param size    parameter to control the height of the feature
   */
  private void drawFeature(Graphics2D g2d, int start, int end, 
                           int ypos, Color colour, int size,
                           boolean selected, float selected_size)
  {
    RoundRectangle2D e = new RoundRectangle2D.Float(start, ypos, 
        end-start,
        getFontHeight()*size, 0, ypos);

    if(colour == null)
      colour = Color.BLACK;
    
    GradientPaint gp = new GradientPaint(start, ypos, 
        colour,
        start, ypos+( (getFontHeight()/2) * size ), 
        Color.white, true);
    g2d.setPaint(gp); 
    g2d.fill(e);
    
    if(selected)
      g2d.setStroke(new BasicStroke(selected_size));
    else
      g2d.setStroke(new BasicStroke(1.f));
    
    // draw boundary
    g2d.setColor(Color.BLACK);
    g2d.draw(e);
  }
  
  /**
   * Get the <code>Color</code> for a feature from its colour attribute.
   * @param feature
   * @return
   */
  private Color getColorFromAttributes(ChadoFeature feature)
  {
    List properties = feature.getFeaturepropList();
    for(int i=0; i<properties.size(); i++)
    {
      ChadoFeatureProp property = (ChadoFeatureProp)properties.get(i);
      
      if(property.getCvterm().getName().equals("colour") ||
         property.getCvterm().getName().equals("color") )
        return Options.getOptions().getColorFromColourNumber(Integer.parseInt(property.getValue()));
    }  
    return Color.CYAN;
  }
  
  /**
   * Get the frame id for a feature segment
   * @param chado_gene  the chado representation of the gene model
   * @param loc         feature location of the feature
   * @param nexon       number of the exon
   * @param exons       List of exons
   * @return frame id
   */
  private int getFrameID(ChadoCanonicalGene chado_gene, 
                         ChadoFeatureLoc loc, 
                         int nexon, List exons)
  {
    final int position_on_strand;
    
    if(loc.getStrand() == -1)
      position_on_strand = chado_gene.getSeqlen()-loc.getFmax();
    else
      position_on_strand = loc.getFmin();
    
    // this will be 0, 1 or 2 depending on which frame the segment is in
    final int start_base_modulo =
      (position_on_strand + getFrameShift(nexon, exons, chado_gene, loc)) % 3;

    if(loc.getStrand() == 1)
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
  }
  
  
  /**
   *  Returns 0, 1 or 2 depending on which translation frame this segment is
   *  in.  A frame shift of zero means that the bases should be translated
   *  starting at the start position of this segment, 1 means start
   *  translating one base ahead of the start position and 2 means start
   *  translating two bases ahead of the start position.
   **/
  private int getFrameShift(int nexon, List exons, 
                            ChadoCanonicalGene chado_gene,
                            ChadoFeatureLoc loc) 
  {
    // find the number of bases in the segments before this one
    int base_count = 0;
    int direction  = 0;
    
    for(int i = 0; i < exons.size(); ++i) 
    {
      ChadoFeature this_feature = (ChadoFeature)exons.get(i);
      ChadoFeatureLoc featureLoc = ChadoFeature.getFeatureLoc(
          this_feature.getFeaturelocsForFeatureId(), chado_gene.getSrcfeature_id());
      
      int this_direction;
      if(featureLoc.getStrand() == 1)
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

        base_count += featureLoc.getFmax()-featureLoc.getFmin();
      }
    }
    
    int codon_start = loc.getPhase();
    int mod_value   = (base_count + 3 - codon_start) % 3;

    if(mod_value == 1) 
      return 2;
    else if(mod_value == 2)
      return 1;
    else 
      return 0;
  }
  
  
  
  private int getFontHeight()
  {
    final FontMetrics fm = this.getFontMetrics(getFont());
    return fm.getHeight();  
  }
  
  /**
   * Display a list of the exons for a transcript.
   * @param embl_transcript
   * @throws InvalidRelationException
   */
  private void showExonsList(final String uniquename)
               throws InvalidRelationException
  {
    final List exons = chado_gene.getExonsOfTranscript(uniquename);
    final DefaultListModel listModel = new DefaultListModel();   
    
    if(exons != null)
    {
      Feature embl_exon = (Feature)exons.get(0);
      final RangeVector ranges = embl_exon.getLocation().getRanges();
      
      for(int k=0; k<ranges.size(); k++)
      {
        Range range = (Range)ranges.get(k);
        listModel.addElement(range.toString());
      }
    } 
    
    final JList list = new JList(listModel);
    JScrollPane listScrollPane = new JScrollPane(list);

    JButton addButt    = new JButton("ADD");
    addButt.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)  
      {
        GFFStreamFeature embl_exon = null;
        
        if(exons != null)
          embl_exon = (GFFStreamFeature)exons.get(0);
        Range range_selected = selection.getSelectionRange();
        
        addExonFeature(embl_exon, range_selected, uniquename);
      }
    });

    JButton deleteButt = new JButton("DELETE");
    deleteButt.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)  
      {
        Object selected[] = list.getSelectedValues();
        Feature embl_exon = (Feature)exons.get(0);
        uk.ac.sanger.artemis.Feature exon =
          (uk.ac.sanger.artemis.Feature)embl_exon.getUserData();
        deleteExonFeatures(exon, selected, listModel, uniquename);       
      }
    });
    
    JButton editButt   = new JButton("EDIT");
    editButt.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)  
      {

      }
    });
    
    //Create a panel that uses BoxLayout.
    JPanel buttonPane = new JPanel();
    buttonPane.setLayout(new BoxLayout(buttonPane,
                                       BoxLayout.LINE_AXIS));
    buttonPane.add(listScrollPane);
    buttonPane.add(Box.createHorizontalStrut(5));
    buttonPane.add(new JSeparator(SwingConstants.VERTICAL));
    buttonPane.add(Box.createHorizontalStrut(5));
    
    Box bdown = Box.createVerticalBox();
    bdown.add(addButt);
    bdown.add(deleteButt);
    bdown.add(editButt);
    buttonPane.add(bdown);
    buttonPane.setBorder(BorderFactory.createEmptyBorder(5,5,5,5));
    buttonPane.add(Box.createHorizontalGlue());
    final JDialog dialog = new JDialog((Frame)null,
                                       "Exon List", false);
    
    dialog.setContentPane(buttonPane);
    dialog.pack();
    dialog.setVisible(true);

  }
 
  /**
   * Create an exons feature
   * @param feature
   * @param range
   * @param transcript_name
   */
  private void addExonFeature(GFFStreamFeature feature, Range range,
                              final String transcript_name)
  {
    try
    {
      if(feature == null)
      {
        QualifierVector qualifiers = new QualifierVector();
        qualifiers.add(new Qualifier("Parent", transcript_name));
      
        String ID = chado_gene.autoGenerateExonName(transcript_name);
        if(ID != null)
          qualifiers.add(new Qualifier("ID", ID));
        
        uk.ac.sanger.artemis.Feature exon = createFeature(
            new Location(range),
            entry_group, new Key("exon"),
            qualifiers);
      
        GFFStreamFeature gff_exon = (GFFStreamFeature)exon.getEmblFeature();
        chado_gene.addExon(transcript_name, gff_exon);
        
        if(ID != null)
        {
          Hashtable id_range_store = new Hashtable();
          id_range_store.put(ID, range);
          gff_exon.setSegmentRangeStore(id_range_store);
        }
      }
      else
      {
        // add new ID
        ((uk.ac.sanger.artemis.Feature)feature.getUserData()).addSegment(range);
      }
    }
    catch(OutOfRangeException e)
    {
      e.printStackTrace();
    }
    catch(InvalidRelationException e)
    {
      e.printStackTrace();
    }
    catch(ReadOnlyException e)
    {
      e.printStackTrace();
    }
  }
  
  private void deleteExonFeatures(uk.ac.sanger.artemis.Feature feature,
                                  Object selected[], 
                                  final DefaultListModel listModel,
                                  final String transcript_uniquename) 
  {
    FeatureSegmentVector segments = feature.getSegments();
    Range range = null;
    FeatureSegment segment = null;
    selection.clear();
    
    try
    {
      int option = JOptionPane.showConfirmDialog(null, 
          "Delete selected feature exons", 
          "Delete selected feature exons", 
          JOptionPane.OK_CANCEL_OPTION);
      
      if(option == JOptionPane.CANCEL_OPTION)
        return;
      
      for(int i = 0; i < segments.size(); i++)
      {
        segment = segments.elementAt(i);
        range   = segment.getRawRange();

        for(int j = 0; j < selected.length; j++)
        {
          if(range.toString().equals((String) selected[j]))
          {
            segment.removeFromFeature();
            listModel.removeElement(range.toString());
          }
        }
      }
    }
    catch(ReadOnlyException e)
    {
      e.printStackTrace();
    }
    catch(LastSegmentException e)
    {
      // delete entire feature here
      Feature transcript = 
         (Feature)chado_gene.getFeatureFromId(transcript_uniquename);
      try
      {
        // delete exon feature and associated transcript
        deleteFeature(segment.getFeature());
        deleteFeature((uk.ac.sanger.artemis.Feature)transcript.getUserData());
        
        listModel.removeElement(range.toString());
        chado_gene.deleteTranscript(transcript_uniquename);
        repaint();
      }
      catch(ReadOnlyException e1)
      {
        e1.printStackTrace();
      }
    }
  }
  
  private void deleteFeature(uk.ac.sanger.artemis.Feature feature) 
          throws ReadOnlyException
  {
    feature.removeFromEntry();
  }
  
  private String getFeatureID(final Feature feature)
  {
    String id = null;
    try
    {
      id = (String)feature.getQualifierByName("ID").getValues().get(0);
    }
    catch(InvalidRelationException e)
    {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }  
    return id;
  }
  
  /**
   * Popup listener
   */
  class PopupListener extends MouseAdapter
  {
    public void mousePressed(MouseEvent e)
    {
      maybeShowPopup(e);
    }

    public void mouseReleased(MouseEvent e)
    {
      maybeShowPopup(e);
    }

    private void maybeShowPopup(MouseEvent e)
    {
      if(e.isPopupTrigger())
      {
        final Feature embl_transcript = getTranscriptAt(e.getPoint());
        if(embl_transcript == null)
          return;
        
        popup.show(e.getComponent(),
                e.getX(), e.getY());
      }
      else if(e.getButton() != e.BUTTON3 &&
              e.getClickCount() == 1)
      {
        if(selection.getMarkerRange() != null &&
           e.isShiftDown())
          return;
        
        click_range = null;
        selection.clear();
        Object feature = getFeatureAt(e.getPoint());
        if(feature == null)
          return;      
        
        if(feature instanceof Feature)
          selection.set((uk.ac.sanger.artemis.Feature)((Feature)feature).getUserData());
        else
          selection.set((FeatureSegment)feature);

        repaint();
      }
    }
  }
 
}