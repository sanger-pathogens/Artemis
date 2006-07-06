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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/genebuilder/GeneViewerPanel.java,v 1.5 2006-07-06 15:10:14 tjc Exp $
 */

package uk.ac.sanger.artemis.components.genebuilder;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.util.List;
import java.util.Collections;
import java.awt.geom.RoundRectangle2D;

import uk.ac.sanger.artemis.FeatureSegment;
import uk.ac.sanger.artemis.FeatureSegmentVector;
import uk.ac.sanger.artemis.chado.ChadoFeature;
import uk.ac.sanger.artemis.chado.ChadoFeatureLoc;
import uk.ac.sanger.artemis.chado.ChadoFeatureProp;
import uk.ac.sanger.artemis.io.Feature;
import uk.ac.sanger.artemis.io.ChadoCanonicalGene;
import uk.ac.sanger.artemis.io.InvalidRelationException;
import uk.ac.sanger.artemis.io.Range;
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

  public GeneViewerPanel(final ChadoCanonicalGene chado_gene)
  {
    this.chado_gene = chado_gene;
    
    Dimension dim = new Dimension(400,400);
    setPreferredSize(dim);
    setBackground(Color.white);
    
//  Popup menu
    addMouseListener(new PopupListener());
    popup = new JPopupMenu();

    JCheckBoxMenuItem menuItem = new JCheckBoxMenuItem("Overlay transcripts");
    menuItem.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)  
      {
        overlay_transcripts = !overlay_transcripts;
        repaint();
        revalidate();
      }
    });

    popup.add(menuItem);
  }

  public void paintComponent(Graphics g)
  {
    super.paintComponent(g);
    Graphics2D g2d = (Graphics2D)g;
    Feature embl_gene = (Feature)chado_gene.getGene();
    uk.ac.sanger.artemis.Feature gene =
      (uk.ac.sanger.artemis.Feature)embl_gene.getUserData();
    
    setFont(uk.ac.sanger.artemis.Options.getOptions().getFont());
    //final FontMetrics fm = this.getFontMetrics(getFont());

    int start = embl_gene.getFirstBase();
    int end   = embl_gene.getLastBase();
    
    BasicStroke stroke = new BasicStroke(4.f);
    Stroke original_stroke = g2d.getStroke();
    g2d.setStroke(stroke);
    g2d.setColor( gene.getColour() );
    
    int ypos = border;
    
    // draw gene
    g2d.drawString(gene.getIDString(), border, ypos);
    drawFeature(g2d, border, getSize().width - border, ypos, gene.getColour());
    
    List transcripts = chado_gene.getTranscripts();
    Feature embl_transcript;
    uk.ac.sanger.artemis.Feature transcript;
    
    float fraction = (float)(getSize().width - (2*border))/
                     (float)(end-start);
    
    ypos += border*2;
    
    for(int i=0; i<transcripts.size(); i++)
    {
      // draw frame lines
      if(!overlay_transcripts || i == 0)
        drawFrameLines(g2d, ypos, stroke,
            start, end, fraction);

      embl_transcript  = (Feature)transcripts.get(i);
      transcript = (uk.ac.sanger.artemis.Feature)embl_transcript.getUserData();
      
      int t_start = border+(int)((embl_transcript.getFirstBase()-start)*fraction);
      int t_end   = border+(int)((embl_transcript.getLastBase()-start)*fraction);
      
      g2d.setColor( transcript.getColour() );
      int nframe;
      if(!embl_transcript.getLocation().isComplement())
        nframe = (3 * getFontHeight() * 2) ;
      else
        nframe = (4 * getFontHeight() * 2) ;
      
      g2d.drawString(transcript.getIDString(), border, ypos+nframe);
      drawFeature(g2d, t_start, t_end, ypos+nframe, transcript.getColour());
      
      try
      {
        List exons = chado_gene.getExonsOfTranscript(
            (String)embl_transcript.getQualifierByName("ID").getValues().get(0));
        
        if(exons == null)
          continue;
        
        int offset = 0;
        
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
            
            drawExons(g2d, ex_start, ex_end, 
                      last_ex_start, last_ex_end, last_ypos,
                      offset, ypos, exon_col,
                      original_stroke, stroke, isForward);
            
            last_ex_end   = ex_end;
            last_ex_start = ex_start;
            last_ypos   = ypos+offset;
          }
          
          if(!overlay_transcripts)
            ypos += 9 * getFontHeight() * 2;
          continue;
        }
        
        
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
                       
            drawExons(g2d, ex_start, ex_end, 
                     last_ex_start, last_ex_end, last_ypos,
                     offset, ypos, exon.getColour(),
                     original_stroke, stroke, segment.isForwardSegment());
            
            last_ex_end   = ex_end;
            last_ex_start = ex_start;
            last_ypos   = ypos+offset;
          }
          
          ypos += (8 * getFontHeight() * 2) + (getFontHeight()/2) + border;
        }
        
      }
      catch(InvalidRelationException e)
      {
        e.printStackTrace();
      }
    }
    setPreferredSize(new Dimension(getSize().width, ypos+border));
  }
  
  
  /**
   * Draw frame lines
   * @param g2d
   * @param ypos
   * @param stroke
   * @param start
   * @param end
   * @param fraction
   */
  private void drawFrameLines(Graphics2D g2d, int ypos,
                              Stroke stroke,
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
    g2d.setStroke(stroke);
  }
  
  private void drawExons(Graphics2D g2d, int ex_start, int ex_end, 
                         int last_ex_start, int last_ex_end, int last_ypos,
                         int offset, int ypos, Color exon_colour,
                         Stroke original_stroke, Stroke stroke, boolean isForward)
  {
    /*
    RoundRectangle2D e = new RoundRectangle2D.Float(ex_start, ypos+offset, 
                                                    ex_end-ex_start,
                                                    getFontHeight(), 0, ypos+offset);
    
    GradientPaint gp = new GradientPaint(ex_start, ypos+offset, 
                                         exon_colour,
                                         ex_start, ypos+offset+(getFontHeight()/2), 
                                         Color.white, true);
    g2d.setPaint(gp); 
    g2d.fill(e);
    */
    
    drawFeature(g2d, ex_start, ex_end, ypos+offset, exon_colour);
    
    // draw connections
    if(last_ex_end != 0 ||
       last_ex_start != 0)
    {
      g2d.setStroke(original_stroke);
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
      g2d.setStroke(stroke);  
    }

  }
  
  private void drawFeature(Graphics2D g2d, int start, int end, 
                           int ypos, Color colour)
  {
    RoundRectangle2D e = new RoundRectangle2D.Float(start, ypos, 
        end-start,
        getFontHeight(), 0, ypos);

    if(colour == null)
      colour = Color.BLACK;
    
    GradientPaint gp = new GradientPaint(start, ypos, 
        colour,
        start, ypos+(getFontHeight()/2), 
        Color.white, true);
    g2d.setPaint(gp); 
    g2d.fill(e);
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
        popup.show(e.getComponent(),
                e.getX(), e.getY());
    }
  }
  
}