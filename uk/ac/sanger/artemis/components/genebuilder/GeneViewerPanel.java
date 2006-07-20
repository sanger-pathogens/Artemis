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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/genebuilder/GeneViewerPanel.java,v 1.9 2006-07-20 10:32:05 tjc Exp $
 */

package uk.ac.sanger.artemis.components.genebuilder;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.util.List;
import java.util.Collections;
import java.awt.geom.RoundRectangle2D;

import uk.ac.sanger.artemis.FeatureChangeEvent;
import uk.ac.sanger.artemis.FeatureChangeListener;
import uk.ac.sanger.artemis.FeatureSegment;
import uk.ac.sanger.artemis.FeatureSegmentVector;
import uk.ac.sanger.artemis.LastSegmentException;
import uk.ac.sanger.artemis.Selection;
import uk.ac.sanger.artemis.chado.ChadoFeature;
import uk.ac.sanger.artemis.chado.ChadoFeatureLoc;
import uk.ac.sanger.artemis.chado.ChadoFeatureProp;
import uk.ac.sanger.artemis.io.Feature;
import uk.ac.sanger.artemis.io.ChadoCanonicalGene;
import uk.ac.sanger.artemis.io.InvalidRelationException;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.util.ReadOnlyException;
import uk.ac.sanger.artemis.Options;

public class GeneViewerPanel extends JPanel
          implements FeatureChangeListener
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

  public GeneViewerPanel(final ChadoCanonicalGene chado_gene,
                         final Selection selection)
  {
    this.chado_gene = chado_gene;
    this.selection  = selection;
    
    try
    {
      addListeners();
    }
    catch(InvalidRelationException e)
    {
      e.printStackTrace();
    } 
    
    
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

  /**
   * Add feature listeners for each artemis feature.
   * @throws InvalidRelationException
   */
  private void addListeners() throws InvalidRelationException
  {
    // add feature listeners
    Feature embl_gene = (Feature)chado_gene.getGene();
    uk.ac.sanger.artemis.Feature gene =
      (uk.ac.sanger.artemis.Feature)embl_gene.getUserData();
    gene.addFeatureChangeListener(this);
    
    List transcripts = chado_gene.getTranscripts();
    for(int i=0; i<transcripts.size(); i++)
    {
       Feature transcript = (Feature)transcripts.get(i);
       uk.ac.sanger.artemis.Feature trans = 
         (uk.ac.sanger.artemis.Feature)transcript.getUserData();
       trans.addFeatureChangeListener(this);
       List exons = chado_gene.getExonsOfTranscript(
           (String)trans.getQualifierByName("ID").getValues().get(0));
       
       if(exons == null)
         continue;
       
       if(exons.get(0) instanceof ChadoFeature)
         return;
       
       for(int j=0; j<exons.size(); j++)
       {
         Feature embl_exon = (Feature)exons.get(j);

         uk.ac.sanger.artemis.Feature exon = 
           (uk.ac.sanger.artemis.Feature)embl_exon.getUserData();
         exon.addFeatureChangeListener(this);
       }
    }
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

    int start = embl_gene.getFirstBase();
    int end   = embl_gene.getLastBase();
    g2d.setColor( gene.getColour() );
    int ypos = border;
    
    // draw gene
    g2d.drawString(gene.getIDString(), border, ypos);
    drawFeature(g2d, border, 
                getSize().width - border, 
                ypos, gene.getColour(), 1, selection.contains(gene));
    
    List transcripts = chado_gene.getTranscripts();   
    float fraction = (float)(getSize().width - (2*border))/
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
    setPreferredSize(new Dimension(getSize().width, ypos+border));
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
    
    int ntranscript = (p.y - (border*2))/getTranscriptSize();
    if(ntranscript < transcripts.size())
      return (Feature)transcripts.get(ntranscript);
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
                ypos+nframe, transcript.getColour(), 1, selection.contains(transcript));
    
    try
    {          
      List exons = chado_gene.getExonsOfTranscript(
          (String)embl_transcript.getQualifierByName("ID").getValues().get(0));
      
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
                    1, isForward, last_segment);
          
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
                   1, segment.isForwardSegment(), last_segment);
          
          last_ex_end   = ex_end;
          last_ex_start = ex_start;
          last_ypos   = ypos+offset; 
        }
      }
      
    }
    catch(InvalidRelationException e)
    {
      e.printStackTrace();
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
                                    int start, int end, int ypos, 
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
                ypos, transcript.getColour(), 1, selection.contains(transcript));

    try
    {          
      List exons = chado_gene.getExonsOfTranscript(
       (String)embl_transcript.getQualifierByName("ID").getValues().get(0));

      if(exons == null)
        return; 

      ypos += border*2;
      
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

          boolean isForward = false;
          if(loc.getStrand() == 1)
            isForward = true;

          if(j == exons.size()-1)
            last_segment = true;
          
          drawExons(g2d, ex_start, ex_end, 
              last_ex_start, last_ex_end, last_ypos,
              0, ypos, exon_col,
              2, isForward, last_segment);

          last_ex_end   = ex_end;
          last_ex_start = ex_start;
          last_ypos   = ypos;
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

          int ex_start = border+(int)((range.getStart()-start)*fraction);
          int ex_end   = border+(int)((range.getEnd()-start)*fraction);

          if(exon.getColour() != null)
            g2d.setColor( exon.getColour() );

          if(k == segments.size()-1)
            last_segment = true;
          
          drawExons(g2d, ex_start, ex_end, 
              last_ex_start, last_ex_end, last_ypos,
              0, ypos, exon.getColour(),
              2, segment.isForwardSegment(),
              last_segment);

          last_ex_end   = ex_end;
          last_ex_start = ex_start;
          last_ypos   = ypos;
        }
      }
    }
    catch(InvalidRelationException e)
    {
      e.printStackTrace();
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
                         boolean isForward, boolean last_segment)
  {   
    drawFeature(g2d, ex_start, ex_end, 
                ypos+offset, exon_colour, size, false);
    
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
                           boolean selected)
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
      g2d.setStroke(new BasicStroke(2.f));
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
    List exons = chado_gene.getExonsOfTranscript(uniquename);
    if(exons == null)
      return;
    

    Feature embl_exon = (Feature)exons.get(0);

    final uk.ac.sanger.artemis.Feature exon = 
          (uk.ac.sanger.artemis.Feature)embl_exon.getUserData();
    
    final FeatureSegmentVector segments = exon.getSegments();
    final DefaultListModel listModel = new DefaultListModel();     
    for(int k=0; k<segments.size(); k++)
    {
      FeatureSegment segment = segments.elementAt(k);
      Range range = segment.getRawRange();
      listModel.addElement(range.toString());
    }
    
    final JList list = new JList(listModel);
    JScrollPane listScrollPane = new JScrollPane(list);

    JButton addButt    = new JButton("ADD");
    addButt.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)  
      {

      }
    });

    JButton deleteButt = new JButton("DELETE");
    deleteButt.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)  
      {
        Object selected[] = list.getSelectedValues();
        final FeatureSegmentVector segments = exon.getSegments();
        Range range = null;
        FeatureSegment segment = null;
        try
        {
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
          // TODO Auto-generated catch block
          e.printStackTrace();
        }
        catch(LastSegmentException e)
        {
          // delete entire feature here
          Feature transcript = (Feature)chado_gene.getFeatureFromId(uniquename);
          try
          {
            ((uk.ac.sanger.artemis.Feature)transcript.getUserData()).removeFromEntry();
            listModel.removeElement(range.toString());
            chado_gene.deleteTranscript(uniquename);
          }
          catch(ReadOnlyException e1)
          {
            e1.printStackTrace();
          }
        }
        
        
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
 
  public void featureChanged(FeatureChangeEvent event)
  {
    Feature feature = event.getFeature().getEmblFeature();
    String uniquename = null;
    try
    {
      uniquename = (String)(feature.getQualifierByName("ID").getValues().get(0));
    }
    catch(InvalidRelationException e)
    {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
    
    repaint();
    //revalidate();
  }
  
  /**
   * Popup listener
   */
  class PopupListener extends MouseAdapter
  {
    private JMenuItem exonMenu = new JMenuItem();
    private ActionListener exonListener;
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
        exonMenu.removeActionListener(exonListener);
        final Feature embl_transcript = getTranscriptAt(e.getPoint());
        if(embl_transcript == null)
          return;
        
        try
        {
          final String uniquename = 
            (String)(embl_transcript.getQualifierByName("ID").getValues().get(0));
            
          exonMenu.setText("Show exon list for "+uniquename);       
          exonListener = new ActionListener()
          {
            public void actionPerformed(ActionEvent event)  
            {
              try
              {
                showExonsList(uniquename);
              }
              catch(InvalidRelationException e)
              {
                e.printStackTrace();
              }
            }
          };
          exonMenu.addActionListener(exonListener);
          popup.add(exonMenu);
        }
        catch(InvalidRelationException e1)
        {
          e1.printStackTrace();
        }
        
        popup.show(e.getComponent(),
                e.getX(), e.getY());
      }
    }
  }
 
}