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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/genebuilder/GeneViewerPanel.java,v 1.2 2006-05-31 15:40:53 tjc Exp $
 */

package uk.ac.sanger.artemis.components.genebuilder;

import javax.swing.*;
import java.awt.*;
import java.util.Vector;
import java.awt.geom.RoundRectangle2D;

import uk.ac.sanger.artemis.io.Feature;
import uk.ac.sanger.artemis.io.ChadoCanonicalGene;
import uk.ac.sanger.artemis.io.InvalidRelationException;

public class GeneViewerPanel extends JPanel
{
  private ChadoCanonicalGene chado_gene;
  private int border = 15;
  
  public GeneViewerPanel(final ChadoCanonicalGene chado_gene)
  {
    this.chado_gene = chado_gene;
    Dimension dim = new Dimension(300,300);
    setPreferredSize(dim);
    setBackground(Color.white);
  }

  public void paintComponent(Graphics g)
  {
    Graphics2D g2d = (Graphics2D)g;
    Feature embl_gene = chado_gene.getGene();
    uk.ac.sanger.artemis.Feature gene = (uk.ac.sanger.artemis.Feature)embl_gene.getUserData();
    
    setFont(uk.ac.sanger.artemis.Options.getOptions().getFont());
    final FontMetrics fm = this.getFontMetrics(getFont());

    int start = embl_gene.getFirstBase();
    int end   = embl_gene.getLastBase();
    boolean complement = embl_gene.getLocation().isComplement();
    
    BasicStroke stroke = new BasicStroke(2.f);
    g2d.setStroke(stroke);
    g2d.setColor( gene.getColour() );
    
    int ypos = border;
    
    ypos += fm.getHeight();
    g2d.drawString(gene.getIDString(), border, ypos);
    ypos += fm.getHeight();
    
    g2d.drawLine(border, ypos, getSize().width - border, ypos);
    
    Vector transcripts = chado_gene.getTranscripts();
    Feature embl_transcript;
    uk.ac.sanger.artemis.Feature transcript;
    
    float fraction = (float)(getSize().width - (2*border))/(float)(end-start);
    for(int i=0; i<transcripts.size(); i++)
    {
      ypos += border;
      embl_transcript  = (Feature)transcripts.get(i);
      transcript = (uk.ac.sanger.artemis.Feature)embl_transcript.getUserData();
      
      int t_start = border+(int)((embl_transcript.getFirstBase()-start)*fraction);
      int t_end   = border+(int)((embl_transcript.getLastBase()-start)*fraction);
      
      ypos += fm.getHeight();
      g2d.drawString(transcript.getIDString(), border, ypos);
      ypos += fm.getHeight();
      
      g2d.setColor( transcript.getColour() );
      g2d.drawLine(t_start, ypos, t_end, ypos);
      
      try
      {
        Vector exons = chado_gene.getExonsOfTranscript(
            (String)embl_transcript.getQualifierByName("ID").getValues().get(0));
        
        ypos += border;
        
        if(exons == null)
          continue;
        
        for(int j=0; j<exons.size(); j++)
        {
          Feature embl_exon = (Feature)exons.get(j);
          uk.ac.sanger.artemis.Feature exon = 
            (uk.ac.sanger.artemis.Feature)embl_exon.getUserData();
          
          int ex_start = border+(int)((embl_exon.getFirstBase()-start)*fraction);
          int ex_end   = border+(int)((embl_exon.getLastBase()-start)*fraction);
          
          if(exon.getColour() != null)
            g2d.setColor( exon.getColour() );
          RoundRectangle2D e = new RoundRectangle2D.Float(ex_start, ypos, ex_end-ex_start,
                                                          border, 0, ypos);
          GradientPaint gp = new GradientPaint(ex_start, ypos, exon.getColour(),
                                               ex_start, ypos+(border/2), Color.white, true);
          g2d.setPaint(gp); 
          g2d.fill(e);

          //g2d.drawLine(ex_start, ypos, ex_end, ypos);
        }
      }
      catch(InvalidRelationException e)
      {
        e.printStackTrace();
      }
    }
  }
}