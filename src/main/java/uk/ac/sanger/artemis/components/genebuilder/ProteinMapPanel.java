/* ProteinMapPanel.java
 *
 * created: 2008
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2008  Genome Research Limited
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
 **/

package uk.ac.sanger.artemis.components.genebuilder;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionListener;
import java.util.Collections;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.List;
import java.util.Vector;

import uk.ac.sanger.artemis.Selection;
import uk.ac.sanger.artemis.components.QualifierTextArea;
import uk.ac.sanger.artemis.io.Feature;
import uk.ac.sanger.artemis.io.ChadoCanonicalGene;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierLazyLoading;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.StringVector;

public class ProteinMapPanel extends MapPanel
{
  private static final long serialVersionUID = 1L;
  public static String[] TMHMM = 
  { 
    "membrane_structure", 
    "non_cytoplasm_location",
    "non_cytoplasmic_polypeptide_region",
    "transmembrane",
    "transmembrane_polypeptide_region",
    "cytoplasm_location",
    "cytoplasmic_polypeptide_region"
  };
  protected static String GPI_ANCHORED   = "GPI_anchor_cleavage_site";
  //private static String PlasmoAP_SCORE = "PlasmoAP_score";
  protected static String[] SIGNALP = 
  { 
    "signal_peptide"
    //"SignalP_prediction", 
    //"signal_peptide_probability",
    //"signal_anchor_probability" 
  };
  
  public static String POLYPEPTIDE_DOMAIN = "polypeptide_domain";
  
  protected static Vector<String> PROTEIN_MAP_ELEMENTS = new Vector<String>();
  static
  {
    Collections.addAll(PROTEIN_MAP_ELEMENTS, TMHMM);
    Collections.addAll(PROTEIN_MAP_ELEMENTS, SIGNALP);
    PROTEIN_MAP_ELEMENTS.add(GPI_ANCHORED);
    //PROTEIN_MAP_ELEMENTS.add(PlasmoAP_SCORE);
    PROTEIN_MAP_ELEMENTS.add(POLYPEPTIDE_DOMAIN);
  }

  protected Hashtable<Rectangle, String> toolTipPositions = new Hashtable<Rectangle, String>();
  protected Feature feature;

  public ProteinMapPanel(final Feature feature, 
                         final ChadoCanonicalGene chado_gene,
                         final Selection selection) 
  {
    this.chado_gene = chado_gene;
    this.selection  = selection;
    this.feature    = feature;
    
    Dimension dim = new Dimension(400,350);
    setPreferredSize(dim);
    setBackground(Color.white);
    setToolTipText("");
    
    addMouseMotionListener(new MouseMotionListener()
    {
      public void mouseDragged(MouseEvent e) {}

      public void mouseMoved(MouseEvent e)
      {
        final String tt = ProteinMapPanel.this.getToolTipText(e);
        if(tt != null && isDatabaseHyperLink(false,tt))
          setCursor(new Cursor(Cursor.HAND_CURSOR));
        else
          setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
      }
      
    });
    
    addMouseListener(new MouseAdapter()
    {
      public void mouseExited(MouseEvent e)
      {
        setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
      }
      
      public void mouseClicked(MouseEvent e)
      {
        final String tt = ProteinMapPanel.this.getToolTipText(e);
        isDatabaseHyperLink(true,tt);
      }
    });
  }

  public void paintComponent(Graphics g)
  {
    super.paintComponent(g);
    
    if(this instanceof BasicProteinMapPanel)
      return;
    
    Graphics2D g2d = (Graphics2D)g;
    
    QualifierVector qualifiers = feature.getQualifiers();
    
    setFont(uk.ac.sanger.artemis.Options.getOptions().getFont());
    
    Feature embl_gene = (Feature)chado_gene.getGene();
    uk.ac.sanger.artemis.Feature gene =
      (uk.ac.sanger.artemis.Feature)embl_gene.getUserData();
    
    setFont(uk.ac.sanger.artemis.Options.getOptions().getFont());

    int ypos = border;
    
    int geneStart = embl_gene.getFirstBase();
    int geneEnd   = embl_gene.getLastBase();
    float fraction = (float)(getSize().width - (2*border))/
                     (float)(geneEnd-geneStart);
    
    final List<Feature> transcripts = chado_gene.getTranscripts();
    for(int i=0; i<transcripts.size(); i++)
    {
      Feature transcript = (Feature)transcripts.get(i);
      String transcriptName = GeneUtils.getUniqueName(transcript);
      
      Feature protein_embl_feature = 
        (Feature)chado_gene.getProteinOfTranscript(transcriptName);
      
      if(protein_embl_feature == null)
        continue;
      
      uk.ac.sanger.artemis.Feature protein = 
        (uk.ac.sanger.artemis.Feature)protein_embl_feature.getUserData();
      
      List<Feature> exons = chado_gene.getSpliceSitesOfTranscript(transcriptName, 
                                            DatabaseDocument.EXONMODEL);
      if(exons == null || exons.size() == 0)
        exons = chado_gene.getSpliceSitesOfTranscript(transcriptName,
                                                 "pseudogenic_exon");
      
      if(exons == null || exons.size() == 0)
        continue;
      
      Feature exon_embl_feature = exons.get(0);
      
      uk.ac.sanger.artemis.Feature exon = 
        (uk.ac.sanger.artemis.Feature)exon_embl_feature.getUserData();
      
      int ppLength = exon.getTranslationBasesLength()/3;
      
      int emblStart = protein_embl_feature.getFirstBase();
      int emblEnd   = protein_embl_feature.getLastBase();
      
      // draw protein
      g2d.drawString(protein.getIDString()+
          ( GeneUtils.isObsolete((GFFStreamFeature)embl_gene) ? " (obsolete)" : "") , border, ypos);
      
      int ppStart = border+(int)((emblStart-geneStart)*fraction);
      int ppEnd   = border+(int)((emblEnd-geneStart)*fraction);
      
      drawFeature(g2d, ppStart, ppEnd, 
                  ypos, protein.getColour(), 1, 
                  selection.contains(gene), 2.f, getFontHeight());
      
      //record position for tooltip
      Rectangle r = new Rectangle(ppStart, ypos, ppEnd-ppStart, 1*getFontHeight());
      toolTipPositions.put(r, protein.getIDString()+" length="+ppLength);
      
      ypos += border*2;
      
      ypos = drawDomain(protein_embl_feature, 
            g2d, ypos, ppStart, ppEnd, ppLength);
      
      if(qualifiers.getQualifierByName(TMHMM[0]) != null)
      {
        g2d.drawString("Transmembrane Domains:", ppStart, ypos);
        drawPrediction(protein_embl_feature, 
                     g2d, ypos, ppStart, ppEnd, ppLength, TMHMM);
        ypos += border * 2;
      }
      
      if(qualifiers.getQualifierByName(SIGNALP[0]) != null)
      {
        g2d.drawString("SignalP:", ppStart, ypos);
        drawPrediction(protein_embl_feature, 
                     g2d, ypos, ppStart, ppEnd, ppLength, 
                     new String[]{ SIGNALP[0] });
        ypos += border * 2;
      }
 
      Qualifier gpiAnchor;
      if((gpiAnchor=qualifiers.getQualifierByName(GPI_ANCHORED)) != null)
      {
        g2d.drawString("GPI anchor cleavage site:", ppStart, ypos);
        drawGPIArrow(g2d, gpiAnchor, ppStart, ppEnd, ppLength, ypos);
      }
      
      if(i != transcripts.size()-1)
        ypos += border * 2;
    }
    setPreferredSize(new Dimension(getSize().width, ypos+border));
  }
  
  /**
   * Draw polypeptide_domain
   * @param feature
   * @param g2d
   * @param ypos
   * @param ppStart
   * @param ppEnd
   * @param ppLength
   * @return
   */
  protected int drawDomain(final Feature feature, 
      final Graphics2D g2d, int ypos,
      final int ppStart, final int ppEnd, final int ppLength)
  {
    float fraction = (float)(ppEnd-ppStart)/
                     (float)(ppLength);

    QualifierVector qualifiers = feature.getQualifiers();
    for(int i = 0; i < qualifiers.size(); i++)
    {
      Qualifier qualifier = (Qualifier) qualifiers.get(i);
      if(!qualifier.getName().equals(POLYPEPTIDE_DOMAIN))
        continue;
      
      if(qualifier instanceof QualifierLazyLoading)
        ((QualifierLazyLoading) qualifier).setForceLoad(true);
      StringVector values = qualifier.getValues();

      Collections.sort(values);

      String lastDomain = null;
      for(int j = 0; j < values.size(); j++)
      {
        StringVector parts = StringVector.getStrings((String) values.get(j),
            ";");

        if(lastDomain == null || !lastDomain.equals((String) parts.get(0)))
          lastDomain = (String) parts.get(0);


        StringBuffer toolTipBuff = new StringBuffer();
        for(int k = 0; k < parts.size(); k++)
          toolTipBuff.append((String) parts.get(k)+"\n");
        
        for(int k = 1; k < parts.size(); k++)
        {
          String part = (String) parts.get(k);

          int domainCoords[] = getCoords(part);
          if(domainCoords == null)
            continue;

          int start = ppStart + (int) (domainCoords[0] * fraction);
          int end   = ppStart + (int) (domainCoords[1] * fraction);

          if(k==1)
            g2d.drawString((String) parts.get(0), start, ypos);
          else
            g2d.drawString((String) parts.get(1), start, ypos);

          final Color col = Color.GRAY;

          drawFeature(g2d, start, end, ypos, col, 1, false, 2.f, getFontHeight());
          
          //record position for tooltip
          Rectangle r = new Rectangle(start, ypos, end-start, 1*getFontHeight());
          toolTipPositions.put(r, toolTipBuff.toString());
          
          ypos += border * 2;
        }
      }
    }
    return ypos;
  }
  
  /**
   * Draw TMHMM
   * @param feature
   * @param g2d
   * @param ypos
   * @param ppStart
   * @param ppEnd
   * @param ppLength
   */
  protected void drawPrediction(final Feature feature, 
                        final Graphics2D g2d, final int ypos,
                        final int ppStart, final int ppEnd, final int ppLength,
                        final String[] prediction)
  {
    
    float fraction = (float)(ppEnd-ppStart)/
                     (float)(ppLength);
    
    final QualifierVector qualifiers = feature.getQualifiers();
    
    for(int i=0;i<qualifiers.size(); i++)
    {
      Qualifier qualifier = (Qualifier) qualifiers.get(i);
      String qualifierName = qualifier.getName();
      
      // predictions
      for(int j = 0; j < prediction.length; j++)
      {
        if(qualifierName.equals(prediction[j]))
        { 
          if(qualifier instanceof QualifierLazyLoading)
            ((QualifierLazyLoading) qualifier).setForceLoad(true);

          StringVector values = qualifier.getValues();
          if(qualifierName.equals(TMHMM[0]))
            continue;
          for(int k = 0; k < values.size(); k++)
          {
            int coords[] = getCoords(((String)values.get(k)));
            if(coords == null)
              continue;
            
            int start = ppStart+(int)(coords[0]*fraction);
            int end   = ppStart+(int)(coords[1]*fraction);
            
            final Color col;
            switch(j)
            {
              case 0:  col = Color.BLUE; break;
              case 1:  col = Color.LIGHT_GRAY; break;
              case 2:  col = Color.GREEN; break;
              default: col = Color.YELLOW; break;
            }
            
            drawFeature(g2d, start, end, 
                  ypos, col, 1, false, 2.f, getFontHeight());
            
            
            StringVector parts = StringVector.getStrings((String)values.get(k), ";");
            String tt = qualifierName+" : "+coords[0]+".."+coords[1];
            for(int l=0; l<parts.size(); l++)
              if( ((String)parts.get(l)).indexOf("query") < 0 )
                tt = tt.concat("\n"+(String)parts.get(l));
            // record position for tooltip
            Rectangle r = new Rectangle(start, ypos, end-start, 1*getFontHeight());
            toolTipPositions.put(r, tt);  
          }
        }
      }
    }
  }
  
  /**
   * Draw GPI anchor site
   * @param g2d
   * @param gpiAnchor
   * @param ppStart
   * @param ppEnd
   * @param ppLength
   * @param ypos
   * @return
   */
  protected int drawGPIArrow(final Graphics2D g2d,
      final Qualifier gpiAnchor,
      final int ppStart, final int ppEnd, final int ppLength,
      int ypos)
  {
    float fraction = (float)(ppEnd-ppStart)/
                     (float)(ppLength);
    ((QualifierLazyLoading)gpiAnchor).setForceLoad(true);
    
    StringVector values = gpiAnchor.getValues();
    for(int j=0; j<values.size(); j++)
    {
      int coords[] = getCoords(((String)values.get(j)));
      int start = ppStart+(int)(coords[0]*fraction);
      
      g2d.setStroke(new BasicStroke(2.f));
      g2d.setColor(Color.RED);
      g2d.drawLine(start-1,ypos,start-1,ypos+15);
      
      
      StringVector parts = StringVector.getStrings((String)values.get(j), ";");
      String tt = gpiAnchor.getName()+" : "+coords[0];
      for(int l=0; l<parts.size(); l++)
        if( ((String)parts.get(l)).indexOf("query") < 0 )
          tt = tt.concat("\n"+(String)parts.get(l));
      // record position for tooltip
      Rectangle r = new Rectangle(start-1, ypos, 2, 15);
      toolTipPositions.put(r, tt);  
    }  
    ypos += border * 2;
    return ypos;
  }
  
  /**
   * Get coordinates from the qualifier value by looking
   * for query in the string
   * @param value
   * @return
   */
  private int[] getCoords(String value)
  {
    int coords[] = new int[2];

    int beginIndex = value.indexOf("query ")+6;
    
    if(beginIndex < 6)
      return null;
    
    int endIndex = value.indexOf(';', beginIndex);
    if(endIndex < 0)
      endIndex = value.length();
    value = value.substring(beginIndex, endIndex);
    
    String splitStr[] = value.split("-");
    
    coords[0] = Integer.parseInt(splitStr[0]);
    coords[1] = Integer.parseInt(splitStr[1]);
    return coords;
  }
  
  /**
   * Determine if the tooltip contains a database link. 
   * @param sendToBrowser if true then the link is opened in a browser
   * @param tt tooltip text
   * @return
   */
  private boolean isDatabaseHyperLink(final boolean sendToBrowser,
                                      final String tt)
  {
    if(tt == null)
      return false;
    
    boolean isLink = false;
    final Vector<String> dbs = QualifierTextArea.DATABASES;
    for(int i=0; i<dbs.size(); i++)
    {
      String db = dbs.get(i);
      int beginIndex;
      if( (beginIndex=tt.indexOf(db+":")) > -1 )
      {
        int endIndex = QualifierTextArea.getEndOfLink(tt, beginIndex);
        isLink = true;
        if(sendToBrowser)
          QualifierTextArea.handleMouseSingleClick(
              tt.substring(beginIndex, endIndex),
              ProteinMapPanel.this);
        else
          return isLink;
      }
    }
    return isLink;
  }
  
  public String getToolTipText(MouseEvent me)
  {
    Point p = me.getPoint();
    Enumeration<Rectangle> rectangles = toolTipPositions.keys();
    while(rectangles.hasMoreElements())
    {
      Rectangle r = rectangles.nextElement();
      if((r.x <= p.x && r.x+r.width >= p.x) &&
         (r.y <= p.y && r.y+r.height >= p.y) )
        return (String) toolTipPositions.get(r);
    }
    return null;
  }

  /**
   * Check if a feature contains a protein map qualifier. Return a list 
   * of the protein features containing protein map qualifiers.
   * @param feature
   * @return
   */
  public static List<Feature> getProteinsWithProteinMapElement(final GFFStreamFeature feature)
  {
    List<Feature> transcripts = feature.getChadoGene().getTranscripts();
    List<Feature> proteins = null;
    if(transcripts != null)
    {
      for(int i=0; i<transcripts.size(); i++)
      {
        Feature transcript = transcripts.get(i);
        String transcriptName = GeneUtils.getUniqueName(transcript);
        Feature protein = feature.getChadoGene().getProteinOfTranscript(transcriptName);
        
        if(protein != null)
        {
          QualifierVector qualifiers = protein.getQualifiers();
          for(int j=0; j<qualifiers.size(); j++)
            if(isProteinMapElement((Qualifier)qualifiers.get(j)))
            {
              if(proteins == null)
                proteins = new Vector<Feature>();
              proteins.add(protein);
            }
        }
      }
    }
    return proteins;
  }
  
  public static boolean isProteinMapElement(final Qualifier this_qualifier)
  {
    final String qualifierName = this_qualifier.getName();
    
    if(PROTEIN_MAP_ELEMENTS.contains(qualifierName))
      return true;
    
    return false;
  }

  public static QualifierVector getProteinMapQualifiers(uk.ac.sanger.artemis.Feature f)
  {
    QualifierVector proteinMapQualifiers = null;
    QualifierVector qualifiers = f.getQualifiers();
    for(int i=0; i<qualifiers.size(); i++)
    {
      if(isProteinMapElement((Qualifier)qualifiers.get(i)))
      {
        if(proteinMapQualifiers == null)
          proteinMapQualifiers = new QualifierVector();
        proteinMapQualifiers.addQualifierValues((Qualifier)qualifiers.get(i));
      }
    }
    return proteinMapQualifiers;
  }
}