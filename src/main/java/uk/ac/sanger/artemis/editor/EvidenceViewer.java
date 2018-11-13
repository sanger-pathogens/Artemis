/*
 *
 * created: Wed Aug 3 2004
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2000  Genome Research Limited
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
 */

package uk.ac.sanger.artemis.editor;

import javax.swing.JPanel;
import javax.swing.JOptionPane;
import javax.swing.JDesktopPane;
import java.awt.*;
import java.awt.event.*;
import java.awt.geom.AffineTransform;
import java.util.StringTokenizer;

import uk.ac.sanger.artemis.io.Location;
import uk.ac.sanger.artemis.*;
import uk.ac.sanger.artemis.sequence.Marker;

public class EvidenceViewer extends JPanel
{

  /** */
  private static final long serialVersionUID = 1L;
  //private String evidenceType;
  //private int evidenceStart;
  //private int evidenceEnd;
  private int featStart;
  private int featEnd;
  /** viewer boundary */
  private int bound = 35;
  private Feature edit_feature;
  private FeatureVector overlapFeature;
  private float unitLength;
  private int YDISPLACEMENT = 20;
  private JDesktopPane desktop;
  
  /** URL of Pfam database. */
  protected String pfamUrl;


  /**
   * Constructor.
   * 
   * @param edit_feature
   * @param overlapFeature
   * @param desktop
   */
  public EvidenceViewer(Feature edit_feature, 
                        FeatureVector overlapFeature, JDesktopPane desktop)
  {
    this.edit_feature   = edit_feature;
    this.overlapFeature = overlapFeature;
    final Location this_loc = edit_feature.getLocation();
    featStart = this_loc.getFirstBase();
    featEnd   = this_loc.getLastBase();

    this.desktop = desktop;
    //this.evidenceStart = evidenceStart;
    //this.evidenceEnd   = evidenceEnd;

    final int hgt = 100+(20*overlapFeature.size());
    final Dimension dim = new Dimension(500,hgt);
    setPreferredSize(dim);
    setMaximumSize(dim);
    
    setPfamUrl();

    addMouseListener(new MouseClickListener());
  }

  /**
  *
  * Override paintComponent
  * @param g    graphics
  *
  */
  public void paintComponent(Graphics g)
  {
// let UI delegate paint first (incl. background filling)
    super.paintComponent(g);

    Font font = new Font("monospaced",Font.PLAIN,10);
    g.setFont(font);
    FontMetrics metrics = g.getFontMetrics();
    int hgtNumber = metrics.getAscent();

    Graphics2D g2   = (Graphics2D)g;
    int resultwidth = (int)getPreferredSize().getWidth()-(2*bound);
    int scale    = 1;
    int npoints  = (int)(scale*10);
    int unit     = (int)(resultwidth/npoints);
    int featUnit = (featEnd-featStart)/npoints;

    unitLength = (float)resultwidth/(float)(featEnd-featStart);

// feature segments
    final FeatureSegmentVector this_feature_segments = edit_feature.getSegments();
    // draw each segment/exon
    for(int i = 0 ; i < this_feature_segments.size() ; ++i)
    {
      final FeatureSegment current_segment =
                           this_feature_segments.elementAt(i);

      drawSegment(g2, current_segment, unitLength, featStart, hgtNumber);

//    if(i + 1 < this_feature_segments.size())
//    {
//        // draw a line between the segments

//      final FeatureSegment next_segment =
//        this_feature_segments.elementAt(i + 1);

//      drawSegmentConnection(g, current_segment, next_segment);
//    }
    }

//
    g2.setColor(Color.black);
    g2.setStroke(new BasicStroke(3.f));
    g2.drawLine(bound,bound+hgtNumber,bound+resultwidth,bound+hgtNumber);

// draw ruler
    g2.setStroke(new BasicStroke(2.f));

    for(int i=0; i<=npoints; i++)
      g2.drawLine(bound+(unit*i),bound+hgtNumber,
                  bound+(unit*i),bound+hgtNumber-6);

    AffineTransform origin = g2.getTransform();
    AffineTransform newOrig = (AffineTransform)(origin.clone());
    newOrig.rotate(Math.toRadians(90.),0,0);
    g2.setTransform(newOrig);

    String strPos = Integer.toString(featStart);
    int strwid = metrics.stringWidth(strPos);
    int hgtNumber2 = hgtNumber/2;
    g2.drawString(strPos,bound-strwid,-bound+hgtNumber2);

    for(int i=1; i<=npoints; i++)
    {
      strPos = Integer.toString((i*featUnit) + featStart);
      strwid = metrics.stringWidth(strPos);
      g2.drawString(strPos,bound-strwid,-bound-(unit*i)+hgtNumber2);
    }

    g2.setTransform(origin);

// 
    int bound2 = bound*2;
    int ydisp  = 0;
    for(int i = 0; i<overlapFeature.size(); ++i)
    {
      Feature this_feature = overlapFeature.elementAt(i);
      String note = this_feature.getNote();

      int indPF   = note.indexOf("PF");  // Pfam ID
      if(indPF == -1)
        continue;

      StringTokenizer tok = new StringTokenizer(note.substring(indPF)," ,;");
      String pfamID = tok.nextToken();
      Location this_loc = this_feature.getLocation();

      int start = this_loc.getFirstBase();
      int end   = this_loc.getLastBase();

      start = (int)((start-featStart)*unitLength);
      end   = (int)((end-featStart)*unitLength);

      g2.setColor(Color.red);
      g2.setStroke(new BasicStroke(15.f));
      g2.drawLine(bound+start, bound2+ydisp,
                  bound+end,   bound2+ydisp);

      strwid = metrics.stringWidth(pfamID);
      g2.setColor(Color.black);
      g2.setStroke(new BasicStroke(3.f));
      g2.drawString(pfamID, bound+start+(end-start-strwid)/2, bound2+ydisp+5);
     
      ydisp += YDISPLACEMENT;
    }
  }


  private void drawSegment(Graphics2D g2, FeatureSegment segment,
                           float unitLength, int featStart, int hgtNumber)
  {
    final Feature segment_feature = segment.getFeature();

    final Marker segment_start_marker = segment.getStart();
    int startPosition = segment_start_marker.getRawPosition();

    final Marker segment_end_marker = segment.getEnd();
    int endPosition = segment_end_marker.getRawPosition();

    final Color feature_colour = segment_feature.getColour();
    g2.setColor(feature_colour);

    g2.setStroke(new BasicStroke(10.f));
    int start = (int)((startPosition-featStart)*unitLength);
    int stop  = (int)((endPosition-featStart)*unitLength);

    g2.drawLine(bound+start, bound+hgtNumber+6,
                bound+stop,  bound+hgtNumber+6);
  }

  public class MouseClickListener implements MouseListener
  {
    public void mousePressed(MouseEvent e) 
    {
    }

    public void mouseReleased(MouseEvent e) 
    {
    }

    public void mouseEntered(MouseEvent e) 
    {
    }

    public void mouseExited(MouseEvent e) 
    {
    }

    public void mouseClicked(MouseEvent e) 
    {
      if(e.getClickCount() < 2)
        return;

      setCursor(new Cursor(Cursor.WAIT_CURSOR));

      int ydisp  = 0;
      int bound2 = bound*2;
      int ydisp2 = YDISPLACEMENT/2;
      Point loc  = e.getPoint();
      for(int i = 0; i<overlapFeature.size(); ++i)
      {
        Feature this_feature = overlapFeature.elementAt(i);
        Location this_loc = this_feature.getLocation();
        int start = this_loc.getFirstBase();
        int end   = this_loc.getLastBase();

        start = (int)((start-featStart)*unitLength)+bound;
        end   = (int)((end-featStart)*unitLength)+bound;

        if(loc.x < end &&
           loc.x > start &&
           loc.y < ydisp+ydisp2+bound2 &&
           loc.y > ydisp-ydisp2+bound2)
        {
          String note = this_feature.getNote();

          int indPF   = note.indexOf("PF");  // Pfam ID
          if(indPF == -1)
            continue;

          StringTokenizer tok = new StringTokenizer(note.substring(indPF)," ,;");
          String pfamID = tok.nextToken();


          String pfam_cmd = pfamUrl+pfamID;
          BrowserControl.displayURL(pfam_cmd);

          if(BigPane.srsTabPane.isSelected())
          {
            try
            {
              DataCollectionPane.setUpSRSFrame(new java.net.URL(pfam_cmd),
                                               pfamID,desktop);
            }
            catch(java.net.ConnectException connect)
            {
              JOptionPane.showMessageDialog(EvidenceViewer.this,
                       "Cannot retrieve "+pfamID+
                       "\nConnection failed to:\n"+pfamID,
                       "Connection Error",
                       JOptionPane.WARNING_MESSAGE);
            }
            catch(Exception exp)
            {
              exp.printStackTrace();
              setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
            }
          }
        }
        ydisp += YDISPLACEMENT;
      }
      setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
    }
  }
  
  /**
   * Determine the PFAM database URL from properties
   * if possible.
   */
  protected void setPfamUrl() 
  {
	  if (pfamUrl == null)
	  {
		  pfamUrl = Options.getOptions().getDatabaseHyperlinkProperty(Options.PFAM_HYPERLINK_PROPERTY_NAME);
	  }
  }

}

