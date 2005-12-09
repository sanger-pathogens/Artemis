/* ContigTool.java
 *
 * created: 2005
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2005  Genome Research Limited
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

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.editor.MultiLineToolTipUI;

import java.awt.BasicStroke;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FontMetrics;
import java.awt.GradientPaint;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Image;
import java.awt.Insets;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.Shape;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.InputEvent;
import java.awt.geom.Ellipse2D;
import java.awt.geom.RoundRectangle2D;
import javax.swing.*;
import javax.swing.border.Border;
import java.awt.datatransfer.StringSelection;
import java.awt.datatransfer.*;
import java.awt.dnd.*;

public class ContigTool extends JPanel
       implements DragGestureListener, DropTargetListener,
             DragSourceListener, Autoscroll

{

  private FeatureVector contig_features;
  private FeatureVector selected_features = new FeatureVector();
  private FeatureDisplay feature_display;

  private int scale  = 1000;
  private int xbound = 50;
  private int length = xbound*2;

  /** pop up menu */
  private JPopupMenu popup = new JPopupMenu();

  private int highlight_drop_base = -1;
  /** AutoScroll margin */
  private static final int AUTOSCROLL_MARGIN = 45;
  /** used by AutoScroll method */
  private Insets autoscrollInsets = new Insets( 0, 0, 0, 0 );

  /** status label */
  final private JLabel status_line = new JLabel("");

  public ContigTool(final FeatureVector contig_features,
                    final FeatureDisplay feature_display,
                    final JScrollPane jsp)
  {
    super();
    this.contig_features = contig_features;
    this.feature_display = feature_display;

    MultiLineToolTipUI.initialize();
    setToolTipText("");   //enable tooltip display

    for(int i=0; i<contig_features.size(); i++)
    {
      final Range this_feature_range = contig_features.elementAt(i).getMaxRawRange();
      length += ((this_feature_range.getEnd() - this_feature_range.getStart())/scale);
    }
   
    Dimension dim = new Dimension(length, 20);
    setPreferredSize(dim);

    DragSource dragSource = DragSource.getDefaultDragSource();

    dragSource.createDefaultDragGestureRecognizer(
       this,                             // component where drag originates
       DnDConstants.ACTION_COPY_OR_MOVE, // actions
       this);                            // drag gesture recognizer

    setDropTarget(new DropTarget(this,this));

    addMouseListener(new MouseAdapter()
    {
      public void mousePressed(MouseEvent event)
      {
        if(event.isPopupTrigger())
        {
          popup.show(event.getComponent(),
                  event.getX(), event.getY());
          return;
        }

        FeatureVector contig_features = ContigTool.this.contig_features;
        if(event.getClickCount() == 1)
        {
          Point p = event.getPoint();
          for(int i=0; i<contig_features.size(); i++)
          {
            final Feature feature = contig_features.elementAt(i);
            final Range this_feature_range = feature.getMaxRawRange();

            int xstart = xbound + this_feature_range.getStart()/scale;
            int xend   = xbound + this_feature_range.getEnd()/scale;

            if(p.x >= xstart && p.x <= xend)
            {
              if(selected_features.contains(feature))
                selected_features.remove(feature);
              else
              {
                selected_features.removeAllElements();

                String tt = this_feature_range.getStart()+".."+
                            this_feature_range.getEnd();
 
                if(feature.getIDString() != null)
                  tt = tt + ", " + feature.getIDString();

                status_line.setText(tt);
                selected_features.add(feature); 
              }
            }
          }
          repaint();
        }
      }
   
      public void mouseReleased(MouseEvent event)
      {
      }
    });

    // set popup menu items
    JMenuItem zoomIn = new JMenuItem("Zoom In - x0.1");
    popup.add(zoomIn);
    zoomIn.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)
      {
        if(scale < 10)
        {
          scale = 1; 
          return;
        }

        scale = scale / 10;  
        adjustSize(jsp);
        repaint();
      }
    });

    JMenuItem zoomOut = new JMenuItem("Zoom Out - x10");
    popup.add(zoomOut);
    zoomOut.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)
      {
        scale = scale * 10;
        adjustSize(jsp);
        repaint();
      }
    });

    // set up status bar
    status_line.setFont(Options.getOptions().getFont());
    final FontMetrics fm =
      this.getFontMetrics(status_line.getFont());

    final int font_height = fm.getHeight()+10;

    status_line.setMinimumSize(new Dimension(100, font_height));
    status_line.setPreferredSize(new Dimension(100, font_height));

    Border loweredbevel = BorderFactory.createLoweredBevelBorder();
    Border raisedbevel = BorderFactory.createRaisedBevelBorder();
    Border compound = BorderFactory.createCompoundBorder(raisedbevel,loweredbevel);
    status_line.setBorder(compound);
  }

  /**
   *
   * Used when changing the scale / magnification.
   * @param jsp  scrollpane to reset viewport size
   *
   */
  private void adjustSize(final JScrollPane jsp)
  {
    length = xbound*2;
    for(int i=0; i<contig_features.size(); i++)
    {
      final Range this_feature_range = contig_features.elementAt(i).getMaxRawRange();
      length += ((this_feature_range.getEnd() - this_feature_range.getStart())/scale);
    }

    Dimension dim = new Dimension(length, 60);
    setPreferredSize(dim);
    jsp.revalidate();
  }

  /**
   *
   * Get the status bar.
   * @return status label
   *
   */
  protected JLabel getStatusBar()
  {
    return status_line;
  }

  /**
   *
   * Override paintComponent()
   *
   */
  protected void paintComponent(Graphics g)
  {
    Graphics2D g2 = (Graphics2D)g;
    super.paintComponent(g2);
    
    setFont(Options.getOptions().getFont());
    final FontMetrics fm =
      this.getFontMetrics(getFont());

    for(int i=0; i<contig_features.size(); i++)
    {
      final Feature feature = contig_features.elementAt(i);
      final Range this_feature_range = feature.getMaxRawRange();     
      Color colour = feature.getColour();
     
      int xstart = xbound + this_feature_range.getStart()/scale;
      int xend   = xbound + this_feature_range.getEnd()/scale;

      RoundRectangle2D e = new RoundRectangle2D.Float(xstart, 10, xend-xstart,
                                                      20, 0, 10);
      GradientPaint gp = new GradientPaint(xstart, 10, colour,
                                           xstart, 10+10, Color.white, true);
      g2.setPaint(gp);
      g2.fill(e);
    }

    BasicStroke stroke  = (BasicStroke)g2.getStroke();
    BasicStroke stroke1 = new BasicStroke(1.f);
    BasicStroke stroke2 = new BasicStroke(2.f);
    g2.setColor(Color.black);

    // draw feature outline
    for(int i=0; i<contig_features.size(); i++)
    {
      final Feature feature = contig_features.elementAt(i);
      final Range this_feature_range = feature.getMaxRawRange();
      int xstart = xbound + this_feature_range.getStart()/scale;
      int xend   = xbound + this_feature_range.getEnd()/scale;

      if(selected_features.contains(feature))
        g2.setStroke(stroke2);
      else
        g2.setStroke(stroke1);
 
      g2.drawRect(xstart, 10, xend-xstart, 20);

      final String label_or_gene = feature.getIDString();
      final int string_width = fm.stringWidth(label_or_gene);

      Shape saved_clip = g.getClip();
      g2.setColor(Color.black);
      g2.setClip(xstart, 10, xend-xstart, 20);
      g2.drawString(label_or_gene, xstart,
                    10 + fm.getMaxAscent() + 1);
      g2.setClip(saved_clip);
    }


    g2.setStroke(stroke);
    if(highlight_drop_base > 0)
    {
      g2.setColor(Color.red);
      final int draw_x_position = xbound + highlight_drop_base/scale;
      int nlines = 16;

      g.drawLine(draw_x_position, 0,
                 draw_x_position, 100);
    }

  }

  /**
   *
   * Determine the tool tip to display
   * @param e    mouse event
   * @return     tool tip
   *
   */
  public String getToolTipText(MouseEvent e)
  {
    Point loc = e.getPoint();
    int pos = loc.x*scale - xbound;
    int first;
    int last;

    for(int i = 0; i < contig_features.size(); i++)
    {
      final Feature this_feature = contig_features.elementAt(i);
      first = this_feature.getRawFirstBase();
      last  = this_feature.getRawLastBase();
      if(pos >= first && pos <=last)
      {
        return first+".."+last;
      }
    }
    return null;
  }


  
////////////////////
// DRAG AND DROP
////////////////////
 
  /**
   *
   * Given a point find the nearest start/stop of a feature and
   * set highlight_drop_base.
   *
   */ 
  private void getNearestFeatureEnd(Point loc)
  {
    final int base_pos = (loc.x-50)*scale;
    int first;
    int last;

    for(int i = 0; i < contig_features.size(); i++)
    {
      final Feature this_feature = contig_features.elementAt(i);

      if(this_feature.getKey().equals("fasta_record"))
      {
        first = this_feature.getRawFirstBase();
        last  = this_feature.getRawLastBase();

        if( Math.abs(first - base_pos) < Math.abs(base_pos - highlight_drop_base) )
          highlight_drop_base = first;
        if( Math.abs(last - base_pos) < Math.abs(base_pos - highlight_drop_base) )
          highlight_drop_base = last+1;
      }
    }
  }

// drop
  public void drop(DropTargetDropEvent e)
  {
    Transferable t = e.getTransferable();
    if(e.isDataFlavorSupported(DataFlavor.stringFlavor))
    {
      feature_display.reorder(highlight_drop_base, 
                              selected_features.elementAt(0));   // rearrange contigs
      repaint();
    }
    highlight_drop_base = -1;
  }

  public void dragExit(DropTargetEvent e)
  {
    highlight_drop_base = -1;
  }

  public void dropActionChanged(DropTargetDragEvent e) {}

  public void dragOver(DropTargetDragEvent e)
  {
    if(e.isDataFlavorSupported(DataFlavor.stringFlavor))
    {
      Point ploc = e.getLocation();
      getNearestFeatureEnd(ploc);
      repaint();
    }
    e.rejectDrag();
  }

  public void dragEnter(DropTargetDragEvent e)
  {
    if(e.isDataFlavorSupported(DataFlavor.stringFlavor))
      e.acceptDrag(DnDConstants.ACTION_COPY_OR_MOVE);
  }

// drag source
  public void dragGestureRecognized(DragGestureEvent e)
  {
    // ignore if mouse popup trigger
    InputEvent ie = e.getTriggerEvent();
    if(ie instanceof MouseEvent)
      if(((MouseEvent)ie).isPopupTrigger())
        return;

    if(selected_features.size() == 1 &&
       ( selected_features.elementAt(0).getKey().equals("source") ||
         selected_features.elementAt(0).getKey().equals("fasta_record") ))
    {
      ClassLoader cl = this.getClass().getClassLoader();
      ImageIcon icon = new ImageIcon(cl.getResource("images/icon.gif"));
      final Image icon_image = icon.getImage();

      TransferableContig tcontig = new TransferableContig(selected_features.elementAt(0));
      StringSelection name = new StringSelection(selected_features.elementAt(0).getGeneName());

      e.startDrag(DragSource.DefaultCopyDrop,     // cursor
                  icon_image, new Point(-1, -1),
                 (Transferable)name,              // transferable data
                                       this);     // drag source listener
    }
  }

  public void dragDropEnd(DragSourceDropEvent e) {}
  public void dragEnter(DragSourceDragEvent e) {}

  public void dragExit(DragSourceEvent e)
  {
    highlight_drop_base = -1;
  }
  public void dragOver(DragSourceDragEvent e) {}
  public void dropActionChanged(DragSourceDragEvent e) {}


////////////////////
// AUTO SCROLLING //
////////////////////
  /**
   *
   * Handles the auto scrolling of the JTree.
   * @param location The location of the mouse.
   *
   */
  public void autoscroll( Point location )
  {
    int top = 0, left = 0, bottom = 0, right = 0;
    Dimension size = getSize();
    Rectangle rect = getVisibleRect();
    int bottomEdge = rect.y + rect.height;
    int rightEdge = rect.x + rect.width;
    if( location.y - rect.y < AUTOSCROLL_MARGIN && rect.y > 0 )
      top = AUTOSCROLL_MARGIN;
    if( location.x - rect.x < AUTOSCROLL_MARGIN && rect.x > 0 )
      left = AUTOSCROLL_MARGIN;
    if( bottomEdge - location.y < AUTOSCROLL_MARGIN && bottomEdge < size.height )
      bottom = AUTOSCROLL_MARGIN;
    if( rightEdge - location.x < AUTOSCROLL_MARGIN && rightEdge < size.width )
      right = AUTOSCROLL_MARGIN;
    rect.x += right - left;
    rect.y += bottom - top;
    scrollRectToVisible( rect );
  }

  /**
   *
   * Gets the insets used for the autoscroll.
   * @return The insets.
   *
   */
  public Insets getAutoscrollInsets()
  {
    Dimension size = getSize();
    Rectangle rect = getVisibleRect();
    autoscrollInsets.top = rect.y + AUTOSCROLL_MARGIN;
    autoscrollInsets.left = rect.x + AUTOSCROLL_MARGIN;
    autoscrollInsets.bottom = size.height - (rect.y+rect.height) + AUTOSCROLL_MARGIN;
    autoscrollInsets.right  = size.width - (rect.x+rect.width) + AUTOSCROLL_MARGIN;
    return autoscrollInsets;
  }

  protected void setScale(int scale)
  {
    this.scale = scale;
  }

  protected int getScale()
  {
    return scale;
  }

  /**
   *
   * Popup listener
   *
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
