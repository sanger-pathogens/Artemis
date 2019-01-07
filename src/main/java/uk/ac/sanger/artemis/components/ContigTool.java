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
import uk.ac.sanger.artemis.SelectionChangeListener;
import uk.ac.sanger.artemis.SelectionChangeEvent;
import uk.ac.sanger.artemis.Selection;

import java.awt.*;
import java.awt.event.*;

import java.util.Vector;

import java.awt.geom.RoundRectangle2D;

import javax.swing.BorderFactory;
import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JScrollPane;
import javax.swing.border.Border;
import java.awt.datatransfer.*;
import java.awt.dnd.*;

public class ContigTool extends JPanel
       implements DragGestureListener, DropTargetListener,
                  DragSourceListener, Autoscroll,
                  SelectionChangeListener
{
  private static final long serialVersionUID = 1L;
  private FeatureVector contig_features;
  private FeatureDisplay feature_display;
  private Selection selection;

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
                    final JScrollPane jsp,
                    final Selection selection)
  {
    super();
    this.contig_features = contig_features;
    this.feature_display = feature_display;
    this.selection       = selection;

    setFocusable(true); // required for KeyEvent to work
    MultiLineToolTipUI.initialize();
    setToolTipText("");   //enable tooltip display

    /*for(int i=0; i<contig_features.size(); i++)
    {
      final Range this_feature_range = contig_features.elementAt(i).getMaxRawRange();
      length += ((this_feature_range.getEnd() - this_feature_range.getStart())/scale);
    }*/
   
    length += contig_features.elementAt(0).getStrand().getSequenceLength() / scale;
    
    
    Dimension dim = new Dimension(length, 20);
    setPreferredSize(dim);

    DragSource dragSource = DragSource.getDefaultDragSource();

    dragSource.createDefaultDragGestureRecognizer(
       this,                             // component where drag originates
       DnDConstants.ACTION_COPY_OR_MOVE, // actions
       this);                            // drag gesture recognizer

    setDropTarget(new DropTarget(this,this));

    getSelection().addSelectionChangeListener(this);

    addKeyListener(new KeyAdapter()
    {
      public void keyPressed(final KeyEvent event)
      {
        // this is done so that menu shortcuts don't cause each action to be
        // performed twice
        if(event.getModifiers() != 0)
          return;

        switch(event.getKeyCode())
        {
          case KeyEvent.VK_UP:
            goToNext(true);
            repaint();
            break;
          case KeyEvent.VK_DOWN:
            goToNext(false);
            repaint();
            break;
          default:
            break;
        }
      }
    });


    addMouseListener(new MouseAdapter()
    {
      public void mouseReleased(MouseEvent event)
      {
        if(event.isPopupTrigger())
        {
          popup.show(event.getComponent(),
                  event.getX(), event.getY());
          return;
        }

        FeatureVector contig_features = ContigTool.this.contig_features;
        if(event.getClickCount() == 1 &&
           event.getID() == MouseEvent.MOUSE_RELEASED)
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
              if(getSelection().contains(feature))
                getSelection().remove(feature);
              else
              {
                clearSelection();

                String tt = this_feature_range.getStart()+".."+
                            this_feature_range.getEnd();
 
                if(feature.getIDString() != null)
                  tt = tt + ", " + feature.getIDString();

                status_line.setText(tt);

                getSelection().add(feature);
              }
            }
          }
          repaint();
        }
      }
   
      public void mousePressed(MouseEvent event)
      {
        if(event.isPopupTrigger())
        {
          popup.show(event.getComponent(),
                  event.getX(), event.getY());
          return;
        }
      }
    });

    JMenu zoomIn = new JMenu("Zoom In");
    popup.add(zoomIn);
    // set popup menu items
    JMenuItem zoomIn5 = new JMenuItem("x 1/5");
    zoomIn.add(zoomIn5);
    zoomIn5.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)
      {
        zoomIn(5,jsp);
      }
    });
    
    JMenuItem zoomIn10 = new JMenuItem("x 1/10");
    zoomIn.add(zoomIn10);
    zoomIn10.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)
      {
        zoomIn(10,jsp);
      }
    });

    
    JMenu zoomOut = new JMenu("Zoom Out");
    popup.add(zoomOut);
    JMenuItem zoomOut5 = new JMenuItem("x5");
    zoomOut.add(zoomOut5);
    zoomOut5.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)
      {
        scale = scale * 5;
        adjustSize(jsp);
        repaint();
      }
    });
    
    JMenuItem zoomOut10 = new JMenuItem("x10");
    zoomOut.add(zoomOut10);
    zoomOut10.addActionListener(new ActionListener()
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

  private void zoomIn(final int factor, final JScrollPane jsp)
  {
    if(scale < factor)
    {
      scale = 1; 
      return;
    }

    scale = scale / factor;  
    adjustSize(jsp);
    repaint();  
  }
  
  /**
   *  
   * Clear all selected features
   *
   **/
  private void clearSelection()
  {
    FeatureVector selected_features = getSelection().getAllFeatures();

    for(int i=0; i<selected_features.size(); i++)
    {
      Feature this_feature = selected_features.elementAt(i);
      getSelection().remove(this_feature);
    }
  }

  private Selection getSelection()
  {
    return selection;
  }

  /**
   *  Implementation of the SelectionChangeListener interface.  We listen to
   *  SelectionChange events so that we can update the list to reflect the
   *  current selection.
   **/
  public void selectionChanged(SelectionChangeEvent event)
  {
    if(!isVisible())
      return;

    // don't bother with events we sent ourself
    if(event.getSource() == this)
      return;

    // if the selected range changes we don't care
    if(getSelection().getMarkerRange() != null &&
       event.getType() == SelectionChangeEvent.OBJECT_CHANGED)
      return;

    repaint();
  }

  /**
   * 
   * Select the next feature in the sequence.
   *
   */
  private void goToNext(boolean up)
  {
    if(getSelection().getSelectedFeatures().size() != 1)
      return;

    final Feature curr_feature = getSelection().getSelectedFeatures().elementAt(0);
    final Range curr_feature_range = curr_feature.getMaxRawRange();
    int start = curr_feature_range.getStart();
    int end   = curr_feature_range.getEnd();

    if(up && start == 1)
      return;

    for(int i=0; i<contig_features.size(); i++)
    {
      final Feature feature = contig_features.elementAt(i);
      final Range this_feature_range = feature.getMaxRawRange();
      if(up && this_feature_range.getEnd() == start-1)
      {
        clearSelection();
        getSelection().add(feature);
        return;
      }
      else if(!up && this_feature_range.getStart() == end+1)
      {
        clearSelection();
        getSelection().add(feature);
        return;
      }
    }
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
    length += contig_features.elementAt(0).getStrand().getSequenceLength() / scale;
    setPreferredSize(new Dimension(length, 60));
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
      if(colour == null)
        colour = Color.white;

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

    FeatureVector selected_features = getSelection().getSelectedFeatures();
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
      final Shape saved_clip = g.getClip();
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
    final Vector contig_keys = FeatureDisplay.getContigKeys();

    for(int i = 0; i < contig_features.size(); i++)
    {
      final Feature this_feature = contig_features.elementAt(i);

      if(contig_keys.contains(this_feature.getKey()))
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
    //Transferable t = e.getTransferable();
    if(e.isDataFlavorSupported(DataFlavor.stringFlavor))
    {
      FeatureVector selected_features = getSelection().getSelectedFeatures();
      feature_display.reorder(highlight_drop_base, 
                              selected_features.elementAt(0));   // rearrange contigs

      // reset status bar
      final Range this_feature_range = selected_features.elementAt(0).getMaxRawRange();
      String tt = this_feature_range.getStart()+".."+
                  this_feature_range.getEnd();

      if(selected_features.elementAt(0).getIDString() != null)
        tt = tt + ", " + selected_features.elementAt(0).getIDString();

      status_line.setText(tt);
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
    else
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

    final Vector contig_keys = FeatureDisplay.getContigKeys();

    FeatureVector selected_features = getSelection().getSelectedFeatures();
    if(selected_features.size() == 1 &&
       contig_keys.contains(selected_features.elementAt(0).getKey()))
    {
      ClassLoader cl = this.getClass().getClassLoader();
      ImageIcon icon = new ImageIcon(cl.getResource("images/icon.gif"));
      final Image icon_image = icon.getImage();

      //TransferableContig tcontig = new TransferableContig(selected_features.elementAt(0));
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
//  highlight_drop_base = -1;
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
