/* AbstractGraphPanel
 *
 * created: 2010
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2010  Genome Research Limited
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
package uk.ac.sanger.artemis.components.alignment;

import java.awt.Color;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionListener;
import java.text.DecimalFormat;

import javax.swing.BorderFactory;
import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JTextField;

import uk.ac.sanger.artemis.io.Range;

public class AbstractGraphPanel extends JPanel
{
  private static final long serialVersionUID = 1L;
  protected int start;
  protected int end;
  protected float pixPerBase;
  protected int nBins;
  
  protected BamView bamView;
  protected int windowSize;
  protected int max;
  protected boolean autoWinSize = true;
  protected int userWinSize = 1;
  protected JPopupMenu popup = new JPopupMenu();
  
  public AbstractGraphPanel()
  {
    super();
    setBorder(BorderFactory.createMatteBorder(0, 0, 1, 0, Color.gray));
    setBackground(Color.white);

    final MouseMotionListener mouseMotionListener = new MouseMotionListener()
    {
      public void mouseDragged(MouseEvent e)
      {
        if(e.getButton() == MouseEvent.BUTTON3) 
          return;

        if(e.getClickCount() > 1)
        {
          bamView.getSelection().clear();
          repaint();
          return;  
        }
        bamView.highlightRange(e, 
            MouseEvent.BUTTON1_DOWN_MASK & MouseEvent.BUTTON2_DOWN_MASK);
      }

      public void mouseMoved(MouseEvent e){}
    };
    addMouseMotionListener(mouseMotionListener);
  }
  
  protected void initPopupMenu(JComponent menu)
  {
    final JMenuItem setScale = new JMenuItem("Set the Window Size...");
    setScale.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event)
      {
        GridBagLayout gridbag = new GridBagLayout();
        GridBagConstraints c = new GridBagConstraints();
        JPanel pane = new JPanel(gridbag);
        final JTextField newWinSize = new JTextField(Integer.toString(userWinSize), 10);
        newWinSize.setEnabled(!autoWinSize);
        final JLabel lab = new JLabel("Window size:");
        c.gridy = 0;
        pane.add(lab, c);
        pane.add(newWinSize, c);

        final JCheckBox autoSet = new JCheckBox("Automatically set window size", autoWinSize);
        autoSet.addActionListener(new ActionListener()
        {
          public void actionPerformed(ActionEvent e)
          {
            lab.setEnabled(!autoSet.isSelected());
            newWinSize.setEnabled(!autoSet.isSelected());
          }
        });
        c.gridy = 1;
        c.gridwidth = GridBagConstraints.REMAINDER;
        pane.add(autoSet, c);

        String window_options[] = { "OK", "Cancel" };
        int select = JOptionPane.showOptionDialog(null, pane, "Window Size",
            JOptionPane.DEFAULT_OPTION, JOptionPane.QUESTION_MESSAGE, null,
            window_options, window_options[0]);

        if (select == 1)
          return;
        autoWinSize = autoSet.isSelected();
        try
        {
          userWinSize = Integer.parseInt(newWinSize.getText().trim());
        }
        catch (NumberFormatException nfe)
        {
          return;
        }
        bamView.repaint();
      }
    });
    menu.add(setScale);
    
    addMouseListener(new PopupListener());
  }
  
  /**
   * Draw maximum average value.
   * @param g2
   */
  protected void drawMax(Graphics2D g2)
  {
    drawMax(g2, (float)max/(float)windowSize);
  }
  
  /**
   * Draw maximum average value.
   * @param g2
   */
  protected void drawMax(Graphics2D g2, float max)
  {
    DecimalFormat df = new DecimalFormat("0.#");
    String maxStr = df.format(max);

    FontMetrics fm = getFontMetrics(getFont());
    g2.setColor(Color.black);
    
    int xpos = bamView.getJspView().getVisibleRect().width - fm.stringWidth(maxStr) - 
               bamView.getJspView().getVerticalScrollBar().getWidth();
    g2.drawString(maxStr, xpos, fm.getHeight());
  }
  
  /**
   * Return the log value if the log scale is selected
   * @param val
   * @return
   */
  protected float getValue(int val, boolean logScale)
  {
    if(val == 0)
      return 0.f;
    return (float) (logScale ? Math.log(val) : val);
  }
  
  protected void drawSelectionRange(final Graphics2D g2,
                                    final float pixPerBase, 
                                    final int start,
                                    final int end,
                                    final int hgt,
                                    final Color c)
  {
    if(bamView != null && bamView.getSelection() != null)
    {
      final Range selectedRange = bamView.getSelection().getSelectionRange();
      if(selectedRange != null)
      {
        int rangeStart = selectedRange.getStart();
        int rangeEnd   = selectedRange.getEnd();

        if(end < rangeStart || start > rangeEnd)
          return;

        int x = (int) (pixPerBase*(rangeStart-bamView.getBaseAtStartOfView()));
        int width = (int) (pixPerBase*(rangeEnd-rangeStart+1));

        g2.setColor(c);
        g2.fillRect(x, 0, width, hgt);
      }
    }
  }

  protected void setStartAndEnd(int start, int end)
  {
    this.start = start;
    this.end = end;
  }

  protected void setPixPerBase(float pixPerBase)
  {
    this.pixPerBase = pixPerBase;
  }
  
  /**
   * Popup menu listener
   */
   class PopupListener extends MouseAdapter
   {
     JMenuItem gotoMateMenuItem;
     JMenuItem showDetails;
     
     public void mouseClicked(MouseEvent e)
     {
       if(e.getClickCount() > 1)
         bamView.getSelection().clear(); 
     }
     
     public void mousePressed(MouseEvent e)
     {
       maybeShowPopup(e);
     }

     public void mouseReleased(MouseEvent e)
     {
       maybeShowPopup(e);
       bamView.dragStart = -1;
     }

     private void maybeShowPopup(MouseEvent e)
     {
       if(e.isPopupTrigger())
       {
         popup.show(e.getComponent(),
                 e.getX(), e.getY());
       }
     }
   }
}