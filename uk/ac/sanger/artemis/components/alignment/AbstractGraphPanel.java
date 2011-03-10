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

public class AbstractGraphPanel extends JPanel
{
  private static final long serialVersionUID = 1L;
  protected int start;
  protected int end;
  protected float pixPerBase;
  
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
  }
  
  protected void initPopupMenu(JComponent menu)
  {
    final JMenuItem setScale = new JMenuItem("Set the Window Size...");
    setScale.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent _)
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
    DecimalFormat df = new DecimalFormat("0.#");
    String maxStr = df.format((float)max/(float)windowSize);

    FontMetrics fm = getFontMetrics(getFont());
    g2.setColor(Color.black);
    
    int xpos = bamView.getJspView().getVisibleRect().width - fm.stringWidth(maxStr) - 
               bamView.getJspView().getVerticalScrollBar().getWidth();
    g2.drawString(maxStr, xpos, fm.getHeight());
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
     }
     
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
         popup.show(e.getComponent(),
                 e.getX(), e.getY());
       }
     }
   }
}