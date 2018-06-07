/* OpenSectionButton
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2009  Genome Research Limited
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

package uk.ac.sanger.artemis.components.genebuilder;

import java.awt.Dimension;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JComponent;

  public class OpenSectionButton extends JButton 
  {
    private static final long serialVersionUID = 1L;
    private JComponent panel;
    
    public OpenSectionButton(final String text, final JComponent panel)
    {
      super(text);
      this.panel = panel;

      Dimension size = new Dimension(55,20);
      setForeground(GeneEditorPanel.STEEL_BLUE);
      Font font = getFont();
      font = font.deriveFont(Font.BOLD);
      setFont(font);
      setBorderPainted(false);
      setOpaque(false);
      setPreferredSize(size);
      setMaximumSize(size);

      addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          updateButton();
        }  
      });
    }
    
    private void updateButton()
    {
      if(getText().equals("-"))
      {
        setText("+");
        panel.setVisible(false);
      }
      else
      {
        setText("-");
        panel.setVisible(true);
      }
    }
    
    public void setOpen(final boolean isOpen)
    {
      if(isOpen && getText().equals("-"))
        return;
      if(!isOpen && getText().equals("+"))
        return;
      updateButton();
    }

  }