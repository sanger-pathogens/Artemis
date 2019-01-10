/* RemoveButton
 * This file is part of Artemis
 *
 * Copyright (C) 2009  Genome Research Limited
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
 */

package uk.ac.sanger.artemis.components.genebuilder.gff;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.event.ActionListener;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.SwingConstants;

class RemoveButton extends JButton
{
  private static final long serialVersionUID = 1L;

  public RemoveButton(ActionListener removeAction)
  {
    super("X");

    setVerticalTextPosition(SwingConstants.TOP);
    setHorizontalTextPosition(SwingConstants.LEFT);

    setForeground(Color.red);
    setFont(getFont().deriveFont(Font.BOLD, 9.f));
    setBorder(BorderFactory.createEmptyBorder());
    setOpaque(false);

    Dimension size = new Dimension(9, 20);
    setPreferredSize(size);
    setMaximumSize(size);

    addActionListener(removeAction);
  }
}