/*
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
 */

package uk.ac.sanger.artemis.circular;

import java.awt.BorderLayout;
import java.awt.Dimension;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JProgressBar;

import uk.ac.sanger.artemis.components.Utilities;


class ProgressFrame extends JFrame
{
  private static final long serialVersionUID = 1L;
  private JProgressBar progress;
  
  public ProgressFrame()
  {
    setUndecorated(true);
    JPanel panel = (JPanel) getContentPane();
    progress = new JProgressBar(1,10);
    progress.setStringPainted(true);
    progress.setValue(2);
    progress.setPreferredSize(new Dimension(450, progress.getPreferredSize().height));
    panel.add(progress, BorderLayout.CENTER);
    pack();
    Utilities.centreFrame(this);
    setVisible(true);
  }
  
  protected void setString(String s)
  {
    toFront();
    progress.setString(s);
  }

  protected void setValue(int i)
  {
    progress.setValue(i);
  }
}