/* PopupMessageFrame
 *
 * created: 2009
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
package uk.ac.sanger.artemis.components.alignment;

import java.awt.Point;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JTextArea;

import uk.ac.sanger.artemis.components.SwingWorker;
import uk.ac.sanger.artemis.components.Utilities;

/**
 * Create undecorated popup frames.
 */
class PopupMessageFrame extends JFrame
{
  private static final long serialVersionUID = 1L;
  private JTextArea textArea = new JTextArea();

  PopupMessageFrame()
  {
    super();
    getRootPane().putClientProperty("Window.alpha", new Float(0.8f));
    textArea.setWrapStyleWord(true);
    textArea.setEditable(false);
    setUndecorated(true);
    getContentPane().add(textArea);
  }

  PopupMessageFrame(String msg)
  {
    this();
    textArea.setFont(textArea.getFont().deriveFont(14.f));
    textArea.setText(msg);
    setTitle(msg);
    pack();
  }

  protected void showWaiting(final String msg, final JPanel mainPanel)
  {
    textArea.setText(msg);
    pack();
    Point p = mainPanel.getLocationOnScreen();
    p.x += (mainPanel.getWidth() - PopupMessageFrame.this.getWidth()) / 2;
    p.y += mainPanel.getHeight() / 2;
    setLocation(p);
    SwingWorker worker = new SwingWorker()
    {
      public Object construct()
      {
        setVisible(true);
        return null;
      }
    };
    worker.start();
  }

  protected void show(final String title, final String msg, final int ypos,
      final long aliveTime)
  {
    SwingWorker worker = new SwingWorker()
    {
      public Object construct()
      {
        setTitle(title);
        textArea.setText(msg);
        pack();
        Utilities.centreJustifyFrame(PopupMessageFrame.this, Math.abs(ypos));
        setVisible(true);
        requestFocus();

        try
        {
          Thread.sleep(aliveTime);
        }
        catch (InterruptedException e)
        {
          e.printStackTrace();
        }
        setVisible(false);
        return null;
      }
    };
    worker.start();
  }
}
