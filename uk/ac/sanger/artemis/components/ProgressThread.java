/* ArtemisMain.java
 *
 * created: Tue May 11 2004
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2004  Genome Research Limited
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

package uk.ac.sanger.artemis.components;

import javax.swing.*;
import java.awt.Dimension;
import java.awt.Color;

public class ProgressThread extends Thread
{
  private JFrame frame;
  private JFrame progress_frame;
  private String msg;
  private JProgressBar progressBar = new JProgressBar();

  public ProgressThread(JFrame frame, String msg)
  {
    this.frame = frame;
    this.msg   = msg;
  }

  public void run()
  {
    try
    {
      progress_frame = new JFrame("Loading...");
      Dimension d = progress_frame.getToolkit().getScreenSize();
      progressBar.setIndeterminate(true);
      progressBar.setBackground(Color.white);
      progress_frame.getContentPane().add(progressBar);
      progress_frame.pack();
      progress_frame.setLocation(
              ((int)d.getWidth()-progress_frame.getWidth())/2,
              ((int)d.getHeight()-progress_frame.getHeight())/2);
      progress_frame.setVisible(true);
    }
    catch(NoSuchMethodError nsme){}
  }

  public void finished()
  {
    if(progress_frame != null)
      progress_frame.setVisible(false);
  }

}
