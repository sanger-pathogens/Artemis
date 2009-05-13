/* ExternalProgram.java
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
 **/
package uk.ac.sanger.artemis.components;

import java.awt.Color;

import javax.swing.JFrame;
import javax.swing.JProgressBar;

public class ProgressBarFrame extends JFrame
{
  private static final long serialVersionUID = 1L;

  public ProgressBarFrame(int seconds, String name)
  {
    super();
    setUndecorated(true);
    final int max = (seconds*1000)/40;
    final JProgressBar progressBar = new JProgressBar(0,max);
    progressBar.setStringPainted(true);
    progressBar.setString("Sending "+name+" process now!");
    progressBar.setBackground(Color.white);

    SwingWorker batchWorker = new SwingWorker()
    {
      public Object construct()
      {
        try
        {
          for(int i=0; i<max; i++)
          {
            Thread.sleep(40);
            progressBar.setValue(i);
          }
          dispose();
        }
        catch(InterruptedException intr){}
        return null;
      }
    };

    getContentPane().add(progressBar);
    pack();
    Utilities.centreFrame(this);
    setVisible(true);
    batchWorker.start();
  }
}