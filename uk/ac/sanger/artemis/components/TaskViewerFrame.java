/* TaskViewerFrame.java
 *
 * created: Mon Aug 18 2003
 *
 * This file is part of Artemis
 * 
 * Copyright(C) 2003  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/TaskViewerFrame.java,v 1.1 2004-06-09 09:47:51 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.jcon.gui.TaskViewer;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

/**
 *  TaskViewerFrame class
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: TaskViewerFrame.java,v 1.1 2004-06-09 09:47:51 tjc Exp $
 **/

public class TaskViewerFrame extends JFrame 
{

  /** TaskViewer that was created in the constructor. */
  private TaskViewer tv = null;


  /**
   *  Create a new TaskViewerFrame to monitor the tasks with the given ids.
   **/
  public TaskViewerFrame(final int [] ids) 
  {
    try 
    {
      tv = new TaskViewer(ids);
    }
    catch(Exception e) 
    {
      dispose();
      new MessageDialog(this, "Exception when viewing tasks: " + 
                         e.getMessage());
      return;
    }

    final JPanel panel = new JPanel();
    getContentPane().add(panel);
    panel.setLayout(new BorderLayout());
    panel.add(tv, "Center");

    final JButton close_button = new JButton("Close");
    close_button.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent e) 
      {
        dispose();
      }
    });

    panel.add(close_button, "South");
    pack();

    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();

    setLocation(new Point((screen.width - getSize().width) / 2,
                           (screen.height - getSize().height) / 2));
  }

}
