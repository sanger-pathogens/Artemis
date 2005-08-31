/* fILEtRANSFERpROGRESSmONITOR.JAVA
 *
 * created: Aug 2005
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
 **/

package uk.ac.sanger.artemis.j2ssh;

import javax.swing.JFrame;
import javax.swing.JProgressBar;
import java.awt.Component;
import java.awt.Point;
import javax.swing.Box;

import com.sshtools.j2ssh.FileTransferProgress;

import uk.ac.sanger.artemis.components.SwingWorker;

public class FileTransferProgressMonitor
{
  private FTProgress ftprogress;
  private JFrame progressFrame;
  private Box yBox = Box.createVerticalBox();

  public FileTransferProgressMonitor(Component parent)
  {
    progressFrame = new JFrame("Transfer");
    progressFrame.getContentPane().add(yBox);

    if(parent != null)
    {
      Point loc = parent.getLocationOnScreen();
      loc.x+=100;
      loc.y+=100;
      progressFrame.setLocation(loc);
    }

    yBox.add(Box.createHorizontalStrut(30));
    progressFrame.pack();
    progressFrame.setVisible(true);
  }

  public FTProgress add(String filename)
  {
    if(filename.length() > 25)
      filename = filename.substring(0,25)+"...";
    JProgressBar progress = new JProgressBar();
    progress.setStringPainted(true);
    progress.setString("    "+filename+" 0 %    ");
    yBox.add(progress);

    progressFrame.pack();
    progressFrame.validate();
    ftprogress = new FTProgress(progress, filename);
    return ftprogress;
  }

  public void close()
  { 
    progressFrame.setVisible(false);
    progressFrame.dispose();
  }
}


