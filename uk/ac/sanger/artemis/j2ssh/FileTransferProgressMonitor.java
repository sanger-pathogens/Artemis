/* FileTransferProgressMonitor.java
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

import com.sshtools.j2ssh.FileTransferProgress;

import uk.ac.sanger.artemis.components.SwingWorker;

public class FileTransferProgressMonitor
              implements FileTransferProgress
{

  private JProgressBar progress;
  private JFrame progressFrame;
  private long bytesTotal;
  private String filename;

  public FileTransferProgressMonitor(Component parent, final String filename)
  {
    this.filename = filename;

    if(filename.length() > 25)
      this.filename = filename.substring(0,25)+"...";

    progress = new JProgressBar();
    progress.setStringPainted(true);
    progress.setString("    "+filename+" 0 %    ");
    progressFrame = new JFrame("Transfer");
    progressFrame.getContentPane().add(progress);
    progressFrame.pack();

    if(parent != null)
    {
      Point loc = parent.getLocationOnScreen();
      loc.x+=100;
      loc.y+=100;
      progressFrame.setLocation(loc);
    }
      
    progressFrame.setVisible(true);
  }

  public void completed()
  { 
    progressFrame.setVisible(false);
    progressFrame.dispose();
  }

  public boolean isCancelled()
  {
    return false;
  }

  public void progressed(long bytesSoFar)
  {
    progress.setValue((new Long(bytesSoFar)).intValue());
    int percent = (int)((bytesSoFar*100)/bytesTotal);
    progress.setString(filename+" "+percent+" %");
  }

  public void started(final long bytesTotal, final String remoteFile)
  {
    this.bytesTotal = bytesTotal;
    progress.setMinimum(0);
    progress.setMaximum((new Long(bytesTotal)).intValue());
  }
}


