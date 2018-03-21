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

import com.sshtools.j2ssh.FileTransferProgress;
import javax.swing.JProgressBar;

public class FTProgress implements FileTransferProgress
{
  private JProgressBar progress;
  private long bytesTotal;
  private String filename;

  public FTProgress(JProgressBar progress, String filename)
  {
    this.progress = progress;
    this.filename = filename;
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

  public void completed()
  {
  }

  public long getProgress()
  {
    return progress.getValue();
  }

}

