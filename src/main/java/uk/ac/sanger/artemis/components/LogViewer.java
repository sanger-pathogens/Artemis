/* LogViewer.java
 *
 * created: Wed Aug 30 2000
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 2000  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/LogViewer.java,v 1.3 2007-03-22 13:41:31 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;
import org.apache.log4j.AppenderSkeleton;
import org.apache.log4j.spi.LoggingEvent;

import java.io.*;

/**
 *  A class for viewing log messages in a FileViewer component.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: LogViewer.java,v 1.3 2007-03-22 13:41:31 tjc Exp $
 **/

public class LogViewer extends AppenderSkeleton implements Logger 
{

  public int maxLogLines;
  /** The FileViewer that is used to show the messages. */
  private FileViewer file_viewer = null;

  /**
   *  Create a new, empty LogViewer component.
   **/
  public LogViewer() 
  {

  }

  /**
   *  Send the given String to the log.
   **/
  public void log(final String message) 
  {
    maybeMakeViewer();
    file_viewer.appendString(message);
  }

  /**
   *  Read from the given Reader and send it to the log.
   **/
  public void log(final Reader reader)
      throws IOException 
  {
    maybeMakeViewer();
    file_viewer.appendFile(reader);
  }

  /**
   *  Call setVisible() on the FileViewer that this object contains.
   **/
  public void setVisible(final boolean flag) 
  {
    maybeMakeViewer();
    file_viewer.setVisible(flag);
  }

  /**
   *  If file_viewer is null make a FileViewer and assign it to file_viewer.
   **/
  private void maybeMakeViewer() 
  {
    if (file_viewer == null) 
    {
      file_viewer = new FileViewer("Log Viewer", false) 
      {
        /** */
        private static final long serialVersionUID = 1L;

        public void dispose()
        {
          // if the FileViewer is deleted we want to know
          file_viewer = null;
          super.dispose();
        }
      };

      file_viewer.pack();
    }
  }

  protected void append(LoggingEvent e)
  {
    String message = this.layout.format(e);
    FileViewer fv = ((LogViewer)Splash.getLogger()).getFileViewer();
    if(fv  != null &&
       maxLogLines < fv.getLineCount())
      fv.setText("");

    Splash.getLogger().log(message);
  }

  public void close()
  {
    
  }

  public boolean requiresLayout()
  {
    return true;
  }

  public int getMaxLogLines()
  {
    return maxLogLines;
  }

  /**
   * Set in log4j.properties by log4j.appender.R.MaxLogLines
   * @param maxLogLines
   */
  public void setMaxLogLines(int maxLogLines)
  {
    this.maxLogLines = maxLogLines;
  }
  
  public FileViewer getFileViewer()
  {
    return this.file_viewer;
  }
 

}
