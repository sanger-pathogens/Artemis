/*
 *
 * created: Wed Aug 3 2004
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2000  Genome Research Limited
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
 *
 * A simple, static class to display a URL in the system browser.
 *
 * Under Unix, the system browser is hard-coded to be 'netscape'.
 * Netscape must be in your PATH for this to work.  
 *
 * Under Windows, this will bring up the default browser under windows,
 * usually either Netscape or Microsoft IE.  The default browser is
 * determined by the OS.  
 *
 * Taken from: http://www.javaworld.com/javaworld/javatips/jw-javatip66_p.html
 *
 */

package uk.ac.sanger.artemis.editor;

import java.io.IOException;

public class BrowserControl
{

  // The default system browser under windows.
  private static final String WIN_PATH = "rundll32";
  // The flag to display a url.
  private static final String WIN_FLAG = "url.dll,FileProtocolHandler";
  // The default browser under unix.
  private static final String UNIX_PATH = "netscape";
  // The flag to display a url.
  private static final String UNIX_FLAG = "-remote openURL";

  /**
   * Display a file in the system browser.  If you want to display a
   * file, you must include the absolute path name.
   *
   * @param url the file's url (the url must start with either "http://"
   *        or "file://").
   */
  public static void displayURL(String url)
  {
    boolean windows = isWindowsPlatform();
    String cmd = null;
    try
    {
      if(windows)
      {
        // cmd = 'rundll32 url.dll,FileProtocolHandler http://...'
        cmd = WIN_PATH + " " + WIN_FLAG + " " + url;
        Process p = Runtime.getRuntime().exec(cmd);
      }
      else
      {
        // Netscape has to be running for the "-remote" command to work
        cmd = UNIX_PATH + " " + UNIX_FLAG + "(" + url + ")";
        Process p = Runtime.getRuntime().exec(cmd);
        try
        {
          // wait for exit code -- if it's 0, command worked,
          // otherwise we need to start the browser up.
          int exitCode = p.waitFor();
          if(exitCode != 0)
          {
            // Command failed, start up the browser
            cmd = UNIX_PATH + " "  + url;
            p = Runtime.getRuntime().exec(cmd);
          }
        }
        catch(InterruptedException x)
        {
          System.err.println("Error bringing up browser, cmd='" +
                             cmd + "'");
          System.err.println("Caught: " + x);
        }
      }
    }
    catch(IOException x)
    {
      // couldn't exec browser
      System.err.println("Could not invoke browser, command=" + cmd);
      System.err.println("Caught: " + x);
    }
  }

  /**
   * Try to determine whether this application is running under Windows
   * or some other platform by examing the "os.name" property.
   *
   * @return true if this application is running under a Windows OS
   */
  public static boolean isWindowsPlatform()
  {
    String os = System.getProperty("os.name");
    if ( os != null && os.startsWith("Windows"))
      return true;
    else
      return false;

  }
  /**
   * Simple example.
   */
  public static void main(String[] args)
  {
    displayURL("http://www.sanger.ac.uk");
  }
}
