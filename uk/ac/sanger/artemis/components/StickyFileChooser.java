/* StickyFileChooser.java
 *
 * created: Mon Sep  1 2003
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 2003  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/StickyFileChooser.java,v 1.1 2004-06-09 09:47:49 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.Options;
import java.io.File;
import javax.swing.JFileChooser;
import javax.swing.JFrame;

/**
 *  A JFileChooser that remembers which directory it is in for next time.
 *  @author Kim Rutherford
 **/

public class StickyFileChooser extends JFileChooser 
{
  private static final long serialVersionUID = 1L;
  /**
   *  Used to remember the directory the JFileChooser was in when the user
   *  pressed OK.  This is used as the starting directory next time.
   **/
  private static File last_directory = null;

  /**
   *  Create a new StickyFileChooser and set the current directory to what is
   *  was after the last call to StickyFileChooser.showOpenDialog() or
   *  StickyFileChooser.showSaveDialog().
   **/
  public StickyFileChooser() 
  {
    super();
    setName("StickyFileChooser");

    if (last_directory == null) 
    {
      if (Options.getOptions ().getProperty ("default_directory") != null) 
      {
        final String default_directory =
          Options.getOptions().getProperty("default_directory");
        setCurrentDirectory (new File (default_directory));
      } 
      else 
      {
        if (last_directory == null) 
          last_directory = new File (System.getProperty ("user.dir"));
      }
    }

    setCurrentDirectory(last_directory);
    setSize(620, 460);
  }

  /**
   *  Calls super.showOpenDialog() then remembers the current directory.
   **/
  public int showOpenDialog(JFrame owner) 
  {
    int status = super.showOpenDialog(owner);
    last_directory = getCurrentDirectory ();
    return status;
  }

  /**
   *  Calls super.showSaveDialog() then remembers the current directory.
   **/
  public int showSaveDialog(JFrame owner) 
  {
    int status = super.showSaveDialog (owner);
    last_directory = getCurrentDirectory ();
    return status;
  }

}
