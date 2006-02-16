/* MacHandler.java
 *
 * created: May 2005
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2005  Genome Research Limited
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

import com.apple.mrj.*;
import javax.swing.SwingUtilities;
import javax.swing.JOptionPane;

class MacHandler
	implements MRJQuitHandler, MRJPrefsHandler, MRJAboutHandler 
{
  private Splash us;
  private String title;

  public MacHandler(Splash theProgram, String title) 
  {
    us = theProgram;
    this.title = title;
//  System.setProperty("com.apple.mrj.application.apple.menu.about.name", title);
    System.setProperty("apple.laf.useScreenMenuBar", "true");
    MRJApplicationUtils.registerAboutHandler(this);
    MRJApplicationUtils.registerPrefsHandler(this);
    MRJApplicationUtils.registerQuitHandler(this);
  }

  public void handleAbout() 
  {
    us.about();
  }

  public void handlePrefs() 
  {
//  us.prefs();
  }

  public void handleQuit() 
  {
    SwingUtilities.invokeLater(new Runnable() 
    {
      public void run() 
      {
        int select = JOptionPane.showConfirmDialog(null,
                          "Quit "+title+"?", "Quit", 
                          JOptionPane.YES_NO_OPTION,
                          JOptionPane.QUESTION_MESSAGE);
        if(select == JOptionPane.YES_OPTION)
          us.exitApp();
      }
    });
    throw new IllegalStateException("Let the quit handler do it");
  }

  

}

