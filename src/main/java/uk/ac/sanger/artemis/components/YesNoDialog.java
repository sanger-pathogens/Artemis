/* YesNoDialog.java
 *
 * created: Sat Dec 12 1998
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 1998,1999,2000  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/YesNoDialog.java,v 1.2 2008-10-15 13:45:58 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import javax.swing.JFrame;
import javax.swing.JOptionPane;

/**
 *  A popup dialog box that displays a message and then waits for the to press
 *  yes or no.
 *
 *  @author Kim Rutherford
 *  @version $Id: YesNoDialog.java,v 1.2 2008-10-15 13:45:58 tjc Exp $
 **/

public class YesNoDialog
{
  private JFrame parent;
  private String title;
  private String message;
  
  /**
   *  Create a new YesNoDialog component.  The constructor does not show ()
   *  the dialog, call getResult () to do that.
   *  @param parent The parent window.
   *  @param title The title of the new dialog JFrame.
   *  @param message The message to display in the JDialog.
   **/
  public YesNoDialog (JFrame parent, String title, String message) 
  {
    this.parent = parent;
    this.title  = title;
    this.message = message;
  }

  /**
   *  Create a new YesNoDialog component.  The constructor does not show ()
   *  the dialog, call getResult () to do that.
   *  @param parent The parent window.
   *  @param message The message to display in the JDialog and to uise as the
   *    title string.
   **/
  public YesNoDialog (JFrame parent, String message) 
  {
    this (parent, null, message);
  }
  
  /**
   *  This method calls show () on this object, and then waits for the user to
   *  press the Yes button or the No button.
   *  @return true is the user pressed Yes, false otherwise.
   **/
  protected boolean getResult () 
  {
    int result = JOptionPane.showConfirmDialog(
        parent, message, title, 
        JOptionPane.YES_NO_OPTION);
    if(result == JOptionPane.YES_OPTION)
      return true;
    return false;
  }
}

