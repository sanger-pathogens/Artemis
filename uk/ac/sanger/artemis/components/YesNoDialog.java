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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/YesNoDialog.java,v 1.1 2004-06-09 09:48:01 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.Options;

import java.awt.*;
import java.awt.event.*;

import javax.swing.*;

/**
 *  A popup dialog box that displays a message and then waits for the to press
 *  yes or no.
 *
 *  @author Kim Rutherford
 *  @version $Id: YesNoDialog.java,v 1.1 2004-06-09 09:48:01 tjc Exp $
 **/

public class YesNoDialog extends JDialog {
  /**
   *  Create a new YesNoDialog component.  The constructor does not show ()
   *  the dialog, call getResult () to do that.
   *  @param parent The parent window.
   *  @param title The title of the new dialog JFrame.
   *  @param message The message to display in the JDialog.
   **/
  public YesNoDialog (JFrame parent, String title, String message) {
    super (parent, message, true);
    
    getContentPane ().add (new JLabel (message), "North");

    final JPanel panel = new JPanel ();

    panel.add (yes_button);
    yes_button.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent e) {
        button_result = true;
        setVisible (false);
        YesNoDialog.this.dispose ();
      }
    });

    panel.add (no_button);
    no_button.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent e) {
        button_result = false;
        setVisible (false);
        YesNoDialog.this.dispose ();
      }
    });

    getContentPane ().add (panel, "South");
    pack ();
        
    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
    setLocation (new Point ((screen.width - getSize ().width) / 2,
                            (screen.height - getSize ().height) / 2));
  }

  /**
   *  Create a new YesNoDialog component.  The constructor does not show ()
   *  the dialog, call getResult () to do that.
   *  @param parent The parent window.
   *  @param message The message to display in the JDialog and to uise as the
   *    title string.
   **/
  public YesNoDialog (JFrame parent, String message) {
    this (parent, message, message);
  }
  
  /**
   *  This method calls show () on this object, and then waits for the user to
   *  press the Yes button or the No button.
   *  @return true is the user pressed Yes, false otherwise.
   **/
  public boolean getResult () {
    setVisible (true);
    
    return button_result;
  }

  final JButton yes_button = new JButton ("Yes");

  final JButton no_button = new JButton ("No");

  /**
   *  Set by the yes and no button action listeners and read by getResult ().
   **/
  boolean button_result;
}

