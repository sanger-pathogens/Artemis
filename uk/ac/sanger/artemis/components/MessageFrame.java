/* MessageFrame.java
 *
 * created: Mon Jan 18 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/MessageFrame.java,v 1.1 2004-06-09 09:47:08 tjc Exp $
 **/

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.Options;

import java.awt.*;
import java.awt.event.*;

import javax.swing.*;

/**
 *  A popup JFrame box that displays a message and has an OK JButton.
 *
 *  @author Kim Rutherford
 *  @version $Id: MessageFrame.java,v 1.1 2004-06-09 09:47:08 tjc Exp $
 **/

public class MessageFrame extends JFrame {
  /**
   *  Create a new MessageFrame component.
   *  @param message The message to display in the JDialog and it's title.
   **/
  public MessageFrame (final String message) {
    this (message, message);

    this.message = new JLabel (message);
  }

  /**
   *  Create a new MessageFrame component.
   *  @param title The title of the new dialog JFrame.
   *  @param message The message to display in the JDialog.
   **/
  public MessageFrame (final String title,
                       final String message) {
    super (title);

    this.message = new JLabel (message);
    
    getContentPane ().add (this.message, "North");

    final JPanel panel = new JPanel ();

    panel.add (ok_button);
    ok_button.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent e) {
        MessageFrame.this.dispose ();
      }
    });

    addWindowListener (new WindowAdapter () {
      public void windowClosing (WindowEvent event) {
        MessageFrame.this.dispose ();
      }
    });

    addKeyListener (new KeyAdapter () {
      public void keyTyped(final KeyEvent e) {
        MessageFrame.this.dispose ();
      }
    });
    
    getContentPane ().add (panel, "South");
    pack ();
    
    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
    setLocation (new Point ((screen.width - getSize ().width) / 2,
                            (screen.height - getSize ().height) / 2));

    setVisible (true);
  }

  final private JButton ok_button = new JButton ("OK");

  /**
   *  This is the message displayed above the OK button.  It can be set with
   *  setMessage ().
   **/
  private JLabel message;
}


