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
 **/

package uk.ac.sanger.artemis.components;

import java.awt.Dimension;
import java.awt.Point;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;

/**
 *  A popup JFrame box that displays a message and has an OK JButton.
 *  @author Kim Rutherford
 **/
public class MessageFrame extends JFrame {
  private static final long serialVersionUID = 1L;

  /**
   *  Create a new MessageFrame component.
   *  @param message The message to display in the JDialog and it's title.
   **/
  public MessageFrame (final String message) {
    this (message, message);
  }

  /**
   *  Create a new MessageFrame component.
   *  @param title The title of the new dialog JFrame.
   *  @param message The message to display in the JDialog.
   **/
  protected MessageFrame (final String title,
                          final String message) {
    super (title);
    setDefaultCloseOperation(DISPOSE_ON_CLOSE);
   
    getContentPane ().add (new JLabel (message), "North");
    
    final JButton ok_button = new JButton ("OK");
    ok_button.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent e) {
        MessageFrame.this.dispose ();
      }
    });

    addKeyListener (new KeyAdapter () {
      public void keyTyped(final KeyEvent e) {
        MessageFrame.this.dispose ();
      }
    });
    
    getContentPane ().add (ok_button, "South");
    pack ();
    
    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
    setLocation (new Point ((screen.width - getSize ().width) / 2,
                            (screen.height - getSize ().height) / 2));
    setVisible (true);
  }
}
