/* MessageDialog.java
 *
 * created: Mon Dec 14 1998
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/MessageDialog.java,v 1.2 2008-10-15 13:56:42 tjc Exp $
 **/

package uk.ac.sanger.artemis.components;

import java.awt.Dimension;
import java.awt.Font;
import java.awt.Point;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import uk.ac.sanger.artemis.Options;

import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;

/**
 *  A popup dialog box that displays a message and has an OK JButton.
 *
 *  @author Kim Rutherford
 *  @version $Id: MessageDialog.java,v 1.2 2008-10-15 13:56:42 tjc Exp $
 **/

public class MessageDialog extends JDialog 
{
  private static final long serialVersionUID = 1L;

  /** Messages longer than this will be put in a TextArea rather than a Label. */
  final private int MESSAGE_SPLIT_SIZE = 100;

  final private JButton ok_button = new JButton ("OK");
  
  /**
   *  Create a blocking MessageDialog component.
   *  @param parent The parent window.
   *  @param message The message to display in the JDialog and to use as the
   *    frame title.
   **/
  public MessageDialog (final JFrame parent, final String message)
  {
    this (parent, message, message, true);
  }
  
  /**
   *  Create a new MessageDialog component.
   *  @param parent The parent window.
   *  @param message The message to display in the JDialog and to use as the
   *    frame title.
   *  @param modal If true, dialog blocks input to the parent window when
   *    shown.
   **/
  public MessageDialog (final JFrame parent, final String message,
                        final boolean modal) 
  {
    this (parent, message, message, modal);
  }
  
  /**
   *  Create a blocking MessageDialog component.
   *  @param parent_frame The parent window.
   *  @param title The title of the new dialog JFrame.
   *  @param message The message to display in the JDialog.
   **/
  public MessageDialog (final JFrame parent_frame,
                        final String title,
                        final String message) 
  {
    this (parent_frame, title, message, true);
  }

  /**
   *  Create a new MessageDialog component.
   *  @param parent_frame The parent window.
   *  @param title The title of the new dialog JFrame.
   *  @param message The message to display in the JDialog.
   *  @param modal If true, dialog blocks input to the parent window when
   *    shown.
   **/
  public MessageDialog (final JFrame parent_frame,
                        final String title,
                        final String message,
                        final boolean modal) 
  {
    super (parent_frame, title, modal);

    final Font font = Options.getOptions ().getFont ();
    setFont (font);

    if (message.length () < MESSAGE_SPLIT_SIZE) 
      getContentPane ().add (new JLabel (message), "North"); 
    else 
    {
      final JTextArea text_area = new JTextArea (18, 90);
      text_area.setText (message);

      getContentPane().add (new JScrollPane (text_area), "North");
      text_area.setEditable (false);
    }

    final JPanel panel = new JPanel ();
    panel.add (ok_button);
    ok_button.addActionListener (new ActionListener () 
    {
      public void actionPerformed (ActionEvent e) 
      {
        MessageDialog.this.dispose ();
      }
    });

    addWindowListener (new WindowAdapter () 
    {
      public void windowClosing (WindowEvent event) 
      {
        MessageDialog.this.dispose ();
      }
    });

    addKeyListener (new KeyAdapter () 
    {
      public void keyTyped(final KeyEvent e) 
      {
        MessageDialog.this.dispose ();
      }
    });
    
    getContentPane ().add (panel, "South");
    pack ();

    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
    setLocation (new Point ((screen.width - getSize ().width) / 2,
                            (screen.height - getSize ().height) / 2));
    setVisible (true);
  }
}

