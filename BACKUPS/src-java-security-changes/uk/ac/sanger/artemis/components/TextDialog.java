/* TextDialog.java
 *
 * created: Mon Jan 11 1999
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/TextDialog.java,v 1.1 2004-06-09 09:47:52 tjc Exp $
 **/

package uk.ac.sanger.artemis.components;

import java.awt.*;
import java.awt.event.*;
import java.util.Vector;

import javax.swing.*;

/**
 *  A popup TextField JDialog with an OK and a Cancel button.
 *
 *  @author Kim Rutherford
 *  @version $Id: TextDialog.java,v 1.1 2004-06-09 09:47:52 tjc Exp $
 **/

public class TextDialog extends JDialog {
  /**
   *  Create a new TextDialog component with the given prompt. Other
   *  components can listen for TextDialogEvent object.
   *  @param parent The parent window.
   *  @param prompt A message that is displayed in the component beside the
   *    TextArea that the user types into.  This String is also used as the
   *    JFrame title.
   *  @param width The width of the new TextField.
   *  @param initial_text The initial text to put in the TextField.
   **/
  public TextDialog (final JFrame parent,
                     final String prompt,
                     final int width,
                     final String initial_text) {
    super (parent, prompt, true);

    getContentPane ().add (new JLabel (prompt), "North");

    final JPanel panel = new JPanel ();

    panel.add (ok_button);
    ok_button.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent e) {
        text = text_field.getText ();
        TextDialog.this.dispose ();
      }
    });

    panel.add (cancel_button);
    cancel_button.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent e) {
        text = null;
        TextDialog.this.dispose ();
      }
    });

    getContentPane ().add (panel, "South");

    text_field = new JTextField (initial_text, width);

    text_field.addKeyListener (new KeyAdapter () {
      public void keyTyped(final KeyEvent e) {
        if (e.getKeyChar () == '\n') {
          text = text_field.getText ();
          TextDialog.this.dispose ();
        }
      }
    });

    getContentPane ().add (text_field, "Center");

    pack ();

    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
    setLocation (new Point ((screen.width - getSize ().width) / 2,
                            (screen.height - getSize ().height) / 2));
  }

  /**
   *  Wait for a user action then return the text if the user hit OK or null
   *  if the user hit Cancel.
   **/
  public String getText () {
    setVisible (true);

    return text;
  }

  private final JButton ok_button = new JButton ("OK");
  private final JButton cancel_button = new JButton ("Cancel");
  private JTextField text_field = null;

  /**
   *  Set to null if and only if the user cancels the dialog, otherwise
   *  contains the text that the user entered.
   **/
  private String text = null;
}
