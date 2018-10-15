/* TextRequester.java
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/TextRequester.java,v 1.1 2004-06-09 09:47:54 tjc Exp $
 **/

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.Options;

import java.awt.*;
import java.awt.event.*;
import java.util.Vector;

import javax.swing.*;

/**
 *  A popup JTextField with an OK and a Cancel button.
 *
 *  @author Kim Rutherford
 *  @version $Id: TextRequester.java,v 1.1 2004-06-09 09:47:54 tjc Exp $
 **/

public class TextRequester extends JFrame {
  /**
   *  Create a new TextRequester component with the given prompt. Other
   *  components can listen for TextRequesterEvent object.
   *  @param prompt A message that is displayed in the component beside the
   *    TextArea that the user types into.  This String is also used as the
   *    JFrame title.
   *  @param width The width of the JTextField in the new requester.
   *  @param initial_text The initial text to put in the JTextField.
   **/
  public TextRequester (final String prompt,
                        final int width,
                        final String initial_text) {
    super (prompt);

    getContentPane ().add (new JLabel (prompt), "North");

    final JPanel panel = new JPanel ();

    panel.add (ok_button);
    ok_button.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent e) {
        performOK ();
      }
    });

    panel.add (cancel_button);
    cancel_button.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent e) {
        performCancel ();
      }
    });

    getContentPane ().add (panel, "South");

    text_field = new JTextField (initial_text, width);

    text_field.addKeyListener (new KeyAdapter () {
      public void keyTyped(final KeyEvent e) {
        if (e.getKeyChar () == '\n') {
          performOK ();
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
   *  Add the given object as a listen for TextRequester events from this
   *  TextRequester.
   **/
  public void addTextRequesterListener (final TextRequesterListener l) {
    listeners.addElement (l);
  }

  /**
   *  Return the text the is currently displayed in this requester.
   **/
  protected String getText () {
    return text_field.getText ();
  }

  /**
   *  Send a TextRequesterEvent of type OK to all the listeners.
   **/
  protected void performOK () {
    final TextRequesterEvent new_event =
      new TextRequesterEvent (this, getText (), TextRequesterEvent.OK);

    sendEvent (new_event);

    TextRequester.this.dispose ();
  }

  /**
   *  Send a TextRequesterEvent of type CANCEL to all the listeners.
   **/
  protected void performCancel () {
    final TextRequesterEvent new_event =
      new TextRequesterEvent (this, getText (), TextRequesterEvent.CANCEL);

    sendEvent (new_event);

    TextRequester.this.dispose ();
  }

  /**
   *  Send the given TextRequesterEvent to all the object that are listening
   *  for the event.
   **/
  private void sendEvent (final TextRequesterEvent event) {
    for (int i = 0 ; i < listeners.size () ; ++i) {
      final TextRequesterListener listener =
        ((TextRequesterListener) listeners.elementAt (i));

      listener.actionPerformed (event);
    }
  }

  private final JButton ok_button = new JButton ("OK");
  private final JButton cancel_button = new JButton ("Cancel");
  private JTextField text_field = null;

  /**
   *  This contains the objects that are listening for TextRequester events
   *  from this TextRequester.
   **/
  private Vector listeners = new Vector ();
}
