/* ChoiceFrame.java
 *
 * created: Tue Aug  6 2002
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 2001  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/ChoiceFrame.java,v 1.1 2004-06-09 09:46:08 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.util.StringVector;

import java.awt.*;
import java.awt.event.*;

import javax.swing.*;

/**
 *  A Choice in a JFrame.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: ChoiceFrame.java,v 1.1 2004-06-09 09:46:08 tjc Exp $
 **/

public class ChoiceFrame extends JFrame {
  /**
   *  Create a new ChoiceFrame component with the given list of Strings.
   **/
  public ChoiceFrame (final String choice_title, final StringVector strings) {
    super (choice_title);

    choice = new JComboBox ();

    for (int i = 0 ; i < strings.size () ; ++i) {
      choice.addItem (strings.elementAt (i));
    }

    final JPanel choice_panel = new JPanel ();
    choice_panel.add (choice);
    
    getContentPane ().add (choice_panel, "Center");

    final JPanel panel = new JPanel ();

    panel.add (ok_button);
    ok_button.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent e) {
        ChoiceFrame.this.dispose ();
      }
    });

    panel.add (close_button);
    close_button.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent e) {
        ChoiceFrame.this.dispose ();
      }
    });

    getContentPane ().add (panel, "South");
    pack ();

    addWindowListener (new WindowAdapter () {
      public void windowClosing (WindowEvent event) {
        ChoiceFrame.this.dispose ();
      }
    });

    pack ();

    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
    setLocation (new Point ((screen.width - getSize ().width) / 2,
                            (screen.height - getSize ().height) / 2));
  }

  /**
   *  Return the Choice component that is displayed in this JFrame.
   **/
  public JComboBox getChoice () {
    return choice;
  }

  /**
   *  Return the reference of the OK button of this Chooser.
   **/
  public JButton getOKButton () {
    return ok_button;
  }

  private JComboBox choice;

  final private JButton ok_button = new JButton ("OK");

  final private JButton close_button = new JButton ("Cancel");
}
