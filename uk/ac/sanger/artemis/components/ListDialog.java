/* ListDialog.java
 *
 * created: Fri Sep  1 2000
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 2000  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/ListDialog.java,v 1.2 2004-12-03 18:11:28 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import java.awt.*;
import java.awt.event.*;

import javax.swing.*;

/**
 *  This component is a JDialog that contains a List.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: ListDialog.java,v 1.2 2004-12-03 18:11:28 tjc Exp $
 **/

public class ListDialog extends JDialog {
  /**
   *  Create a new ListDialog component.
   **/
  public ListDialog (final JFrame parent, final String title) {
    super (parent, title, true);

    JScrollPane scroll_pane = new JScrollPane (getList ());

    getContentPane ().add (scroll_pane, "Center");

    getList ().setBackground (Color.white);
    
    button_panel = new JPanel ();

    button_panel.add (ok_button);

    getContentPane ().add (button_panel, "South");

    ok_button.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent e) {
        selected_item = getList ().getSelectedValue ();
        ListDialog.this.dispose ();
      }
    });

    button_panel.add (cancel_button);
    cancel_button.addActionListener (new ActionListener () {
      public void actionPerformed (ActionEvent e) {
        selected_item = null;
        ListDialog.this.dispose ();
      }
    });
    packme ();
  }
    
  /**
   *  This method will call pack () and then move the JDialog to the centre of
   *  the screen.
   **/
  private void packme () {
    pack ();

    setSize (750, 400);

    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();

    final int x_position = (screen.width - getSize ().width) / 2;
    int y_position = (screen.height - getSize ().height) / 2;

    if (y_position < 10) {
      y_position = 10;
    }

    setLocation (new Point (x_position, y_position));
  }

  /**
   *  Show the dialog and then return the selected item, or null if the user
   *  hits cancel.
   **/
  public Object getSelectedValue () {
    setVisible(true);
    return selected_item;
  }

  /**
   *  Return the reference of the List.
   **/
  public JList getList () {
    return list;
  }

  /**
   *  The List.
   **/
  final private JList list = new JList ();

  private final JButton ok_button = new JButton ("OK");
  private final JButton cancel_button = new JButton ("Cancel");

  private JPanel button_panel = null;

  /**
   *  Set to null if and only if the user cancels the dialog, otherwise
   *  contains the selected item.
   **/
  private Object selected_item = null;
}
