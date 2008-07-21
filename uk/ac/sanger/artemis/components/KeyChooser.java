/* KeyChooser.java
 *
 * created: Mon Sep  6 1999
 *
 * This file is part of Artemis
 * 
 * Copyright (C) 1999  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/KeyChooser.java,v 1.2 2008-07-21 16:03:18 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.io.EntryInformation;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JPanel;


/**
 *  This component is a KeyChoice component in a JFrame (see addItemListener ()
 *  and removeItemListener () for details).
 *
 *  @author Kim Rutherford
 *  @version $Id: KeyChooser.java,v 1.2 2008-07-21 16:03:18 tjc Exp $
 **/
public class KeyChooser extends JFrame 
{
  private static final long serialVersionUID = 1L;

  private KeyChoice key_choice;

  final private JButton ok_button = new JButton ("OK");
  
  /**
   *  Create a new KeyChooser component with CDS as the default key.
   *  @param entry_information The object to get the list of possible
   *    keys from.
   **/
  public KeyChooser (final EntryInformation entry_information) 
  {
    this (entry_information, Key.CDS);
  }

  /**
   *  Create a new KeyChooser component with the given key as the default.
   *  @param entry_information The object to get the list of possible
   *    keys from.
   **/
  public KeyChooser (final EntryInformation entry_information,
                     final Key default_key) 
  {
    key_choice = new KeyChoice (entry_information, default_key);

    getContentPane ().add (key_choice, "Center");

    final JPanel panel = new JPanel ();

    panel.add (ok_button);
    ok_button.addActionListener (new ActionListener () 
    {
      public void actionPerformed (ActionEvent e) 
      {
        KeyChooser.this.dispose ();
      }
    });

    final JButton close_button = new JButton ("Cancel");
    panel.add (close_button);
    close_button.addActionListener (new ActionListener ()
    {
      public void actionPerformed (ActionEvent e) 
      {
        KeyChooser.this.dispose ();
      }
    });

    getContentPane ().add (panel, "South");

    addWindowListener (new WindowAdapter () 
    {
      public void windowClosing (WindowEvent event)
      {
        KeyChooser.this.dispose ();
      }
    });

    pack();

    Utilities.centreFrame(this);
  }

  /**
   *  Return the KeyChoice component that is displayed in this JFrame.
   **/
  protected KeyChoice getKeyChoice () 
  {
    return key_choice;
  }

  /**
   *  Return the reference of the OK button of this KeyChooser.
   **/
  protected JButton getOKButton ()
  {
    return ok_button;
  }
}

