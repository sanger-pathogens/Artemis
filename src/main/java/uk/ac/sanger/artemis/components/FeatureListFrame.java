/* FeatureListFrame.java
 *
 * created: Fri Sep  3 1999
 *
 * This file is part of Artemis
 *
 * Copyright (C) 1999,2000,2001  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/FeatureListFrame.java,v 1.3 2005-12-09 16:17:13 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;

import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

/**
 *  A JFrame that contains a FeatureList component and an Close button.
 *
 *  @author Kim Rutherford
 *  @version $Id: FeatureListFrame.java,v 1.3 2005-12-09 16:17:13 tjc Exp $
 **/

public class FeatureListFrame extends JFrame
    implements EntryGroupChangeListener 
{
  /**
   *  Create a new FeatureListFrame component.  The constructor does not call
   *  setVisible (true).
   *  @param title The title to use for the new JFrame.
   *  @param feature_list The FeatureList to show.
   *  @param selection The Selection that the commands in the menus will
   *    operate on.
   *  @param entry_group The EntryGroup object where new features/entries will
   *    be added.
   *  @param goto_event_source The object that the menu item will call
   *    makeBaseVisible() on.
   **/
  public FeatureListFrame(final String title,
                          final Selection selection,
                          final GotoEventSource goto_event_source,
                          final EntryGroup entry_group,
                          final BasePlotGroup base_plot_group) 
  {
    super(title);

    this.entry_group = entry_group;

    feature_list = new FeatureList(entry_group, selection, goto_event_source,
                                   base_plot_group);
    
    final Font default_font = Options.getOptions().getFont();
    setFont(default_font);

    final JMenuBar menu_bar = new JMenuBar();
    menu_bar.setFont(default_font);
    setJMenuBar(menu_bar);

    final JMenu file_menu = new JMenu("File");
    final JMenuItem close = new JMenuItem("Close");
    close.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent event) 
      {
        setVisible(false);
        FeatureListFrame.this.dispose();
        feature_list.stopListening();        
      }
    });

    file_menu.add(close);
    menu_bar.add(file_menu);
    
    final JMenu select_menu =
      new SelectMenu(this, selection, goto_event_source, entry_group,
                     base_plot_group);
    menu_bar.add(select_menu);

    final JMenu view_menu =
      new ViewMenu(this, selection, goto_event_source, entry_group,
                   base_plot_group);
    menu_bar.add(view_menu);

    final JMenu goto_menu =
      new GotoMenu(this, selection, goto_event_source, entry_group);
    menu_bar.add(goto_menu);

    if(Options.readWritePossible())
    {
      final JMenu edit_menu =
        new EditMenu(this, selection, goto_event_source, entry_group,
                     base_plot_group, null);
      menu_bar.add(edit_menu);

      final JMenu write_menu = new WriteMenu(this, selection, entry_group);
      menu_bar.add(write_menu);

      final JMenu run_menu = new RunMenu(this, selection);
      menu_bar.add(run_menu);
    }

    final JScrollPane jsp_feature_list = new JScrollPane(feature_list);
    getContentPane().add(jsp_feature_list, "Center");

    final JPanel panel = new JPanel();

    final JButton close_button = new JButton("Close");

    panel.add(close_button);
    close_button.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent e) 
      {
        setVisible(false);
        FeatureListFrame.this.dispose();
        feature_list.stopListening();
      }
    });

    getContentPane().add(panel, "South");
    pack();

    addWindowListener(new WindowAdapter()
    {
      public void windowClosing(WindowEvent event) 
      {
        setVisible(false);
        entry_group.removeEntryGroupChangeListener(FeatureListFrame.this);
        FeatureListFrame.this.dispose();
        feature_list.stopListening();
      }
    });

    entry_group.addEntryGroupChangeListener(this);

    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();

    int screen_height = screen.height;
    int screen_width = screen.width;

    if(screen_width <= 700 || screen_height <= 400) 
      setSize(screen_width * 9 / 10, screen_height * 9 / 10);
    else 
      setSize(700, 400);

    final int hgt = entry_group.getAllFeaturesCount() *
                    feature_list.getLineHeight();
    feature_list.setPreferredSize(new Dimension(getSize().width*4,hgt));
    jsp_feature_list.getVerticalScrollBar().setUnitIncrement(feature_list.getLineHeight());

    setLocation(new Point((screen.width - getSize().width) / 2,
                          (screen.height - getSize().height) / 2));
  }

  /**
   *  Implementation of the EntryGroupChangeListener interface.  We listen to
   *  EntryGroupChange events so that we can get rid of the Navigator when the
   *  EntryGroup is no longer in use (for example when the EntryEdit is
   *  closed).
   **/
  public void entryGroupChanged(final EntryGroupChangeEvent event) 
  {
    switch(event.getType())
    {
      case EntryGroupChangeEvent.DONE_GONE:
        entry_group.removeEntryGroupChangeListener(this);
        dispose();
        break;
    }
  }

  /**
   *  Return the FeatureList that this JFrame is displaying.
   **/
  public FeatureList getFeatureList() 
  {
    return feature_list;
  }

  /**
   *  The FeatureList that this JFrame is displaying.
   **/
  private FeatureList feature_list;

  /**
   *  The EntryGroup that was passed to the constructor.
   **/
  private EntryGroup entry_group;
}
