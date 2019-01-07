/* AlignMatchViewer.java
 *
 * created: Tue Feb 13 2001
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/AlignMatchViewer.java,v 1.3 2008-05-29 14:18:48 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;

import java.awt.*;
import java.awt.event.*;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import javax.swing.*;
import javax.swing.event.ListSelectionListener;
import javax.swing.event.ListSelectionEvent;
import java.util.Comparator;
import java.util.Collections;

/**
 *  A component for viewing AlignMatchVectors selected in an AlignmentViewer.
 *
 *  @author Kim Rutherford
 *  @version $Id: AlignMatchViewer.java,v 1.3 2008-05-29 14:18:48 tjc Exp $
 **/

public class AlignMatchViewer extends JFrame 
{

  /**
   *  The AlignmentViewer to call alignAt () and setSelection () on when the
   *  user clicks on a match.
   **/
  private AlignmentViewer alignment_viewer = null;

  /**
   *  The Vector of AlignMatch objects that was passed to the constructor.
   **/
  private AlignMatchVector matches = null;

  /**
   *  The list that contains the matches.
   **/
  private JList list;

  /**
   *  If selected the list items will be sorted by score.
   **/
  private JCheckBoxMenuItem sort_by_score_menu_item =
    new JCheckBoxMenuItem ("Sort by Score");

  /**
   *  If selected the list items will be sorted by percent identity.
   **/
  private JCheckBoxMenuItem sort_by_percent_id_menu_item =
    new JCheckBoxMenuItem ("Sort by Percent Identity");

  /**
   *  If selected the list items will be sorted by the start position of the
   *  query.
   **/
  private JCheckBoxMenuItem sort_by_query_start =
    new JCheckBoxMenuItem ("Sort by Hit Query Start");

  /**
   *  If selected the list items will be sorted by the start position of the
   *  subject.
   **/
  private JCheckBoxMenuItem sort_by_subject_start =
    new JCheckBoxMenuItem ("Sort by Hit Subject start");

  /**
   *  Create a new AlignMatchViewer which shows the given matches.
   *  @param alignment_viewer The AlignmentViewer to call alignAt () and
   *    setSelection () on when the user clicks on a match.
   *  @param matches The Vector of AlignMatch objects to show.
   **/
  public AlignMatchViewer(final AlignmentViewer alignment_viewer,
                          final AlignMatchVector matches) 
  {
    this.alignment_viewer = alignment_viewer;
    this.matches = matches;

    final JMenuBar menu_bar = new JMenuBar();
    final JMenu file_menu = new JMenu("File");
    
    final JMenuItem save = new JMenuItem("Save List to File...");
    save.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent arg0)
      {
        saveMatchList();
      }
    });
    
    final JMenuItem close = new JMenuItem("Close");
    close.addActionListener(new ActionListener () 
    {
      public void actionPerformed (ActionEvent event) 
      {
        setVisible (false);
        AlignMatchViewer.this.dispose ();
      }
    });

    file_menu.add (save);
    file_menu.add (close);
    menu_bar.add (file_menu);

    final JMenu sort_menu = new JMenu("Sort");
    sort_menu.add(sort_by_score_menu_item);

    sort_by_score_menu_item.addItemListener(new ItemListener () 
    {
      public void itemStateChanged(ItemEvent e) 
      {
        if (sort_by_score_menu_item.getState()) 
        {
          sort_by_percent_id_menu_item.setState (false);
          sort_by_query_start.setState (false);
          sort_by_subject_start.setState (false);
        }
        setList ();
      }
    });

    sort_menu.add (sort_by_percent_id_menu_item);

    sort_by_percent_id_menu_item.addItemListener(new ItemListener()
    {
      public void itemStateChanged(ItemEvent e) 
      {
        if(sort_by_percent_id_menu_item.getState())
        {
          sort_by_score_menu_item.setState(false);
          sort_by_query_start.setState(false);
          sort_by_subject_start.setState (false);
        }
        setList ();
      }
    });

    sort_menu.add (sort_by_query_start);

    sort_by_query_start.addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent e) 
      {
        if(sort_by_query_start.getState()) 
        {
          sort_by_percent_id_menu_item.setState(false);
          sort_by_score_menu_item.setState(false);
          sort_by_subject_start.setState(false);
        }

        setList();
      }
    });

    sort_menu.add(sort_by_subject_start);

    sort_by_subject_start.addItemListener(new ItemListener() 
    {
      public void itemStateChanged(ItemEvent e) 
      {
        if(sort_by_subject_start.getState()) 
        {
          sort_by_percent_id_menu_item.setState(false);
          sort_by_score_menu_item.setState(false);
          sort_by_query_start.setState(false);
        }
        setList ();
      }
    });

    menu_bar.add(sort_menu);

    setJMenuBar(menu_bar);

    list = new JList();

    list.setBackground(Color.white);

    getContentPane().add(new JScrollPane(list), "Center");

    final JPanel panel = new JPanel();

    final JButton close_button = new JButton("Close");

    panel.add(close_button);
    close_button.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent e)
      {
        AlignMatchViewer.this.dispose ();
      }
    });

    close_button.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent e) 
      {
        AlignMatchViewer.this.dispose();
      }
    });

    getContentPane().add(panel, "South");

    addWindowListener(new WindowAdapter() 
    {
      public void windowClosing(WindowEvent event) 
      {
        AlignMatchViewer.this.dispose();
      }
    });

    pack();

    setList();

    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();

    int screen_height = screen.height;
    int screen_width  = screen.width;

    if(screen_height <= 600) 
      setSize (350, screen_height * 9 / 10);
    else 
      setSize (350, screen_height - 200);

    setLocation(new Point (screen.width - getSize ().width - 5,
                           (screen.height - getSize ().height) / 2));
  }

  /**
   *  Sort the matches depending on the setting of sort_by_score_menu_item and
   *  sort_by_percent_id_menu_item.
   **/
  private AlignMatchVector getSortedMatches() 
  {
    final AlignMatchVector matches_copy = (AlignMatchVector)matches.clone();

    final Comparator comparator;

    if(sort_by_score_menu_item.getState())
    {
      comparator =
        new Comparator () 
        {
          public int compare (Object fst, Object snd) 
          {
            final int fst_score = ((AlignMatch)fst).getScore();
            final int snd_score = ((AlignMatch)snd).getScore();

            if (fst_score < snd_score) 
              return 1;
            else
            {
              if (fst_score > snd_score) 
                return -1;
              else 
                return 0;
            }
          }
        };
    }
    else 
    {
      if(sort_by_percent_id_menu_item.getState()) 
      {
        comparator =
          new Comparator() 
          {
            public int compare (Object fst, Object snd) 
            {
              final int fst_value = ((AlignMatch)fst).getPercentID();
              final int snd_value = ((AlignMatch)snd).getPercentID();
              
              if (fst_value < snd_value) 
                return 1;
              else 
              {
                if(fst_value > snd_value) 
                  return -1;
                else 
                  return 0;
              }
            }
          };
      } 
      else
      {
        if(sort_by_query_start.getState())
        {
          comparator =
            new Comparator() 
            {
              public int compare(Object fst, Object snd) 
              {
                final int fst_value =
                  ((AlignMatch)fst).getQuerySequenceStart();
                final int snd_value =
                  ((AlignMatch)snd).getQuerySequenceStart();
                
                if(fst_value > snd_value) 
                  return 1;
                else
                {
                  if (fst_value < snd_value) 
                    return -1;
                  else 
                    return 0;
                }
              }
            };
        } 
        else
        {
          if(sort_by_subject_start.getState()) 
          {
            comparator =
              new Comparator() 
              {
                public int compare (Object fst, Object snd) 
                {
                  final int fst_value =
                    ((AlignMatch)fst).getSubjectSequenceStart();
                  final int snd_value =
                    ((AlignMatch)snd).getSubjectSequenceStart();
                  
                  if (fst_value > snd_value) 
                    return 1;
                  else 
                  {
                    if (fst_value < snd_value) 
                      return -1;
                    else 
                      return 0;
                  }
                }
              };
          } 
          else 
            return matches;
        }
      }
    }

    matches_copy.sort(comparator);
    return matches_copy;
  }

  /**
   *  Clear the List and then fill it with the matches in the order
   **/
  private void setList () 
  {
    final AlignMatchVector sorted_matches = getSortedMatches();

    list.setEnabled(false);
    list.setVisible(false);
    list.removeAll();

    String listItems[] = new String[sorted_matches.size()];

    for(int i = 0 ; i<sorted_matches.size() ; ++i)  
    {
      final AlignMatch this_align_match = sorted_matches.elementAt(i);

      listItems[i] = new String(
        this_align_match.getQuerySequenceStart() + ".." +
        this_align_match.getQuerySequenceEnd() + " -> " +
        this_align_match.getSubjectSequenceStart() + ".." +
        this_align_match.getSubjectSequenceEnd() + " " +
        (this_align_match.isRevMatch() ? "rev " : "") +
        this_align_match.getPercentID() + "% id score " +
        this_align_match.getScore() );

    }

    list.setListData(listItems);

    list.addListSelectionListener(new ListSelectionListener() 
    {
      public void valueChanged(ListSelectionEvent e) 
      {
        final int item_number = list.getSelectedIndex();

        final AlignMatch selected_match = matches.elementAt(item_number);

        alignment_viewer.setSelection(selected_match);
        alignment_viewer.alignAt(selected_match);
      }
    });

    list.setEnabled(true);
    list.setVisible(true);
  }
  
  
  /**
   *  Save the text of the match list to a file.
   **/
  private void saveMatchList()
  {
    final StickyFileChooser file_dialog = new StickyFileChooser();

    file_dialog.setDialogTitle("Choose save file ...");
    file_dialog.setDialogType(JFileChooser.SAVE_DIALOG);
    final int status = file_dialog.showSaveDialog(this);

    if(status != JFileChooser.APPROVE_OPTION ||
       file_dialog.getSelectedFile() == null)
      return;

    final File write_file =
      new File(file_dialog.getCurrentDirectory(),
                file_dialog.getSelectedFile().getName());

    if(write_file.exists())
    {
      final YesNoDialog yes_no_dialog =
        new YesNoDialog(this,
                         "this file exists: " + write_file +
                         " overwrite it?");
      if(yes_no_dialog.getResult())
      {
        // yes - continue
      }
      else
        return;
    }

    try
    {
      final PrintWriter writer =
        new PrintWriter(new FileWriter(write_file));

      for(int i = 0 ; i < list.getModel().getSize() ; ++i)
        writer.println(list.getModel().getElementAt(i));
      writer.close();
    }
    catch(IOException e)
    {
      new MessageDialog(this, "error while writing: " + e.getMessage());
    }
  }


}
