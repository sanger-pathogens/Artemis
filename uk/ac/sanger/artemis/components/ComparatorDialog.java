/* ComparatorDialog.java
 *
 * created: Wed May 10 2000
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2000  Genome Research Limited
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or(at your option) any later version.
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/ComparatorDialog.java,v 1.3 2006-10-18 14:25:23 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;

import uk.ac.sanger.artemis.util.InputStreamProgressListener;

import java.awt.*;
import java.awt.event.*;
import java.io.*;
import java.util.Vector;

import javax.swing.*;

/**
 *  ComparatorDialog a dialog that allows the user to the choose two files to
 *  compare with a Comparator component and a file containing comparison data.
 *
 *  @author Kim Rutherford <kmr@sanger.ac.uk>
 *  @version $Id: ComparatorDialog.java,v 1.3 2006-10-18 14:25:23 tjc Exp $
 **/

public class ComparatorDialog extends JFrame 
{
  /**
   *  Create a dialog that allow the user to the choose two files to compare
   *  and a file containing comparison data.
   *  @param ActMain the object to call makeMultiComparator() on
   **/
  public ComparatorDialog(final ActMain act_main)
  {
    final JPanel top_panel = new JPanel();

    getContentPane().add(top_panel, "Center");

    final GridBagLayout gridbag = new GridBagLayout();
    top_panel.setLayout(gridbag);

    final Vector text_field_vector = new Vector();

    for(int i = 0; i < 3; ++i) 
    {
      final String label;
      switch(i) 
      {
        case 0:
          label = "Sequence file 1";
          break;
        case 1:
          label = "Comparison file 1";
          break;
        case 2:
          label = "Sequence file 2";
          break;
        default:
          throw new Error("internal error");
      }

      final JTextField text_field =
        makeFileNamePanel(label, top_panel, gridbag);

      text_field_vector.addElement(text_field);
    }


    final GridBagConstraints c = new GridBagConstraints();
    c.fill = GridBagConstraints.HORIZONTAL;
    c.anchor = GridBagConstraints.NORTH;
    c.weighty = 0;
    c.gridwidth = 1;

    final JButton more_button = new JButton("more files ...");

    final JPanel more_button_panel = new JPanel();
    more_button_panel.setLayout(new FlowLayout(FlowLayout.LEFT));
    more_button_panel.add(more_button);

    c.gridwidth = 1;
    gridbag.setConstraints(more_button_panel, c);
    top_panel.add(more_button_panel);

    more_button.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent e) 
      {
        top_panel.remove(more_button_panel);

        final boolean add_sequence_flag;

        if(text_field_vector.size() % 2 == 0 ||
            Options.getOptions().isNoddyMode()) 
          add_sequence_flag = true;
        else 
          add_sequence_flag = false;

        if(text_field_vector.size() % 2 == 1 ||
            Options.getOptions().isNoddyMode()) 
        {
          final String seq_label =
            "Comparison file " +(text_field_vector.size() / 2 + 1);
          final JTextField seq_text_field =
            makeFileNamePanel(seq_label, top_panel, gridbag);

          text_field_vector.addElement(seq_text_field);
        }

        if(add_sequence_flag) 
        {
          final String comp_label =
            "Sequence file " +(text_field_vector.size() / 2 + 1);
          final JTextField comp_text_field =
            makeFileNamePanel(comp_label, top_panel, gridbag);

          text_field_vector.addElement(comp_text_field);
        }

        top_panel.add(more_button_panel);

        packAndCentre();
      }
    });


    final JButton apply_button = new JButton("Apply");

    apply_button.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent e) 
      {
        doApply(act_main, text_field_vector);
      }
    });


    final JButton close_button = new JButton("Close");

    close_button.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent e) 
      {
        ComparatorDialog.this.dispose();
      }
    });


    final FlowLayout flow_layout =
      new FlowLayout(FlowLayout.CENTER, 15, 5);

    final JPanel close_and_apply_panel = new JPanel(flow_layout);

    close_and_apply_panel.add(apply_button);
    close_and_apply_panel.add(close_button);

    getContentPane().add(close_and_apply_panel, "South");


    addWindowListener(new WindowAdapter() 
    {
      public void windowClosing(WindowEvent event) 
      {
        ComparatorDialog.this.dispose();
      }
    });

    packAndCentre();
  }

  /**
   *  Call pack() and then centre the JFrame in the middle of the screen.
   **/
  private void packAndCentre() 
  {
    pack();

    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
    setLocation(new Point((screen.width - getSize().width) / 2,
                           (screen.height - getSize().height) / 2));
  }

  /**
   *  Make a panel for choosing one file name panel.
   **/
  private JTextField makeFileNamePanel(final String label_string,
                                       final JPanel parent_panel,
                                       final GridBagLayout gridbag) 
  {
    GridBagConstraints c = new GridBagConstraints();

    c.fill = GridBagConstraints.HORIZONTAL;
    c.anchor = GridBagConstraints.NORTH;
    c.weighty = 0;

    final JLabel label = new JLabel(label_string);

    final JPanel panel = new JPanel();
    panel.setLayout(new FlowLayout(FlowLayout.LEFT));
    panel.add(label);

    c.gridwidth = 1;
    gridbag.setConstraints(panel, c);
    parent_panel.add(panel);

    final TextFieldSink text_field = new TextFieldSink("", 28);
    gridbag.setConstraints(text_field, c);
    parent_panel.add(text_field);

    final JButton choose_button = new JButton("Choose ...");

    choose_button.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent actionEvent) 
      {
        final StickyFileChooser file_dialog =
          new StickyFileChooser();

        file_dialog.setDialogTitle("Choose first sequence ...");
        file_dialog.setFileSelectionMode(JFileChooser.FILES_ONLY);
        file_dialog.setDialogType(JFileChooser.OPEN_DIALOG);

        file_dialog.showOpenDialog(ComparatorDialog.this);

        if(file_dialog.getSelectedFile() != null) 
        {
          final File selected_file =
            new File(file_dialog.getCurrentDirectory(),
                      file_dialog.getSelectedFile().getName());

          text_field.setText(selected_file.toString());
        }
      }
    });

    c.gridwidth = GridBagConstraints.REMAINDER;
    gridbag.setConstraints(choose_button, c);
    parent_panel.add(choose_button);

    return text_field;
  }

  /**
   *  Attempt to call ActMain.makeComparator() for the given file names.
   *  @param ActMain the object to call makeMultiComparator() on
   **/
  private void doApply(final ActMain act_main,
                        final Vector text_field_vector) 
  {
    if(text_field_vector.size() < 3) 
      throw new Error("internal error - not enough file names given to " +
                       "ComparatorDialog.doApply()");

    Object[] file_names = new Object[text_field_vector.size()];

    for(int i = 0; i < text_field_vector.size(); ++i) 
    {
      TextFieldSink tfield = (TextFieldSink)text_field_vector.elementAt(i);
      
      if(tfield.getDbNode() != null)
        file_names[i] = tfield.getDbNode();
      else
        file_names[i] = tfield.getText().trim();

      if(file_names[i] instanceof String && ((String)file_names[i]).length() == 0) 
      {
        // ignore the problem if there are at least 3 files listed(ie. seq1
        // comp1v2 seq2) and the remaining files lengths are zero
        if(i > 3)
        {
          // set to true if a text field is found that isn't zero length
          boolean found_file = false;

          for(int sub_index = i; sub_index < text_field_vector.size();
              ++sub_index) 
          {
            final JTextField text_field =
             (JTextField) text_field_vector.elementAt(sub_index);

            final String this_text = text_field.getText().trim();

            if(this_text.length() != 0) 
              found_file = true;
          }

          if(!found_file)
          {
            // truncate file_names
            final String [] new_file_names = new String [i];
            System.arraycopy(file_names, 0, new_file_names, 0, i);
            file_names = new_file_names;
            break;
          }
        }

        new MessageDialog(this, "one of the file names is missing");
        return;
      }
    }

    final MessageFrame reading_message = new MessageFrame("reading ...");

    final InputStreamProgressListener progress_listener =
      act_main.getInputStreamProgressListener();

    if(ActMain.makeMultiComparator(act_main, progress_listener,
                                   file_names)) 
      ComparatorDialog.this.dispose();
    
    reading_message.dispose();
  }
}
