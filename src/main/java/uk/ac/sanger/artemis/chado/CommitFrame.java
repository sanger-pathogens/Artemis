/* CommitFrame
 *
 * created: July 2006
 *
 * This file is part of Artemis
 *
 * Copyright (C) 2006  Genome Research Limited
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
 */

package uk.ac.sanger.artemis.chado;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Vector;

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.ListSelectionModel;

import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.GotoEventSource;
import uk.ac.sanger.artemis.Selection;
import uk.ac.sanger.artemis.components.BasePlotGroup;
import uk.ac.sanger.artemis.components.EntryEdit;
import uk.ac.sanger.artemis.components.Utilities;


public class CommitFrame extends JFrame
{
  private static final long serialVersionUID = 1L;
  private JList commitList = new JList();
  private Color HIGHLIGHT = Color.red;
  private static Cursor CBUSY = new Cursor(Cursor.WAIT_CURSOR);
  private static Cursor CDONE = new Cursor(Cursor.DEFAULT_CURSOR);
  private JButton testCommitButton = new JButton("Test Commit");
  
  public CommitFrame(final ChadoTransactionManager ctm,
                     final EntryGroup entry_group,
                     final EntryEdit entry_edit,
                     final Selection selection,
                     final GotoEventSource goto_event_source,
                     final BasePlotGroup base_plot_group)
  {
    setList(ctm);
    
    final JPanel panel = (JPanel) getContentPane();
    panel.setLayout(new BorderLayout());
    
    final JScrollPane jsp =  new JScrollPane(commitList);
    
    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
    jsp.setPreferredSize(new Dimension(screen.width/2, jsp.getPreferredSize().height));
    panel.add(jsp, BorderLayout.CENTER);
    
    commitList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
    commitList.setSelectionBackground(HIGHLIGHT);
    commitList.setVisibleRowCount(10);
    
    final Box xBox = Box.createHorizontalBox();
    xBox.add(Box.createHorizontalGlue());
    xBox.add(testCommitButton);
    
    testCommitButton.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent arg0)
      {
        System.setProperty("nocommit","true");
        setCursor(CBUSY);
        EntryEdit.commitToDatabase(entry_group, ctm, entry_edit, 
            selection, goto_event_source, base_plot_group);
        setCursor(CDONE);
        if(ctm.getCommitReturnValue() < ctm.getTransactionCount())
        {
          commitList.setSelectedIndex(ctm.getCommitReturnValue());
          commitList.ensureIndexIsVisible(ctm.getCommitReturnValue());
          
          JOptionPane.showMessageDialog(CommitFrame.this, 
              "Test commit failed!",
              "Commit Test Result", 
              JOptionPane.INFORMATION_MESSAGE);
        }
        else
          JOptionPane.showMessageDialog(CommitFrame.this, 
              "Test commit (of "+ctm.getCommitReturnValue()+
              " changes) succeeded!",
              "Commit Test Result", 
              JOptionPane.INFORMATION_MESSAGE);

        System.setProperty("nocommit","false");
      }
    });
    
    final JButton closeButton = new JButton("Close");
    xBox.add(closeButton);
    closeButton.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        CommitFrame.this.dispose();
      }
    });
    
    panel.add(xBox, BorderLayout.NORTH);
    pack();
    
    Utilities.centreFrame(this);
    setVisible(true);
  }
  
  /**
   * Set the data for the JList
   * @param ctm
   */
  public void setList(final ChadoTransactionManager ctm)
  {
    setTitle("Commit List :: "+(ctm.getTransactionCount()>0 ? 
        ctm.getTransactionCount()+" commit(s):" : "Nothing to commit"));
    
    final Vector transactions = new Vector(ctm.getTransactionCount());
    for(int i=0; i<ctm.getTransactionCount(); i++)
    {
      uk.ac.sanger.artemis.chado.ChadoTransaction tsn = 
        ctm.getTransactionAt(i);
      
      transactions.add(tsn.getLogComment());
    }
    
    commitList.setListData(transactions);
    setCommitButtonColour(ctm);
    
    validate();
  }

  /**
   * Set the colour of the button depending on whether there is
   * anything to commit 
   * @param ctm
   */
  private void setCommitButtonColour(final ChadoTransactionManager ctm)
  {
    if(ctm.hasTransactions())
    {
      testCommitButton.setFont(testCommitButton.getFont().deriveFont(Font.BOLD));
      testCommitButton.setForeground(Color.red);
    }
    else
    {
      testCommitButton.setFont(testCommitButton.getFont().deriveFont(Font.PLAIN));
      testCommitButton.setForeground(Color.black);
    }
  }
}