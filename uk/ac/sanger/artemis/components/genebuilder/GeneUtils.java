/* GeneUtils.java
 *
 * This file is part of Artemis
 * Copyright (C) 2007  Genome Research Limited
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
 **/
package uk.ac.sanger.artemis.components.genebuilder;

import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Vector;

import javax.swing.Box;
import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;

import uk.ac.sanger.artemis.io.DatabaseDocumentEntry;
import uk.ac.sanger.artemis.io.Feature;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.EntryVector;
import uk.ac.sanger.artemis.FeatureVector;


public class GeneUtils
{
  private static final long serialVersionUID = 1L;
  private static Vector showFeatures = new Vector();
  private static Vector hideFeatures = new Vector();
  
  static 
  {
    showFeatures.add("gene");
    showFeatures.add("pseudogene");
    showFeatures.add("exon-model");
    showFeatures.add("pseudogenic_exon");
  }
  
  static
  {
    hideFeatures.add("polypeptide");
    hideFeatures.add("mRNA");
    hideFeatures.add("pseudogenic_transcript");
  }
  
  /**
   * Given a collection of features, determine if these should be
   * shown or hidden in the Artemis display
   * @param features
   */
  public static void defineShowHideGeneFeatures(final FeatureVector features)
  {
    final DefaultListModel showListModel = new DefaultListModel();
    for(int i=0; i<showFeatures.size(); i++)
      showListModel.addElement(showFeatures.get(i));
    final JList displayList = new JList(showListModel);
    
    
    final DefaultListModel hideListModel = new DefaultListModel();
    for(int i=0; i<hideFeatures.size(); i++)
      hideListModel.addElement(hideFeatures.get(i));
    final JList hideList = new JList(hideListModel);
    
    
    final JButton hide_butt = new JButton("HIDE");
    hide_butt.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        while(!displayList.isSelectionEmpty())
        {
          final String hideKey = (String)displayList.getSelectedValue();
          hideListModel.addElement(hideKey);
          showListModel.removeElement(hideKey);
          
          hideFeatures.add(hideKey);
          if(showFeatures.contains(hideKey))
            showFeatures.remove(hideKey);
        }
      }
    });

    Box bdown = Box.createVerticalBox();
    bdown.add(new JLabel("Features Displayed:"));
    bdown.add(new JScrollPane(displayList));
    bdown.add(hide_butt);
    
    final JPanel hideShowPanel = new JPanel(new BorderLayout());
    hideShowPanel.add(bdown, BorderLayout.CENTER);

    
    final JButton show_butt = new JButton("SHOW");
    show_butt.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        while(!hideList.isSelectionEmpty())
        {
          final String showKey = (String)hideList.getSelectedValue();
          showListModel.addElement(showKey);
          hideListModel.removeElement(showKey);
          
          if(hideFeatures.contains(showKey))
            hideFeatures.remove(showKey);
          showFeatures.add(showKey);
        }
      }
    });

    bdown = Box.createVerticalBox();
    bdown.add(Box.createVerticalGlue());
    bdown.add(new JLabel("Features Hidden:"));
    bdown.add(new JScrollPane(hideList));
    bdown.add(show_butt);
    hideShowPanel.add(bdown, BorderLayout.EAST);

    int select = JOptionPane.showConfirmDialog(null, hideShowPanel,
                            "Gene Model Features Displayed...",
                             JOptionPane.OK_CANCEL_OPTION,
                             JOptionPane.QUESTION_MESSAGE);

    if(select == JOptionPane.CANCEL_OPTION)
      return;
    
    showHideGeneFeatures(features);
  }
  
  public static void showHideGeneFeatures(final FeatureVector features)
  {
    for(int i=0; i<features.size(); i++)
    {
      final Feature feature = features.elementAt(i).getEmblFeature();
      
      if(feature instanceof GFFStreamFeature)
      {
        final String key = feature.getKey().getKeyString();
        if(hideFeatures.contains(key))
          ((GFFStreamFeature)feature).setVisible(false);
        else if(showFeatures.contains(key))
          ((GFFStreamFeature)feature).setVisible(true);
      }
    }
  }
  
  public static boolean isHiddenFeature(final String key)
  {
    return hideFeatures.contains(key);
  }
  
  /**
   * Given an group of entries determine if they contain a database entry
   * @param entryGroup
   * @return
   */
  public static boolean isDatabaseEntry(final EntryGroup entryGroup)
  {
    final EntryVector entries = entryGroup.getActiveEntries();
    
    for(int i=0; i<entries.size(); i++)
    {
      if( entries.elementAt(i).getEMBLEntry() instanceof DatabaseDocumentEntry )
        return true;
    }
    return false;
  }
  
  public static void main(String args[])
  {
    GeneUtils.defineShowHideGeneFeatures(new FeatureVector());
  }
}
