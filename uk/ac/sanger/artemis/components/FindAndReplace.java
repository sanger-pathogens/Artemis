/* 
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2008  Genome Research Limited
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
 */

package uk.ac.sanger.artemis.components;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.regex.Pattern;

import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JRadioButton;
import javax.swing.JSeparator;
import javax.swing.JTextField;

import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.EntryGroupChangeEvent;
import uk.ac.sanger.artemis.EntryGroupChangeListener;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureKeyQualifierPredicate;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.FilteredEntryGroup;
import uk.ac.sanger.artemis.GotoEventSource;
import uk.ac.sanger.artemis.Selection;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.EntryInformationException;
import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.util.ReadOnlyException;
import uk.ac.sanger.artemis.util.StringVector;

public class FindAndReplace extends JFrame
       implements EntryGroupChangeListener 
{
  private static final long serialVersionUID = 1L;
  private EntryGroup entry_group;
  
  public FindAndReplace(final Selection selection,
      final GotoEventSource goto_event_source,
      final EntryGroup entry_group,
      final BasePlotGroup base_plot_group)
  {
    super("Find/Replace Qualifier Text");

    GridBagLayout gridbag = new GridBagLayout();
    getContentPane().setLayout(gridbag);

    GridBagConstraints c = new GridBagConstraints();

    c.fill = GridBagConstraints.HORIZONTAL;
    c.anchor = GridBagConstraints.WEST;
    c.ipadx = 5;
    c.ipady = 5;


    /*final FeatureVector features = entry_group.getAllFeatures ();
    
    for(int i=0; i<features.size(); i++)
      findOrDeleteDuplicate(features.elementAt(i), true);*/
    
    Entry default_entry = entry_group.getDefaultEntry();
    if(default_entry == null)
      default_entry = entry_group.elementAt(0);
  
    final EntryInformation default_entry_information =
                        default_entry.getEntryInformation();

    final KeyChoice key_selector = new KeyChoice(default_entry_information);
    key_selector.setEnabled(false);
    
    //final StringVector qualifierNames =
    //  Feature.getAllQualifierNames (selected_features);
    
    int ypos = 0;
    
    final JTextField findTextField = new JTextField(15);
    c.gridx = 0;
    c.gridy = ++ypos;
    getContentPane().add(new JLabel("Find:"),c);
    c.gridx = 1;
    getContentPane().add(findTextField, c);
    
    final JTextField replaceTextField = new JTextField(15);
    c.gridx = 0;
    c.gridy = ++ypos;
    getContentPane().add(new JLabel("Replace with:"), c);
    c.gridx = 1;
    getContentPane().add(replaceTextField, c);
    
    c.gridx = 0;
    c.gridy = ++ypos;
    getContentPane().add(key_selector,c);
    
    final JCheckBox selectedKeyButton = new JCheckBox("Restrict to Selected Key", false);
    selectedKeyButton.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        if(selectedKeyButton.isSelected())
          key_selector.setEnabled(true);
        else
          key_selector.setEnabled(false);
      }
    });
    c.gridx = 1;
    getContentPane().add(selectedKeyButton,c);
    
    

    final JCheckBox caseSensitive = new JCheckBox("Case sensitive", true);
    c.gridx = 0;
    c.gridy = ++ypos;
    getContentPane().add(caseSensitive, c);
    
    final JCheckBox qualifierValueSubString = new JCheckBox("Match substring", true);
    c.gridy = ++ypos;
    getContentPane().add(qualifierValueSubString, c);
    
    final JButton findButton = new JButton("Find");
    findButton.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        if(findTextField.getText().equals(""))
          JOptionPane.showMessageDialog(FindAndReplace.this,
              "No text entered!", "No Text", JOptionPane.WARNING_MESSAGE);
        
        Key key = null;
        if(selectedKeyButton.isSelected())
          key = key_selector.getSelectedItem();
        
        FeatureKeyQualifierPredicate predicate = 
                         new FeatureKeyQualifierPredicate(key, null,
                                                          findTextField.getText(), 
                                                          qualifierValueSubString.isSelected(), 
                                                          !caseSensitive.isSelected());
        
        final FilteredEntryGroup filtered_entry_group =
          new FilteredEntryGroup(entry_group, predicate, findTextField.getText());

        final FeatureListFrame feature_list_frame =
          new FeatureListFrame("Features Found",
                               selection,
                               goto_event_source, filtered_entry_group,
                               base_plot_group);

        feature_list_frame.setVisible(true);
      }
    });
    c.gridx = 0;
    c.gridy = ++ypos;
    getContentPane().add(findButton, c);
    
    final JButton replaceButton = new JButton("Replace");
    replaceButton.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        Key key = null;
        if(selectedKeyButton.isSelected())
          key = key_selector.getSelectedItem();
        
        FeatureKeyQualifierPredicate predicate = 
          new FeatureKeyQualifierPredicate(key, null,
                                           findTextField.getText(), 
                                           qualifierValueSubString.isSelected(), 
                                           !caseSensitive.isSelected());

        final FilteredEntryGroup filtered_entry_group =
          new FilteredEntryGroup(entry_group, predicate, findTextField.getText());

        entry_group.getActionController().startAction();
        
        FeatureVector features = filtered_entry_group.getAllFeatures();
        for(int i=0; i<features.size(); i++)
        {
          final Feature feature = features.elementAt(i);
          feature.findOrReplaceText(findTextField.getText(),
              !caseSensitive.isSelected(), qualifierValueSubString.isSelected(), 
              null, replaceTextField.getText());
        }
        
        entry_group.getActionController().endAction();
        
/*        final FeatureListFrame feature_list_frame =
          new FeatureListFrame("Features Found and Text Replaced",
                               selection,
                               goto_event_source, filtered_entry_group,
                               base_plot_group);

        feature_list_frame.setVisible(true);*/
      }
    });
    c.gridx = 1;
    getContentPane().add(replaceButton, c);
    
    addWindowListener(new WindowAdapter() 
    {
      public void windowClosing(WindowEvent event) 
      {
        entry_group.removeEntryGroupChangeListener(FindAndReplace.this);
        FindAndReplace.this.dispose();
      }
    });
    pack();
    
    Utilities.centreFrame(this);
    setVisible(true);
  }
  
  /**
   *  Implementation of the EntryGroupChangeListener interface.  We listen to
   *  EntryGroupChange events so that we can get rid of the FinAndReplace when the
   *  EntryGroup is no longer in use(for example when the EntryEdit is
   *  closed).
   **/
  public void entryGroupChanged(EntryGroupChangeEvent event)
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
   *  Return true if duplicate qualifier values found
   *  @param feature The text to search for.
   *  @param delete
   **/
  private boolean findOrDeleteDuplicate(final Feature feature, final boolean delete) 
  {
    final QualifierVector qualifiers = feature.getQualifiers();
    final QualifierVector newQualifiers = new QualifierVector();
 
    for(int i = 0; i  < qualifiers.size(); ++i) 
    {
      final Qualifier this_qualifier = (Qualifier)qualifiers.elementAt(i);
      final StringVector values = this_qualifier.getValues();

      if(values != null)
      {
        StringVector newValues = new StringVector();
        int val_size = values.size();
        
        for(int values_index = 0; values_index < val_size; 
             ++values_index) 
        {
          String this_value_string = (String)values.elementAt(values_index);
          if(this_value_string == null) 
            continue;
            
          if(!newValues.contains(this_value_string))
            newValues.add(this_value_string);
          else if(!delete)
            return true;
        }

        newQualifiers.setQualifier(
              new Qualifier(this_qualifier.getName(),newValues));
      }
    }
    
    if(!delete)
      return false;
    
    for(int i = 0; i < newQualifiers.size(); i++)
    {
      try
      {
        feature.setQualifier((Qualifier) newQualifiers.elementAt(i));
      }
      catch(ReadOnlyException e)
      {
        e.printStackTrace();
      }
      catch(EntryInformationException e)
      {
        e.printStackTrace();
      }
    }
    
    return false;
  }
  
}