/* ReadCountDialog.java
 *
 * created: 2012
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2012  Genome Research Limited
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

package uk.ac.sanger.artemis.components.alignment;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JPanel;

import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.components.FeatureDisplay;
import uk.ac.sanger.artemis.components.ListSelectionPanel;
import uk.ac.sanger.artemis.io.InvalidRelationException;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.util.StringVector;


  class ReadCountDialog extends JDialog 
  {
    private static final long serialVersionUID = 1L;
    private static ListSelectionPanel systematicListSelectionPanel;
    private static StringVector qNames;
    private int status = 0;
    
    ReadCountDialog(final JFrame frame, final String title, 
                    final FeatureDisplay feature_display, final Box yBox) 
    {
      super(frame, title, true);

      qNames = null;
      final Object systematic_names[] =
          Options.getOptions().getSystematicQualifierNames().toArray();
      final String[] description = 
            { "Qualifiers to search when getting the feature(s) name." } ;
      
      if(systematicListSelectionPanel == null)
        systematicListSelectionPanel =
              new ListSelectionPanel(feature_display.getEntryGroup(), systematic_names,
                     description, false);
      systematicListSelectionPanel.setVisible(false);
      final JButton advanced = new JButton("Name Options >>");
      advanced.setToolTipText("Define the qualifier list to search for the feature name");
      advanced.addActionListener(new ActionListener(){

        public void actionPerformed(ActionEvent arg0)
        {
          systematicListSelectionPanel.setVisible(
              !systematicListSelectionPanel.isVisible());
          if(systematicListSelectionPanel.isVisible())
            advanced.setText("Name Options <<");
          else
            advanced.setText("Name Options >>");
          pack();
        }
      });
      
      final JPanel pane = new JPanel(new GridBagLayout());
      final GridBagConstraints c = new GridBagConstraints();
      c.gridx = 0;
      c.gridy = 0;
      c.ipady = 10;
      c.anchor = GridBagConstraints.NORTHWEST;
      pane.add(yBox, c);
      
      c.gridy = 1;
      pane.add(advanced, c);
      
      c.gridy = 2;
      c.fill = GridBagConstraints.HORIZONTAL;
      pane.add(systematicListSelectionPanel, c);
      
      Box xBox = Box.createHorizontalBox();
      final JButton okButton = new JButton("OK");
      okButton.addActionListener(new ActionListener(){
        public void actionPerformed(ActionEvent arg0)
        {
          setVisible(false); 
          dispose();
          status = 0;
          
          Options.getOptions().setSystematicQualifierNames(
              systematicListSelectionPanel.getResultString());
        }
      });
      xBox.add(okButton);

      final JButton cancelButton = new JButton("Cancel");
      cancelButton.addActionListener(new ActionListener(){
        public void actionPerformed(ActionEvent arg0)
        {
          setVisible(false); 
          dispose();
          status = -1;
        }
      });
      xBox.add(cancelButton);
      c.fill = GridBagConstraints.NONE;
      c.gridy = 3;
      pane.add(xBox, c);
      
      getContentPane().add(pane);
      pack();
      setLocationRelativeTo(frame);
      setVisible(true);
    }
    
    protected static String getFeatureName(Feature f)
    {
      if(qNames == null)
      {
        qNames = new StringVector();
        Object objs[];
        if(systematicListSelectionPanel != null)
          objs = systematicListSelectionPanel.getResultArray();
        else
          objs = Options.getOptions().getSystematicQualifierNames().toArray();
        for(Object o: objs)
          qNames.add((String)o);
      }
      return pickName(f, qNames);
    }
    
    protected int getStatus()
    {
      return status;
    }
    
    /**
     *  Look at the qualifier_names one-by-one and return the first value of the
     *  first qualifier found.
     **/
    private static String pickName(final Feature f, final StringVector qualifier_names) 
    {
      int qn_size = qualifier_names.size();
      for(int i = 0; i < qn_size; ++i)
      {
        try
        {
          final Qualifier qualifier =
            f.getQualifierByName((String)qualifier_names.elementAt(i));

          if(qualifier != null)
          {
            final StringVector values = qualifier.getValues();

            if(values != null && values.size() > 0)
            {
              for(int j=0; j<values.size(); j++)
              {
                final String value = (String)values.elementAt(j);
                if(value != null && !value.endsWith("current=false") && !value.equals(""))
                  return value;
              }
            }
          }
        } 
        catch(InvalidRelationException e){}
      }

      return f.getIDString();
    }
  }