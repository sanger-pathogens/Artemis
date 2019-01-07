/* 
 * This file is part of Artemis
 *
 * Copyright(C) 2013  Genome Research Limited
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

import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSeparator;
import javax.swing.JSlider;
import javax.swing.event.ChangeEvent;

import uk.ac.sanger.artemis.circular.TextFieldFloat;

/**
 * Enable fine tuning of the heights of the individual panels
 * (plots, BAM view, VCF view, comparisons).
 */
class ActPanelResizer extends JFrame
{
  private static final long serialVersionUID = 1L;
  private GridBagConstraints c = new GridBagConstraints();
  
  ActPanelResizer(final MultiComparator comparator, final GridBagLayout layout)
  {
    setTitle("Panel Height Weights");
    final JPanel panel = (JPanel) getContentPane();
    panel.setLayout(new GridBagLayout());
    c.anchor = GridBagConstraints.WEST;

    for(int i = 0 ; i < comparator.getEntryGroupArray().length ; i++)
      showWeighty(comparator, layout, i);
    
    c.gridx = 0;
    c.gridy+=1;
    final JButton close = new JButton("CLOSE");
    close.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent arg0)
      {
        dispose();
      }
    });
    getContentPane().add(close, c);
    
    pack();
    Utilities.centreFrame(this);
    setVisible(true);
    setDefaultCloseOperation(DISPOSE_ON_CLOSE);
  }
  
  /**
   * Display the panel height weighting for an entry
   * @param comp
   * @param layout
   * @param idx - entry index
   */
  private void showWeighty(final MultiComparator comp, final GridBagLayout layout, int idx)
  {
    if(idx > 0)
      addSeparator();
    
    c.gridx = 0;
    c.gridy+=1;
    final JLabel l = new JLabel("SEQUENCE "+(idx+1));
    l.setFont(l.getFont().deriveFont(Font.BOLD));
    getContentPane().add(l, c);

    // plots
    if(comp.getBasePlotGroupArray()[idx].getVisibleCount() > 0)
      showComponentWgtY(comp, layout, comp.getBasePlotGroupArray()[idx], "Plots");
    
    // bams
    if(comp.getBamPanelArray()[idx].isVisible())
      showComponentWgtY(comp, layout, comp.getBamPanelArray()[idx], "BAM view");
    
    // vcfs
    if(comp.getVcfPanelArray()[idx].isVisible())
      showComponentWgtY(comp, layout, comp.getVcfPanelArray()[idx], "VCF view");

    if (idx < comp.getAlignmentViewerArray().length )
    {
      addSeparator();
      // comparison panel
      showComponentWgtY(comp, layout, comp.getAlignmentViewerArray()[idx], "Comparison");
    }
  }
  

  /**
   * Display the panel height weight for a component
   * @param comparator
   * @param layout
   * @param component
   * @param label
   */
  private void showComponentWgtY(final MultiComparator comparator, 
                                 final GridBagLayout layout, 
                                 final JComponent component, 
                                 final String label)
  {
    c.gridx = 0;
    c.gridy+=1;

    getContentPane().add(new JLabel(label), c);
    final Dimension d = new Dimension(80,25);
    final GridBagConstraints cc = layout.getConstraints(component);
    final TextFieldFloat wt = new TextFieldFloat();
    final JSlider slider = new JSlider(0, 10, (int) (cc.weighty * 10));
    slider.setToolTipText("set the y-weight to between 0 - 1");
    wt.setPreferredSize(d);
    wt.setMaximumSize(d);
    wt.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent arg0)
      {
        slider.setValue((int) (wt.getValue()*10));
        cc.weighty = wt.getValue();
        layout.setConstraints(component, cc);
        component.revalidate();
        comparator.validate();
        comparator.repaint();
      }
    });

    slider.addChangeListener(new javax.swing.event.ChangeListener()
    {
      public void stateChanged(ChangeEvent e)
      {
        double y = slider.getValue() / 10.d;
        wt.setValue(y);
        cc.weighty = y;
        layout.setConstraints(component, cc);
        component.revalidate();
        comparator.validate();
        comparator.repaint();
      }
    });
    
    c.gridx+=1;
    getContentPane().add(slider, c);
    wt.setValue(cc.weighty);
    c.gridx+=1;
    getContentPane().add(wt, c);
  }
  
  /**
   * Add a separator between rows
   */
  private void addSeparator()
  {
    c.gridx = 1;
    c.gridy+=1;
    c.fill = GridBagConstraints.HORIZONTAL;
    getContentPane().add(new JSeparator(), c);
  }
}