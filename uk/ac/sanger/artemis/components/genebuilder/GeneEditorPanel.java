/* GeneEditorPanel.java
 *
 * created: 2007
 *
 * This file is part of Artemis
 *
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

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;

import uk.ac.sanger.artemis.components.QualifierTextArea;
import uk.ac.sanger.artemis.components.genebuilder.cv.CVPanel;
import uk.ac.sanger.artemis.components.genebuilder.gff.GffPanel;
import uk.ac.sanger.artemis.components.genebuilder.ortholog.MatchPanel;


/**
 * Panel for display controlled vocabulary terms for Chado
 */
public class GeneEditorPanel extends JPanel
{

  /** */
  private static final long serialVersionUID = 1L;
  private static Color STEEL_BLUE = new Color(25, 25, 112);
  private static Color LIGHT_STEEL_BLUE = new Color(176, 196, 222);

  /**
   * Gene editor panel - showing annotation in a single panel.
   * @param qualifier_text_area
   * @param cvForm
   * @param matchForm
   * @param gffPanel
   */
  public GeneEditorPanel(final QualifierTextArea qualifier_text_area,
                         final CVPanel cvForm, 
                         final MatchPanel matchForm,
                         final GffPanel gffPanel)
  {
    setLayout( new BoxLayout(this, BoxLayout.PAGE_AXIS) );
    setBackground(Color.WHITE);
    JScrollPane jspCore = new JScrollPane(qualifier_text_area);
    jspCore.setPreferredSize(new Dimension(jspCore.getPreferredSize().width, 100));
    
    addDarkSeparator(this);
    addOpenClosePanel("Core",jspCore);
    add(jspCore);
    
    addDarkSeparator(this);
    addOpenClosePanel("Controlled Vocabulary", cvForm);
    add(cvForm);
    
    addDarkSeparator(this);
    addOpenClosePanel("Match", matchForm);
    add(matchForm);
    
    addDarkSeparator(this);
    addOpenClosePanel("GFF", gffPanel);
    add(gffPanel);
    add(Box.createVerticalGlue());
  }

  /**
   * Add a Separator to a given component
   * @param comp
   * @param parent
   * @param useLightColor
   */
  private static void addSeparator(final JComponent comp,
                                   final boolean useLightColor)
  {
    final JSeparator separator = new JSeparator();

    if(useLightColor)
      separator.setForeground(LIGHT_STEEL_BLUE);
    else
      separator.setForeground(STEEL_BLUE);
    
    //separator.setPreferredSize(new Dimension(width,10));
    // add glue to expand the separator horizontally
    comp.add(Box.createHorizontalGlue());
    comp.add(Box.createVerticalStrut(5));
    separator.setPreferredSize(
        new Dimension(comp.getPreferredSize().width,
                      separator.getPreferredSize().height));
    separator.setMaximumSize(new Dimension(1500,10));
    comp.add(separator);
  }
  
  /**
   * Add a light Separator to a given component
   * @param comp
   * @param parent
   */
  public static void addLightSeparator(final JComponent comp)
  {
    addSeparator(comp, true);
  }
  
  /**
   * Add a dark Separator to a given component
   * @param comp
   * @param parent
   */
  public static void addDarkSeparator(final JComponent comp)
  {
    addSeparator(comp, false);
  }
  
  private void addOpenClosePanel(final String name,
                                 final JComponent panel)
  {
    final JPanel bannerPanel = new JPanel();
    bannerPanel.setLayout(new BoxLayout(bannerPanel, BoxLayout.LINE_AXIS));
    bannerPanel.setBackground(LIGHT_STEEL_BLUE);
    
    final JLabel nameLabel = new JLabel(name);
    nameLabel.setForeground(STEEL_BLUE);
    Font font = nameLabel.getFont();
    font = font.deriveFont(Font.BOLD);
    nameLabel.setFont(font);
    
    Dimension size = new Dimension(35,20);

    final JButton openButton = new JButton("-");
    openButton.setForeground(STEEL_BLUE);
    openButton.setFont(font);
    openButton.setBorderPainted(false);
    openButton.setOpaque(false);
    openButton.setPreferredSize(size);
    openButton.setMaximumSize(size);
    openButton.addActionListener(new ActionListener()
    {

      public void actionPerformed(ActionEvent e)
      {
        if(openButton.getText().equals("-"))
        {
          openButton.setText("+");
          panel.setVisible(false);
        }
        else
        {
          openButton.setText("-");
          panel.setVisible(true);
        }
      }
      
    });
    
    bannerPanel.add(nameLabel);
    bannerPanel.add(Box.createHorizontalGlue());
    bannerPanel.add(openButton);
    bannerPanel.setPreferredSize(
       new Dimension(bannerPanel.getPreferredSize().width, 20));
    
    add(bannerPanel);
  }
}