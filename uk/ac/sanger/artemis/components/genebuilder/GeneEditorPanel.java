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
  
  private OpenSectionButton coreButton;
  private OpenSectionButton cvButton;
  private OpenSectionButton matchButton;
  private OpenSectionButton propertiesButton;
  
  private QualifierTextArea qualifier_text_area;
  private CVPanel cvForm;
  private MatchPanel matchForm;
  private GffPanel propertiesPanel;

  /**
   * Gene editor panel - showing annotation in a single panel.
   * @param qualifier_text_area
   * @param cvForm
   * @param matchForm
   * @param propertiesPanel
   */
  public GeneEditorPanel(final QualifierTextArea qualifier_text_area,
                         final CVPanel cvForm, 
                         final MatchPanel matchForm,
                         final GffPanel propertiesPanel)
  {
    this.qualifier_text_area = qualifier_text_area;
    this.cvForm = cvForm;
    this.matchForm = matchForm;
    this.propertiesPanel = propertiesPanel;
    
    setLayout( new BoxLayout(this, BoxLayout.PAGE_AXIS) );
    setBackground(Color.WHITE);
    JScrollPane jspCore = new JScrollPane(qualifier_text_area);
    jspCore.setPreferredSize(new Dimension(jspCore.getPreferredSize().width, 100));
    
    addDarkSeparator(this);
    propertiesButton = addOpenClosePanel("Properties", propertiesPanel);
    add(propertiesPanel);
    
    addDarkSeparator(this);
    coreButton = addOpenClosePanel("Core",jspCore);
    add(jspCore);
    
    addDarkSeparator(this);
    cvButton = addOpenClosePanel("Controlled Vocabulary", cvForm);
    add(cvForm);
    
    addDarkSeparator(this);
    matchButton = addOpenClosePanel("Match", matchForm);
    add(matchForm);
    
    add(Box.createVerticalGlue());
  }

  /**
   * Open/close the sections if they contain elements or
   * are empty.
   */
  public void updatePanelState()
  {
    if(qualifier_text_area.getText().equals(""))
      coreButton.setOpen(false);
    else
      coreButton.setOpen(true);
    
    if(cvForm.isEmpty())
      cvButton.setOpen(false);
    else
      cvButton.setOpen(true);
    
    if(matchForm.isEmpty())
      matchButton.setOpen(false);
    else
      matchButton.setOpen(true);
    
    if(propertiesPanel.isEmpty())
      propertiesButton.setOpen(false);
    else
      propertiesButton.setOpen(true);
  }
  
  /**
   * Add a Separator to a given component
   * @param comp
   * @param parent
   * @param useLightColor
   */
  public static JSeparator getSeparator(final JComponent comp,
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

    return separator;
  }
  
  /**
   * Add a light Separator to a given component
   * @param comp
   * @param parent
   */
  public static void addLightSeparator(final JComponent comp)
  {
    JSeparator separator = getSeparator(comp, true);
    comp.add(separator);
  }
  
  /**
   * Add a dark Separator to a given component
   * @param comp
   * @param parent
   */
  public static void addDarkSeparator(final JComponent comp)
  {
    JSeparator separator = getSeparator(comp, false);
    comp.add(separator);
  }
  
  private OpenSectionButton addOpenClosePanel(final String name,
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
    
    final OpenSectionButton openButton = new OpenSectionButton("-", panel);

    bannerPanel.add(nameLabel);
    bannerPanel.add(Box.createHorizontalGlue());
    bannerPanel.add(openButton);
    bannerPanel.setPreferredSize(
       new Dimension(bannerPanel.getPreferredSize().width, 20));
    
    add(bannerPanel);
    return openButton;
  }
  
  public class OpenSectionButton extends JButton 
  {
    private static final long serialVersionUID = 1L;
    private JComponent panel;
    
    public OpenSectionButton(final String text, final JComponent panel)
    {
      super(text);
      this.panel = panel;

      Dimension size = new Dimension(35,20);
      setForeground(STEEL_BLUE);
      Font font = getFont();
      font = font.deriveFont(Font.BOLD);
      setFont(font);
      setBorderPainted(false);
      setOpaque(false);
      setPreferredSize(size);
      setMaximumSize(size);

      addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          updateButton();
        }  
      });
    }
    
    private void updateButton()
    {
      if(getText().equals("-"))
      {
        setText("+");
        panel.setVisible(false);
      }
      else
      {
        setText("-");
        panel.setVisible(true);
      }
    }
    
    public void setOpen(final boolean isOpen)
    {
      if(isOpen && getText().equals("-"))
        return;
      if(!isOpen && getText().equals("+"))
        return;
      updateButton();
    }

  }
}