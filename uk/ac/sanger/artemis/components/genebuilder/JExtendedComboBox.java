/* JExtendedComboBox.java
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

import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.util.Vector;

import javax.swing.Box;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JList;
import javax.swing.JOptionPane;
import javax.swing.JPopupMenu;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.ListCellRenderer;
import javax.swing.border.EmptyBorder;
import javax.swing.event.PopupMenuEvent;
import javax.swing.event.PopupMenuListener;
import javax.swing.text.JTextComponent;

import org.gmod.schema.cv.CvTerm;


/**
 * JComboBox with horizontal scrollbar
 */
public class JExtendedComboBox extends JComboBox
{ 
  private static final long serialVersionUID = 1L;
  public static String SEPARATOR = "SEPARATOR";
  private int current;
  private boolean highLightCurrent = false;
  
  /**
   * @param str
   * @param useAutoComplete  set true to use auto-complete
   */
  public JExtendedComboBox(String str[], final boolean useAutoComplete)
  { 
    super(str);
    init(useAutoComplete);
  } 
  
  /**
   * @param str
   * @param useAutoComplete  set true to use auto-complete
   */
  public JExtendedComboBox(Vector<?> vector, final boolean useAutoComplete)
  { 
    super(vector);
    init(useAutoComplete);
  } 
  
  /**
   * @param str
   * @param useAutoComplete  set true to use auto-complete
   */
  public JExtendedComboBox(final boolean useAutoComplete)
  { 
    super();
    init(useAutoComplete);
  } 

  public JExtendedComboBox(String str[])
  { 
    this(str, false);
  } 
  
  public JExtendedComboBox(Vector<?> vector)
  { 
    this(vector, false);
  }

  /**
   * Set up renderer, auto-complete and horizontal scroll bar
   * @param useAutoComplete
   */
  private void init(final boolean useAutoComplete)
  {
    setRenderer(new ComboBoxRenderer());
    
    if(useAutoComplete)
    {
      setEditable(true);
      JTextComponent editor = (JTextComponent) getEditor().getEditorComponent();
      editor.setDocument(new AutoCompleteComboDocument(this));
    }
    
    addPopupMenuListener(new ComboPopupMenuLister());
    //setUI(new ComboUI());
  }
  
  class ComboPopupMenuLister implements PopupMenuListener
  {
    public void popupMenuCanceled(PopupMenuEvent e){}
    public void popupMenuWillBecomeInvisible(PopupMenuEvent e){}
    
    public void popupMenuWillBecomeVisible(PopupMenuEvent e)
    {
      Object comp = getUI().getAccessibleChild(JExtendedComboBox.this, 0);
      if (!(comp instanceof JPopupMenu)) 
        return;

      JComponent scrollPane = (JComponent) ((JPopupMenu) comp).getComponent(0);
      if (scrollPane instanceof JScrollPane) 
      {
        JScrollPane sp = (JScrollPane) scrollPane;
        sp.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_AS_NEEDED);
      }
    }
  }
  
  class ComboBoxRenderer extends JLabel implements ListCellRenderer 
  {
    private static final long serialVersionUID = 1L;
    private JSeparator separator;

    public ComboBoxRenderer() 
    {
      setOpaque(true);
      setBorder(new EmptyBorder(1, 1, 1, 1));
      separator = new JSeparator(JSeparator.HORIZONTAL);
    }

    public Component getListCellRendererComponent(
        final JList list, final Object value,
        final int index,  final boolean isSelected,
        final boolean cellHasFocus) 
    {
      String str;
      if(value instanceof String || value == null)
        str = (value == null) ? "" : value.toString();
      else
        str = ((CvTerm)value).getName();

      if (SEPARATOR.equals(str))
        return separator;
 
      if (isSelected) 
      {
        setBackground(list.getSelectionBackground());
        setForeground(list.getSelectionForeground());
      } 
      else 
      {
        setBackground(list.getBackground());
        setForeground(list.getForeground());
      }
      
      Font f = list.getFont();
      if(isHighLightCurrent() && index == getCurrent())
        f = f.deriveFont(Font.BOLD);

      setFont(f);
      setText(str);
      return this;
    }
  }
  
  
  /*public class ComboUI extends BasicComboBoxUI
  {
    protected ComboPopup createPopup()
    {
      return new BasicComboPopup(comboBox)
      {
        private static final long serialVersionUID = 1L;

        protected JScrollPane createScroller() 
        {
          return new JScrollPane( list, ScrollPaneConstants.VERTICAL_SCROLLBAR_AS_NEEDED,
                  ScrollPaneConstants.HORIZONTAL_SCROLLBAR_AS_NEEDED );
        }
      };
    }
  }*/

  public int getCurrent()
  {
    return current;
  }

  public void setCurrent(int current)
  {
    this.current = current;
  }

  public boolean isHighLightCurrent()
  {
    return highLightCurrent;
  }

  public void setHighLightCurrent(boolean highLightCurrent)
  {
    this.highLightCurrent = highLightCurrent;
  }
  
  public static void main(String args[])
  {
    /*
    uk.ac.sanger.artemis.components.database.DatabaseEntrySource entry_source = 
      new uk.ac.sanger.artemis.components.database.DatabaseEntrySource();
    entry_source.setLocation(true);
    uk.ac.sanger.artemis.util.DatabaseDocument doc = entry_source.getDatabaseDocument();
    Vector terms = doc.getCvTermsByCvName(
        uk.ac.sanger.artemis.util.DatabaseDocument.PRODUCTS_TAG_CVNAME);
    Collections.sort(terms, 
        new uk.ac.sanger.artemis.components.genebuilder.cv.CvTermsComparator());
    */
    final String options[] = { "<PREV", "CANCEL", "NEXT>"};   
    
    String[] terms = { "aaaa", "bbbb", "cccc" };
    JExtendedComboBox term_list = new JExtendedComboBox(terms);
    term_list.setCurrent(0);
    term_list.setHighLightCurrent(true);
    
    Box xbox = Box.createHorizontalBox();
    xbox.add(term_list);
    
    Dimension d = new Dimension(70,term_list.getPreferredSize().height);
    term_list.setPreferredSize(d);
    term_list.setMaximumSize(d);
   
    JOptionPane.showOptionDialog(null, xbox,
        "CV term selection",
         JOptionPane.YES_NO_CANCEL_OPTION,
         JOptionPane.QUESTION_MESSAGE,
         null,
         options,
         options[2]);
  }
}