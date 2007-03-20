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
import java.awt.FontMetrics;
import java.util.Vector;

import javax.swing.JComboBox;
import javax.swing.JScrollPane;
import javax.swing.plaf.basic.BasicComboPopup;

import org.gmod.schema.cv.CvTerm;

/**
 * JComboBox with horizontal scrollbar
 */
public class JExtendedComboBox extends JComboBox
{ 
  /** */
  private static final long serialVersionUID = 1L;

  public JExtendedComboBox(String str[])
  { 
    super(str);
    setHorizontalScrollBar();
  } 
  
  public JExtendedComboBox(Vector vector)
  { 
    super(vector);
    setHorizontalScrollBar();
  } 

  private void setHorizontalScrollBar()
  { 
    FontMetrics fm = getFontMetrics(getFont()); 
    BasicComboPopup popup = (BasicComboPopup)getUI().getAccessibleChild(this,0);//Popup 
    if(popup==null)
      return; 
   
    int size = (int)getPreferredSize().getWidth(); 
    
    for(int i=0; i<getItemCount(); i++)
    { 
      String str;
      if(getItemAt(i) instanceof String)
        str = (String)getItemAt(i); 
      else
        str = ((CvTerm)getItemAt(i)).getName();
      
      if(size<fm.stringWidth(str)) 
        size = fm.stringWidth(str); 
    } 
    
    Component comp = popup.getComponent(0); //JScrollPane 
    JScrollPane scrollpane = null; 
    if(comp instanceof JScrollPane) 
    { 
      scrollpane = (JScrollPane)comp; 
      scrollpane.setHorizontalScrollBarPolicy(JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);
    } 
  }
  
}