/* AutoCompleteComboDocument.java
 * from http://www.orbital-computer.de/JComboBox/
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

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.FocusAdapter;
import java.awt.event.FocusEvent;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;

import javax.swing.JComboBox;
import javax.swing.text.AttributeSet;
import javax.swing.text.BadLocationException;
import javax.swing.text.JTextComponent;
import javax.swing.text.PlainDocument;

import org.gmod.schema.cv.CvTerm;

  public class AutoCompleteComboDocument extends PlainDocument 
  {
    private static final long serialVersionUID = 1L;
    private JComboBox comboBox;

    private JTextComponent editor;
    // flag to indicate if setSelectedItem has been called
    // subsequent calls to remove/insertString should be ignored
    boolean selecting=false;
    boolean hidePopupOnFocusLoss;
    boolean hitBackspace=false;
    boolean hitBackspaceOnSelection;

    
    public AutoCompleteComboDocument(final JComboBox comboBox)
    {
      this.comboBox = comboBox;

      editor = (JTextComponent) comboBox.getEditor().getEditorComponent();
      editor.setDocument(this);
      
      if(comboBox.getModel().getSize() <= comboBox.getMaximumRowCount())
        comboBox.setMaximumRowCount(comboBox.getModel().getSize()-1);

      comboBox.addActionListener(new ActionListener() 
      {
          public void actionPerformed(ActionEvent e) 
          {
              if (!selecting) highlightCompletedText(0);
          }
      });
      
      editor.addKeyListener(new KeyAdapter() 
      {
        public void keyPressed(KeyEvent e) 
        {
          if(comboBox.isDisplayable()) 
            comboBox.setPopupVisible(true);
          hitBackspace=false;
          switch (e.getKeyCode()) 
          {
            // determine if the pressed key is backspace (needed by the remove method)
            case KeyEvent.VK_BACK_SPACE : 
              hitBackspace=true;
              hitBackspaceOnSelection=editor.getSelectionStart()!=editor.getSelectionEnd();
              break;
          }
        }
      });
      
      // Bug 5100422 on Java 1.5: Editable JComboBox won't hide popup when tabbing out
      hidePopupOnFocusLoss=System.getProperty("java.version").startsWith("1.5");
      // Highlight whole text when gaining focus
      editor.addFocusListener(new FocusAdapter() 
      {
          public void focusGained(FocusEvent e) 
          {
              highlightCompletedText(0);
          }
          public void focusLost(FocusEvent e) 
          {
            // Workaround for Bug 5100422 - Hide Popup on focus loss
            if (hidePopupOnFocusLoss) comboBox.setPopupVisible(false);
          }
      });
      // Handle initially selected object
      Object selected = comboBox.getSelectedItem();
      if(selected!=null)
        setText(getStringValue(selected));
      highlightCompletedText(0);
    }

    public void remove(int offs, int len) throws BadLocationException 
    {
      // return immediately when selecting an item
      if (selecting) 
        return;
      if (hitBackspace)
      {
        // user hit backspace => move the selection backwards
        // old item keeps being selected
        if(offs>0)
          if(hitBackspaceOnSelection) 
            offs--;
        else 
        {
          // User hit backspace with the cursor positioned on the start => beep
          comboBox.getToolkit().beep(); // when available use: UIManager.getLookAndFeel().provideErrorFeedback(combo
        }
        highlightCompletedText(offs);
      }
      else
        super.remove(offs, len);
    }

    public void insertString(int offs, String str, AttributeSet a)
        throws BadLocationException
    {
      // return immediately when selecting an item
      if(selecting)
        return;
      super.insertString(offs, str, a);
      
      // lookup and select a matching item
      Object item = lookupItem(getText(0, getLength()));

      boolean listContainsSelectedItem = true;
      if(item == null)
      {
        item = comboBox.getModel().getSelectedItem();
        listContainsSelectedItem = false;
      }
      setSelectedItem(item);
      setText(getStringValue(item));
      // select the completed part
      if(listContainsSelectedItem)
        highlightCompletedText(offs + str.length());
    }

    private void setText(String text)
    {
      try
      {
        super.remove(0, getLength());
        super.insertString(0, text, null);
      }
      catch(BadLocationException e)
      {
        throw new RuntimeException(e.toString());
      }
    }

    private void highlightCompletedText(int start)
    {
      editor.setCaretPosition(getLength());
      editor.moveCaretPosition(start);
    }

    private void setSelectedItem(Object item)
    {
      selecting = true;
      comboBox.getModel().setSelectedItem(item);
      selecting = false;
    }

    private Object lookupItem(String pattern)
    {
      Object selectedItem = comboBox.getModel().getSelectedItem();
      String selectedItemStr = getStringValue(selectedItem);

      // only search if the currently selected does not match
      if(selectedItemStr != null
          && startsWithIgnoreCase(selectedItemStr, pattern))
        return selectedItem;
      else
      {
        // iterate over all items
        for(int i = 0, n = comboBox.getModel().getSize(); i < n; i++)
        {
          Object currentItem = comboBox.getModel().getElementAt(i);
          String currentItemStr = getStringValue(currentItem);

          // current item starts with the pattern?
          if(startsWithIgnoreCase(currentItemStr, pattern))
            return currentItem;
        }
      }

      return null;
    }
    
    private String getStringValue(Object item)
    {
      String itemStr = null;
      if(item != null)
      {
        if(item instanceof String)
          itemStr = item.toString();
        else
          itemStr = ((CvTerm)item).getName();
      }
      return itemStr;
    }

    private boolean startsWithIgnoreCase(String str1, String str2)
    {
      return str1.toUpperCase().startsWith(str2.toUpperCase());
    }
  }