/* SelectionMenu.java
 *
 * created: Sun Jan 10 1999
 *
 * This file is part of Artemis
 * 
 * Copyright(C) 1998,1999,2000,2001,2002  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/SelectionMenu.java,v 1.9 2008-06-11 15:17:43 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.sequence.MarkerRange;
import uk.ac.sanger.artemis.*;

import java.awt.event.*;
import javax.swing.*;

import java.awt.*;
import java.util.Vector;

/**
 *  This a super class for EditMenu, ViewMenu and GotoMenu.  It is a JMenu that
 *  knows how to get hold of the Selection.  It also has the method
 *  getParentFrame() to find the owning JFrame of the menu.
 *
 *  @author Kim Rutherford
 *  @version $Id: SelectionMenu.java,v 1.9 2008-06-11 15:17:43 tjc Exp $
 **/

public class SelectionMenu extends JMenu 
{

  /** */
  private static final long serialVersionUID = 1L;

  /** The Selection that was passed to the constructor. */
  private Selection selection;

  /** The JFrame reference that was passed to the constructor. */
  private JFrame frame = null;

  /**
   *  Create a new SelectionMenu object.
   *  @param frame The JFrame that owns this JMenu.
   *  @param name The string to use as the menu title.
   *  @param selection The Selection that the commands in the menu will
   *    operate on.
   **/
  public SelectionMenu(final JFrame frame,
                        final String menu_name,
                        final Selection selection) 
  {
    super(menu_name);
    this.frame = frame;
    this.selection = selection;
  }

  /**
   *  Return the JFrame that was passed to the constructor.
   **/
  public JFrame getParentFrame() 
  {
    return frame;
  }

  /**
   *  Check that the are some Features in the current selection.  If there are
   *  some features then return true.  If there are no features selected then
   *  popup a message saying so and return false.
   **/
  protected boolean checkForSelectionFeatures() 
  {
    return checkForSelectionFeatures(getParentFrame(), getSelection());
  }

  /**
   *  Check that the are some Features in the given selection.  If there are
   *  some features then return true.  If there are no features selected then
   *  popup a message saying so and return false.
   *  @param frame The JFrame to use for MessageDialog components.
   *  @param selection The selected features to check.
   **/
  static boolean checkForSelectionFeatures(final JFrame frame,
                                           final Selection selection) 
  {
    final FeatureVector features_to_check = selection.getAllFeatures();

    if(features_to_check.size() == 0) 
    {
      new MessageDialog(frame, "No features selected");
      return false;
    } 
    else 
      return true;
  }

  /**
   *  Check that the are some Features in the current selection and that the
   *  selection isn't too big.  If there is more than one feature in the
   *  selection and less than the given maximum then return true.  If there
   *  are no features selected then popup a message saying so and return
   *  false.  If there are more selected features than the given maximum
   *  than display the message in a YesNoDialog component and return the
   *  result of the dialog.
   **/
  protected boolean checkForSelectionFeatures(final int maximum_size,
                                              final String message) 
  {
    return checkForSelectionFeatures(getParentFrame(), getSelection(),
                                      maximum_size, message);
  }

  /**
   *  Check that the are some Features in the given selection and that the
   *  selection isn't too big.  If there is more than one feature in the
   *  selection and less than the given maximum then return true.  If there
   *  are no features selected then popup a message saying so and return
   *  false.  If there are more selected features than the given maximum
   *  than display the message in a YesNoDialog component and return the
   *  result of the dialog.
   *  @param frame The JFrame to use for MessageDialog components.
   *  @param selection The selected features to check.
   **/
  static boolean checkForSelectionFeatures(final JFrame frame,
                                           final Selection selection,
                                           final int maximum_size,
                                           final String message) 
  {
    final FeatureVector features_to_check = selection.getAllFeatures();

    if(features_to_check.size() == 0) 
    {
      new MessageDialog(frame, "No features selected");
      return false;
    }
    else 
    {
      if(features_to_check.size() > maximum_size) 
      {
        final YesNoDialog dialog =
          new YesNoDialog(frame, message);

        return dialog.getResult();
      }
      else 
        return true;
    }
  }

  /**

   *  Check that there are only CDS features in the given selection, that there
   *  are some features selected and that the selection isn't too big.  If
   *  there is more than one feature in the selection and less than the given
   *  maximum then return true.  If there are no features selected then popup
   *  a message saying so and return false.  If there are more selected
   *  features than the given maximum than display the message in a
   *  YesNoDialog component and return the result of the dialog.
   *  @param frame The Frame to use for MessageDialog components.
   *  @param selection The selected features to check.
   **/
  static boolean checkForSelectionCDSFeatures(final JFrame frame,
                                               final Selection selection,
                                               final int maximum_size,
                                               final String max_message) 
  {
    final FeatureVector features_to_check = selection.getAllFeatures();

    if(features_to_check.size() == 0) 
    {
      new MessageDialog(frame, "No CDS features selected");
      return false;
    }
    else
    {
      for(int i = 0 ; i < features_to_check.size() ; ++i) 
      {
        final Feature this_feature = features_to_check.elementAt(i);
        if(!this_feature.isCDS()) 
        {
          final String message =
            "a non-CDS feature(" + this_feature.getIDString() +
            ") is selected - can't continue";
          final MessageDialog dialog =
            new MessageDialog(frame, message);
          return false;
        }
      }

      if(features_to_check.size() > maximum_size) 
      {
        final YesNoDialog dialog =
          new YesNoDialog(frame, max_message);

        return dialog.getResult();
      }
      else 
        return true;
    }
  }

  /**
   *  Check that the are some FeatureSegments in the current selection and
   *  that the selection isn't too big.  If there is more than one
   *  FeatureSegments in the selection and less than the given maximum then
   *  return true.  If there are no FeatureSegments selected then popup a
   *  message saying so and return false.  If there are more selected
   *  FeatureSegments than the given maximum than display the message in a
   *  YesNoDialog component and return the result of the dialog.
   **/
  protected boolean checkForSelectionFeatureSegments(final int maximum_size,
                                                     final String message) 
  {
    final FeatureSegmentVector segments_to_check =
      getSelection().getAllSegments();

    if(segments_to_check.size() == 0) 
    {
      new MessageDialog(getParentFrame(), "No exons selected");
      return false;
    } 
    else 
    {
      if(segments_to_check.size() > maximum_size) 
      {
        final YesNoDialog dialog =
          new YesNoDialog(getParentFrame(), message);

        return dialog.getResult();
      } 
      else 
        return true;
    }
  }

  /**
   *  Check that the current selection contains a MarkerRange.  If it does
   *  then return true.  If not then popup a message saying so and return
   *  false.
   **/
  protected boolean checkForSelectionRange() 
  {
    return checkForSelectionRange(getParentFrame(), getSelection());
  }

  /**
   *  Check that the current selection contains a MarkerRange.  If it does
   *  then return true.  If not then popup a message saying so and return
   *  false.
   *  @param frame The JFrame to use for MessageDialog components.
   *  @param selection The selected features to check.
   **/
  public static boolean checkForSelectionRange(final JFrame frame,
                                        final Selection selection) 
  {
    final MarkerRange marker_range = selection.getMarkerRange();

    if(marker_range == null) 
    {
      new MessageDialog(frame, "No bases selected");
      return false;
    } 
    else 
      return true;
  }

  /**
   *  Return the Selection object that was passed to the constructor.
   **/
  protected Selection getSelection() 
  {
    return selection;
  }
  
  /**
   *
   **/
  protected static KeyStroke makeMenuKeyStroke(final int key_code)
  {
    return KeyStroke.getKeyStroke(key_code, 
                  Toolkit.getDefaultToolkit().getMenuShortcutKeyMask()); //InputEvent.CTRL_MASK);
  }

  private void getJMenuItems(JMenu menu, Vector<JMenuItem> menu_items)
  {
    final Component menus[] = menu.getMenuComponents();
    for(int i=0; i<menus.length; i++)
    {
      if(menus[i] instanceof JMenuItem)
        menu_items.add((JMenuItem) menus[i]);
    }
  }
  
  private Vector<JMenuItem> getSelectionMenuItems()
  {
    final Vector<JMenuItem> menu_items = new Vector<JMenuItem>();
    final Component menus[] = getMenuComponents();
    for(int i=0; i<menus.length; i++)
    {
      if(menus[i] instanceof JMenu)
      {
        menu_items.add((JMenu)menus[i]);
        getJMenuItems((JMenu)menus[i], menu_items);
      }
      else if(menus[i] instanceof JMenuItem)
        menu_items.add((JMenuItem)menus[i]);
    }
    return menu_items;
  }

  protected Vector<String> getUsedShortCutKeys()
  {
    final Vector<String> sc = new Vector<String>();
    final Vector<JMenuItem> menu_items = getSelectionMenuItems();
    for(JMenuItem mi : menu_items)
      if(mi.getAccelerator() != null)
        sc.add(getKeyText(mi.getAccelerator().getKeyCode()));
    return sc;
  }
  
  protected JScrollPane getShortCuts(Vector<String> usedShortCutKeys)
  {
    usedShortCutKeys.add("C"); // Create menu
    usedShortCutKeys.add("G"); // Goto menu
    final Vector<JMenuItem> menu_items = getSelectionMenuItems();
    final Box bdown = Box.createVerticalBox();
    final String alist[] = 
             { "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K",
               "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V",
               "W", "X", "Y", "Z", "--"};
    final String mod_list[] = { "Default", "Alt", "Ctrl", "Shift" };
    
    for(final JMenuItem mi : menu_items)
    {
      final Box bacross = Box.createHorizontalBox();
      if(mi instanceof JMenu)
      {
        bacross.add(new JLabel( mi.getText() ));
        bacross.add(Box.createHorizontalGlue());
        bdown.add(bacross);
        continue;
      }

      bacross.add(new JLabel(mi.getText()));

// short cut
      final JComboBox combo = new JComboBox(alist);
      combo.setRenderer(new ComboBoxRenderer(usedShortCutKeys));
      final JComboBox mod_combo = new JComboBox(mod_list);
      final KeyStroke ks = mi.getAccelerator();
      if(ks != null)
      {
        combo.setSelectedItem( getKeyText(ks.getKeyCode()) );
        mod_combo.setSelectedItem( getModifierFromInt(ks.getModifiers()) );
      }
      else
      {
        combo.setSelectedItem("--");
        mod_combo.setSelectedItem("Default");
      }

      Dimension dim = combo.getPreferredSize();
      dim = new Dimension(100, dim.height);
      combo.setPreferredSize(dim);
      combo.setMaximumSize(dim);
      combo.addItemListener(new ItemListener()
      {
        public void itemStateChanged(ItemEvent e) 
        {
          if(e.getStateChange() == ItemEvent.SELECTED)
            setAccelerator(combo,mod_combo,mi);
        }
      });

// modifier
      mod_combo.setPreferredSize(dim);
      mod_combo.setMaximumSize(dim);
      mod_combo.addItemListener(new ItemListener()
      {
        public void itemStateChanged(ItemEvent e)
        {
          if(e.getStateChange() == ItemEvent.SELECTED)
            setAccelerator(combo,mod_combo,mi);
        }
      });

      bacross.add(Box.createHorizontalGlue());
      bacross.add(combo);
      bacross.add(mod_combo);
      bdown.add(bacross);
    }

    bdown.add(Box.createVerticalGlue());

    final JScrollPane jsp = new JScrollPane(bdown);
    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
    jsp.setPreferredSize(new Dimension(jsp.getPreferredSize().width,
                         screen.height/2));
    return jsp;
  }

  private void setAccelerator(final JComboBox combo, final JComboBox mod_combo, 
                              final JMenuItem mi)
  {
    if(combo.getSelectedItem() == "--")
      mi.setAccelerator(null);
    else
    {
      char c[] = ((String)combo.getSelectedItem()).toCharArray();
      int modifier = getModifierFromString((String)mod_combo.getSelectedItem()); 
      mi.setAccelerator(KeyStroke.getKeyStroke(c[0], modifier));
    }

    if(ShortCut.usingCache())
      new ShortCut(getText(), mi.getText(), mi.getAccelerator());
  }
  
  private int getModifierFromString(String modStr)
  {
    int modifier = Toolkit.getDefaultToolkit().getMenuShortcutKeyMask();
    if( modStr.equals("Alt") )
      modifier = InputEvent.ALT_MASK;
    else if( modStr.equals("Ctrl") )
      modifier = InputEvent.CTRL_MASK;
    else if( modStr.equals("Shift") )
      modifier = InputEvent.SHIFT_MASK;
    return modifier;
  }
  
  private String getModifierFromInt(int mod)
  {
    if( (InputEvent.ALT_MASK & mod) == InputEvent.ALT_MASK)
      return "Alt";
    else if( (InputEvent.CTRL_MASK & mod) == InputEvent.CTRL_MASK)
      return "Ctrl";
    else if( (InputEvent.SHIFT_MASK & mod) == InputEvent.SHIFT_MASK)
      return "Shift";
    return "Default";
  }
  
  public static String getKeyText(int keyCode)
  {
    if(keyCode >= KeyEvent.VK_0 && keyCode <= KeyEvent.VK_9 ||
       keyCode >= KeyEvent.VK_A && keyCode <= KeyEvent.VK_Z) 
      return String.valueOf((char)keyCode);
    
    switch(keyCode) 
    {
      case KeyEvent.VK_BACK_SPACE: return "BACK_SPACE";
      case KeyEvent.VK_CONTROL: return "CONTROL";
      case KeyEvent.VK_LEFT: return "LEFT";
      case KeyEvent.VK_UP: return "UP";
      case KeyEvent.VK_RIGHT: return "RIGHT";
      case KeyEvent.VK_DOWN: return "DOWN";
    
      case KeyEvent.VK_DELETE: return "DELETE";
    
      case KeyEvent.VK_BACK_QUOTE: return "BACK_QUOTE";
    
      case KeyEvent.VK_KP_UP: return "KP_UP";
      case KeyEvent.VK_KP_DOWN: return "KP_DOWN";
      case KeyEvent.VK_KP_LEFT: return "KP_LEFT";
      case KeyEvent.VK_KP_RIGHT: return "KP_RIGHT";
    }
    
    return "unknown(0x" + Integer.toString(keyCode, 16) + ")";
  }
  
  /**
   * True if this menu has editable shortcuts
   * @return
   */
  protected boolean isEditableShortCutMenu()
  {
    if(this instanceof SelectMenu ||
       this instanceof EditMenu ||
       this instanceof ViewMenu )
      return true;
    return false;
  }
  
  /**
   * Override to apply any cached shortcuts
   */
  public JMenuItem add(JMenuItem menuItem)
  {
    JMenuItem mi = super.add(menuItem);
    if(isEditableShortCutMenu() && ShortCut.usingCache())
      ShortCut.applyShortCutFromCache(getText(), menuItem);
    return mi;
  }
  
}

class ComboBoxRenderer implements ListCellRenderer 
{
  protected DefaultListCellRenderer defaultRenderer = new DefaultListCellRenderer();
  private Vector<String> usedShortCutKeys;
  public ComboBoxRenderer(Vector<String> usedShortCutKeys)
  {
    this.usedShortCutKeys = usedShortCutKeys;
  }
  
  public Component getListCellRendererComponent(JList list, Object value, int index,
      boolean isSelected, boolean cellHasFocus) 
  {
    JLabel renderer = (JLabel) defaultRenderer.getListCellRendererComponent(list, value, index,
        isSelected, cellHasFocus);
    if (usedShortCutKeys.contains(value)) 
      renderer.setBackground(Color.lightGray);
    return renderer;
  }
}

/**
 * Editable shortcut sub-menu
 */
class SelectionSubMenu extends JMenu
{
  private static final long serialVersionUID = 1L;
  private String parentMenuStr;
  public SelectionSubMenu(SelectionMenu parentMenu, String str)
  {
    super(str);
    parentMenuStr = parentMenu.getText();
  }

  public JMenuItem add(JMenuItem menuItem)
  {
    JMenuItem mi = super.add(menuItem);
    if(ShortCut.usingCache())
      ShortCut.applyShortCutFromCache(parentMenuStr, menuItem);
    return mi;
  }
}
