/*
 *
 * created: Wed Aug 3 2004
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2000  Genome Research Limited
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

package uk.ac.sanger.artemis.editor;


import java.awt.*;
import java.awt.event.*;
import java.util.*;

import javax.swing.*;
import javax.swing.plaf.ToolTipUI;
import javax.swing.plaf.ComponentUI;

/**
*
* UI for multiple line tooltips
*
*/
public class MultiLineToolTipUI extends ToolTipUI 
{
  static MultiLineToolTipUI SINGLETON = new MultiLineToolTipUI();
  static boolean DISPLAY_ACCELERATOR=true;
  static boolean displayAccelerator;
  int accelerator_offset = 15;
  int inset = 3;
  Graphics g;

  private MultiLineToolTipUI() {}

  public static void initialize() 
  {
    // don't hardcode class name, this way we can obfuscate.
    String key = "ToolTipUI";
    Class cls = SINGLETON.getClass();
    String name = cls.getName();
    UIManager.put(key,name);
    UIManager.put(name,cls);
  }

  public static ComponentUI createUI(JComponent c) 
  {
    return SINGLETON;
  }

  public void installUI(JComponent c) 
  {
    LookAndFeel.installColorsAndFont(c, "ToolTip.background",
      "ToolTip.foreground", "ToolTip.font");
    LookAndFeel.installBorder(c, "ToolTip.border");
  }

  public void uninstallUI(JComponent c) 
  {
    LookAndFeel.uninstallBorder(c);
  }

  public static void setDisplayAcceleratorKey(boolean val) 
  {
    displayAccelerator=val;
  }

  public Dimension getPreferredSize(JComponent c) 
  {
    Font font = c.getFont();
    String tipText = ((JToolTip)c).getTipText();
    MyToolTip mtt = new MyToolTip();
    FontMetrics fontMetrics = mtt.toolTipFontMetrics(font);
    int fontHeight = fontMetrics.getHeight();

    if (tipText == null) 
      tipText = "";
    String lines[] = PlafMacros.breakupLines(tipText);
    int num_lines = lines.length;

    Dimension dimension;
    int width, height, onewidth;
    height = num_lines * fontHeight;
    width = 0;
    for (int i=0; i<num_lines; i++) 
    {
      onewidth = fontMetrics.stringWidth(lines[i]);
      if (displayAccelerator && i == num_lines - 1) 
      {
        String keyText = getAcceleratorString((JToolTip)c);
        if (!keyText.equals(""))
           onewidth += fontMetrics.stringWidth(keyText) 
            + accelerator_offset;
      }
      width = Math.max(width,onewidth);
    }
    return new Dimension(width+inset*2,height+inset*2);
  }

  public Dimension getMinimumSize(JComponent c) 
  {
    return getPreferredSize(c);
  }

  public Dimension getMaximumSize(JComponent c) 
  {
    return getPreferredSize(c);
  }

  public void paint(Graphics g, JComponent c) 
  {
    Font font = c.getFont();
    MyToolTip mtt = new MyToolTip();
    FontMetrics fontMetrics = mtt.toolTipFontMetrics(font);
    Dimension dimension = c.getSize();
    int fontHeight = fontMetrics.getHeight();
    int fontAscent = fontMetrics.getAscent();
    String tipText = ((JToolTip)c).getTipText();
    String lines[] = PlafMacros.breakupLines(tipText);
    int num_lines = lines.length;
    int height;
    int i;

    g.setColor(c.getBackground());
    g.fillRect(0, 0, dimension.width, dimension.height);
    g.setColor(c.getForeground());
    for (i=0, height=2+fontAscent; 
         i<num_lines; i++, height+=fontHeight) 
    {
      g.drawString(lines[i], inset, height);
      if (displayAccelerator && i == num_lines - 1) 
      {
        String keyText = getAcceleratorString((JToolTip)c);
        if (!keyText.equals("")) 
        {
          Font smallFont = new Font(font.getName(), 
            font.getStyle(), font.getSize()-2);
          g.setFont(smallFont);
          g.drawString(keyText, fontMetrics.stringWidth(lines[i]) 
            + accelerator_offset, height);
        }
      }
    }
  }

  public String getAcceleratorString(JToolTip tip) 
  {
    JComponent comp = tip.getComponent();
    if (comp == null)
      return "";
    KeyStroke[] keys =comp.getRegisteredKeyStrokes();
    String controlKeyStr = "";
    KeyStroke postTip=KeyStroke.getKeyStroke(
      KeyEvent.VK_F1,Event.CTRL_MASK);

    for (int i = 0; i < keys.length; i++) 
    {
      // Ignore ToolTipManager postTip action, 
      // in swing1.1beta3 and onward
      if (postTip.equals(keys[i])) 
        continue;
      char c = (char)keys[i].getKeyCode();
      int mod = keys[i].getModifiers();
      if ( mod == InputEvent.CTRL_MASK ) 
      {
        controlKeyStr = "Ctrl+"+(char)keys[i].getKeyCode();
        break;
      } 
      else if (mod == InputEvent.ALT_MASK) 
      {
        controlKeyStr = "Alt+"+(char)keys[i].getKeyCode();
        break;
      } 
    }
    return controlKeyStr;
  }

  /**
  *
  * Use this to getFontMetrics
  *
  */
  private class MyToolTip extends JToolTip
  {
    protected FontMetrics toolTipFontMetrics(Font currentFont)
    {
      return getFontMetrics(currentFont);
    }
  }

}

  
