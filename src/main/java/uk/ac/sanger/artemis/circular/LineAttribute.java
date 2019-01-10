/*
 * Copyright (C) 2008  Genome Research Limited
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
 *  @author: Tim Carver
 */


package uk.ac.sanger.artemis.circular;

import javax.swing.Box;
import javax.swing.ButtonGroup;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JSlider;
import javax.swing.KeyStroke;
import javax.swing.event.*;
import java.awt.event.*;
import java.awt.Dimension;
import java.util.Hashtable;

public class LineAttribute extends JPanel
{
  private static final long serialVersionUID = 1L;
  private Hashtable lineAttr = new Hashtable();
  private TextFieldInt start;
  private TextFieldInt end;
  private TextFieldInt lineSize;

  public LineAttribute(final DNADraw draw)
  {
    super();
 
    Dimension d = new Dimension(100,25);
    Box bdown = Box.createVerticalBox();
    bdown.add(Box.createVerticalStrut(4));
    Box bacross = Box.createHorizontalBox();

    // circular or linear dna
    ButtonGroup group = new ButtonGroup();
    final JRadioButton jcirc = new JRadioButton("Circular");
    jcirc.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        lineAttr.put("circular",new Boolean(jcirc.isSelected()));
        if(draw != null)
          draw.repaint();
      }
    });
    group.add(jcirc);
    JRadioButton jline = new JRadioButton("Linear");
    jline.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        lineAttr.put("circular",new Boolean(jcirc.isSelected()));
        if(draw != null)
          draw.repaint();
      }
    });
    group.add(jline);
    jcirc.setSelected( draw.isCircular() );
    jline.setSelected( !draw.isCircular() );
    lineAttr.put("circular",new Boolean(draw.isCircular()));
    bacross.add(jcirc); 
    bacross.add(jline);
    bacross.add(Box.createHorizontalGlue());
    bdown.add(bacross);

    // start position
    bdown.add(Box.createVerticalStrut(4));
    bacross = Box.createHorizontalBox();
    start = new TextFieldInt();
    start.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        int s = getStart();
        lineAttr.put("start",new Integer(s));
        draw.setStart(s);
        if(draw != null)
          draw.repaint();
      }
    });

    if(draw != null)
      start.setValue(draw.getStart());
    else
      start.setValue(1);

    start.setPreferredSize(d);
    start.setMaximumSize(d);
    bacross.add(start);
    bacross.add(new JLabel(" start"));
    bacross.add(Box.createHorizontalGlue());
    bdown.add(bacross);
 
    // end position
    bacross = Box.createHorizontalBox();
    end = new TextFieldInt();
    end.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent ev)
      {
        int e = getEnd();
        lineAttr.put("end",new Integer(e));
        draw.setEnd(e);
        if(draw != null)
          draw.repaint();
      }
    });
    end.setPreferredSize(d);
    end.setMaximumSize(d);

    if(draw != null)
      end.setValue(draw.getEnd());
    else
      end.setValue(1000);

    bacross.add(end);
    bacross.add(new JLabel(" stop"));
    bacross.add(Box.createHorizontalGlue());
    bdown.add(bacross);

    bacross = Box.createHorizontalBox();
    lineSize = new TextFieldInt();

    // line size
    int lsize = 5;
    if(draw != null)
      lsize = draw.getLineSize();
    lineAttr.put("lsize",new Integer(lsize));

    final JSlider slider = new JSlider(1,20,lsize);
    lineSize.setPreferredSize(d);
    lineSize.setMaximumSize(d);
    lineSize.setValue(lsize);
    bacross.add(lineSize);
    bacross.add(new JLabel(" line width"));
    bacross.add(Box.createHorizontalGlue());
    bdown.add(bacross);
    // change line size on carriage return 
    lineSize.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        int size = lineSize.getValue();
        slider.setValue(size);
        lineAttr.put("lsize",new Integer(size));
        if(draw != null)
        {
          draw.setLineSize(size);
          draw.repaint();
        }
      }
    });

    bacross = Box.createHorizontalBox();
    slider.addChangeListener(new ChangeListener()
    {
      public void stateChanged(ChangeEvent e)
      {
        int size = slider.getValue();
        lineSize.setValue(size);
        lineAttr.put("lsize",new Integer(size));
        if(draw != null)
        {
          draw.setLineSize(size);
          draw.repaint();
        }
      }
    });

    bacross.add(slider);
    bacross.add(Box.createHorizontalGlue());
    bdown.add(bacross);
   
    add(bdown);

  }

  protected JMenuBar createMenuBar(final JFrame f)
  {
    JMenuBar menuBar = new JMenuBar();

    JMenu fileMenu = new JMenu("File");
    fileMenu.setMnemonic(KeyEvent.VK_F);
    menuBar.add(fileMenu);

    JMenuItem closeMenu = new JMenuItem("Close");
    closeMenu.setAccelerator(KeyStroke.getKeyStroke(
              KeyEvent.VK_E, ActionEvent.CTRL_MASK));

    closeMenu.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        f.dispose();
      }
    });
    fileMenu.add(closeMenu);
    return menuBar;
  }

  protected Hashtable getLineAttr()
  {
    lineAttr.put("start",new Integer(start.getValue()));
    lineAttr.put("end",new Integer(end.getValue()));
    lineAttr.put("lsize",new Integer(lineSize.getValue()));
    
    return lineAttr;
  }

  protected int getStart()
  {
    return start.getValue();
  }


  protected int getEnd()
  {
    return end.getValue();
  }

}

