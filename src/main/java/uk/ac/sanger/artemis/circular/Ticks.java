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

import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.KeyStroke;


public class Ticks extends JPanel
{

  private TextFieldInt interval;
  private TextFieldInt minor;
  private TextFieldInt startTick;


  public Ticks(final DNADraw draw, boolean setTicks)
  {
    super();
   
    Box bdown = Box.createVerticalBox();
    Dimension d = new Dimension(100,25);
    bdown.add(Box.createVerticalStrut(4));
    Box bacross = Box.createHorizontalBox();

//tick interval
    bacross = Box.createHorizontalBox();
    interval = new TextFieldInt();
    if(draw != null)
      interval.setValue(draw.getTickInterval());
    else
      interval.setValue(200);
    interval.setMaximumSize(d);
    interval.setPreferredSize(d);
    bacross.add(interval);
    bacross.add(new JLabel(" Tick Interval"));
    bacross.add(Box.createHorizontalGlue());
    bdown.add(bacross);

//minor ticks
    bacross = Box.createHorizontalBox();
    minor = new TextFieldInt();
    if(draw != null)
      minor.setValue(draw.getMinorTickInterval());
    else
      minor.setValue(100);
    minor.setMaximumSize(d);
    minor.setPreferredSize(d);
    bacross.add(minor);
    bacross.add(new JLabel(" Minor Tick Interval"));
    bacross.add(Box.createHorizontalGlue());
    bdown.add(bacross);

//start ticks
    bacross = Box.createHorizontalBox();
    startTick = new TextFieldInt();
    if(draw != null)
      startTick.setValue(draw.getStartTick());
    else
      startTick.setValue(0);
    startTick.setMaximumSize(d);
    startTick.setPreferredSize(d);
    bacross.add(startTick);
    bacross.add(new JLabel(" Start Tick"));
    bacross.add(Box.createHorizontalGlue());
    bdown.add(bacross);
    
    if(setTicks)
    {
      bacross = Box.createHorizontalBox();
      JButton setTick = new JButton("Set");
      setTick.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          if(draw != null)
          {
            if(!draw.setStartTick(startTick.getValue()))
            {
              JOptionPane.showMessageDialog(null,
                    "Start tick is out of range for this sequence",
                    "Out of Range",JOptionPane.ERROR_MESSAGE);
            }
          
            if(!draw.setTickInterval(interval.getValue()))
            {
               JOptionPane.showMessageDialog(null, 
                    "Tick interval is out of range", 
                    "Out of Range",JOptionPane.ERROR_MESSAGE);
            }

            if(!draw.setMinorTickInterval(minor.getValue()))
            {
               JOptionPane.showMessageDialog(null,
                    "Minor Tick interval is out of range",
                    "Out of Range",JOptionPane.ERROR_MESSAGE);
            }

            draw.calculateTickPosistions();
            draw.repaint();
          }
        }
      });
      bacross.add(setTick);
      bdown.add(bacross);
    }

    add(bdown);
  }

  protected JMenuBar createMenuBar(final JFrame f)
  {
    final JMenuBar menuBar = new JMenuBar();

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


  /**
  *
  * Return the position tick marks start at
  *
  */
  protected int getStartTick()
  {
    return startTick.getValue();
  }


  /**
  *
  * Return the interval for the tick marks
  *
  */
  protected int getTickInterval()
  {
    return interval.getValue();
  }


  /**
  *
  * Return the interval for the minor tick marks
  *
  */
  protected int getMinorTickInterval()
  {
    return minor.getValue();
  }



}

