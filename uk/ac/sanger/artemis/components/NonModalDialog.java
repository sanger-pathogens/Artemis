/* NonModalDialog
 *
 * created: 2011
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2011  Genome Research Limited
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

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.FontMetrics;
import java.awt.GraphicsConfiguration;
import java.awt.GraphicsDevice;
import java.awt.GraphicsEnvironment;
import java.awt.Insets;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;

public class NonModalDialog extends JDialog
{
  private static final long serialVersionUID = 1L;

  /**
   * Create a non-modal dialog with an array of strings to
   * display on separate lines.
   * @param f
   * @param labelStr
   */
  public NonModalDialog(JFrame f, String labelStr[])
  {
    super(f, "Check the reference is correctly selected", false);
    
    GraphicsEnvironment ge = GraphicsEnvironment.getLocalGraphicsEnvironment();
    GraphicsDevice gd = ge.getDefaultScreenDevice();
    GraphicsConfiguration gc = gd.getDefaultConfiguration();
    Insets ins = Toolkit.getDefaultToolkit().getScreenInsets(gc);
    int sw = gc.getBounds().width - ins.left - ins.right;
    int sh = gc.getBounds().height - ins.top - ins.bottom;

    setDefaultCloseOperation(JDialog.DISPOSE_ON_CLOSE);
    Container cp = getContentPane();
    Box yBox = Box.createVerticalBox();
    cp.setLayout(new BorderLayout());
        
    int width = 10;
    for(int i=0; i<labelStr.length; i++)
    {
      JLabel label = new JLabel(labelStr[i]);
      yBox.add(label);
      
      FontMetrics fm = label.getFontMetrics(label.getFont());
      int thisWidth = fm.stringWidth(label.getText());
      if(thisWidth > width)
        width = thisWidth;
    }
    cp.add(yBox, BorderLayout.NORTH);
    
    JButton okButton = new JButton("OK");
    cp.add(okButton, BorderLayout.SOUTH);
    okButton.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent arg0)
      {
        dispose();
      }
    });

    setBounds((sw - width)/2, sh/2, width, 100);
    setVisible(true);
  }
}