/* FileViewer.java
 *
 * created: Thu Nov 19 1998
 *
 * This file is part of Artemis
 *
 * Copyright(C) 1998,1999,2000,2001  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/FileViewer.java,v 1.3 2005-01-11 15:54:32 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import uk.ac.sanger.artemis.*;

import java.awt.*;
import java.awt.event.*;
import java.io.Reader;
import java.io.BufferedReader;
import java.io.IOException;

import javax.swing.*;

/**
 *  This class implements a simple file viewer.  In fact any Reader object can
 *  be viewed.
 *
 *  @author Kim Rutherford
 *  @version $Id: FileViewer.java,v 1.3 2005-01-11 15:54:32 tjc Exp $
 *
 **/

public class FileViewer extends JFrame
{
  /** A JPanel to hold the close button. */
  private JPanel button_panel;

  /** The main component we use for displaying the file. */
  private JTextArea text_area = null;

  /**
   *  The size of the last FileViewer JFrame to be resized.  When a new
   *  FileViewer component is created it is given the current value of this
   *  variable.  This is updated when any FileViewer frame is resized.
   **/
  private static Dimension saved_size = null;

  /**
   *  The position of the last FileViewer JFrame to be moved.  When a new
   *  FileViewer component is created it is given the current value of this
   *  variable.  This is updated when any FileViewer frame is moved.
   **/
  private static Point saved_position = null;

  /**
   *  Create a new FileViewer component and make it visible.
   *  @param label The name to attach to the new JFrame.
   **/
  public FileViewer(final String label) 
  {
    this(label, true);
  }

  /**
   *  Create a new FileViewer component.
   *  @param label The name to attach to the new JFrame.
   *  @param visible The new FileViewer will be made visible if and only if
   *    this argument is true.
   **/
  public FileViewer(final String label, final boolean visible) 
  {
    super(label);

    getContentPane().setLayout(new BorderLayout());
    final Font font = Options.getOptions().getFont();
    setFont(font);

    text_area = new JTextArea();
    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
    
    JScrollPane scroller = new JScrollPane(text_area);
    scroller.setPreferredSize(new Dimension((int)screen.getWidth()/2,
                                            (int)screen.getHeight()/2));
    text_area.setEditable(false);
    text_area.setFont(font);
    text_area.setBackground(Color.white);

    getContentPane().add(scroller, "Center");

    button_panel = new JPanel();
    getContentPane().add(button_panel, "South");

    final JButton close_button = new JButton("Close");
    close_button.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent e) 
      {
        dispose();
      }
    });
    button_panel.add(close_button);

    addWindowListener(new WindowAdapter() 
    {
      public void windowClosing(WindowEvent event) 
      {
        dispose();
      }
    });

    addComponentListener(new ComponentAdapter() 
    {
      public void componentResized(ComponentEvent e) 
      {
        saved_size = FileViewer.this.getSize();
        saved_position = FileViewer.this.getLocation();
      }
      public void componentMoved(ComponentEvent e) 
      {
        saved_size = FileViewer.this.getSize();
        saved_position = FileViewer.this.getLocation();
      }
    });

    pack();

//  final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();

    if(saved_position == null) 
    {
      int screen_height = screen.height;
      int screen_width = screen.width;

//    if(screen_width <= 800 || screen_height <= 700) 
//      setSize(screen_width * 9 / 10, screen_height * 9 / 10);
//    else 
//      setSize(800, 700);

      setLocation(new Point((screen.width - getSize().width) / 2,
                            (screen.height - getSize().height) / 2));
    } 
    else
    {
      if(saved_position.x < 0 || saved_position.x + 20 > screen.width) 
        saved_position.x = 20;

      if(saved_position.y < 0 || saved_position.y + 20 > screen.height) 
        saved_position.y = 20;

      if(saved_size.width < 50) 
        saved_size.width = 50;

      if(saved_size.height < 50) 
        saved_size.height = 50;
      
      setLocation(saved_position);
      setSize(saved_size);
    }
    
    setVisible(visible);
  }

  /**
   *  Clear the viewer window.
   **/
  protected void clear()  
  {
    text_area.setText("");
  }

  /**
   *  Read from the given Reader and append the text to this FileViewer.
   *  @param read_stream The stream to read the contents of the viewer from.
   **/
  protected void appendFile(Reader read_stream)
      throws IOException 
  {
    final BufferedReader buffered_reader = new BufferedReader(read_stream);
    String line;

    while((line = buffered_reader.readLine()) != null) 
    {
      appendString(line + "\n");
      Thread.yield();
    }

    buffered_reader.close();
  }

  /**
   *  Clear the viewer window and display the lines from the given Reader.
   *  @param read_stream The stream to read the contents of the viewer from.
   **/
  public void readFile(Reader read_stream)
      throws IOException 
  {
    final BufferedReader buffered_reader = new BufferedReader(read_stream);

    String line;

    final StringBuffer line_buffer = new StringBuffer();

    while((line = buffered_reader.readLine()) != null) 
      line_buffer.append(line).append('\n');

    buffered_reader.close();

    final String new_text = line_buffer.toString();
    text_area.setText(new_text);
    text_area.setCaretPosition(0);
  }

  /**
   *  Clear the viewer window and display the lines from the given String.
   *  @param read_string The string to read the contents of the viewer from.
   **/
  protected void setText(String read_string) 
  {
    if(!read_string.equals(text_area.getText())) 
    {
      text_area.setText(read_string);
      text_area.setCaretPosition(0);
    }
  }


  /**
   *  Clear the viewer window and display the lines from the given String.
   *  @param read_string The string to read the contents of the viewer from.
   **/
  protected void appendString(String read_string) 
  {
    text_area.append(read_string);
    text_area.getCaret().setDot(0);
  }

  /**
   *  Return a String containing the text that this component is displaying.
   **/
  protected String getText() 
  {
    return getTextArea().getText();
  }

  /**
   *  Destroy this component.
   **/
  public void dispose() 
  {
    setVisible(false);
    
    saved_size     = getSize();
    saved_position = getLocation();

    super.dispose();
  }

  /**
   *  return the TextArea component from this FileViewer.
   **/
  protected JTextArea getTextArea() 
  {
    return text_area;
  }

  /**
   *  Return the JPanel containing the close button of this FileViewer.
   **/
  protected JPanel getButtonPanel() 
  {
    return button_panel;
  }

}
