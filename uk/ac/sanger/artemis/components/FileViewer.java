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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/FileViewer.java,v 1.19 2008-11-12 16:50:37 tjc Exp $
 */

package uk.ac.sanger.artemis.components;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.Point;
import java.awt.Toolkit;
import java.awt.event.*;
import java.io.File;
import java.io.FileWriter;
import java.io.Reader;
import java.io.BufferedReader;
import java.io.IOException;
import java.util.Hashtable;

import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextPane;
import javax.swing.text.BadLocationException;
import javax.swing.text.Document;
import javax.swing.text.Element;
import javax.swing.text.MutableAttributeSet;
import javax.swing.text.SimpleAttributeSet;
import javax.swing.text.StyleConstants;

import org.apache.log4j.Level;

import uk.ac.sanger.artemis.Options;


/**
 *  This class implements a simple file viewer.  In fact any Reader object can
 *  be viewed.
 *
 *  @author Kim Rutherford
 *  @version $Id: FileViewer.java,v 1.19 2008-11-12 16:50:37 tjc Exp $
 *
 **/

public class FileViewer extends JFrame
{
  /** */
  private static final long serialVersionUID = 1L;

  /** A JPanel to hold the close button. */
  protected JPanel button_panel;

  /** The main component we use for displaying the file. */
  private JTextPane textPane = null;

  private Hashtable<Level, MutableAttributeSet> fontAttributes;

  private final static Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
  
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
  
  private boolean isHideOnClose = false;
  
  /**
   *  Create a new FileViewer component and make it visible.
   *  @param label The name to attach to the new JFrame.
   **/
  public FileViewer(final String label) 
  {
    this(label, true);
  }

  public FileViewer(final String label, final boolean visible)
  {
    this(label, visible, true, false);  
  }
  
  public FileViewer(final String label, final boolean visible, 
                    final boolean showClearButton,
                    final boolean showSaveButton) 
  {
    this(label, visible, showClearButton, showSaveButton, null); 
  }
  
  /**
   *  Create a new FileViewer component.
   *  @param label The name to attach to the new JFrame.
   *  @param visible The new FileViewer will be made visible if and only if
   *    this argument is true.
   **/
  public FileViewer(final String label, final boolean visible, 
                    final boolean showClearButton,
                    final boolean showSaveButton,
                    final ActionListener saveListener) 
  {
    super(label);

    final JMenuBar mBar = new JMenuBar();
    setJMenuBar(mBar);
    final JMenu fileMenu = new JMenu("File");
    final JMenuItem close = new JMenuItem("Close");
    close.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent arg0)
      {
        if(isHideOnClose())
          setVisible(false);
        else
          dispose();
      }
    });
    fileMenu.add(close);
    mBar.add(fileMenu);
    
    getContentPane().setLayout(new BorderLayout());
    final Font font = Options.getOptions().getFont();
    setFont(font);

    // ensure wrapping is turned off
    textPane = new JTextPane()
    {
      /** */
      private static final long serialVersionUID = 1L;

      public boolean getScrollableTracksViewportWidth()
      {
        return false;
      }
    };
    
    final JScrollPane scroller = new JScrollPane(textPane);
    scroller.getViewport().setBackground(Color.white);

    textPane.setEditable(false);
    textPane.setFont(font);
    textPane.setBackground(Color.white);

    getContentPane().add(scroller, "Center");

    button_panel = new JPanel(new FlowLayout());
    getContentPane().add(button_panel, "South");

    final JButton close_button = new JButton("Close");
    close_button.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent e) 
      {
        if(isHideOnClose())
          setVisible(false);
        else
          dispose();
      }
    });
    button_panel.add(close_button);
    
    if(showClearButton)
    {
      final JButton clearbutton = new JButton("Clear");
      clearbutton.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          clear();
        }
      });
      button_panel.add(clearbutton);
    }
    
    if(showSaveButton)
    {
      final JButton saveToFile = new JButton("Save");
      if(saveListener != null)
        saveToFile.addActionListener(saveListener);
      else
        saveToFile.addActionListener(new ActionListener()
        {
          public void actionPerformed(ActionEvent e)
          {
            writeToFile(getText());
           }
        });
      button_panel.add(saveToFile);
    }
    
    addWindowListener(new WindowAdapter() 
    {
      public void windowClosing(WindowEvent event) 
      {
        if(isHideOnClose())
          setVisible(false);
        else
          dispose();
      }
    });

    if(saved_position == null) 
    {
      Dimension d = new Dimension((int)screen.getWidth()/2,
          (int)screen.getHeight()/2);
      scroller.setPreferredSize(d);
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
      
      if(this instanceof ValidateViewer)
        scroller.setPreferredSize(new Dimension((int)screen.getWidth()/3,
            (int)screen.getHeight()/3));
      else
        scroller.setPreferredSize(saved_size);
    }
    
    pack();
    if(saved_position == null) 
      Utilities.centreFrame(this);
    
    if(visible)
      setVisible(true);
    createDefaultFontAttributes();
    
    addComponentListener(new ComponentAdapter() 
    {
      public void componentResized(ComponentEvent e) 
      {
        saved_size = scroller.getSize();
        saved_position = FileViewer.this.getLocation();
      }
      public void componentMoved(ComponentEvent e) 
      {
        saved_size = scroller.getSize();
        saved_position = FileViewer.this.getLocation();
      }
    });
  }

  
  /**
   *  Clear the viewer window.
   **/
  protected void clear()  
  {
    textPane.setText("");
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
  protected void readFile(Reader read_stream)
      throws IOException 
  {
    final BufferedReader buffered_reader = new BufferedReader(read_stream);

    String line;
    final StringBuffer line_buffer = new StringBuffer();
    while((line = buffered_reader.readLine()) != null) 
      line_buffer.append(line).append('\n');

    textPane.setText(line_buffer.toString());
    textPane.setCaretPosition(0);
    buffered_reader.close();
  }

  /**
   *  Clear the viewer window and display the lines from the given String.
   *  @param read_string The string to read the contents of the viewer from.
   **/
  protected void setText(String read_string) 
  {
    textPane.setText(read_string);
  }


  /**
   *  Clear the viewer window and display the lines from the given String.
   *  @param read_string The string to read the contents of the viewer from.
   **/
  protected void appendString(String read_string) 
  { 
    final Level level;
    if(read_string.indexOf("FATAL") > -1)
      level = Level.FATAL;
    else if(read_string.indexOf("ERROR") > -1)
      level = Level.ERROR;
    else if(read_string.indexOf("WARN") > -1)
      level = Level.WARN;
    else if(read_string.indexOf("INFO") > -1)
      level = Level.INFO;
    else
      level = Level.DEBUG;
    appendString(read_string, level);
  }
  
  public void appendString(final String read_string, final Level level) 
  {
    final Document doc = textPane.getStyledDocument();
    try
    {
      doc.insertString(doc.getLength(), read_string, 
          (MutableAttributeSet)fontAttributes.get(level));
      textPane.setCaretPosition(doc.getLength());
    }
    catch(BadLocationException e)
    {
      e.printStackTrace();
    }
  }

  private void createDefaultFontAttributes() 
  {
    Level[] prio = { Level.FATAL, Level.ERROR, 
           Level.WARN, Level.INFO, Level.DEBUG };

    fontAttributes = new Hashtable<Level, MutableAttributeSet>();
    for (int i=0; i<prio.length;i++) 
    {
      MutableAttributeSet att = new SimpleAttributeSet();
      fontAttributes.put(prio[i], att);
    }

    setTextColor(Level.FATAL, Color.red);
    setTextColor(Level.ERROR, Color.magenta);
    setTextColor(Level.WARN, Color.orange);
    setTextColor(Level.INFO, Color.black);
    setTextColor(Level.DEBUG, Color.blue);
  }
  
  void setTextColor(Level p, Color c) 
  {
    StyleConstants.setForeground(
          (MutableAttributeSet)fontAttributes.get(p),c);
  }
  
  /**
   *  Return a String containing the text that this component is displaying.
   **/
  protected String getText() 
  {
    return textPane.getText();
  }


  /**
   *  Return the JPanel containing the close button of this FileViewer.
   **/
  protected JPanel getButtonPanel() 
  {
    return button_panel;
  }

  /**
   * 
   * @return
   */
  protected int getLineCount()
  {
    int caretPosition = textPane.getCaretPosition();
    Element root = textPane.getDocument().getDefaultRootElement();
 
    return root.getElementIndex( caretPosition ) + 1;
    
    /*try
    {
      int offset     = textPane.getDocument().getLength();
      Rectangle r    = textPane.modelToView( offset );
      
      return (int) r.y/lineHeight;
    }
    catch(Exception e)
    {
      return 0;
    }*/

  }

  public JTextPane getTextPane()
  {
    return textPane;
  }

  public boolean isHideOnClose()
  {
    return isHideOnClose;
  }

  public void setHideOnClose(boolean isHideOnClose)
  {
    this.isHideOnClose = isHideOnClose;
  }
  
  public static void writeToFile(final String txt)
  {
    StickyFileChooser fc = new StickyFileChooser();
    fc.showSaveDialog(null);
    
    File f = fc.getSelectedFile();
    if(f.exists() && f.canWrite())
    {
      int status = JOptionPane.showConfirmDialog(null, 
          f.getName()+" exists overwrite?", 
          "Overwrite", JOptionPane.OK_CANCEL_OPTION);
      if(status == JOptionPane.CANCEL_OPTION)
        return;
    }
    try
    {
      FileWriter writer = new FileWriter(f);
      writer.write(txt);
      writer.close();
    }
    catch (IOException e1)
    {
      JOptionPane.showMessageDialog(null,
          e1.getMessage(), "Problem Writing", JOptionPane.WARNING_MESSAGE);
    }
  }
}
