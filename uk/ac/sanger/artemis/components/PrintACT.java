/* PrintACT.java
 *
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

package uk.ac.sanger.artemis.components;

import java.awt.*;
import java.awt.print.Paper;
import java.awt.print.PageFormat;
import java.awt.print.PrinterJob;
import java.awt.image.BufferedImage;
import java.awt.image.RenderedImage;
import java.awt.event.*;
import javax.swing.*;
import java.util.*;
import java.io.*;
import javax.swing.border.*;
import javax.imageio.*;
import javax.imageio.stream.*;

import uk.ac.sanger.artemis.editor.ScrollPanel;

public class PrintACT extends ScrollPanel
{

  private MultiComparator mc;

  public PrintACT(MultiComparator mc)
  {
    this.mc = mc;
  }


  public void paintComponent(Graphics g)
  {
// let UI delegate paint first (incl. background filling)
    super.paintComponent(g);
    Graphics2D g2d = (Graphics2D)g.create();

    for(int i = 0; i < mc.getFeatureDisplayArray().length; ++i)
    {
      if(!(i == mc.getEntryGroupArray().length - 1 &&
                mc.getEntryGroupArray().length == 2))
      {
        Component c[] = mc.getBasePlotGroupArray()[i].getComponents();
        for(int j=0; j<c.length; j++)
        {
          if(c[j] instanceof BasePlot)
          {
            ((BasePlot)c[j]).paintCanvas(g2d);
            g2d.translate(0,((BasePlot)c[j]).getHeight());
          }
        }
      }


      mc.getFeatureDisplayArray()[i].paintComponent(g2d);
      g2d.translate(0,mc.getFeatureDisplayArray()[i].getHeight());

      if(i == mc.getEntryGroupArray().length - 1 &&
              mc.getEntryGroupArray().length == 2)
      {
        Component c[] = mc.getBasePlotGroupArray()[i].getComponents();
        for(int j=0; j<c.length; j++)
        {
          if(c[j] instanceof BasePlot)
          {
            ((BasePlot)c[j]).paintCanvas(g2d);
            g2d.translate(0,((BasePlot)c[j]).getHeight());
          }
        }
      }

      if(i < mc.getAlignmentViewerArray().length)
      {
        mc.getAlignmentViewerArray()[i].paintComponent(g2d);
        g2d.translate(0,mc.getAlignmentViewerArray()[i].getHeight());
      }
    }

  }

  /**
  *
  * Display a print preview page
  *
  */
  protected void printPreview()
  {
    int width  = 999999;
    int height = 0;
    for(int i = 0; i < mc.getFeatureDisplayArray().length; ++i)
    {
      if(!(i == mc.getEntryGroupArray().length - 1 &&
                mc.getEntryGroupArray().length == 2))
      {
        Component c[] = mc.getBasePlotGroupArray()[i].getComponents();
        for(int j=0; j<c.length; j++)
          if(c[j] instanceof BasePlot)
          {
            height += ((BasePlot)c[j]).getHeight();
            if(((BasePlot)c[j]).getSize().width < width &&
               ((BasePlot)c[j]).getSize().width  > 0)
              width = ((BasePlot)c[j]).getCanvas().getSize().width;
          }
      }

      height += mc.getFeatureDisplayArray()[i].getHeight();
      if(mc.getFeatureDisplayArray()[i].getWidth() < width)
        width = mc.getFeatureDisplayArray()[i].getWidth();

      if(i == mc.getEntryGroupArray().length - 1 &&
              mc.getEntryGroupArray().length == 2)
      {
        Component c[] = mc.getBasePlotGroupArray()[i].getComponents();
        for(int j=0; j<c.length; j++)
          if(c[j] instanceof BasePlot)
          {
            height += ((BasePlot)c[j]).getHeight();
            if(((BasePlot)c[j]).getSize().width < width &&
               ((BasePlot)c[j]).getSize().width  > 0)
              width = ((BasePlot)c[j]).getCanvas().getSize().width;
          }
      }

      if(i < mc.getAlignmentViewerArray().length)
      {
        height += mc.getAlignmentViewerArray()[i].getHeight();
        if(mc.getAlignmentViewerArray()[i].getWidth() < width)
          width = mc.getAlignmentViewerArray()[i].getWidth();
      }
    }

    final JFrame f = new JFrame("Print Preview");
    JPanel jpane = (JPanel)f.getContentPane();
    JScrollPane scrollPane = new JScrollPane(this);
    jpane.setLayout(new BorderLayout());
    jpane.add(scrollPane,BorderLayout.CENTER);
    
    final Dimension dScreen = f.getToolkit().getScreenSize();
    Dimension d = new Dimension((int)(dScreen.getWidth()*3/4),
                                (int)(dScreen.getHeight()/2));
    f.setSize(d);
    
    setPreferredSize(new Dimension(width,height));

    JMenuBar menuBar = new JMenuBar();
    JMenu filemenu = new JMenu("File");
    menuBar.add(filemenu);

// print png/jpeg
    JMenuItem printImage = new JMenuItem("Print Image Files (png/jpeg)...");
    printImage.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        print();
      }
    });
    filemenu.add(printImage);

// close
    filemenu.add(new JSeparator());
    JMenuItem menuClose = new JMenuItem("Close");
    menuClose.setAccelerator(KeyStroke.getKeyStroke(
              KeyEvent.VK_E, ActionEvent.CTRL_MASK));

    filemenu.add(menuClose);
    menuClose.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        f.dispose();
      }
    });

    f.setJMenuBar(menuBar);
    f.setVisible(true);
  }

  /**
  *
  * Print to a jpeg or png file
  *
  */
  public void print()
  {

    // file chooser
    String cwd = System.getProperty("user.dir");
    JFileChooser fc = new JFileChooser(cwd);
    File fselect = new File(cwd+
                            System.getProperty("file.separator")+
                            "act.jpeg");
    fc.setSelectedFile(fselect);
     
    // file name prefix
    Box YBox = Box.createVerticalBox();
    JLabel labFormat = new JLabel("Select Format:");
    Font font = labFormat.getFont();
    labFormat.setFont(font.deriveFont(Font.BOLD));
    YBox.add(labFormat);

    Box bacross = Box.createHorizontalBox();
    JComboBox formatSelect =
       new JComboBox(javax.imageio.ImageIO.getWriterFormatNames());
    Dimension d = formatSelect.getPreferredSize();
    formatSelect.setMaximumSize(d);
    bacross.add(Box.createHorizontalGlue());
    bacross.add(formatSelect);
    YBox.add(bacross);

    // file prefix & format options
    fc.setAccessory(YBox);
    int n = fc.showSaveDialog(null);
    if(n == JFileChooser.CANCEL_OPTION)
      return;

    // remove file extension
    String fsave = fc.getSelectedFile().getAbsolutePath().toLowerCase();
    if(fsave.endsWith(".png") ||
       fsave.endsWith(".jpg") ||
       fsave.endsWith(".jpeg") )
    {
      int ind = fsave.lastIndexOf(".");
      fsave = fc.getSelectedFile().getAbsolutePath();
      fsave = fsave.substring(0,ind);
    }
    else
      fsave = fc.getSelectedFile().getAbsolutePath();

    // image type
    String ftype = (String)formatSelect.getSelectedItem();
    try
    {
      RenderedImage rendImage = createImage();
      writeImageToFile(rendImage, new File(fsave+"."+ftype),
                       ftype);
    }
    catch(NoClassDefFoundError ex)
    {
      JOptionPane.showMessageDialog(this,
            "This option requires Java 1.4 or higher.");
    }
  }

  /**
  *
  *  Returns a generated image
  *  @param pageIndex   page number
  *  @return            image
  *
  */
  private RenderedImage createImage()
  {
    int width  = 999999;
    int height = 0;
    for(int i = 0; i < mc.getFeatureDisplayArray().length; ++i)
    {
      if(!(i == mc.getEntryGroupArray().length - 1 &&
                mc.getEntryGroupArray().length == 2))
      {
      }

      height += mc.getFeatureDisplayArray()[i].getHeight();
      if(mc.getFeatureDisplayArray()[i].getWidth() < width)
        width = mc.getFeatureDisplayArray()[i].getWidth();

      if(i == mc.getEntryGroupArray().length - 1 &&
              mc.getEntryGroupArray().length == 2)
      {
      }
 
      if(i < mc.getAlignmentViewerArray().length)
      {
        height += mc.getAlignmentViewerArray()[i].getHeight();
        if(mc.getAlignmentViewerArray()[i].getWidth() < width)
          width = mc.getAlignmentViewerArray()[i].getWidth();
      }
    }

    // Create a buffered image in which to draw
    BufferedImage bufferedImage = new BufferedImage(
                                  width,height,
                                  BufferedImage.TYPE_INT_RGB);
    // Create a graphics contents on the buffered image
    Graphics2D g2d = bufferedImage.createGraphics();
 
    for(int i = 0; i < mc.getFeatureDisplayArray().length; ++i)
    {
      if(!(i == mc.getEntryGroupArray().length - 1 &&
                mc.getEntryGroupArray().length == 2))
      {
//      getContentPane().add(base_plot_group_array[i]);
      }

 
      mc.getFeatureDisplayArray()[i].paintComponent(g2d);
      g2d.translate(0,mc.getFeatureDisplayArray()[i].getHeight());

      if(i == mc.getEntryGroupArray().length - 1 &&
              mc.getEntryGroupArray().length == 2)
      {
//      getContentPane().add(base_plot_group_array[i]);
      }

      if(i < mc.getAlignmentViewerArray().length)
      {
        mc.getAlignmentViewerArray()[i].paintComponent(g2d);
        g2d.translate(0,mc.getAlignmentViewerArray()[i].getHeight());
//      getContentPane().add(getAlignmentViewerArray()[i]);
      }
    }
 
    return bufferedImage;
  }


  /**
  *
  * Write out the image
  * @param image        image
  * @param file         file to write image to
  * @param type         type of image
  *
  */
  private void writeImageToFile(RenderedImage image,
                               File file, String type)
  {
    try
    {
      javax.imageio.ImageIO.write(image,type,file);
    }
    catch ( IOException e )
    {
      System.out.println("Java 1.4+ is required");
      e.printStackTrace();
    }
  }

}
