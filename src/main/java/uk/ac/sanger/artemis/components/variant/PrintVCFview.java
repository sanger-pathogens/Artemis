/* PrintVCFview
 *
 * Copyright(C) 2010  Genome Research Limited
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

package uk.ac.sanger.artemis.components.variant;

import java.awt.*;
import java.awt.image.BufferedImage;
import java.awt.image.RenderedImage;
import java.awt.print.PageFormat;
import java.awt.print.Printable;
import java.awt.print.PrinterException;
import java.awt.print.PrinterJob;
import java.awt.event.*;
import javax.swing.*;
import java.io.*;

import uk.ac.sanger.artemis.editor.ScrollPanel;

/**
* Use to print images from VCFview
*/
public class PrintVCFview extends ScrollPanel implements Printable 
{
  private static final long serialVersionUID = 1L;
  private VCFview vcfView;

  public PrintVCFview(VCFview vcfView)
  {
    super();
    setBackground(Color.white);
    this.vcfView = vcfView;
  }

  /**
  *
  * Override paintComponent to draw entry
  *
  */
  public void paintComponent(Graphics g)
  {
// let UI delegate paint first (incl. background filling)
    super.paintComponent(g);
    vcfView.paintComponent(g);
  }


  /**
  *
  * Set the size of the image
  *
  */
  private void setImageSize()
  {
    setPreferredSize(vcfView.getPreferredSize());
  }

  /**
  *
  * Display a print preview page
  *
  */
  protected void printPreview()
  {
    final JFrame f = new JFrame("Print Preview");
    JPanel jpane   = (JPanel)f.getContentPane();

    JScrollPane scrollPane = new JScrollPane(this);

    jpane.setLayout(new BorderLayout());
    jpane.add(scrollPane,BorderLayout.CENTER);
    
    final Dimension dScreen = f.getToolkit().getScreenSize();
    Dimension d = new Dimension((int)(3*dScreen.getWidth()/4),
                                (int)(dScreen.getHeight()/2));
    f.setSize(d);
    setImageSize();

    JMenuBar menuBar = new JMenuBar();
    JMenu filemenu = new JMenu("File");
    menuBar.add(filemenu);

// print png/jpeg
    JMenuItem printImage = new JMenuItem("Save As Image Files (png/jpeg)...");
    printImage.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        print();
      }
    });
    filemenu.add(printImage);
    
//  print PostScript
    JMenuItem printPS = new JMenuItem("Print...");
    printPS.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        doPrintActions();
      }
    });
    filemenu.add(printPS);

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
                            "vcfView.png");
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
    formatSelect.setSelectedItem("png");

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
            "This option requires Java 1.9 or higher.");
    }
  }
  
  protected void doPrintActions()
  {
    final PrinterJob pj=PrinterJob.getPrinterJob();
    pj.setPrintable(PrintVCFview.this);
    pj.printDialog();
    try
    {
      pj.print();
    }
    catch (Exception PrintException) {}
  }
  
  /**
  *  Returns a generated image
  *  @param pageIndex   page number
  *  @return            image
  */
  private RenderedImage createImage()
  {
    setImageSize();
    // Create a buffered image in which to draw
    BufferedImage bufferedImage = new BufferedImage(
                                  vcfView.getWidth(), vcfView.getHeight(),
                                  BufferedImage.TYPE_INT_RGB);
    // Create a graphics contents on the buffered image
    Graphics2D g2d = bufferedImage.createGraphics();
    paintComponent(g2d);

    return bufferedImage;
  }


  /**
  * Write out the image
  * @param image        image
  * @param file         file to write image to
  * @param type         type of image
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
      System.out.println("An error occurred: " + e.getMessage());
      e.printStackTrace();
    }
  }

  /**
  *
  * The method @print@ must be implemented for @Printable@ interface.
  * Parameters are supplied by system.
  *
  */
  public int print(Graphics g, PageFormat pf, int pageIndex) throws PrinterException
  {
    setImageSize();
    Graphics2D g2 = (Graphics2D) g;

//  RepaintManager.currentManager(this).setDoubleBufferingEnabled(false);
    Dimension d = this.getSize();    //get size of document
    double panelWidth  = d.width;    //width in pixels
    double panelHeight = d.height;   //height in pixels
    
    if(panelWidth == 0)
    {
      d = this.getPreferredSize();
      panelWidth  = d.width;    
      panelHeight = d.height;  
    }
    double pageHeight = pf.getImageableHeight();   //height of printer page
    double pageWidth  = pf.getImageableWidth();    //width of printer page
    double scale = pageWidth/panelWidth;
    int totalNumPages = (int)Math.ceil(scale * panelHeight / pageHeight);
    // Make sure not print empty pages
    if(pageIndex >= totalNumPages)
     return Printable.NO_SUCH_PAGE;

    // Shift Graphic to line up with beginning of print-imageable region
    g2.translate(pf.getImageableX(), pf.getImageableY());
    // Shift Graphic to line up with beginning of next page to print
    g2.translate(0f, -pageIndex*pageHeight);
    // Scale the page so the width fits...
    g2.scale(scale, scale);
    paintComponent(g2);
    return Printable.PAGE_EXISTS;
  }
}
