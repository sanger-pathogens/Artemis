/* PrintACT.java
 *
 *
 * Copyright(C) 2004  Genome Research Limited
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
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.image.BufferedImage;
import java.awt.image.RenderedImage;
import java.awt.print.PageFormat;
import java.awt.print.Printable;
import java.awt.print.PrinterException;
import java.awt.print.PrinterJob;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import java.io.Writer;

import javax.swing.Box;
import javax.swing.JCheckBox;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JComboBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.KeyStroke;

import org.apache.batik.dom.GenericDOMImplementation;
import org.apache.batik.svggen.SVGGeneratorContext;
import org.apache.batik.svggen.SVGGraphics2D;
import org.apache.batik.svggen.SVGGraphics2DIOException;
import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Document;

import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.editor.ScrollPanel;

/**
* Use to print images from ACT
*/
public class PrintACT extends ScrollPanel implements Printable
{
  private static final long serialVersionUID = 1L;
  /** act display to create image from */
  private MultiComparator mc;
  private JCheckBox drawLabels = new JCheckBox("Show labels on alignment");

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
          if(c[j] instanceof BasePlot && c[j].isVisible())
          {
            ((BasePlot)c[j]).paintComponent(g2d);
            g2d.translate(0,((BasePlot)c[j]).getHeight());
          }
        }
        
        if(mc.getBamPanelArray()[i].isVisible())
        {
          mc.getBamPanelArray()[i].paintComponents(g2d);
          g2d.translate(0,mc.getBamPanelArray()[i].getHeight());
        }
        
        if(mc.getVcfPanelArray()[i].isVisible())
        {
          mc.getVcfPanelArray()[i].paintComponents(g2d);
          g2d.translate(0,mc.getVcfPanelArray()[i].getHeight());
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
          if(c[j] instanceof BasePlot && c[j].isVisible())
          {
            ((BasePlot)c[j]).paintComponent(g2d);
            g2d.translate(0,((BasePlot)c[j]).getHeight());
          }
        }
        
        if(mc.getBamPanelArray()[i].isVisible())
        {
          mc.getBamPanelArray()[i].paintComponents(g2d);
          g2d.translate(0,mc.getBamPanelArray()[i].getHeight());
        }
        
        if(mc.getVcfPanelArray()[i].isVisible())
        {
          mc.getVcfPanelArray()[i].paintComponents(g2d);
          g2d.translate(0,mc.getVcfPanelArray()[i].getHeight());
        }
      }

      if(i < mc.getAlignmentViewerArray().length)
      {
        mc.getAlignmentViewerArray()[i].paintComponentForPrint(g2d,drawLabels.isSelected());
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
          if(c[j] instanceof BasePlot && c[j].isVisible())
          {
            height += ((BasePlot)c[j]).getHeight();
            if(((BasePlot)c[j]).getSize().width < width &&
               ((BasePlot)c[j]).getSize().width  > 0)
              width = ((BasePlot)c[j]).getSize().width;
          }
      }
      
      if(mc.getBamPanelArray()[i].isVisible())
        height += mc.getBamPanelArray()[i].getHeight();

      if(mc.getVcfPanelArray()[i].isVisible())
        height += mc.getVcfPanelArray()[i].getHeight();
      
      height += mc.getFeatureDisplayArray()[i].getHeight();
      if(mc.getFeatureDisplayArray()[i].getWidth() < width)
        width = mc.getFeatureDisplayArray()[i].getWidth();

      if(i == mc.getEntryGroupArray().length - 1 &&
              mc.getEntryGroupArray().length == 2)
      {
        Component c[] = mc.getBasePlotGroupArray()[i].getComponents();
        for(int j=0; j<c.length; j++)
          if(c[j] instanceof BasePlot && c[j].isVisible())
          {
            height += ((BasePlot)c[j]).getHeight();
            if(((BasePlot)c[j]).getSize().width < width &&
               ((BasePlot)c[j]).getSize().width  > 0)
              width = ((BasePlot)c[j]).getSize().width;
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

    JMenu optionsmenu = new JMenu("Options");
    menuBar.add(optionsmenu);

// draw labels
    JCheckBoxMenuItem showLabels = new JCheckBoxMenuItem("Display Labels",
                                                         drawLabels.isSelected());
    showLabels.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        drawLabels.setSelected(!drawLabels.isSelected());
        repaint();
      }
    });
    optionsmenu.add(showLabels);

    f.setJMenuBar(menuBar);
    f.setVisible(true);
  }

  /**
  * Print to a jpeg or png file
  */
  public void print()
  {
    // file chooser
    final StickyFileChooser fc = new StickyFileChooser();
    File fselect = new File(fc.getCurrentDirectory()+
        System.getProperty("file.separator")+
        "act.png");
    fc.setSelectedFile(fselect);
     
    // file name prefix
    Box YBox = Box.createVerticalBox();
    JLabel labFormat = new JLabel("Select Format:");
    Font font = labFormat.getFont();
    labFormat.setFont(font.deriveFont(Font.BOLD));
    YBox.add(labFormat);

    Box bacross = Box.createHorizontalBox();
    final JComboBox formatSelect = new JComboBox(PrintArtemis.getImageFormats());
    formatSelect.setSelectedItem("png");
    formatSelect.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent arg0)
      {
        String selected;
        if(fc.getSelectedFile() != null)
        {
          selected = fc.getSelectedFile().getAbsolutePath();
          String fmts[] = PrintArtemis.getImageFormats();
          for(int i=0; i<fmts.length; i++)
            selected = selected.replaceAll("."+fmts[i]+"$", "");
        }
        else
          selected = "act";
        
        fc.setSelectedFile(new File(selected+"."+
               formatSelect.getSelectedItem()));
      }
    });
    formatSelect.setSelectedItem("png");

    Dimension d = formatSelect.getPreferredSize();
    formatSelect.setMaximumSize(d);
    bacross.add(Box.createHorizontalGlue());
    bacross.add(formatSelect);
    YBox.add(bacross);

    bacross = Box.createHorizontalBox();
    bacross.add(Box.createHorizontalGlue());
    bacross.add(drawLabels);
    YBox.add(bacross);
    
    // file prefix & format options
    fc.setAccessory(YBox);
    int n = fc.showSaveDialog(null);
    if(n == JFileChooser.CANCEL_OPTION)
      return;

    // remove file extension
    String fsave = fc.getSelectedFile().getAbsolutePath().toLowerCase();
    if(fsave.endsWith(".svg"))
    {
      createSVG(fc.getSelectedFile());
      return;
    }
    
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
            "This option requires Java 1.8 or higher.");
    }
  }

  private void createSVG(final File fout)
  {
    final DOMImplementation domImpl =
        GenericDOMImplementation.getDOMImplementation();
    final Document doc = domImpl.createDocument(
        "http://www.w3.org/2000/svg", "svg", null);

    SVGGeneratorContext ctx = SVGGeneratorContext.createDefault(doc);
    ctx.setComment("Generated by ACT with Batik SVG Generator");
    final SVGGraphics2D svgG = new SVGGraphics2D(ctx, true);
    svgG.setFont(Options.getOptions().getFont());
    svgG.setSVGCanvasSize( getImageSize() );
    paintComponent(svgG);

    try
    {
      final Writer out = new OutputStreamWriter(
          new FileOutputStream(fout), "UTF-8");
      svgG.stream(out, true);
    }
    catch (UnsupportedEncodingException e)
    {
      e.printStackTrace();
    }
    catch (SVGGraphics2DIOException e)
    {
      e.printStackTrace();
    }
    catch (FileNotFoundException e)
    {
      e.printStackTrace();
    }

    return;
  }
  
  /**
  *  Returns a generated image
  *  @param pageIndex   page number
  *  @return            image
  */
  private RenderedImage createImage()
  {
    Dimension d = getImageSize();

    // Create a buffered image in which to draw
    BufferedImage bufferedImage = new BufferedImage(
                                  d.width,d.height,
                                  BufferedImage.TYPE_INT_RGB);
    // Create a graphics contents on the buffered image
    Graphics2D g2d = bufferedImage.createGraphics();
    paintComponent(g2d);

    return bufferedImage;
  }

  /**
   * Get the size of the image
   * @return
   */
  private Dimension getImageSize()
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
              width = ((BasePlot)c[j]).getSize().width;
          }
      }

      if(mc.getBamPanelArray()[i].isVisible())
        height += mc.getBamPanelArray()[i].getHeight();
      
      if(mc.getVcfPanelArray()[i].isVisible())
        height += mc.getVcfPanelArray()[i].getHeight();
      
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
              width = ((BasePlot)c[j]).getSize().width;
          }
      }

      if(i < mc.getAlignmentViewerArray().length)
      {
        height += mc.getAlignmentViewerArray()[i].getHeight();
        if(mc.getAlignmentViewerArray()[i].getWidth() < width)
          width = mc.getAlignmentViewerArray()[i].getWidth();
      }
    }
    return new Dimension(width, height);
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
      System.out.println("Java 1.8+ is required");
      e.printStackTrace();
    }
  }
  
  protected void doPrintActions()
  {
    final PrinterJob pj=PrinterJob.getPrinterJob();
    pj.setPrintable(PrintACT.this);
    pj.printDialog();
    try
    {
      pj.print();
    }
    catch (Exception PrintException) {}
  }
  
  /**
  *
  * The method @print@ must be implemented for @Printable@ interface.
  * Parameters are supplied by system.
  *
  */
  public int print(Graphics g, PageFormat pf, int pageIndex) throws PrinterException
  {
    Graphics2D g2 = (Graphics2D) g;

//  RepaintManager.currentManager(this).setDoubleBufferingEnabled(false);
    Dimension d = getImageSize();    //get size of document
    double panelWidth  = d.width;    //width in pixels
    double panelHeight = d.height;   //height in pixels
    
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
