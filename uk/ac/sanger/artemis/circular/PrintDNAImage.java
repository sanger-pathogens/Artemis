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

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.print.PageFormat;
import java.awt.print.Paper;
import java.awt.print.PrinterJob;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.image.BufferedImage;
import java.awt.image.RenderedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import java.io.Writer;

import javax.swing.BorderFactory;
import javax.swing.Box;
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
import javax.swing.JTextField;
import javax.swing.KeyStroke;
import javax.swing.border.Border;

import org.apache.batik.dom.GenericDOMImplementation;
import org.apache.batik.svggen.SVGGeneratorContext;
import org.apache.batik.svggen.SVGGraphics2D;
import org.apache.batik.svggen.SVGGraphics2DIOException;
import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Document;

import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.components.PrintArtemis;
import uk.ac.sanger.artemis.components.StickyFileChooser;




/**
*
* Print png/jpeg image and print preview.
* Java 1.4 or higher is required for the imageio package
* which is used here to create jpeg and png images of the
* DNA diagram.
*
*/
public class PrintDNAImage extends ScrollPanel
{
  private static final long serialVersionUID = 1L;

  /** page format */
  private PageFormat format = null;

  /** alignment sequence panel */
  private DNADraw dna;
  /** status field for print preview */
  private JTextField statusField = new JTextField("");
  /** type (jpeg/png) */
  private String type;

  /**
  *
  * @param dna   dna panel
  *
  */
  public PrintDNAImage(DNADraw dna)
  {
    super();
    this.dna = dna;

    setBackground(Color.white);
  }


  /**
  *
  * Override this method to draw the sequences
  * @return Graphics g
  *
  */
  public void paintComponent(Graphics g)
  {
// let UI delegate paint first (incl. background filling)
    super.paintComponent(g);
    Graphics2D g2d = (Graphics2D) g.create();
    dna.drawAll(g2d,true);
  }

  /**
  *
  * Get a default page format
  * @return     page format
  *
  */
  protected PageFormat getFormatDialog()
  {
    PrinterJob printerJob = PrinterJob.getPrinterJob();
    format = new PageFormat();
    format = printerJob.pageDialog(format);
    return format;
  }


  /**
  *
  *  Returns a generated image
  *  @param pageIndex   page number
  *  @return            image
  *
  */
  private RenderedImage createDNAImage(int pageIndex)
  {
    int width  = (int)format.getWidth();
    int height = (int)format.getHeight();
    // Create a buffered image in which to draw
    BufferedImage bufferedImage = new BufferedImage(
                                  width,height,
                                  BufferedImage.TYPE_INT_RGB);

    // Create a graphics contents on the buffered image
    Graphics2D g2d = bufferedImage.createGraphics();
    g2d.setColor(Color.white);
    g2d.fillRect(0,0,width,height);
    // Draw graphics
    dna.drawAll(g2d,true);

    return bufferedImage;
  }


  /**
  *
  * Display a print preview page
  *
  */
  protected void printPreview()
  {
    Border loweredbevel = BorderFactory.createLoweredBevelBorder();
    Border raisedbevel = BorderFactory.createRaisedBevelBorder();
    Border compound = BorderFactory.createCompoundBorder(raisedbevel,loweredbevel);
    statusField.setBorder(compound);
    statusField.setEditable(false);

    if(format == null)
      format = getFormatDialog();

    statusField.setText("DNA map");
    final JFrame f = new JFrame("Print Preview");
    JPanel jpane = (JPanel)f.getContentPane();
    JScrollPane scrollPane = new JScrollPane(this);
    jpane.setLayout(new BorderLayout());
    jpane.add(scrollPane,BorderLayout.CENTER);
    jpane.add(statusField,BorderLayout.SOUTH);

    Dimension d = new Dimension((int)format.getWidth(),
                                (int)format.getHeight());
    setPreferredSize(d);
    f.setSize(d);

    JMenuBar menuBar = new JMenuBar();
    JMenu filemenu = new JMenu("File");
    menuBar.add(filemenu);

// print postscript
    JMenu printMenu = new JMenu("Print");
    filemenu.add(printMenu);

    JMenuItem print = new JMenuItem("Print Postscript...");
    print.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        dna.doPrintActions();
      }
    });
    printMenu.add(print);

// print png/jpeg
    JMenuItem printImage = new JMenuItem("Print png/jpeg Image...");
    printImage.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        printAsSinglePage();
      }
    });
    printMenu.add(printImage);

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
  * Provide some options for the image created
  */
  private File showOptions()
  {
    final StickyFileChooser fc = new StickyFileChooser();
    File fselect = new File(fc.getCurrentDirectory()+
        System.getProperty("file.separator")+
        "dnaplotter.png");
    fc.setSelectedFile(fselect);

// file name prefix
    Box bdown = Box.createVerticalBox();
    bdown.add(Box.createVerticalGlue());

    JLabel labFormat = new JLabel("Select Format:");
    Font font = labFormat.getFont();
    labFormat.setFont(font.deriveFont(Font.BOLD));

    bdown.add(labFormat);

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
          selected = "dnaplotter";
        
        fc.setSelectedFile(new File(selected+"."+
               formatSelect.getSelectedItem()));
      }
    });

    Dimension d = formatSelect.getPreferredSize();
    formatSelect.setMaximumSize(d);
    bacross.add(Box.createHorizontalGlue());
    bacross.add(formatSelect);
    bdown.add(bacross);

// file prefix & format options
    fc.setAccessory(bdown);
    int n = fc.showSaveDialog(null);
    if(n == JFileChooser.CANCEL_OPTION)
      return null;

    type = (String)formatSelect.getSelectedItem();
    return fc.getSelectedFile();
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
    catch (IOException e)
    {
      System.out.println("Java 1.4+ is required");
    }
  }

  /**
  * Print to one jpeg or png file
  */
  public void printAsSinglePage()
  {
    //PrinterJob printerJob = PrinterJob.getPrinterJob();
    format = new PageFormat();

    File file = showOptions();
    if(file == null)
      return;
    
    Dimension d = dna.getSize();
    double imageWidth  = d.getWidth();
    double imageHeight = d.getHeight();

    if(type.equals("svg"))
    {
      createSVG(file, d);
      return;
    }
    
    Paper paper  = format.getPaper();
    paper.setSize(imageWidth,imageHeight);
    paper.setImageableArea(0,0,
                           imageWidth,imageHeight+imageHeight);
    format.setPaper(paper);

    try
    {
      RenderedImage rendImage = createDNAImage(0);
      writeImageToFile(rendImage,file,type);
    }
    catch(NoClassDefFoundError ex)
    {
      JOptionPane.showMessageDialog(this,
            "This option requires Java 1.4 or higher.");
    }
  }
  
  private void createSVG(final File fout, Dimension d)
  {
    final DOMImplementation domImpl =
        GenericDOMImplementation.getDOMImplementation();
    final Document doc = domImpl.createDocument(
        "http://www.w3.org/2000/svg", "svg", null);

    SVGGeneratorContext ctx = SVGGeneratorContext.createDefault(doc);
    ctx.setComment("Generated by DNAPlotter with Batik SVG Generator");
    final SVGGraphics2D svgG = new SVGGraphics2D(ctx, true);
    svgG.setFont(Options.getOptions().getFont());
    svgG.setSVGCanvasSize(d);
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

}


