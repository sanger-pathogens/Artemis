/* PrintArtemis.java
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
import java.awt.image.BufferedImage;
import java.awt.image.RenderedImage;
import java.awt.print.PageFormat;
import java.awt.print.Printable;
import java.awt.print.PrinterException;
import java.awt.print.PrinterJob;
import java.awt.event.*;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.UnsupportedEncodingException;
import java.io.Writer;
import java.util.Arrays;
import java.util.HashSet;

import javax.swing.*;

import org.apache.batik.dom.GenericDOMImplementation;
import org.apache.batik.svggen.SVGGeneratorContext;
import org.apache.batik.svggen.SVGGraphics2D;
import org.apache.batik.svggen.SVGGraphics2DIOException;
import org.w3c.dom.DOMImplementation;
import org.w3c.dom.Document;

import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.editor.ScrollPanel;

/**
*
* Use to print images from Artemis
*
*/
public class PrintArtemis extends ScrollPanel implements Printable 
{
  private static final long serialVersionUID = 1L;

  /** entry to create image from */
  private EntryEdit entry;

  private JCheckBox selectDisplay   = new JCheckBox("Show Selection Header",true);
  private JCheckBox featDisplay     = new JCheckBox("Show Feature Display",true);
  private JCheckBox groupsDisplay   = new JCheckBox("Show Entries Loaded",true);
  private JCheckBox plotsDisplay    = new JCheckBox("Show Graphs",true);
  private JCheckBox jamDisplay      = new JCheckBox("Show Read Alignment",true);
  private JCheckBox vcfDisplay      = new JCheckBox("Show VCF",true);
  private JCheckBox onelineDisplay  = new JCheckBox("Show One Line Display",true);
  private JCheckBox baseDisplay     = new JCheckBox("Show Bases Display",true);
  private JCheckBox featListDisplay = new JCheckBox("Show Feature List",true);
  private int width;
  private int height;

  public PrintArtemis(EntryEdit entry)
  {
    super();
    this.entry = entry;
    setBackground(Color.white);
  }

  /**
  * Override paintComponent to draw entry
  */
  public void paintComponent(Graphics g2d)
  {
// let UI delegate paint first (incl. background filling)
    super.paintComponent(g2d);

    // feature list
    if(featListDisplay.isSelected())
    {
      FeatureList flist = entry.getFeatureList();
      Point ploc = flist.getViewport().getViewPosition(); 
      try
      {
        int translateX = 0;
        if(selectDisplay.isSelected())
          translateX += entry.getSelectionInfoDisplay().getHeight();
        if(groupsDisplay.isSelected())
          translateX += entry.getEntryGroupDisplay().getHeight();
        if(plotsDisplay.isSelected())
          translateX += entry.getBasePlotGroup().getHeight();
        if(jamDisplay.isSelected() && entry.getBamPanel() != null && entry.getBamPanel().isVisible())
          translateX += entry.getBamPanel().getHeight()-1;
        if(vcfDisplay.isSelected() && entry.getVcfView() != null && entry.getVcfView().isVisible())
          translateX += entry.getVcfPanel().getHeight();
        if(onelineDisplay.isSelected())
          translateX += entry.getOneLinePerEntryDisplay().getHeight();
        if(featDisplay.isSelected())
          translateX += entry.getFeatureDisplay().getHeight();
        if(baseDisplay.isSelected())
          translateX += entry.getBaseDisplay().getHeight();

        translateX-=2+ploc.y;
        g2d.translate(0,translateX);
        flist.paintComponent(g2d);
        g2d.translate(0,-translateX);
      }
      catch(IllegalArgumentException e){} // thrown if the list is not visible
    }

    // selection info
    if(selectDisplay.isSelected())
    {
      entry.getSelectionInfoDisplay().paintComponent(g2d);
      g2d.translate(0,entry.getSelectionInfoDisplay().getHeight());
    }

    // entry groups
    if(groupsDisplay.isSelected())
    {
      entry.getEntryGroupDisplay().printComponent(g2d);
      g2d.translate(0,entry.getEntryGroupDisplay().getHeight());
    }

    // plots
    if(plotsDisplay.isSelected())
      entry.getBasePlotGroup().printComponent(g2d);
//  g2d.translate(0,entry.getBasePlotGroup().getHeight());

    if(jamDisplay.isSelected() && entry.getBamPanel() != null && entry.getBamPanel().isVisible())
    {
      entry.getBamPanel().paintComponents(g2d);
      g2d.translate(0,entry.getBamPanel().getHeight()-1);
    }
    
    if(vcfDisplay.isSelected() && entry.getVcfView() != null && entry.getVcfView().isVisible())
    {
      entry.getVcfPanel().paintComponents(g2d);
      g2d.translate(0,entry.getVcfPanel().getHeight());
    }
    
    // one line per entry
    if(onelineDisplay.isSelected())
    {
      entry.getOneLinePerEntryDisplay().paintComponent(g2d);
      g2d.translate(0,entry.getOneLinePerEntryDisplay().getHeight());
    }

    // feature display
    if(featDisplay.isSelected())
    {
      FeatureDisplay fd = entry.getFeatureDisplay();
      fd.paintComponent(g2d);
      g2d.translate(0,entry.getFeatureDisplay().getHeight());
    }

    // base display
    if(baseDisplay.isSelected())
    {
      entry.getBaseDisplay().paintComponent(g2d);
      g2d.translate(0,entry.getBaseDisplay().getHeight());
    }
  }
  
  


  /**
  *
  * Set the size of the image
  *
  */
  private Dimension getImageSize()
  {
    height = 0;
    width  = entry.getFeatureDisplay().getDisplayWidth();
    if(selectDisplay.isSelected())
      height += entry.getSelectionInfoDisplay().getHeight();

    if(groupsDisplay.isSelected())
      height += entry.getEntryGroupDisplay().getHeight();

    if(jamDisplay.isSelected() && 
       entry.getJamView() != null && entry.getJamView().isVisible())
      height += entry.getBamPanel().getHeight();

    if(vcfDisplay.isSelected() && 
        entry.getVcfView() != null && entry.getVcfView().isVisible())
       height += entry.getVcfPanel().getHeight();
    
    if(plotsDisplay.isSelected())
      height += entry.getBasePlotGroup().getHeight();

    if(onelineDisplay.isSelected())
      height += entry.getOneLinePerEntryDisplay().getHeight();

    if(baseDisplay.isSelected())
      height += entry.getBaseDisplay().getHeight();

    if(featDisplay.isSelected())
      height += entry.getFeatureDisplay().getHeight();

    if(featListDisplay.isSelected())
      height += entry.getFeatureList().getViewport().getExtentSize().height;
    return new Dimension(width,height);
  }

  private void setImageSize()
  {
    setPreferredSize(getImageSize());
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

    JMenu optionsmenu = new JMenu("Options");
    menuBar.add(optionsmenu);

// draw selection info
    JCheckBoxMenuItem showSelection = new JCheckBoxMenuItem("Show Selection Header",
                                                            selectDisplay.isSelected());
    showSelection.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        selectDisplay.setSelected(!selectDisplay.isSelected());
        repaint();
      }
    });
    optionsmenu.add(showSelection);

// draw entry groups
    JCheckBoxMenuItem showEntryGroups = new JCheckBoxMenuItem("Show Entries Loaded",
                                                            groupsDisplay.isSelected());
    showEntryGroups.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        groupsDisplay.setSelected(!groupsDisplay.isSelected());
        repaint();
      }
    });
    optionsmenu.add(showEntryGroups);

// draw graphs
    JCheckBoxMenuItem showPlots = new JCheckBoxMenuItem("Show Graphs",
                                                         plotsDisplay.isSelected());

// only enable if graphs displayed
    if(entry.getBasePlotGroup().getNumberBasePlots() == 0)
      showPlots.setEnabled(false);

    showPlots.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        plotsDisplay.setSelected(!plotsDisplay.isSelected());
        repaint();
      }
    });
    optionsmenu.add(showPlots);
    
    
// draw read alignment viewer
    JCheckBoxMenuItem showJam = new JCheckBoxMenuItem("Show Read Alignment View",
                                                      jamDisplay.isSelected());
    
    if(entry.getJamView() == null || !entry.getJamView().isVisible())
      showJam.setEnabled(false);
    
    showJam.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        jamDisplay.setSelected(!jamDisplay.isSelected());
        repaint();
      }
    });
    optionsmenu.add(showJam);
    
 // draw vcf viewer
    JCheckBoxMenuItem showVcf = new JCheckBoxMenuItem("Show VCF View",
                                                      vcfDisplay.isSelected());
    
    if(entry.getVcfPanel() == null || !entry.getVcfPanel().isVisible())
      showJam.setEnabled(false);
    
    showVcf.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        vcfDisplay.setSelected(!vcfDisplay.isSelected());
        repaint();
      }
    });
    optionsmenu.add(showVcf);
    
// draw one line 
    JCheckBoxMenuItem showOneLine = new JCheckBoxMenuItem("Show One Line Display",
                                                          onelineDisplay.isSelected());
    if(!entry.getOneLinePerEntryDisplay().isVisible())
      showOneLine.setEnabled(false);
    
    showOneLine.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        onelineDisplay.setSelected(!onelineDisplay.isSelected());
        repaint();
      }
    });
    optionsmenu.add(showOneLine);

// draw features
    JCheckBoxMenuItem showFeatures = new JCheckBoxMenuItem("Show Feature Display",
                                                            featDisplay.isSelected());
    showFeatures.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        featDisplay.setSelected(!featDisplay.isSelected());
        repaint();
      }
    });
    optionsmenu.add(showFeatures);

// draw base display
    JCheckBoxMenuItem showBases = new JCheckBoxMenuItem("Show Bases Display",
                                                         baseDisplay.isSelected());
    showBases.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        baseDisplay.setSelected(!baseDisplay.isSelected());
        repaint();
      }
    });
    optionsmenu.add(showBases);

// draw feature list
    JCheckBoxMenuItem showFeatureList = new JCheckBoxMenuItem("Show Feature List",
                                                               featListDisplay.isSelected());
    showFeatureList.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        featListDisplay.setSelected(!featListDisplay.isSelected());
        repaint();
      }
    });
    optionsmenu.add(showFeatureList);

    f.setJMenuBar(menuBar);
    f.setVisible(true);
  }

  public static String[] getImageFormats()
  {
    final String fmts[] = javax.imageio.ImageIO.getWriterFormatNames();
    final HashSet<String> list = new HashSet<String>();
    for(int i=0; i<fmts.length; i++)
      list.add(fmts[i].toLowerCase());

    final String tmpFmts[] = new String[list.size()+1];
    System.arraycopy(list.toArray(), 0, tmpFmts, 0, list.size());
    tmpFmts[tmpFmts.length-1] = "svg";
    Arrays.sort(tmpFmts);

    return tmpFmts;
  }
  
  /**
  *
  * Print to a jpeg or png file
  *
  */
  public void print()
  {
    // file chooser
    final StickyFileChooser fc = new StickyFileChooser();
    File fselect = new File(fc.getCurrentDirectory()+
                            System.getProperty("file.separator")+
                            "artemis.png");
    fc.setSelectedFile(fselect);
     
    // file name prefix
    Box YBox = Box.createVerticalBox();
    JLabel labFormat = new JLabel("Select Format:");
    Font font = labFormat.getFont();
    labFormat.setFont(font.deriveFont(Font.BOLD));
    YBox.add(labFormat);

    Box bacross = Box.createHorizontalBox();
    final JComboBox formatSelect = new JComboBox(getImageFormats());
    formatSelect.setSelectedItem("png");
    formatSelect.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent arg0)
      {
        String selected;
        if(fc.getSelectedFile() != null)
        {
          selected = fc.getSelectedFile().getAbsolutePath();
          String fmts[] = getImageFormats();
          for(int i=0; i<fmts.length; i++)
            selected = selected.replaceAll("."+fmts[i]+"$", "");
        }
        else
          selected = "artemis";
        
        fc.setSelectedFile(new File(selected+"."+
               formatSelect.getSelectedItem()));
      }
    });

    Dimension d = formatSelect.getPreferredSize();
    formatSelect.setMaximumSize(d);
    bacross.add(Box.createHorizontalGlue());
    bacross.add(formatSelect);
    YBox.add(bacross);

    bacross = Box.createHorizontalBox();
    bacross.add(selectDisplay);
    bacross.add(Box.createHorizontalGlue());
    YBox.add(bacross);

    bacross = Box.createHorizontalBox();
    bacross.add(groupsDisplay);
    bacross.add(Box.createHorizontalGlue());
    YBox.add(bacross);

    if(entry.getBasePlotGroup().getNumberBasePlots() > 0)
    {
      bacross = Box.createHorizontalBox();
      bacross.add(plotsDisplay);
      bacross.add(Box.createHorizontalGlue());
      YBox.add(bacross);
    }

    if(entry.getBamPanel() != null && entry.getBamPanel().isVisible())
    {
      bacross = Box.createHorizontalBox();
      bacross.add(jamDisplay);
      bacross.add(Box.createHorizontalGlue());
      YBox.add(bacross);
    }
    
    if(entry.getVcfView() != null && entry.getVcfView().isVisible())
    {
      bacross = Box.createHorizontalBox();
      bacross.add(vcfDisplay);
      bacross.add(Box.createHorizontalGlue());
      YBox.add(bacross);
    }

    if(!entry.getOneLinePerEntryDisplay().isVisible())
    {
      bacross = Box.createHorizontalBox();
      bacross.add(onelineDisplay);
      bacross.add(Box.createHorizontalGlue());
      YBox.add(bacross);
    }

    bacross = Box.createHorizontalBox();
    bacross.add(featDisplay);
    bacross.add(Box.createHorizontalGlue());
    YBox.add(bacross);

    bacross = Box.createHorizontalBox();
    bacross.add(baseDisplay);
    bacross.add(Box.createHorizontalGlue());
    YBox.add(bacross);

    bacross = Box.createHorizontalBox();
    bacross.add(featListDisplay);
    bacross.add(Box.createHorizontalGlue());
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
    ctx.setComment("Generated by Artemis with Batik SVG Generator");
    final SVGGraphics2D svgG = new SVGGraphics2D(ctx, true);
    svgG.setFont(Options.getOptions().getFont());
    final FontMetrics fm = svgG.getFontMetrics();
    final Dimension d = getImageSize();
    svgG.setSVGCanvasSize( new Dimension(
        d.width+fm.stringWidth(" "), d.height+fm.getHeight()) );
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
  
  protected void doPrintActions()
  {
    final PrinterJob pj=PrinterJob.getPrinterJob();
    pj.setPrintable(PrintArtemis.this);
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
                                  width,height,
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
      System.out.println("Java 1.8+ is required");
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
    Graphics2D g2 = (Graphics2D)g.create();

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
