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
import java.awt.event.*;
import javax.swing.*;
import java.io.*;
import javax.imageio.*;
import javax.imageio.stream.*;

import uk.ac.sanger.artemis.editor.ScrollPanel;

/**
*
* Use to print images from Artemis
*
*/
public class PrintArtemis extends ScrollPanel
{

  /** entry to create image from */
  private EntryEdit entry;

  private JCheckBox selectDisplay   = new JCheckBox("Show Selection Header",true);
  private JCheckBox featDisplay     = new JCheckBox("Show Feature Display",true);
  private JCheckBox groupsDisplay   = new JCheckBox("Show Entries Loaded",true);
  private JCheckBox plotsDisplay    = new JCheckBox("Show Graphs",true);
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
  *
  * Override paintComponent to draw entry
  *
  */
  public void paintComponent(Graphics g)
  {
// let UI delegate paint first (incl. background filling)
    super.paintComponent(g);
    Graphics2D g2d = (Graphics2D)g.create();

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

    // one line per entry
    if(onelineDisplay.isSelected())
    {
      entry.getOneLinePerEntryDisplay().paintComponent(g2d);
      g2d.translate(0,entry.getOneLinePerEntryDisplay().getHeight());
    }

    // feature display
    if(featDisplay.isSelected())
    {
      entry.getFeatureDisplay().paintComponent(g2d);
      g2d.translate(0,entry.getFeatureDisplay().getHeight());
    }

    // base display
    if(baseDisplay.isSelected())
    {
      entry.getBaseDisplay().paintComponent(g2d);
      g2d.translate(0,entry.getBaseDisplay().getHeight());
    }

    // feature list
    if(featListDisplay.isSelected())
    {
      FeatureList flist = entry.getFeatureList();
      Point ploc = flist.getViewport().getViewPosition();
//    flist.setOpaque(false);
      g2d.translate(0,-ploc.y);
      flist.paintComponent(g2d);
//    flist.setOpaque(true);
    }
  }


  /**
  *
  * Set the size of the image
  *
  */
  private void setImageSize()
  {
    height = 0;
    width  = entry.getFeatureDisplay().getDisplayWidth();
    if(selectDisplay.isSelected())
      height += entry.getSelectionInfoDisplay().getHeight();

    if(groupsDisplay.isSelected())
      height += entry.getEntryGroupDisplay().getHeight();

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
    setPreferredSize(new Dimension(width,height));
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
                            "artemis.png");
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
