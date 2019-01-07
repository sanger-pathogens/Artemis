/* VCFview
 *
 * created: July 2010
 *
 * This file is part of Artemis
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
import java.awt.AlphaComposite;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Composite;
import java.awt.Container;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.FontMetrics;
import java.awt.Frame;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.Point;
import java.awt.Rectangle;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionListener;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Vector;
import java.util.regex.Pattern;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.ButtonGroup;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JRadioButton;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.JScrollBar;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.border.Border;
import javax.swing.border.EmptyBorder;

import htsjdk.samtools.util.BlockCompressedInputStream;

import org.apache.log4j.Level;

import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.EntryGroup;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureKeyPredicate;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.Selection;
import uk.ac.sanger.artemis.SelectionChangeEvent;
import uk.ac.sanger.artemis.SelectionChangeListener;
import uk.ac.sanger.artemis.SimpleEntryGroup;
import uk.ac.sanger.artemis.components.DisplayAdjustmentEvent;
import uk.ac.sanger.artemis.components.DisplayAdjustmentListener;
import uk.ac.sanger.artemis.components.EntryEdit;
import uk.ac.sanger.artemis.components.EntryFileDialog;
import uk.ac.sanger.artemis.components.FeatureDisplay;
import uk.ac.sanger.artemis.components.FileViewer;
import uk.ac.sanger.artemis.components.IndexReferenceEvent;
import uk.ac.sanger.artemis.components.MessageDialog;
import uk.ac.sanger.artemis.components.MultiComparator;
import uk.ac.sanger.artemis.components.SequenceComboBox;
import uk.ac.sanger.artemis.components.Utilities;
import uk.ac.sanger.artemis.components.alignment.FileSelectionDialog;
import uk.ac.sanger.artemis.components.alignment.LineAttributes;
import uk.ac.sanger.artemis.editor.MultiLineToolTipUI;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.io.RangeVector;
import uk.ac.sanger.artemis.sequence.Bases;
import uk.ac.sanger.artemis.sequence.MarkerRange;
import uk.ac.sanger.artemis.sequence.NoSequenceException;
import uk.ac.sanger.artemis.util.Document;
import uk.ac.sanger.artemis.util.DocumentFactory;
import uk.ac.sanger.artemis.util.FTPSeekableStream;
import uk.ac.sanger.artemis.util.OutOfRangeException;


public class VCFview extends JPanel
             implements DisplayAdjustmentListener, SelectionChangeListener
{
  
  private static final long serialVersionUID = 1L;
  
  private JScrollPane jspView;
  
  private JScrollBar scrollBar;
  private JPanel vcfPanel;
  private AbstractVCFReader vcfReaders[];
  private List<String> vcfFiles;
  private List<Integer> hideVcfList = new Vector<Integer>();
  
  private FeatureDisplay feature_display;
  private Selection selection;
  private int nbasesInView;
  protected int seqLength;
  private EntryGroup entryGroup;
  private String chr;
  private Point lastMousePoint;
  private VCFRecord mouseVCF;
  private int mouseOverIndex = -1;
  private int mouseOverSampleIndex = -1;
  
  private GraphPanel graphPanel;

//record of where a mouse drag starts
  private int dragStart = -1;
  private JPopupMenu popup;
  private JMenu vcfFilesMenu = new JMenu("VCF files");
  private int LINE_HEIGHT = 14;
  
  protected boolean showSynonymous = true;
  protected boolean showNonSynonymous = true;
  protected boolean showDeletions = true;
  protected boolean showInsertions = true;
  protected boolean showMultiAlleles = true;
  protected boolean showHomozygous = true;
  // show variants that do not overlap CDS
  protected boolean showNonOverlappings = true;
  protected boolean showNonVariants = false;
  
  private Map<String, Boolean> manualHash = new HashMap<String, Boolean>();
  
  //private boolean markAsNewStop = false;
  
  private boolean showLabels = false;
  
  private JCheckBoxMenuItem markNewStops =
    new JCheckBoxMenuItem("Mark new stops within CDS features", true);
  
  private static int VARIANT_COLOUR_SCHEME = 0;
  private static int SYN_COLOUR_SCHEME     = 1;
  private static int QUAL_COLOUR_SCHEME    = 2;
  
  private int colourScheme = 0;
  private Color colMap[] = makeColours(Color.RED, 255);
  private Color lighterGrey = new Color(220,220,220);

  private VCFFilter filter;
  Hashtable<String, Integer> offsetLengths = null;
  private boolean concatSequences = false;
  private boolean splitSamples = true;
  
  protected static Pattern tabPattern = Pattern.compile("\t");
  
  public static String VCFFILE_SUFFIX = ".*\\.[bv]{1}cf(\\.gz)*$";
  private static String FILE_SUFFIX = "\\.[bv]{1}cf(\\.gz)*$";

  private List<Integer> cacheVariantLines;
  private SequenceComboBox combo;

  public static org.apache.log4j.Logger logger4j = 
    org.apache.log4j.Logger.getLogger(VCFview.class);

  public VCFview(final JFrame frame,
                 final JPanel vcfPanel,
                 final List<String> vcfFiles, 
                 final int nbasesInView,
                 final int seqLength,
                 final String chr,
                 final String reference,
                 final EntryEdit entry_edit,
                 final FeatureDisplay feature_display)
  {
    super();
    
    this.nbasesInView = nbasesInView;
    this.seqLength = seqLength;
    this.chr = chr;

    this.feature_display = feature_display;
    this.vcfPanel = vcfPanel;
    this.vcfFiles = vcfFiles;
 
    jspView = new JScrollPane(this, 
        JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
        JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
    
    setBackground(Color.white);
    MultiLineToolTipUI.initialize();
    setToolTipText("");
    vcfPanel.setPreferredSize(new Dimension(900, 
        (vcfFiles.size()+1)*(LINE_HEIGHT+5)));
    
    if(feature_display != null)
      this.entryGroup = feature_display.getEntryGroup();
    else if(reference != null)
      this.entryGroup = getReference(reference);
    if(entryGroup != null)
      this.seqLength = entryGroup.getSequenceEntry().getBases().getLength();
    
    try
    {
      vcfReaders = new AbstractVCFReader[vcfFiles.size()];
      
      for(int i=0; i<vcfFiles.size(); i++)
      {
        readHeader(vcfFiles.get(i), i);
      }
    }
    catch(java.lang.UnsupportedClassVersionError err)
    {
      JOptionPane.showMessageDialog(null, 
          "This requires Java 1.9 or higher.", 
          "Check Java Version", JOptionPane.WARNING_MESSAGE);
    }

    vcfPanel.setLayout(new BorderLayout());
    vcfPanel.add(jspView, BorderLayout.CENTER);
    
    JPanel bottomPanel = new JPanel(new BorderLayout());
    graphPanel = new GraphPanel(this);
    graphPanel.setBorder(BorderFactory.createMatteBorder(1, 0, 1, 0, Color.gray));
    graphPanel.setVisible(false);
    
    bottomPanel.add(graphPanel, BorderLayout.CENTER);
    vcfPanel.add(bottomPanel, BorderLayout.SOUTH);
    
    if(this.nbasesInView > this.seqLength)
      this.nbasesInView = this.seqLength/2;

    scrollBar = new JScrollBar(JScrollBar.HORIZONTAL, 1, this.nbasesInView, 1, this.seqLength);
    scrollBar.setUnitIncrement(nbasesInView/20);
    scrollBar.addAdjustmentListener(new AdjustmentListener()
    {
      public void adjustmentValueChanged(AdjustmentEvent e)
      {
        repaint();
      }
    });
    
    //
    //
    addMouseListener(new PopupListener());
    
    //
    createTopPanel(frame, entry_edit);
    createMenus();
    setDisplay();
    
    if(feature_display == null)
    {
      bottomPanel.add(scrollBar, BorderLayout.SOUTH);
      selection = new Selection(null);
    }
    else
    {
      Border empty = new EmptyBorder(0,0,0,0);
      jspView.setBorder(empty);
      selection = feature_display.getSelection();
    }
  }
  
  private void createMenus()
  { 
    // popup menu
    popup = new JPopupMenu();
    
    JMenuItem addVCFMenu = new JMenuItem("Add VCF ...");
    addVCFMenu.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent e)
      {
        FileSelectionDialog fileSelection = new FileSelectionDialog(
            null, true, "VCFview", "VCF");
        List<String> vcfFileList = fileSelection.getFiles(VCFFILE_SUFFIX);
        vcfFiles.addAll(vcfFileList);

        int count = vcfFileList.size();
        int oldSize = vcfReaders.length;
        
        AbstractVCFReader[] trTmp = new AbstractVCFReader[count + vcfReaders.length];
        System.arraycopy(vcfReaders, 0, trTmp, 0, vcfReaders.length);
        vcfReaders = trTmp;
        
        for (int i = 0; i < vcfFileList.size(); i++)
          readHeader(vcfFileList.get(i), i+oldSize);
        
        for(int i=0; i<vcfFileList.size(); i++)
          addToViewMenu(i+oldSize);

        setDisplay();
        repaint();
        jspView.revalidate();
      }
    });
    popup.add(addVCFMenu);
    popup.add(vcfFilesMenu);
    
    for(int i=0; i<vcfFiles.size(); i++)
      addToViewMenu(i);
    
    final JMenu lineHgt = new JMenu("Row Height");
    popup.add(lineHgt);
    final ButtonGroup groupLnHgt = new ButtonGroup();
    for(int i=8; i<37; i+=2)
    {
      final int ii = i;
      final JCheckBoxMenuItem hgtMenu = new JCheckBoxMenuItem(
          Integer.toString(i), (i == LINE_HEIGHT));
      groupLnHgt.add(hgtMenu);
      hgtMenu.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          if(hgtMenu.isSelected())
          {
            LINE_HEIGHT = ii;
            setDisplay();
            revalidate();
          }
        }
      });
      lineHgt.add(hgtMenu);
    }
    
    popup.addSeparator();

    markNewStops.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        repaint();
      }
    });
    popup.add(markNewStops);
    
    
    final JCheckBoxMenuItem split = new JCheckBoxMenuItem("Separate out into samples", splitSamples);
    split.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        splitSamples = split.isSelected();
        setDisplay();
        revalidate();
      }
    });
    popup.add(split);
    
    
    final JMenuItem byQuality = new JMenuItem("Filter ...");
    if(!Options.getOptions().getPropertyTruthValue("java.awt.headless"))
      filter = new VCFFilter(VCFview.this);
    
    byQuality.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        filter.setVisible(true);
      }
    });
    popup.add(byQuality);
    
    final JMenu colourBy = new JMenu("Colour By");
    popup.add(colourBy);
    ButtonGroup group = new ButtonGroup();
    final JRadioButtonMenuItem colByAlt   = new JRadioButtonMenuItem("Variant");
    colByAlt.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        if(colByAlt.isSelected())
          colourScheme = 0;
        repaint();
      }
    });
    colourBy.add(colByAlt);
    
    final JRadioButtonMenuItem colBySyn   = new JRadioButtonMenuItem("Synonymous/Non-synonymous");
    colBySyn.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        if(colBySyn.isSelected())
          colourScheme = 1;
        repaint();
      }
    });
    colourBy.add(colBySyn);
    
    final JRadioButtonMenuItem colByScore = new JRadioButtonMenuItem("Score");
    colByScore.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        if(colByScore.isSelected())
          colourScheme = 2;
        repaint();
      }
    });
    colourBy.add(colByScore);
    
    group.add(colByAlt);
    group.add(colBySyn);
    group.add(colByScore);
    colByAlt.setSelected(true);
    
    popup.addSeparator();

    if (feature_display != null)
    {
      final JMenu create = new JMenu("Create");
      final JMenuItem createTab = new JMenuItem(
          "Features from variants");
      popup.add(create);
      create.add(createTab);
      createTab.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          Container f = getVcfContainer();
          try
          {
            f.setCursor(new Cursor(Cursor.WAIT_CURSOR));
            IOUtils.createFeatures(VCFview.this, entryGroup);
          }
          finally
          {
            f.setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
          }
        }
      });
    }

    final JMenu export = new JMenu("Write");
    popup.add(export);
    
    final JMenuItem exportVCF = new JMenuItem("Filtered VCF");
    exportVCF.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        Container f = getVcfContainer();
        try
        {
          f.setCursor(new Cursor(Cursor.WAIT_CURSOR));
          IOUtils.export(manualHash, vcfFiles, VCFview.this);
        }
        finally
        {
          f.setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
        }
      }
    });
    export.add(exportVCF);
    
    export.add(new JSeparator());
    
    final JMenuItem exportFastaSelected = new JMenuItem("FASTA of selected feature(s) ...");
    exportFastaSelected.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        Container f = getVcfContainer();
        try
        {
          f.setCursor(new Cursor(Cursor.WAIT_CURSOR));
          IOUtils.exportFasta(VCFview.this, selection.getAllFeatures(), false, null);
        }
        finally
        {
          f.setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
        }
      }
    });
    export.add(exportFastaSelected);
    
    final JMenuItem exportFasta = new JMenuItem("FASTA of selected base range ...");
    exportFasta.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        Container f = getVcfContainer();
        try
        {
          f.setCursor(new Cursor(Cursor.WAIT_CURSOR));
          IOUtils.exportFastaByRange(VCFview.this, selection, false, null);
        }
        finally
        {
          f.setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
        }
      }
    });
    export.add(exportFasta);
    
    final JMenuItem viewMinimalFasta = new JMenuItem("FASTA of variant sites only ...");
    viewMinimalFasta.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        Container f = getVcfContainer();
        try
        {
          f.setCursor(new Cursor(Cursor.WAIT_CURSOR));
          IOUtils.exportVariantFasta(VCFview.this);
        }
        finally
        {
          f.setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
        }
      }
    });
    export.addSeparator();
    export.add(viewMinimalFasta);
    
    
    final JMenu view = new JMenu("View");
    popup.add(view);
    final JMenuItem viewFastaSelected = new JMenuItem("FASTA of selected feature(s) ...");
    viewFastaSelected.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        Container f = getVcfContainer();
        try
        {
          f.setCursor(new Cursor(Cursor.WAIT_CURSOR));
          IOUtils.exportFasta(VCFview.this, selection.getAllFeatures(), true, null);
        }
        finally
        {
          f.setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
        }
      }
    });
    view.add(viewFastaSelected);
    
    final JMenuItem viewFasta = new JMenuItem("FASTA of selected base range ...");
    viewFasta.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        Container f = getVcfContainer();
        try
        {
          f.setCursor(new Cursor(Cursor.WAIT_CURSOR));
          IOUtils.exportFastaByRange(VCFview.this, selection, true, null);
        }
        finally
        {
          f.setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
        }
      }
    });
    view.add(viewFasta);


    JMenu graph = new JMenu("Graph");
    popup.add(graph);
    
    final JCheckBoxMenuItem graphSNP = new JCheckBoxMenuItem("SNP");
    final JCheckBoxMenuItem graphDP = new JCheckBoxMenuItem("Depth (DP)");
    final JCheckBoxMenuItem graphSim = new JCheckBoxMenuItem("Base Similarity (%)");
    graph.add(graphSNP);
    graphSNP.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        graphPanel.setVisible(graphSNP.isSelected());
        graphDP.setSelected(false);
        graphSim.setSelected(false);
        graphPanel.setType(0);
        setGraphSize();
      }
    });
    
    graph.add(graphDP);
    graphDP.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        graphPanel.setVisible(graphDP.isSelected());
        graphSNP.setSelected(false);
        graphSim.setSelected(false);
        graphPanel.setType(1);
        setGraphSize();
      }
    });
    
    graph.add(graphSim);
    graphSim.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        graphPanel.setVisible(graphSim.isSelected());
        graphSNP.setSelected(false);
        graphDP.setSelected(false);
        graphPanel.setType(2);
        setGraphSize();
      }
    });
    
    final JMenuItem snpOverview = new JMenuItem("Overview for selected features");
    popup.add(snpOverview);
    snpOverview.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        Container f = getVcfContainer();
        try
        {
          f.setCursor(new Cursor(Cursor.WAIT_CURSOR));
          IOUtils.countVariants(VCFview.this, selection.getAllFeatures());
        }
        catch (IOException e1)
        {
          JOptionPane.showMessageDialog(null, e1.getMessage());
        }
        finally
        {
          f.setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
        }
      }
    });
    
    final JCheckBoxMenuItem labels = new JCheckBoxMenuItem("Show Labels", showLabels);
    labels.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        showLabels = labels.isSelected();
        repaint();
      }
    });
    popup.add(new JSeparator());
    popup.add(labels);
  }
  
  private void setGraphSize()
  {
    repaint();
    if(graphPanel.isVisible())
      graphPanel.setPreferredSize(new Dimension(900, 70));
    else
      graphPanel.setPreferredSize(new Dimension(900, 1));

    vcfPanel.setPreferredSize(new Dimension(900, 
        graphPanel.getPreferredSize().height+getPreferredSize().height));
    vcfPanel.revalidate();
  }
  
  private void createTopPanel(final JFrame frame, final EntryEdit entry_edit)
  {
    final JComponent topPanel;
    if(feature_display != null)
      topPanel = new JPanel(new FlowLayout(FlowLayout.LEADING, 0, 0));
    else
    {
      markNewStops.setSelected(false);
      markNewStops.setEnabled(false);
      topPanel = new JMenuBar();
      if(frame != null)
        frame.setJMenuBar((JMenuBar)topPanel);
      
      JMenu fileMenu = new JMenu("File");
      topPanel.add(fileMenu);
    
      JMenuItem printImage = new JMenuItem("Save As Image Files (png/jpeg)...");
      printImage.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          PrintVCFview part = new PrintVCFview(VCFview.this);
          part.print();
        }
      });
      fileMenu.add(printImage);
      
      JMenuItem printPS = new JMenuItem("Print...");
      printPS.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          PrintVCFview part = new PrintVCFview(VCFview.this);
          part.validate();
          part.doPrintActions();
        }
      });
      fileMenu.add(printPS);
      

      JMenuItem close = new JMenuItem("Close");
      fileMenu.add(close);
      close.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          VCFview.this.setVisible(false);
          Component comp = VCFview.this;
          
          while( !(comp instanceof JFrame) )
            comp = comp.getParent();
          ((JFrame)comp).dispose();
        } 
      });
      
      JButton zoomIn = new JButton("-");
      Insets ins = new Insets(1,1,1,1);
      zoomIn.setMargin(ins);
      zoomIn.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          setZoomLevel((int) (VCFview.this.nbasesInView * 1.1));
        }
      });
      topPanel.add(zoomIn);

      JButton zoomOut = new JButton("+");
      zoomOut.setMargin(ins);
      zoomOut.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          setZoomLevel((int) (VCFview.this.nbasesInView * .9));
        }
      });
      topPanel.add(zoomOut);
    }
    
    combo = new SequenceComboBox(vcfReaders[0].getSeqNames()){
      private static final long serialVersionUID = 1L;
      public void update(IndexReferenceEvent event)
      {
        if(combo.getSelectedItem().equals("Combine References"))
          concatSequences = true;
        else 
        {
          VCFview.this.chr = (String) combo.getSelectedItem();
          concatSequences = false;
        }
        repaint();
      }
    };
    
    if(vcfReaders[0].getSeqNames().length > 1)
      combo.addItem("Combine References");
    
    if(chr == null)
      this.chr = vcfReaders[0].getSeqNames()[0];
    combo.setSelectedItem(this.chr);

    topPanel.add(combo);
    if(topPanel instanceof JPanel)
      vcfPanel.add(topPanel, BorderLayout.NORTH);
    
    // auto hide top panel
    final JCheckBox buttonAutoHide = new JCheckBox("Hide", true);
    buttonAutoHide.setToolTipText("Auto-Hide");
    topPanel.add(buttonAutoHide);
    final MouseMotionListener mouseMotionListener = new MouseMotionListener()
    {
      public void mouseDragged(MouseEvent event)
      {
        handleCanvasMouseDrag(event);
      }
      
      public void mouseMoved(MouseEvent e)
      {
        lastMousePoint = e.getPoint();

        int thisHgt = HEIGHT;
        if (thisHgt < 5)
          thisHgt = 15;

        int y = (int) (e.getY() - jspView.getViewport().getViewRect().getY());
        if (y < thisHgt)
          topPanel.setVisible(true);
        else
        {
          if (buttonAutoHide.isSelected())
            topPanel.setVisible(false);
        }

      }
    };
    addMouseMotionListener(mouseMotionListener);
    
    
    if(feature_display != null)
    {
      JButton close = new JButton("Close");
      topPanel.add(close);
      close.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          setVisible(false);
          if(entry_edit != null)
            entry_edit.setNGDivider();
          else
            vcfPanel.setVisible(false);
        }
      });
    }
  }
  
  /**
   * Create an icon of a box using the given colour.
   * @param c
   * @return
   */
  protected ImageIcon getImageIcon(Color c)
  {
    BufferedImage image = (BufferedImage)this.createImage(10, 10);
    if(image == null)
      return null;
    Graphics2D g2 = image.createGraphics();
    g2.setColor(c);
    g2.fillRect(0, 0, 10, 10);
    return new ImageIcon(image);
  }
  
  private void addToViewMenu(final int thisBamIndex)
  {
    LineAttributes ln[] = GraphPanel.getLineAttributes(vcfReaders.length);
    final JCheckBoxMenuItem cbBam = new JCheckBoxMenuItem(
        getLabel(thisBamIndex), 
        getImageIcon(ln[thisBamIndex].getLineColour()), 
        true);

    vcfFilesMenu.add(cbBam);
    cbBam.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        if(cbBam.isSelected())
          hideVcfList.remove(new Integer(thisBamIndex));
        else
          hideVcfList.add(new Integer(thisBamIndex));
        repaint();
      } 
    });
  }

  /**
   * Test and download if on a http server
   * @param fileName
   * @return
   */
  private String testForURL(String fileName, boolean isBCF)
  {
    if(!fileName.startsWith("http:") && !fileName.startsWith("ftp:"))
      return fileName;

    return download(fileName+".tbi", ".tbi");
  }
  
  protected String download(String f, String suffix)
  {
    try
    {
      final URL urlFile = new URL(f);
      InputStream is = urlFile.openStream();

      // Create temp file.
      String name = urlFile.getFile();
      int ind = name.lastIndexOf('/');
      if(ind > -1)
        name = name.substring(ind+1);
      File bcfFile = File.createTempFile(name.replaceAll("[\\/\\s]", "_"), suffix);
      
      bcfFile.deleteOnExit();

      FileOutputStream out = new FileOutputStream(bcfFile);
      int c;
      while ((c = is.read()) != -1)
        out.write(c);
      out.flush();
      out.close();
      is.close();
      return bcfFile.getAbsolutePath();
    }
    catch(IOException ioe)
    {
      JOptionPane.showMessageDialog(null, 
          "Problem downloading\n"+f, 
          "Problem", JOptionPane.WARNING_MESSAGE);
    }
    return null;
  }
  
  private static EntryGroup getReference(String reference)
  {
    EntryGroup entryGroup = new SimpleEntryGroup();
    final Document entry_document = DocumentFactory.makeDocument(reference);
    final EntryInformation artemis_entry_information =
      Options.getArtemisEntryInformation();
  
    final uk.ac.sanger.artemis.io.Entry new_embl_entry =
      EntryFileDialog.getEntryFromFile(null, entry_document,
                                       artemis_entry_information,
                                       false);
    if(new_embl_entry != null) // the read failed
    {
      Entry entry = null;
      Bases bases = null;
      try
      {
        if (entryGroup.getSequenceEntry() != null)
          bases = entryGroup.getSequenceEntry().getBases();
        if (bases == null)
        {
          entry = new Entry(new_embl_entry);
          bases = entry.getBases();
        }
        else
          entry = new Entry(bases, new_embl_entry);
        entryGroup.add(entry);
      }
      catch (OutOfRangeException e)
      {
        new MessageDialog(null, "read failed: one of the features in "
            + reference + " has an out of range " + "location: "
            + e.getMessage());
      }
      catch (NoSequenceException e)
      {
        // TODO Auto-generated catch block
        e.printStackTrace();
      }
    }
    return entryGroup;
  }

  /**
   * Read the vcf header
   * @param fileName
   * @return
   */
  private void readHeader(String fileName, int index)
  {
    StringBuffer buff = new StringBuffer();
    //buff.append(fileName+"\n");
    try
    {
      if(IOUtils.isBCF(fileName))
      {
        vcfReaders[index] = new BCFReader(fileName);
        String hdr = ((BCFReader)vcfReaders[index]).headerToString();
        if(hdr.indexOf("VCFv4") > -1)
          vcfReaders[index].setVcf_v4(true);
        
        vcfReaders[index].setHeader(hdr);
        return;
      }

      String indexfileName = testForURL(fileName, false);
      BlockCompressedInputStream is;
      if(fileName.startsWith("http")|| fileName.startsWith("ftp"))
      {
        URL url = new URL(fileName);
        if(fileName.startsWith("ftp"))
          vcfReaders[index] = new TabixReader(indexfileName.substring(0, indexfileName.length()-4), new FTPSeekableStream(url));
        else
          vcfReaders[index] = new TabixReader(indexfileName.substring(0, indexfileName.length()-4), url);
        is = new BlockCompressedInputStream(url);
      }
      else
      {
        vcfReaders[index] = new TabixReader(fileName);
        is = new BlockCompressedInputStream(new FileInputStream(fileName));
      }

      String line;
      while( (line = TabixReader.readLine(is) ) != null )
      {
        if(!line.startsWith("#"))
          break;
        
        if(line.indexOf("VCFv4") > -1)
          vcfReaders[index].setVcf_v4(true);
        buff.append(line+"\n");
      }
      is.close();
    }
    catch (IOException e)
    {
      e.printStackTrace();
    }
    
    vcfReaders[index].setHeader(buff.toString());
    return;
  }
  
  /**
   * Set the number of bases being displayed
   * @param nbasesInView
   */
  private void setZoomLevel(final int nbasesInView)
  {
    int startValue = scrollBar.getValue();
    this.nbasesInView = nbasesInView;
    float pixPerBase = getPixPerBaseByWidth(); 
    this.nbasesInView = (int)(getWidth()/pixPerBase);
    
    if(scrollBar != null)
    {
      scrollBar.setValues(startValue, nbasesInView, 1, seqLength);
      scrollBar.setUnitIncrement(nbasesInView/20);
      scrollBar.setBlockIncrement(nbasesInView);
    }
  }
  
  public String getToolTipText()
  {
    if(vcfReaders == null)
      return null;
    
    mouseVCF = null;
    findVariantAtPoint(lastMousePoint);
    if(mouseVCF == null)
      return null;

    String msg = 
           "Seq: "+mouseVCF.getChrom()+"\n";
    msg += "Pos: "+mouseVCF.getPos()+"\n";
    msg += "ID:  "+mouseVCF.getID()+"\n";
    msg += "Variant: "+mouseVCF.getRef()+" -> "+mouseVCF.getAlt().toString()+"\n";
    msg += "Qual: "+mouseVCF.getQuality()+"\n";
    String dp;
    
    if(splitSamples && mouseOverSampleIndex >= 0)
    {
      msg += "Genotype ";
      msg += mouseVCF.getFormat();
      msg += "\n";
      msg += mouseVCF.getFormatValueForSample(mouseOverSampleIndex);
    }
    else if((dp = mouseVCF.getInfoValue("DP")) != null)
    {
      msg += "DP:"+dp;
    }
    return msg;
  }

  
  /**
   * For VCF files with multiple references sequences, calculate
   * the offset from the start of the concatenated sequence for 
   * a given reference.
   * @param refName
   * @return
   */
  protected int getSequenceOffset(String refName)
  {
    if(!concatSequences)
      return 0;
    
    if(offsetLengths == null)
    {   
      String[] contigs = vcfReaders[0].getSeqNames();
      FeatureVector features = entryGroup.getAllFeatures();
      offsetLengths = new Hashtable<String, Integer>(contigs.length);
      for(int i=0; i<contigs.length; i++)
      {
        FeatureContigPredicate predicate = new FeatureContigPredicate(contigs[i]);
        for(int j=0; j<features.size(); j++)
        {
          if(predicate.testPredicate(features.elementAt(j)))
          {
            offsetLengths.put(contigs[i], features.elementAt(j).getFirstBase()-1);
            break;
          }
        }
      }
      
      if(offsetLengths.size() != contigs.length)
      {
        System.err.println("Number of contigs found : "+offsetLengths.size() +
                         "\nNumber of contigs in VCF: "+contigs.length);
        JOptionPane.showMessageDialog(this, 
            "There is a problem matching the reference sequences\n"+
            "to the names in the VCF file. This may mean the labels\n"+
            "on the reference features do not match those in the in\n"+
            "the VCF file.", 
            "Problem Found", JOptionPane.WARNING_MESSAGE);
        concatSequences = false;
        return 0;
      }
    }
    return offsetLengths.get(refName);
  }
  
  public void repaint()
  {
    super.repaint();
    if(graphPanel != null && graphPanel.isVisible())
      graphPanel.repaint();
  }
  
  protected void paintComponent(Graphics g)
  {
    super.paintComponent(g);
    
    Graphics2D g2d = (Graphics2D)g;
    mouseVCF = null;

    float pixPerBase = getPixPerBaseByWidth();
    int start = getBaseAtStartOfView();
    int end   = start+nbasesInView;
    
    drawSelectionRange((Graphics2D)g, pixPerBase, start, end);

    int sumSamples = 0;
    FeatureVector features = getCDSFeaturesInRange(start, end);

    for (int i = 0; i < vcfReaders.length; i++)
    {
      if(hideVcfList.contains(i))
        continue;      
      
      if(concatSequences) 
      {
        String[] contigs = vcfReaders[0].getSeqNames();
        for(int j=0; j<contigs.length; j++)
        {
          int offset = getSequenceOffset(contigs[j]);
          int nextOffset;
          if(j<contigs.length-1)
            nextOffset = getSequenceOffset(contigs[j+1]);
          else
            nextOffset = seqLength;
          
          if( (offset >= start && offset < end) ||
              (offset < start && start < nextOffset) )
          {
            int thisStart = start - offset;
            if(thisStart < 1)
              thisStart = 1;
            int thisEnd   = end - offset;
            
            drawRegion(g2d, contigs[j], thisStart, thisEnd, i, sumSamples, start, pixPerBase, features);
            
          }
        }

      } 
      else
      {
        int thisStart = start;
        if(thisStart < 1)
          thisStart = 1;
        drawRegion(g2d, chr, thisStart, end, i, sumSamples, start, pixPerBase, features); 
      }
      sumSamples += vcfReaders[i].getNumberOfSamples();
    }

    if(feature_display == null)
      drawScale(g2d, start, end, pixPerBase, getHeight());
    
    // show labels for each VCF
    if(showLabels)
       showLabels(g2d);
  }
  
  private void showLabels(Graphics2D g2d)
  {
    int max = 0;
    final FontMetrics fm = getFontMetrics(getFont());
    
    for (int i = 0; i < vcfReaders.length; i++)
    {
      if(hideVcfList.contains(i))
        continue;
      
      if(splitSamples)
      {
        for(int sampleIdx = 0; sampleIdx < vcfReaders[i].getNumberOfSamples(); sampleIdx++)
        {
          if(vcfReaders[i].sampleNames == null)
          {
            int width = fm.stringWidth(getLabel(i));
            if (max < width)
              max = width;
            break;
          }
         
          String labStr = vcfReaders[i].sampleNames[sampleIdx];
          int width = fm.stringWidth(labStr);
          if (max < width)
            max = width;
        }
      }
      else
      {
        String labStr = getLabel(i);
        int width = fm.stringWidth(labStr);
        if (max < width)
          max = width;
      }
    }
    
    Rectangle square = new Rectangle(0, 0, max, getHeight());
    Composite originalComposite = g2d.getComposite();
    g2d.setPaint(Color.lightGray);
    g2d.setComposite(makeComposite(0.75f));
    g2d.fill(square);
    g2d.setComposite(originalComposite);

    g2d.setColor(Color.black);
    g2d.drawLine(max+1, 0, max+1, getHeight());
    
    int sumSample = 0;
    for (int i = 0; i < vcfReaders.length; i++)
    {
      if(hideVcfList.contains(i))
        continue;
      
      if(splitSamples)
      {
        for(int sampleIdx=0; sampleIdx < vcfReaders[i].getNumberOfSamples(); sampleIdx++)
        {
          if(vcfReaders[i].sampleNames == null)
            g2d.drawString(getLabel(i), 1, getYPostion(sumSample+sampleIdx));
          else
            g2d.drawString(vcfReaders[i].sampleNames[sampleIdx], 1, getYPostion(sumSample+sampleIdx));
        }
        sumSample += vcfReaders[i].getNumberOfSamples();
      }
      else
        g2d.drawString(getLabel(i), 1, getYPostion(i));
    }
  }
  
  private String getLabel(int index)
  {
    return vcfReaders[index].getName().replaceAll(FILE_SUFFIX, "");
  }
  
  private AlphaComposite makeComposite(float alpha) 
  {
    return(AlphaComposite.getInstance(AlphaComposite.SRC_OVER, alpha));
  }

  
  private void drawRegion(final Graphics2D g,
                          final String chr,
                          final int sbeg,
                          final int send,
                          final int vcfFileIndex,
                          final int sumSamples,
                          final int start, 
                          final float pixPerBase, 
                          final FeatureVector features) 
  {
    cacheVariantLines = new Vector<Integer>(5);
    try
    {
      VCFRecord record;

      // viewport position and height
      int viewIndex = getHeight()/(LINE_HEIGHT+5) - jspView.getViewport().getViewPosition().y/(LINE_HEIGHT+5);
      int viewHgt = jspView.getViewport().getExtentSize().height/(LINE_HEIGHT+5);

      while((record = vcfReaders[vcfFileIndex].getNextRecord(chr, sbeg, send)) != null)
      {
        int basePosition = record.getPos() + getSequenceOffset(record.getChrom());
        if(!splitSamples)
        {
          drawVariantCall(g, record, start, vcfFileIndex, -1, -1, pixPerBase, features, 
              vcfReaders[vcfFileIndex], basePosition);
          continue;
        }
        
        for(int sampleIndex = 0; sampleIndex < vcfReaders[vcfFileIndex].getNumberOfSamples(); sampleIndex++)
        {
          if(sampleIndex+sumSamples <= viewIndex+2 && sampleIndex+sumSamples >= viewIndex-viewHgt-2)
          {
            drawVariantCall(g, record, start, vcfFileIndex, sampleIndex, sumSamples, pixPerBase, features, 
              vcfReaders[vcfFileIndex], basePosition);
          }
        }
      }
    }
    catch (IOException e)
    {
      logger4j.warn(chr+":"+sbeg+"-"+send+"\n"+e.getMessage());
      e.printStackTrace();
    }
  }
  
  protected FeatureVector getCDSFeaturesInRange(int start, int end)
  {
    if(entryGroup == null)
      return null;
    try
    {
      Range range = new Range(start, end);
      FeatureVector features = entryGroup.getFeaturesInRange(range);
           
      FeatureKeyPredicate predicate = new FeatureKeyPredicate(Key.CDS);
      final FeatureVector cdsFeatures = new FeatureVector();

      for(int i=0; i<features.size(); i++)
      {
        final Feature this_feature = features.elementAt(i);
        if(predicate.testPredicate(this_feature))
          cdsFeatures.add(this_feature);
      }
      return cdsFeatures;
    }
    catch (OutOfRangeException e1)
    {
      e1.printStackTrace();
    }
    return null;
  }
  
  /**
   * Highlight a selected range
   * @param g2
   * @param pixPerBase
   * @param start
   * @param end
   */
  private void drawSelectionRange(Graphics2D g2, float pixPerBase, int start, int end)
  {
    if(getSelection() != null)
    {
      Range selectedRange = getSelection().getSelectionRange();

      if(selectedRange != null)
      {
        int rangeStart = selectedRange.getStart();
        int rangeEnd   = selectedRange.getEnd();
        
        if(end < rangeStart || start > rangeEnd)
          return;
        
        int x = (int) (pixPerBase*(rangeStart-getBaseAtStartOfView()));
        int width = (int) (pixPerBase*(rangeEnd-rangeStart+1));

        g2.setColor(Color.pink);
        g2.fillRect(x, 0, width, getHeight());
      }
    }
  }
  
  private Selection getSelection()
  {
    return selection;
  }
  
  protected int getBaseAtStartOfView()
  {
    if(feature_display != null)
      return feature_display.getForwardBaseAtLeftEdge();
    else
      return scrollBar.getValue();
  }
  
  protected boolean showVariant(final VCFRecord record, final FeatureVector features, final int basePosition, 
      final AbstractVCFReader vcfReader, final int nsample, final int vcfIndex)
  { 
    return VCFFilter.passFilter(manualHash, record, vcfReader, features, basePosition, nsample, vcfIndex);
  }
  
  private void setAsStop(VCFRecord record, FeatureVector features, int basePosition, AbstractVCFReader vcfReader)
  {
    if(!record.isMarkAsNewStop() &&
        markNewStops.isSelected() &&
       !record.getAlt().isDeletion(vcfReader.isVcf_v4()) && 
       !record.getAlt().isInsertion(vcfReader.isVcf_v4()) && 
        record.getAlt().length() == 1 && 
        record.getRef().length() == 1)
    {
      short isSyn = record.getSynFlag(features, basePosition);
      if(isSyn == 2)
        record.setMarkAsNewStop(true);
    }
  }
  
  protected static boolean isOverlappingFeature(FeatureVector features, int basePosition)
  {
    for(int i = 0; i<features.size(); i++)
    {
      Feature feature = features.elementAt(i);
      if(feature.getRawFirstBase() < basePosition && feature.getRawLastBase() > basePosition)
      {
        RangeVector ranges = feature.getLocation().getRanges();
        for(int j=0; j< ranges.size(); j++)
        {
          Range range = (Range) ranges.get(j);
          if(range.getStart() < basePosition && range.getEnd() > basePosition)
            return true;
        }
      }
    }
    return false;
  }
  
  protected static boolean isOverlappingFeature(List<CDSFeature> cdsFeatures, int basePosition) {
      for (CDSFeature cdsFeature : cdsFeatures) {
          if (cdsFeature.firstBase < basePosition && cdsFeature.lastBase > basePosition) 
          {
              for(int i = 0 ; i < cdsFeature.ranges.size()  ; ++i) 
              {
                   Range range = (Range)cdsFeature.ranges.elementAt(i);
                   if (range.getStart() < basePosition && range.getEnd() > basePosition) 
                       return true;
              }
          }
      }
      return false;
  }
  
  /**
   * Draw the VCF record
   * @param g
   * @param record
   * @param start
   * @param vcfIndex
   * @param sampleIndex
   * @param sumSamples
   * @param pixPerBase
   * @param features
   * @param vcfReader
   * @param basePosition
   */
  private void drawVariantCall(final Graphics2D g, 
                               final VCFRecord record, 
                               final int start, 
                               final int vcfIndex,
                               final int sampleIndex,
                               final int sumSamples,
                               final float pixPerBase, 
                               final FeatureVector features, 
                               final AbstractVCFReader vcfReader,
                               final int basePosition)
  {
    boolean show = showVariant(record, features, basePosition, vcfReader, sampleIndex, vcfIndex);
    if( !show )
      return;
    
    setAsStop(record, features, basePosition, vcfReader);
    final int pos[];
    if(sampleIndex < 0)
      pos = getScreenPosition(basePosition, pixPerBase, start, vcfIndex);
    else
      pos = getScreenPosition(basePosition, pixPerBase, start, sampleIndex+sumSamples);

    if (colourScheme == QUAL_COLOUR_SCHEME)
      g.setColor(getQualityColour(record));
    else if (record.getAlt().isDeletion(vcfReader.isVcf_v4()))
      g.setColor(Color.gray);
    else if (record.getAlt().isInsertion(vcfReader.isVcf_v4()))
      g.setColor(Color.magenta);
    else if (record.getAlt().isMultiAllele(sampleIndex))
    {
      g.setColor(Color.orange);
      g.fillArc(pos[0] - 3, pos[1] - LINE_HEIGHT - 3, 6, 6, 0, 360);
    }
    else if (record.getAlt().length() == 1 && record.getRef().length() == 1)
    {
      g.setColor(getColourForSNP(record, features, basePosition));
      if (record.getAlt().isNonVariant())
      {
        // use the cache to avoid drawing over a variant with a non-variant
        if (!cacheVariantLines.contains(pos[0]))
          g.drawLine(pos[0], pos[1], pos[0], pos[1] - LINE_HEIGHT + 6);
        return;
      }
    }
    else
      g.setColor(Color.pink);

    if (record.isMarkAsNewStop())
      g.fillArc(pos[0] - 3, pos[1] - (LINE_HEIGHT / 2) - 3, 6, 6, 0, 360);

    if (cacheVariantLines.size() == 5)
      cacheVariantLines.clear();
    cacheVariantLines.add(pos[0]);

    g.drawLine(pos[0], pos[1], pos[0], pos[1] - LINE_HEIGHT);
  }
  
  /**
   * Determine the colour depending on the colour scheme in use.
   * @param record
   * @param features
   * @param basePosition
   * @return
   */
  private Color getColourForSNP(VCFRecord record, FeatureVector features, int basePosition)
  {
    if(colourScheme == VARIANT_COLOUR_SCHEME)
      return getVariantColour(record.getAlt().toString());
    else if(colourScheme == SYN_COLOUR_SCHEME)  // synonymous / non-synonymous
    {
      if(!record.getAlt().isNonVariant())
      {
        short synFlag = record.getSynFlag(features, basePosition);
        if(synFlag == 1)
          return Color.red;
        else if(synFlag == 0 || synFlag == 2)
          return Color.blue;
      }
      return getVariantColour(record.getAlt().toString());
    }
    else // score
      return getQualityColour(record);
  }
  
  private Color getQualityColour(VCFRecord record)
  {
    if(colMap == null)
      colMap = makeColours(Color.RED, 255);
    int idx = (int) record.getQuality()-1;
    if(idx > colMap.length-1)
      idx = colMap.length-1;
    else if(idx < 0)
      idx = 0;
    return colMap[idx];
  }
  
  private Color getVariantColour(String variant)
  {
    if(variant.equals("C"))
      return Color.red;
    else if(variant.equals("A"))
      return Color.green;
    else if(variant.equals("G"))
      return Color.blue;
    else if(variant.equals("T"))
      return Color.black;
    else
      return lighterGrey; // non-variant
  }
  
  /**
   * Generate the colours for heat map plots.
   * @param col
   * @param NUMBER_OF_SHADES
   * @return
   */
  private Color[] makeColours(Color col, int NUMBER_OF_SHADES)
  {
    Color definedColour[] = new Color[NUMBER_OF_SHADES];
    for(int i = 0; i < NUMBER_OF_SHADES; ++i)
    {
      int R = col.getRed();
      int G = col.getGreen();
      int B = col.getBlue();

      float scale = ((float)(NUMBER_OF_SHADES-i) * (float)(255 / NUMBER_OF_SHADES )) ;
      
      if((R+scale) <= 255)
        R += scale;
      else
        R = 254;
      if((G+scale) <= 255)
        G += scale;
      else
        G = 254;
      if((B+scale) <= 255)
        B += scale;
      else
        B = 254;

      definedColour[i] = new Color(R,G,B);
    }
    return definedColour;
  }
  
  private int[] getScreenPosition(int base, float pixPerBase, int start, int vcfFileIndex)
  {
    int pos[] = new int[2];
    pos[0] = Math.round((base - start)*pixPerBase);
    pos[1] = getYPostion(vcfFileIndex); 
    return pos;
  }
  
  private int getYPostion(int vcfFileIndex)
  {
    int pos = 0;
    if(hideVcfList.size() == 0 || splitSamples)
      pos = vcfFileIndex;
    else
    {
      for(int i=0; i<vcfFileIndex; i++)
        if(!hideVcfList.contains(i))
          pos++; 
    }

    return getHeight() - 15 - (pos*(LINE_HEIGHT+5));
  }
  
  private void findVariantAtPoint(Point mousePoint)
  {
    float pixPerBase = getPixPerBaseByWidth();
    int startBase = getBaseAtStartOfView();
    int start = startBase + (int)(mousePoint.getX()/pixPerBase) - 20;
    int end   = start+20;
    FeatureVector features = getCDSFeaturesInRange(start, end);
    int sumSamples = 0;
    
    for (int i = 0; i < vcfReaders.length; i++)
    {

        if(concatSequences) 
        {
          String[] contigs = vcfReaders[0].getSeqNames();
          for(int k=0; k<contigs.length; k++)
          {
            int offset = getSequenceOffset(contigs[k]);
            int nextOffset;
            if(k<contigs.length-1)
              nextOffset = getSequenceOffset(contigs[k+1]);
            else
              nextOffset = seqLength;
            
            if( (offset >= start && offset < end) ||
                (offset < start && start < nextOffset) )
            {
              int thisStart = start - offset;
              if(thisStart < 1)
                thisStart = 1;
              int thisEnd   = end - offset;
              searchRegion(contigs[k], thisStart, thisEnd, i, sumSamples, mousePoint, features, startBase, pixPerBase);
            }
          }
        } 
        else
        {
          int thisStart = start;
          if(thisStart < 1)
            thisStart = 1;
          searchRegion(chr, thisStart, end, i, sumSamples, mousePoint, features, startBase, pixPerBase);
        }


      sumSamples += vcfReaders[i].getNumberOfSamples();
    }
  }
  
  private void searchRegion(final String chr, 
                            final int sbeg, final int send, 
                            final int fileIndex, 
                            final int sumSamples,
                            final Point mousePoint, FeatureVector features,
                            int start, float pixPerBase) 
  {
    try
    {
      VCFRecord bcfRecord;
      while((bcfRecord = vcfReaders[fileIndex].getNextRecord(chr, sbeg, send)) != null)
      {
        
        if(splitSamples)
        {
          for(int sampleIndex=0; sampleIndex<vcfReaders[fileIndex].getNumberOfSamples(); sampleIndex++)
          {
            int ypos = getYPostion(sampleIndex+sumSamples);
            if(mousePoint.getY() > ypos &&
               mousePoint.getY() < ypos-LINE_HEIGHT)
              continue;
            
            isMouseOver(mousePoint, bcfRecord, features, fileIndex, sampleIndex, sumSamples, start, pixPerBase, vcfReaders[fileIndex]); 
          }
        }
        else
        {
          int ypos = getYPostion(fileIndex);
          if(mousePoint.getY() > ypos &&
             mousePoint.getY() < ypos-LINE_HEIGHT)
            continue;
          
          isMouseOver(mousePoint, bcfRecord, features, fileIndex, -1, sumSamples, start, pixPerBase, vcfReaders[fileIndex]);
        }
        
        
      }
    }
    catch (IOException e)
    {
      logger4j.warn(chr+":"+sbeg+"-"+send+"\n"+e.getMessage());
      e.printStackTrace();
    }
  }
  
  private void isMouseOver(final Point mousePoint, 
                           final VCFRecord record, 
                           final FeatureVector features, 
                           final int vcfFileIndex, 
                           final int sampleIndex,
                           final int sumSamples,
                           final int start, final float pixPerBase, 
                           final AbstractVCFReader vcfReader)
  {
    int basePosition = record.getPos() + getSequenceOffset(record.getChrom());
    if( !showVariant(record, features, basePosition, vcfReader, sampleIndex, vcfFileIndex) )
      return;

    int pos[];
    if(!splitSamples)
      pos = getScreenPosition(basePosition, pixPerBase, start, vcfFileIndex);
    else
      pos = getScreenPosition(basePosition, pixPerBase, start, sampleIndex+sumSamples);

    if(mousePoint != null &&
       mousePoint.getY() < pos[1] &&
       mousePoint.getY() > pos[1]-LINE_HEIGHT &&
       mousePoint.getX() > pos[0]-3 &&
       mousePoint.getX() < pos[0]+3)
     {
       mouseVCF = record;
       mouseOverIndex = vcfFileIndex;
       mouseOverSampleIndex = sampleIndex;
     }
  }
  
  private void drawScale(Graphics2D g2, int start, int end, float pixPerBase, int ypos)
  {
    g2.setColor(Color.black);
    g2.drawLine( 0, ypos-14,
                 (int)((end - start)*pixPerBase),   ypos-14);
    int interval = end-start;

    if(interval > 20000000)
      drawTicks(g2, start, end, pixPerBase, 10000000, ypos);
    else if(interval > 4000000)
      drawTicks(g2, start, end, pixPerBase, 2000000, ypos);
    else if(interval > 800000)
      drawTicks(g2, start, end, pixPerBase, 400000, ypos);
    else if(interval > 160000)
      drawTicks(g2, start, end, pixPerBase, 80000, ypos);
    else if(interval > 50000)
      drawTicks(g2, start, end, pixPerBase, 25000, ypos);
    else if(interval > 16000)
      drawTicks(g2, start, end, pixPerBase, 8000, ypos);
    else if(interval > 4000)
      drawTicks(g2, start, end, pixPerBase, 2000, ypos);
    else if(interval > 1000)
      drawTicks(g2, start, end, pixPerBase, 500, ypos);
    else
      drawTicks(g2, start, end, pixPerBase, 100, ypos);
  }
  
  /**
   *  Handle a mouse drag event on the drawing canvas.
   **/
  private void handleCanvasMouseDrag(final MouseEvent event)
  {
    if(event.getButton() == MouseEvent.BUTTON3) 
      return;

    if(event.getClickCount() > 1)
    {
      getSelection().clear();
      repaint();
      return;  
    }

    highlightRange(event, 
        MouseEvent.BUTTON1_DOWN_MASK & MouseEvent.BUTTON2_DOWN_MASK);
  }
  
  /**
   * 
   * @param event
   * @param onmask
   */
  private void highlightRange(MouseEvent event, int onmask)
  {
    if(entryGroup == null)
      return;
    float pixPerBase = getPixPerBaseByWidth();
    int start = (int) ( ( (event.getPoint().getX())/pixPerBase) + getBaseAtStartOfView() );
    
    if(start < 1)
      start = 1;
    if(start > seqLength)
      start = seqLength;
    
    if (dragStart < 0 && (event.getModifiersEx() & onmask) == onmask)
      dragStart = start;
    else if((event.getModifiersEx() & onmask) != onmask)
      dragStart = -1;
    
    MarkerRange drag_range;
    try
    {
      if(dragStart < 0)
        drag_range = new MarkerRange (entryGroup.getSequenceEntry().getBases().getForwardStrand(), start, start);
      else
        drag_range = new MarkerRange (entryGroup.getSequenceEntry().getBases().getForwardStrand(), dragStart, start);
      getSelection().setMarkerRange(drag_range);
      repaint();
    }
    catch (OutOfRangeException e)
    {
      e.printStackTrace();
    }
  }
 
  private void drawTicks(Graphics2D g2, int start, int end, 
                         float pixPerBase, int division, int ypos)
  {
    int markStart = (Math.round(start/division)*division);
    
    if(markStart < 1)
      markStart = 1;
    
    int sm = markStart-(division/2);
    float x;
    if(sm > start)
    {
      x = (sm-start)*pixPerBase;
      g2.drawLine((int)x, ypos-14,(int)x, ypos-12);
    }
    
    for(int m=markStart; m<end; m+=division)
    {
      x = (m-start)*pixPerBase;
      g2.drawString(Integer.toString(m), x, ypos-1);
      g2.drawLine((int)x, ypos-14,(int)x, ypos-11);
      
      sm = m+(division/2);
      
      if(sm < end)
      {
        x = (sm-start)*pixPerBase;
        g2.drawLine((int)x, ypos-14,(int)x, ypos-12);
      }
      
      if(m == 1)
        m = 0;
    }
  }

  protected float getPixPerBaseByWidth()
  {
    return (float)vcfPanel.getWidth() / (float)nbasesInView;
  }
  
  protected int getBasesInView()
  {
    return nbasesInView;
  }
  
  
  protected EntryGroup getEntryGroup()
  {
    return entryGroup;
  }
  
  protected String getChr()
  {
    return chr;
  }
  
  protected boolean isConcatenate()
  {
    return concatSequences;
  }
  
  /**
   * @return the combo
   */
  public SequenceComboBox getCombo()
  {
    return combo;
  }
  
  /**
   * @return the vcfReaders
   */
  protected AbstractVCFReader[] getVcfReaders()
  {
    return vcfReaders;
  }
  
  
  private Container getVcfContainer()
  {
    try
    {
      Frame fs[] = JFrame.getFrames();
      for(Frame f: fs)
      {
        if( f instanceof JFrame && 
           ((JFrame)f) instanceof EntryEdit ||
           ((JFrame)f) instanceof MultiComparator)
          return ((JFrame)f).getContentPane();
      }
    }
    catch(Exception e){}
    return VCFview.this;
  }
  
  /**
   * Popup menu listener
   */
   class PopupListener extends MouseAdapter
   {
     JMenuItem showDetails;
     JMenuItem annotate;
     
     public void mouseClicked(MouseEvent e)
     {
       if(e.isPopupTrigger() ||
          e.getButton() == MouseEvent.BUTTON3)
         return;
       
       if(e.getClickCount() > 1)
       {
         getSelection().clear();
         repaint();
       }
     }
     
     public void mousePressed(MouseEvent e)
     {
       maybeShowPopup(e);
     }

     public void mouseReleased(MouseEvent e)
     {
       dragStart = -1;
       maybeShowPopup(e);
     }

     private void maybeShowPopup(MouseEvent e)
     {
       if(e.isPopupTrigger())
       {
         if(showDetails != null)
         {
           popup.remove(showDetails);
           popup.remove(annotate);
         }
         
         mouseVCF = null;
         findVariantAtPoint(e.getPoint());
         final VCFRecord thisMouseVCF = mouseVCF;
         final int thisMouseOverIndex = mouseOverIndex;
         if( thisMouseVCF != null )
         {
           showDetails = new JMenuItem("Show details of : "+
               thisMouseVCF.getChrom()+":"+thisMouseVCF.getPos()+" "+thisMouseVCF.getID());
           showDetails.addActionListener(new ActionListener()
           {
             public void actionPerformed(ActionEvent e) 
             {
               FileViewer viewDetail = new FileViewer(
                   thisMouseVCF.getChrom()+":"+thisMouseVCF.getPos()+" "+thisMouseVCF.getID(), true, false, true);
               
               viewDetail.appendString(vcfReaders[mouseOverIndex].getHeader()+"\n", Level.INFO);
               
               viewDetail.appendString("Seq   : "+thisMouseVCF.getChrom()+"\n", Level.DEBUG);
               viewDetail.appendString("Pos   : "+thisMouseVCF.getPos()+"\n", Level.DEBUG);
               viewDetail.appendString("ID    : "+thisMouseVCF.getID()+"\n", Level.DEBUG);
               viewDetail.appendString("Ref   : "+thisMouseVCF.getRef()+"\n", Level.DEBUG);
               viewDetail.appendString("Alt   : "+thisMouseVCF.getAlt().toString()+"\n", Level.DEBUG);
               viewDetail.appendString("Qual  : "+thisMouseVCF.getQuality()+"\n", Level.DEBUG);
               viewDetail.appendString("Filter: "+thisMouseVCF.getFilter()+"\n", Level.DEBUG);
               viewDetail.appendString("Info  : "+thisMouseVCF.getInfo()+"\n", Level.DEBUG);
               
               if(thisMouseVCF.getFormat() != null)
               {
                 viewDetail.appendString("\nGenotype information:\n", Level.INFO);
                 viewDetail.appendString("Format: "+thisMouseVCF.getFormat()+"\n", Level.DEBUG);
                 viewDetail.appendString(thisMouseVCF.getSampleDataString()+"\n", Level.DEBUG);
               }
               
               viewDetail.getTextPane().setCaretPosition(0);
             }  
           });
           popup.add(showDetails);
           
           annotate = new JMenuItem("Manual PASS / FAIL");
           annotate.addActionListener(new ActionListener()
           {
             public void actionPerformed(ActionEvent e) 
             {
               new ManualAnnotation(thisMouseVCF, thisMouseOverIndex);
             }
           });
          popup.add(annotate);
         }
         popup.show(e.getComponent(),
             e.getX(), e.getY());
       }
     }
   }
   
   private void setDisplay()
   {
     Dimension d = new Dimension();
     
     if(splitSamples)
     {
       int count = 0;
       for(int i=0; i<vcfReaders.length; i++)
         count += vcfReaders[i].getNumberOfSamples();
       d.setSize(nbasesInView*getPixPerBaseByWidth(), (count+1)*(LINE_HEIGHT+5));
     }
     else
       d.setSize(nbasesInView*getPixPerBaseByWidth(), (vcfReaders.length+1)*(LINE_HEIGHT+5));
     setPreferredSize(d);
   }
   
   public void displayAdjustmentValueChanged(DisplayAdjustmentEvent event)
   {
     nbasesInView = feature_display.getMaxVisibleBases();
     setDisplay();
     repaint();
   }

   public void selectionChanged(SelectionChangeEvent event)
   {
     repaint();
   }
   
   /**
    * Manual annotation of a variant record.
    */
   class ManualAnnotation extends JFrame
   {
     private static final long serialVersionUID = 1L;

    ManualAnnotation(final VCFRecord thisMouseVCF, final int thisMouseOverIndex)
     {
       super("Manual PASS / FAIL");
       final JPanel pane = (JPanel) getContentPane();
       pane.setLayout(new BorderLayout());
       final JPanel panel = new JPanel(new GridBagLayout());
       final JScrollPane jsp = new JScrollPane(panel);
       jsp.setPreferredSize(new Dimension(350, 180));
       pane.add(jsp, BorderLayout.NORTH);
       
       final GridBagConstraints c = new GridBagConstraints();
       c.gridx = 0;
       c.gridy = 0;
       c.gridwidth = 3;
       c.anchor = GridBagConstraints.NORTHWEST;
       panel.add(new JLabel("Seq   : "+thisMouseVCF.getChrom()), c);
       c.gridy++;
       panel.add(new JLabel("Pos   : "+thisMouseVCF.getPos()), c);
       c.gridy++;
       panel.add(new JLabel("ID    : "+thisMouseVCF.getPos()), c);
       c.gridy++;
       panel.add(new JLabel("Ref   : "+thisMouseVCF.getRef()), c);
       c.gridy++;
       panel.add(new JLabel("Alt   : "+thisMouseVCF.getAlt().toString()), c);
       c.gridy++;
       panel.add(new JLabel("Qual  : "+thisMouseVCF.getQuality()), c);
       c.gridy++;
       panel.add(new JLabel("Filter: "+thisMouseVCF.getFilter()), c);
       c.gridy++;
       panel.add(new JLabel("Info  : "+thisMouseVCF.getInfo()), c);
       
       final JRadioButton passB = new JRadioButton("PASS");
       final JRadioButton failB = new JRadioButton("FAIL");
       final JRadioButton noManualB = new JRadioButton("NO MANUAL FILTER", true);

       final ButtonGroup group = new ButtonGroup();
       group.add(passB);
       group.add(failB);
       group.add(noManualB);

       int res = VCFFilter.checkManualHash(manualHash, thisMouseVCF, thisMouseOverIndex, false);
       switch(res)
       {
         case 1:  passB.setSelected(true);
                  break;
         case 2:  failB.setSelected(true);
                  break;
         default: noManualB.setSelected(true);
                  break;
       }
       
       c.gridy++;
       panel.add(Box.createVerticalStrut(10), c);
       
       c.gridwidth = 1;
       c.gridy++;
       panel.add(passB, c);
       c.gridx++;
       panel.add(failB, c);
       c.gridx++;
       panel.add(noManualB, c);

       final JButton apply = new JButton("Apply");
       apply.addActionListener(new ActionListener()
       {
         public void actionPerformed(ActionEvent arg0)
         {
           setVisible(false);
           if(passB.isSelected())
             manualHash.put(thisMouseVCF.getPos()+":"+thisMouseVCF.getChrom()+":"+thisMouseOverIndex, true);
           else if(failB.isSelected())
             manualHash.put(thisMouseVCF.getPos()+":"+thisMouseVCF.getChrom()+":"+thisMouseOverIndex, false);
           else
             manualHash.remove(thisMouseVCF.getPos()+":"+thisMouseVCF.getChrom()+":"+thisMouseOverIndex);
           VCFview.this.repaint();
           ManualAnnotation.this.dispose();
         }
       });
       pane.add(apply, BorderLayout.SOUTH);
       pack();
       Utilities.centreFrame(this);
       setVisible(true);
     }
   }
   
  public static void main(String args[])
  {
    List<String> vcfFileList = new Vector<String>();
    String reference = null;
    if(args.length == 0)
    {
      System.setProperty("default_directory", System.getProperty("user.dir"));
      FileSelectionDialog fileSelection = new FileSelectionDialog(
          null, true, "VCFview", "VCF");
      vcfFileList = fileSelection.getFiles(VCFFILE_SUFFIX);
      reference = fileSelection.getReferenceFile();
      if(reference.equals(""))
        reference = null;
      
      if(vcfFileList == null || vcfFileList.size() < 1)
        System.exit(0);
    }
    else if(!args[0].startsWith("-"))
    {
      for(int i=0; i< args.length; i++)
        vcfFileList.add(args[i]);
    }
    
    int nbasesInView = 5000000;
    
    for(int i=0;i<args.length; i++)
    {
      if(args[i].equals("-f"))
      {
        while(i < args.length-1 && !args[++i].startsWith("-"))
        {
          String filename = args[i];
          if(FileSelectionDialog.isListOfFiles(filename))
            vcfFileList.addAll(FileSelectionDialog.getListOfFiles(filename));
          else
            vcfFileList.add(filename);
        }
        --i;
      }
      else if(args[i].equals("-r"))
        reference = args[++i];
      else if(args[i].equals("-v"))
        nbasesInView = Integer.parseInt(args[++i]);
      else if(args[i].startsWith("-h"))
      { 
        System.out.println("-h\t show help");
        
        System.out.println("-f\t VCF file to display");
        System.out.println("-r\t reference file (optional)");
        System.out.println("-v\t number of bases to display in the view (optional)");
        //System.out.println("-t\t chr:start-end - this writes out the given region");

        System.exit(0);
      }
    }
    
    JFrame f = new JFrame();
    new VCFview(f, (JPanel) f.getContentPane(), vcfFileList, 
          nbasesInView, 100000000, null, reference, null, null);
    
    f.pack();
    f.setVisible(true);
  }
}