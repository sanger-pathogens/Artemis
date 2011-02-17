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
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Insets;
import java.awt.Point;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.AdjustmentEvent;
import java.awt.event.AdjustmentListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionListener;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.Hashtable;
import java.util.List;
import java.util.Vector;
import java.util.regex.Pattern;

import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.JScrollBar;
import javax.swing.JScrollPane;
import javax.swing.border.Border;
import javax.swing.border.EmptyBorder;

import net.sf.samtools.util.BlockCompressedInputStream;

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
import uk.ac.sanger.artemis.components.EntryFileDialog;
import uk.ac.sanger.artemis.components.FeatureDisplay;
import uk.ac.sanger.artemis.components.FileViewer;
import uk.ac.sanger.artemis.components.MessageDialog;
import uk.ac.sanger.artemis.components.alignment.FileSelectionDialog;
import uk.ac.sanger.artemis.components.variant.BCFReader.BCFReaderIterator;
import uk.ac.sanger.artemis.editor.MultiLineToolTipUI;
import uk.ac.sanger.artemis.io.EmblStreamFeature;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.InvalidRelationException;
import uk.ac.sanger.artemis.io.Key;
import uk.ac.sanger.artemis.io.Location;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.io.RangeVector;
import uk.ac.sanger.artemis.sequence.Bases;
import uk.ac.sanger.artemis.sequence.MarkerRange;
import uk.ac.sanger.artemis.sequence.NoSequenceException;
import uk.ac.sanger.artemis.util.Document;
import uk.ac.sanger.artemis.util.DocumentFactory;
import uk.ac.sanger.artemis.util.OutOfRangeException;


public class VCFview extends JPanel
             implements DisplayAdjustmentListener, SelectionChangeListener
{
  private static final long serialVersionUID = 1L;
  private JScrollBar scrollBar;
  private JPanel vcfPanel;
  private AbstractVCFReader vcfReaders[];
  private List<String> vcfFiles;
  private String header[];
  private FeatureDisplay feature_display;
  private Selection selection;
  private int nbasesInView;
  private int seqLength;
  private EntryGroup entryGroup;
  private String chr;
  private VCFRecord mouseVCF;
  private int mouseOverIndex = -1;
  
  private boolean vcf_v4 = false;
//record of where a mouse drag starts
  private int dragStart = -1;
  private JPopupMenu popup;
  private int LINE_HEIGHT = 15;
  
  protected boolean showSynonymous = true;
  protected boolean showNonSynonymous = true;
  protected boolean showDeletions = true;
  protected boolean showInsertions = true;
  protected boolean showMultiAlleles = true;
  // show variants that do not overlap CDS
  protected boolean showNonOverlappings = true;
  protected boolean showNonVariants = false;
  
  private boolean markAsNewStop = false;
  private JCheckBoxMenuItem markNewStops =
    new JCheckBoxMenuItem("Mark new stops within CDS features", true);
  
  private static int VARIANT_COLOUR_SCHEME = 0;
  private static int SYN_COLOUR_SCHEME     = 1;
  private static int QUAL_COLOUR_SCHEME    = 2;
  
  private int colourScheme = 0;
  private Color colMap[] = makeColours(Color.RED, 255);
  
  Hashtable<String, Integer> offsetLengths = null;
  private boolean concatSequences = false;
  protected static Pattern tabPattern = Pattern.compile("\t");
  
  public static String VCFFILE_SUFFIX = ".*\\.[bv]{1}cf(\\.gz)*$";
  
  public static org.apache.log4j.Logger logger4j = 
    org.apache.log4j.Logger.getLogger(VCFview.class);

  public VCFview(final JFrame frame,
                 final JPanel vcfPanel,
                 final List<String> vcfFiles, 
                 final int nbasesInView,
                 final int seqLength,
                 final String chr,
                 final String reference,
                 final FeatureDisplay feature_display)
  {
    super();
    
    this.nbasesInView = nbasesInView;
    this.seqLength = seqLength;
    this.chr = chr;
    this.feature_display = feature_display;
    this.vcfPanel = vcfPanel;
    this.vcfFiles = vcfFiles;
 
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
      header = new String[vcfFiles.size()];
      
      for(int i=0; i<vcfFiles.size(); i++)
      {
        header[i] = readHeader(vcfFiles.get(i), i);
      }
    }
    catch(java.lang.UnsupportedClassVersionError err)
    {
      JOptionPane.showMessageDialog(null, 
          "This requires Java 1.6 or higher.", 
          "Check Java Version", JOptionPane.WARNING_MESSAGE);
    }
    
    final JScrollPane jspView = new JScrollPane(this, 
        JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED,
        JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
    
    vcfPanel.setLayout(new BorderLayout());
    vcfPanel.add(jspView, BorderLayout.CENTER);
    
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
    createMenus(frame, jspView);
    setDisplay();
    
    if(feature_display == null)
    {
      vcfPanel.add(scrollBar, BorderLayout.SOUTH);
      frame.pack();
      frame.setVisible(true);
      selection = new Selection(null);
    }
    else
    {
      Border empty = new EmptyBorder(0,0,0,0);
      jspView.setBorder(empty);
      selection = feature_display.getSelection();
    }
  }
  
  private void createMenus(JFrame frame, final JScrollPane jspView)
  {
    final JComponent topPanel;
    if(feature_display != null)
      topPanel = new JPanel(new FlowLayout(FlowLayout.LEADING, 0, 0));
    else
    {
      markNewStops.setSelected(false);
      markNewStops.setEnabled(false);
      topPanel = new JMenuBar();
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
    
    final JComboBox combo = new JComboBox(vcfReaders[0].getSeqNames());
    
    if(vcfReaders[0].getSeqNames().length > 1)
      combo.addItem("Combine References");
    
    if(chr == null)
      this.chr = vcfReaders[0].getSeqNames()[0];

    combo.setSelectedItem(this.chr);
    combo.setEditable(false);
    combo.setMaximumRowCount(20);
    
    combo.addItemListener(new ItemListener()
    {
      public void itemStateChanged(ItemEvent e)
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
    });
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
        findVariantAtPoint(e.getPoint());

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
        
        String[] hdTmp = new String[count + vcfReaders.length];
        System.arraycopy(header, 0, hdTmp, 0, header.length);
        header = hdTmp;
        
        for (int i = 0; i < vcfFileList.size(); i++)
          header[i+oldSize] = readHeader(vcfFileList.get(i), i+oldSize);

        setDisplay();
        repaint();
        jspView.revalidate();
      }
    });
    popup.add(addVCFMenu);
    popup.addSeparator();
    

    markNewStops.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        if(!markNewStops.isSelected())
          markAsNewStop = false;
        repaint();
      }
    });
    popup.add(markNewStops);
    
    
    final JMenuItem byQuality = new JMenuItem("Filter ...");
    byQuality.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        new VCFFilter(VCFview.this);
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
    
    final JMenu export = new JMenu("Export / View");
    popup.add(export);
    
    final JMenuItem exportVCF = new JMenuItem("Write Filtered VCF");
    exportVCF.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        VCFview.this.setCursor(new Cursor(Cursor.WAIT_CURSOR));
        IOUtils.export(entryGroup, vcfFiles, VCFview.this);
        VCFview.this.setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
      }
    });
    export.add(exportVCF);
    
    final JMenuItem exportFastaSelected = new JMenuItem("Write FASTA of selected feature(s)");
    exportFastaSelected.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        VCFview.this.setCursor(new Cursor(Cursor.WAIT_CURSOR));
        IOUtils.exportFasta(entryGroup, vcfReaders, chr, VCFview.this, vcf_v4,
            selection.getAllFeatures(), false);
        VCFview.this.setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
      }
    });
    export.add(exportFastaSelected);
    
    final JMenuItem viewFastaSelected = new JMenuItem("View FASTA of selected feature(s)");
    viewFastaSelected.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        VCFview.this.setCursor(new Cursor(Cursor.WAIT_CURSOR));
        IOUtils.exportFasta(entryGroup, vcfReaders, chr, VCFview.this, vcf_v4,
            selection.getAllFeatures(), true);
        VCFview.this.setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
      }
    });
    export.add(viewFastaSelected);
    
    final JMenuItem exportFasta = new JMenuItem("FASTA");
    exportFasta.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        VCFview.this.setCursor(new Cursor(Cursor.WAIT_CURSOR));
        IOUtils.exportFasta(entryGroup, vcfReaders, chr, VCFview.this, vcf_v4);
        VCFview.this.setCursor(new Cursor(Cursor.DEFAULT_CURSOR));
      }
    });
    export.add(exportFasta);
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
   * Test and download if on a http server
   * @param fileName
   * @return
   */
  private String testForURL(String fileName)
  {
    if(!fileName.startsWith("http:"))
      return fileName;

    String newFileName = download(fileName, null, ".bcf");
    download(fileName+".bci", newFileName, ".bci");
    return newFileName;
  }
  
  protected String download(String f, String newName, String suffix)
  {
    try
    {
      final URL urlFile = new URL(f);
      InputStream is = urlFile.openStream();

      // Create temp file.
      File bcfFile;
      if(newName != null)
        bcfFile = new File(newName+suffix);
      else
        bcfFile = File.createTempFile(urlFile.getFile().replaceAll(
        "[\\/\\s]", "_"), suffix);
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
  
  /**
   * Read the vcf header
   * @param fileName
   * @return
   */
  private String readHeader(String fileName, int index)
  {
    StringBuffer buff = new StringBuffer();
    buff.append(fileName+"\n");
    try
    {
      if(IOUtils.isBCF(fileName))
      {
        fileName = testForURL(fileName);
        vcfReaders[index] = new BCFReader(new File(fileName));
  
        String hdr = ((BCFReader)vcfReaders[index]).headerToString();
        if(hdr.indexOf("VCFv4") > -1)
          vcf_v4 = true;
        return hdr;
      }

      BlockCompressedInputStream is = 
        new BlockCompressedInputStream(new FileInputStream(fileName));
      String line;
      while( (line = TabixReader.readLine(is) ) != null )
      {
        if(!line.startsWith("##"))
          break;
        
        if(line.indexOf("VCFv4") > -1)
          vcf_v4 = true;
        
        buff.append(line+"\n");
      }
      
      vcfReaders[index] = new TabixReader(fileName);
    }
    catch (IOException e)
    {
      e.printStackTrace();
    }
    
    return buff.toString();
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
    if(mouseVCF == null)
      return null;

    String msg = 
           "Seq: "+mouseVCF.getChrom()+"\n";
    msg += "Pos: "+mouseVCF.getPos()+"\n";
    msg += "ID:  "+mouseVCF.getID()+"\n";
    msg += "Variant: "+mouseVCF.getRef()+" -> "+mouseVCF.getAlt().toString()+"\n";
    msg += "Qual: "+mouseVCF.getQuality()+"\n";
    String pl;
    if((pl = mouseVCF.getFormatValue("PL")) != null && pl.split(",").length > 1)
      msg += "Genotype likelihood (PL): "+pl+"\n";
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
        JOptionPane.showMessageDialog(this, 
            "There is a problem matching the reference sequences\n"+
            "to the names in the VCF file. This may mean the labels\n"+
            "on the reference features do not match those in the in\n"+
            "the VCF file.", 
            "Problem Found", JOptionPane.WARNING_MESSAGE);
    }
    return offsetLengths.get(refName);
  }
  
  protected void paintComponent(Graphics g)
  {
    super.paintComponent(g);
    mouseVCF = null;

    float pixPerBase = getPixPerBaseByWidth();
    int start = getBaseAtStartOfView();
    int end   = start+nbasesInView;
    
    drawSelectionRange((Graphics2D)g, pixPerBase, start, end);

    FeatureVector features = getCDSFeaturesInRange(start, end);
    for (int i = 0; i < vcfReaders.length; i++)
    {
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
            
            drawRegion(g, contigs[j], thisStart, thisEnd, i, start, pixPerBase, features); 
          }
        }

      } 
      else
      {
        int thisStart = start;
        if(thisStart < 1)
          thisStart = 1;
        drawRegion(g, chr, thisStart, end, i, start, pixPerBase, features); 
      }
    }

    if(feature_display == null)
      drawScale((Graphics2D)g, start, end, pixPerBase, getHeight());
  }
  
  private void drawRegion(Graphics g,
                          String chr,
                          int sbeg,
                          int send,
                          int i, 
                          int start, 
                          float pixPerBase, 
                          FeatureVector features) 
  {
    String s;
    
    if(vcfReaders[i] instanceof BCFReader)
    {
      try
      {
        BCFReader bcfReader = (BCFReader)vcfReaders[i];
        BCFReaderIterator it = bcfReader.query(chr, sbeg, send);
        VCFRecord bcfRecord;
        while((bcfRecord = it.next()) != null)
          drawVariantCall(g, bcfRecord, start, i, pixPerBase, features);
      }
      catch (IOException e)
      {
        logger4j.warn(e.getMessage());
        e.printStackTrace();
      }
      
    }
    else
    {
      TabixReader.Iterator iter = 
        ((TabixReader)vcfReaders[i]).query(chr+":"+sbeg+"-"+send); // get the iterator
      if (iter == null)
        return;
      try
      {
        while ((s = iter.next()) != null)
        {
          VCFRecord vcfRecord = VCFRecord.parse(s);
          drawVariantCall(g, vcfRecord, start, i, pixPerBase, features);
        }
      }
      catch (IOException e)
      {
        logger4j.warn(e.getMessage());
        e.printStackTrace();
      }
    }
  }
  
  private FeatureVector getCDSFeaturesInRange(int start, int end)
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
  

  
  protected boolean showVariant(VCFRecord record, FeatureVector features, int basePosition)
  {  
    if(!showDeletions && record.getAlt().isDeletion(vcf_v4))
      return false;
    
    if(!showInsertions && record.getAlt().isInsertion(vcf_v4))
      return false;

    if(!VCFFilter.passFilter(record))
      return false;
    
    if(!showNonOverlappings && !isOverlappingFeature(features, basePosition))
      return false;
    
    if(!showNonVariants && record.getAlt().isNonVariant())
      return false;
    
    short isSyn = -1;
    markAsNewStop = false;
    if(markNewStops.isSelected() &&
       !record.getAlt().isDeletion(vcf_v4) && 
       !record.getAlt().isInsertion(vcf_v4) && 
        record.getAlt().length() == 1 && 
        record.getRef().length() == 1)
    {
      isSyn = record.getSynFlag(features, basePosition);
      if(isSyn == 2)
        markAsNewStop = true;
    }
    
    if( (!showSynonymous || !showNonSynonymous) &&
         !record.getAlt().isDeletion(vcf_v4) && 
         !record.getAlt().isInsertion(vcf_v4) && 
         record.getAlt().length() == 1 && 
         record.getRef().length() == 1)
    {
      if(isSyn == -1)
        isSyn = record.getSynFlag(features, basePosition);
      
      if( (!showSynonymous && isSyn == 1) ||
          (!showNonSynonymous && (isSyn == 0 || isSyn == 2) ) )
        return false;
    }
    
    if(!showMultiAlleles && record.getAlt().isMultiAllele())
      return false;
    
    return true;
  }
  
  private boolean isOverlappingFeature(FeatureVector features, int basePosition)
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
  
  
  private void drawVariantCall(Graphics g, VCFRecord record, int start, int index, float pixPerBase, FeatureVector features)
  {
    int basePosition = record.getPos() + getSequenceOffset(record.getChrom());
   
    if( !showVariant(record, features, basePosition) )
      return;
    
    int pos[] = getScreenPosition(basePosition, pixPerBase, start, index);

    if(colourScheme == QUAL_COLOUR_SCHEME)
      g.setColor(getQualityColour(record));
    else if(record.getAlt().isMultiAllele())
    {
      g.setColor(Color.orange);
      g.fillArc(pos[0]-3, pos[1]-LINE_HEIGHT-3, 6, 6, 0, 360);
    }
    else if(record.getAlt().isDeletion(vcf_v4))
      g.setColor(Color.gray);
    else if(record.getAlt().isInsertion(vcf_v4))
      g.setColor(Color.yellow);
    else if(record.getAlt().length() == 1 && record.getRef().length() == 1)
    {
      g.setColor(getColourForSNP(record, features, basePosition));
      if(record.getAlt().isNonVariant())
      {
        g.drawLine(pos[0], pos[1]-2, pos[0], pos[1]-LINE_HEIGHT+4);
        return;
      }
    }
    else
      g.setColor(Color.pink);

    if(markAsNewStop)
      g.fillArc(pos[0]-3, pos[1]-(LINE_HEIGHT/2)-3, 6, 6, 0, 360);
    
    g.drawLine(pos[0], pos[1], pos[0], pos[1]-LINE_HEIGHT);
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
      short synFlag = record.getSynFlag(features, basePosition);
      if(synFlag == 1)
        return Color.red;
      else if(synFlag == 0 || synFlag == 2)
        return Color.blue;
      else
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
      return Color.magenta; // non-variant
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
    pos[1] = getHeight() - 15 - (vcfFileIndex*(LINE_HEIGHT+5)); 
    return pos;
  }
  
  private void findVariantAtPoint(Point mousePoint)
  {
    float pixPerBase = getPixPerBaseByWidth();
    int start = getBaseAtStartOfView();
    int end   = start+nbasesInView;
    FeatureVector features = getCDSFeaturesInRange(start, end);
    
    for (int i = 0; i < vcfReaders.length; i++)
    {
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
            searchRegion(contigs[j], thisStart, thisEnd, i, mousePoint, features, start, pixPerBase);
          }
        }
      } 
      else
      {
        int thisStart = start;
        if(thisStart < 1)
          thisStart = 1;
        searchRegion(chr, thisStart, end, i, mousePoint, features, start, pixPerBase);
      }
    }
  }
  
  private void searchRegion(String chr, int sbeg, int send, int i, 
                            Point mousePoint, FeatureVector features,
                            int start, float pixPerBase) 
  {
    if(vcfReaders[i] instanceof BCFReader)
    {
      try
      {
        BCFReader bcfReader = (BCFReader)vcfReaders[i];
        BCFReaderIterator it = bcfReader.query(chr, sbeg, send);
        VCFRecord bcfRecord;
        while((bcfRecord = it.next()) != null)
          isMouseOver(mousePoint, bcfRecord, features, i, start, pixPerBase);
      }
      catch (IOException e)
      {
        // TODO Auto-generated catch block
        e.printStackTrace();
      }
    }
    else
    {
      TabixReader.Iterator iter = 
        ((TabixReader)vcfReaders[i]).query(chr+":"+sbeg+"-"+send); // get the iterator
      if (iter == null)
        return;
      try
      { 
        String s;
        while ((s = iter.next()) != null)
        {
          VCFRecord vcfRecord = VCFRecord.parse(s);
          isMouseOver(mousePoint, vcfRecord, features, i, start, pixPerBase);
        }
      }
      catch (IOException e)
      {
        e.printStackTrace();
      }
    }
  }
  
  private void isMouseOver(Point mousePoint, 
                           VCFRecord record, 
                           FeatureVector features, 
                           int i, 
                           int start, float pixPerBase)
  {
    int basePosition = record.getPos() + getSequenceOffset(record.getChrom());

    if( !showVariant(record, features, basePosition) )
      return;
    
    int pos[] = getScreenPosition(basePosition, pixPerBase, start, i);
    if(mousePoint != null &&
       mousePoint.getY() < pos[1] &&
       mousePoint.getY() > pos[1]-LINE_HEIGHT &&
       mousePoint.getX() > pos[0]-3 &&
       mousePoint.getX() < pos[0]+3)
     {
       mouseVCF = record;
       mouseOverIndex = i;
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
  
  /**
   * Write an Artemis tab file
   * @param vcfFiles
   * @param seq
   */
  private static void writeGffFiles(final List<String> vcfFiles,
                                    final String seq)
  {
    try
    {         
      TabixReader tr[] = new TabixReader[vcfFiles.size()];      
      String s;
      int step = 50000;
      int regions[] = null;
      String chr = seq.split(":")[0];
      
      for(int i = 0; i < tr.length; i++)
      {
        tr[i] = new TabixReader(vcfFiles.get(i));
        
        System.out.println(vcfFiles.get(i)+".tab");
        FileWriter writer = new FileWriter(new File(vcfFiles.get(i)+".tab"));
        
        if(regions == null)
          regions = tr[i].parseReg(seq);
        int start = regions[1];
        int end = regions[2];
        
        for(int j = start; j < end + step; j += step)
        {
          String region = chr+":"+j+"-"+(j+step);
          
          TabixReader.Iterator iter = tr[i].query(region); // get the iterator
          if (iter == null)
            continue;
          try
          {
            
            while ((s = iter.next()) != null)
            {
              String parts[] = tabPattern.split(s, 0);;

              Location location = new Location(parts[1]+".."+parts[1]);
              QualifierVector qualifiers = new QualifierVector();           
              Qualifier note = new Qualifier("note", parts[3]+" "+parts[4]);
              qualifiers.addQualifierValues(note);
              Qualifier score = new Qualifier("score", parts[5]);
              qualifiers.addQualifierValues(score);
              
              try
              {
                EmblStreamFeature emblFeature = new EmblStreamFeature(
                    new Key("SNP"), location, qualifiers);
                emblFeature.writeToStream(writer);
              }
              catch (InvalidRelationException e)
              {
                e.printStackTrace();
              }
            }
          }
          catch (IOException e)
          {
            e.printStackTrace();
          }
        }
        writer.close();
      }
    }
    catch (IOException e)
    {
      // TODO Auto-generated catch block
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

  private float getPixPerBaseByWidth()
  {
    return (float)vcfPanel.getWidth() / (float)nbasesInView;
  }
  
  protected EntryGroup getEntryGroup()
  {
    return entryGroup;
  }
  
  /**
   * Popup menu listener
   */
   class PopupListener extends MouseAdapter
   {
     JMenuItem showDetails;
     
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
           popup.remove(showDetails);
         
         if( mouseVCF != null )
         {
           showDetails = new JMenuItem("Show details of : "+
               mouseVCF.getChrom()+":"+mouseVCF.getPos()+" "+mouseVCF.getID());
           showDetails.addActionListener(new ActionListener()
           {
             public void actionPerformed(ActionEvent e) 
             {
               FileViewer viewDetail = new FileViewer(
                   mouseVCF.getChrom()+":"+mouseVCF.getPos()+" "+mouseVCF.getID(), true, false);
               
               viewDetail.appendString(header[mouseOverIndex]+"\n", Level.INFO);
               
               viewDetail.appendString("Seq   : "+mouseVCF.getChrom()+"\n", Level.DEBUG);
               viewDetail.appendString("Pos   : "+mouseVCF.getPos()+"\n", Level.DEBUG);
               viewDetail.appendString("ID    : "+mouseVCF.getID()+"\n", Level.DEBUG);
               viewDetail.appendString("Ref   : "+mouseVCF.getRef()+"\n", Level.DEBUG);
               viewDetail.appendString("Alt   : "+mouseVCF.getAlt().toString()+"\n", Level.DEBUG);
               viewDetail.appendString("Qual  : "+mouseVCF.getQuality()+"\n", Level.DEBUG);
               viewDetail.appendString("Filter: "+mouseVCF.getFilter()+"\n", Level.DEBUG);
               viewDetail.appendString("Info  : "+mouseVCF.getInfo()+"\n", Level.DEBUG);
               
               if(mouseVCF.getFormat() != null)
               {
                 viewDetail.appendString("\nGenotype information:\n", Level.INFO);
                 viewDetail.appendString("Format: "+mouseVCF.getFormat()+"\n", Level.DEBUG);
                 viewDetail.appendString(mouseVCF.getSampleDataString()+"\n", Level.DEBUG);
               }
               
               viewDetail.getTextPane().setCaretPosition(0);
             }  
           });
           popup.add(showDetails);
         }
         popup.show(e.getComponent(),
             e.getX(), e.getY());
       }
     }
   }
   
   private void setDisplay()
   {
     Dimension d = new Dimension();
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
    boolean writeTab = false;
    String seq = null;
    
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
      else if(args[i].equals("-t"))
      {
        writeTab = true;
        seq = args[++i];
      }
      else if(args[i].startsWith("-h"))
      { 
        System.out.println("-h\t show help");
        
        System.out.println("-f\t VCF file to display");
        System.out.println("-r\t reference file (optional)");
        System.out.println("-v\t number of bases to display in the view (optional)");
        System.out.println("-t\t chr:start-end - this writes out the given region");

        System.exit(0);
      }
    }
    
    if(writeTab)
      writeGffFiles(vcfFileList, seq);
    else
    {
      JFrame f = new JFrame();
      new VCFview(f, (JPanel) f.getContentPane(), vcfFileList, 
          nbasesInView, 100000000, null, reference, null);
    }
  }
}