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
import java.io.FileWriter;
import java.io.IOException;
import java.util.List;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.swing.JButton;
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
import uk.ac.sanger.artemis.sequence.AminoAcidSequence;
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
  private TabixReader tr[];
  private String header[];
  private FeatureDisplay feature_display;
  private Selection selection;
  private int nbasesInView;
  private int seqLength;
  private EntryGroup entryGroup;
  private String chr;
  private String mouseOverVCFline;
  private int mouseOverIndex = -1;
//record of where a mouse drag starts
  private int dragStart = -1;
  private JPopupMenu popup;
  private int LINE_HEIGHT = 15;
  
  private boolean showSynonymous = true;
  private boolean showNonSynonymous = true;
  private boolean showDeletions = true;
  private boolean showInsertions = true;
  private boolean showMultiAlleles = true;
  private float MIN_QUALITY = -10;
  
  private Pattern multiAllelePattern = Pattern.compile("^[AGCT],[AGCT,]+$");

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
      tr = new TabixReader[vcfFiles.size()];
      header = new String[vcfFiles.size()];
      
      for(int i=0; i<vcfFiles.size(); i++)
      {
        header[i] = readHeader(vcfFiles.get(i));
        tr[i] = new TabixReader(vcfFiles.get(i));
      }
    }
    catch (IOException e)
    {
      e.printStackTrace();
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
    final MouseMotionListener mouseMotionListener = new MouseMotionListener()
    {
      public void mouseDragged(MouseEvent event)
      {
        handleCanvasMouseDrag(event);
      }
      
      public void mouseMoved(MouseEvent e)
      {
        findVariantAtPoint(e.getPoint());
      }
    };
    addMouseMotionListener(mouseMotionListener);
    addMouseListener(new PopupListener());
    
    //
    createMenus(frame);
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
  
  private void createMenus(JFrame frame)
  {
    JComponent topPanel;
    if(feature_display != null)
      topPanel = new JPanel(new FlowLayout(FlowLayout.LEADING, 0, 0));
    else
    { 
      topPanel = new JMenuBar();
      frame.setJMenuBar((JMenuBar)topPanel);
    }
    
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
    
    final JComboBox combo = new JComboBox(tr[0].getmSeq());
    if(chr == null)
      this.chr = tr[0].getmSeq()[0];
    combo.setSelectedItem(this.chr);
    combo.setEditable(false);
    combo.setMaximumRowCount(20);
    
    combo.addItemListener(new ItemListener()
    {
      public void itemStateChanged(ItemEvent e)
      {
        VCFview.this.chr = (String) combo.getSelectedItem();
        repaint();
      }
    });
    topPanel.add(combo);
    
    // popup menu
    popup = new JPopupMenu();
    final JCheckBoxMenuItem showSyn = new JCheckBoxMenuItem(
        "Show synonymous", showSynonymous);
    showSyn.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        showSynonymous = showSyn.isSelected();
        repaint();
      }
    });
    popup.add(showSyn);
    
    final JCheckBoxMenuItem showNonSyn = new JCheckBoxMenuItem(
        "Show non-synonymous", showNonSynonymous);
    showNonSyn.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        showNonSynonymous = showNonSyn.isSelected();
        repaint();
      }
    });
    popup.add(showNonSyn);
    
    
    final JCheckBoxMenuItem showDeletionsMenu = new JCheckBoxMenuItem(
        "Show deletions", showDeletions);
    showDeletionsMenu.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        showDeletions = showDeletionsMenu.isSelected();
        repaint();
      }
    });
    popup.add(showDeletionsMenu);
    
    final JCheckBoxMenuItem showInsertionsMenu = new JCheckBoxMenuItem(
        "Show insertions", showInsertions);
    showInsertionsMenu.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        showInsertions = showInsertionsMenu.isSelected();
        repaint();
      }
    });
    popup.add(showInsertionsMenu);
    
    final JCheckBoxMenuItem showMultiAllelesMenu = new JCheckBoxMenuItem(
        "Show multiple alleles", showMultiAlleles);
    showMultiAllelesMenu.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        showMultiAlleles = showMultiAllelesMenu.isSelected();
        repaint();
      }
    });
    popup.add(showMultiAllelesMenu);
    popup.addSeparator();
    
    final JMenuItem filterByQuality = new JMenuItem("Filter by quality");
    filterByQuality.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        //
        String inputValue = JOptionPane.showInputDialog(null, 
            "Enter a minimum quality score:", MIN_QUALITY);
        
        if(inputValue == null)
          return;
        try
        {
          MIN_QUALITY = Float.parseFloat(inputValue);
          repaint();
        }
        catch(NumberFormatException ex)
        {
          JOptionPane.showMessageDialog(null, 
              "Number "+inputValue+" not recognised.", 
              "Format Error", JOptionPane.ERROR_MESSAGE);
        }
      }
    });
    popup.add(filterByQuality);
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
  private String readHeader(final String fileName)
  {
    StringBuffer buff = new StringBuffer();
    buff.append(fileName+"\n");
    try
    {
      FileInputStream fileStream = new FileInputStream(fileName);
      BlockCompressedInputStream inputStrean = new BlockCompressedInputStream(fileStream);
      
      String line;
      while( (line = TabixReader.readLine(inputStrean) ) != null )
      {
        if(!line.startsWith("##"))
          break;
        buff.append(line+"\n");
      }
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
    if(mouseOverVCFline == null)
      return null;
    
    String parts[] = mouseOverVCFline.split("\t");
    String msg = 
           "Seq: "+parts[0]+"\n";
    msg += "Pos: "+parts[1]+"\n";
    msg += "ID:  "+parts[2]+"\n";
    msg += "Variant: "+parts[3]+" -> "+parts[4]+"\n";
    msg += "Qual: "+parts[5]+"\n";
    
    return msg;
  }
  
  protected void paintComponent(Graphics g)
  {
    super.paintComponent(g);
    mouseOverVCFline = null;

    float pixPerBase = getPixPerBaseByWidth();
    String s;
    
    int start = getBaseAtStartOfView();
    int end   = start+nbasesInView;
    String region = chr+":"+start+"-"+end;
       
    drawSelectionRange((Graphics2D)g, pixPerBase, start, end);
 
    FeatureVector features = getCDSFeaturesInRange(start, end);
    
    // a region is specified; random access
    for (int i = 0; i < tr.length; i++)
    {
      TabixReader.Iterator iter = tr[i].query(region); // get the iterator
      if (iter == null)
        return;
      try
      {
        while ((s = iter.next()) != null)
        {
          drawVariantCall(g, s, start, i, pixPerBase, features);
        }
      }
      catch (IOException e)
      {
        // TODO Auto-generated catch block
        e.printStackTrace();
      }
    }

    if(feature_display == null)
      drawScale((Graphics2D)g, start, end, pixPerBase, getHeight());
  }
  
  private FeatureVector getCDSFeaturesInRange(int start, int end)
  {
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
  
  /**
   * Is this a deletion type.
   * @param variant
   * @return
   */
  private boolean isDeletion(String variant)
  {
    if(variant.indexOf("D")>-1)
      return true;
    return false;
  }
  
  /**
   * Is this an insertion type.
   * @param variant
   * @return
   */
  private boolean isInsertion(String variant)
  {
    if(variant.indexOf("I")>-1)
      return true;
    return false;
  }
  
  private boolean showVariant(String variant, FeatureVector features, int basePosition, float quality)
  {
    if(quality < MIN_QUALITY)
      return false;
    
    if(!showDeletions && isDeletion(variant))
      return false;
    
    if(!showInsertions && isInsertion(variant))
      return false;
    
    if( (!showSynonymous || !showNonSynonymous) &&
         !isDeletion(variant) && !isInsertion(variant) && variant.length() == 1)
    {
      boolean isSyn = isSynonymous(features, basePosition, variant.toLowerCase().charAt(0));
      
      if( (!showSynonymous && isSyn) ||
          (!showNonSynonymous && !isSyn) )
        return false;
    }
    
    if(!showMultiAlleles && multiAllelePattern.matcher(variant).matches())
      return false;
    
    return true;
  }
  
  private void drawVariantCall(Graphics g, String line, int start, int index, float pixPerBase, FeatureVector features)
  {
    String parts[] = line.split("\\t");
    int basePosition = Integer.parseInt(parts[1]);
    float quality = Float.parseFloat(parts[5]);
    
    if( !showVariant(parts[4], features, basePosition, quality) )
      return;
    
    int pos[] = getScreenPosition(basePosition, pixPerBase, start, index);
    
    if(parts[4].equals("C"))
      g.setColor(Color.red);
    else if(parts[4].equals("A"))
      g.setColor(Color.green);
    else if(parts[4].equals("G"))
      g.setColor(Color.blue);
    else if(parts[4].equals("T"))
      g.setColor(Color.black);
    else if(isDeletion(parts[4]))
      g.setColor(Color.gray);
    else if(isInsertion(parts[4]))
      g.setColor(Color.yellow);
    else
    {
      Matcher m = multiAllelePattern.matcher(parts[4]);
      if(m.matches())
      {
        g.setColor(Color.orange);
        g.fillArc(pos[0]-3, pos[1]-LINE_HEIGHT-3, 6, 6, 0, 360);
      }
      else
        g.setColor(Color.pink);
    }

    g.drawLine(pos[0], pos[1], pos[0], pos[1]-LINE_HEIGHT);
  }
  
  private boolean isSynonymous(FeatureVector features, int basePosition, char variant)
  {
    int intronlength = 0;
    Range lastRange = null;
    
    for(int i = 0; i<features.size(); i++)
    {
      Feature feature = features.elementAt(i);
      
      if(feature.getRawFirstBase() < basePosition && feature.getRawLastBase() > basePosition)
      {
        RangeVector ranges = feature.getLocation().getRanges();
        intronlength = 0;

        for(int j=0; j< ranges.size(); j++)
        {
          Range range = (Range) ranges.get(j);
          
          if(j > 0)
          {
            if(feature.isForwardFeature())
              intronlength+=range.getStart()-lastRange.getEnd()-1;
            else
              intronlength+=lastRange.getStart()-range.getEnd()-1;
            
            if(intronlength < 0)
              intronlength = 0;
          }
          
          if(range.getStart() < basePosition && range.getEnd() > basePosition)
          {
            int mod;
            int codonStart;
            
            if(feature.isForwardFeature())
            {
              mod = (basePosition-feature.getRawFirstBase())%3;
              codonStart = basePosition-feature.getRawFirstBase()-mod;
            }
            else
            {
              mod = (feature.getRawLastBase()-basePosition)%3;
              codonStart = feature.getRawLastBase()-basePosition-mod;
            }

            codonStart-=intronlength;
            
            try
            {
              char codon[] = feature.getBases().substring(codonStart,
                  codonStart + 3).toLowerCase().toCharArray();

              // String oldBase = new String(codon);
              char aaRef = AminoAcidSequence.getCodonTranslation(codon[0],
                  codon[1], codon[2]);

              codon[mod] = variant;
              char aaNew = AminoAcidSequence.getCodonTranslation(codon[0],
                  codon[1], codon[2]);

              if (aaNew == aaRef)
                return true;
              else
                return false;
            }
            catch(Exception e)
            {
              for(int k=0; k<ranges.size(); k++)
                System.out.println(k+" "+ ((Range)ranges.get(k)).getStart() );
              
              System.out.println(feature.getIDString()+"  "+codonStart+" "+intronlength+" basePosition="+basePosition+" segment="+range.getStart()+".."+range.getEnd()+" mod="+mod);
              throw new RuntimeException(e);
            }
          }

          lastRange = range;
        }
      }
    }
    
    return false;
  }
  
  private int[] getScreenPosition(int base, float pixPerBase, int start, int vcfFileIndex)
  {
    int pos[] = new int[2];
    pos[0] = Math.round((base - start)*pixPerBase);
    pos[1] = getHeight() - 15 - (vcfFileIndex*(LINE_HEIGHT+5)); 
    return pos;
  }
  
  private void findVariantAtPoint(Point lastMousePoint)
  {
    float pixPerBase = getPixPerBaseByWidth();
    String s;
    int start = getBaseAtStartOfView();
    int end   = start+nbasesInView;
    String region = chr+":"+start+"-"+end;
    FeatureVector features = getCDSFeaturesInRange(start, end);
    
    for (int i = 0; i < tr.length; i++)
    {
      TabixReader.Iterator iter = tr[i].query(region); // get the iterator
      if (iter == null)
        return;
      try
      {
        while ((s = iter.next()) != null)
        {
          String parts[] = s.split("\\t");
          int basePosition = Integer.parseInt(parts[1]);
          float quality = Float.parseFloat(parts[5]);
          
          if( !showVariant(parts[4], features, basePosition, quality) )
            continue;
          
          int pos[] = getScreenPosition(basePosition, pixPerBase, start, i);
          if(lastMousePoint != null &&
             lastMousePoint.getY() < pos[1] &&
             lastMousePoint.getY() > pos[1]-LINE_HEIGHT &&
             lastMousePoint.getX() > pos[0]-3 &&
             lastMousePoint.getX() < pos[0]+3)
           {
             mouseOverVCFline = s;
             mouseOverIndex = i;
           }
        }
      }
      catch (IOException e)
      {
        e.printStackTrace();
      }
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
              String parts[] = s.split("\\t");

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
         
         if( mouseOverVCFline != null )
         {
           final String parts[] = mouseOverVCFline.split("\\t");
      
           showDetails = new JMenuItem("Show details of : "+parts[0]+":"+parts[1]+" "+parts[2]);
           showDetails.addActionListener(new ActionListener()
           {
             public void actionPerformed(ActionEvent e) 
             {
               FileViewer viewDetail = new FileViewer(parts[0]+":"+parts[1]+" "+parts[2], true, false);
               
               viewDetail.appendString(header[mouseOverIndex]+"\n", Level.INFO);
               
               viewDetail.appendString("Seq   : "+parts[0]+"\n", Level.DEBUG);
               viewDetail.appendString("Pos   : "+parts[1]+"\n", Level.DEBUG);
               viewDetail.appendString("ID    : "+parts[2]+"\n", Level.DEBUG);
               viewDetail.appendString("Ref   : "+parts[3]+"\n", Level.DEBUG);
               viewDetail.appendString("Alt   : "+parts[4]+"\n", Level.DEBUG);
               viewDetail.appendString("Qual  : "+parts[5]+"\n", Level.DEBUG);
               viewDetail.appendString("Filter: "+parts[6]+"\n", Level.DEBUG);
               viewDetail.appendString("Info  : "+parts[7]+"\n", Level.DEBUG);
               
               if(parts.length > 8)
               {
                 viewDetail.appendString("\nGenotype information:\n", Level.INFO);
                 viewDetail.appendString("Format: "+parts[8]+"\n", Level.DEBUG);
                 for(int i=9; i<parts.length; i++)
                   viewDetail.appendString(parts[i]+"\n", Level.DEBUG);
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
     d.setSize(nbasesInView*getPixPerBaseByWidth(), (tr.length+1)*(LINE_HEIGHT+5));
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
      vcfFileList = fileSelection.getFiles(".vcf");
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