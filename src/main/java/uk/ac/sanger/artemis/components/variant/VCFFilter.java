/*
 * created: 2011
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2011  Genome Research Limited
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
import java.awt.Dimension;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.regex.Pattern;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;

import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.circular.TextFieldFloat;
import uk.ac.sanger.artemis.circular.TextFieldInt;
import uk.ac.sanger.artemis.components.Utilities;

public class VCFFilter extends JFrame
{
  private static final long serialVersionUID = 1L;

  //private static boolean useHeader = true;
/*  private static int MIN_DP = 0;
  private static float MIN_MQ = 0;
  private static float MIN_AF1 = 0;
  private static float MAX_CI95 = 10;
  private static Pattern COMMA_PATTERN = Pattern.compile(",");*/
  private static Pattern SEMICOLON_PATTERN = Pattern.compile(";");
  
  private static Pattern HOMOZYGOUS_PATTERN = Pattern.compile("^0(/0)*+|(\\|0)*+$");
  private FilteredPanel filterPanel;
  protected static boolean manualFilter = false; // show manual filtering

  /**
   * Filter VCF records by the variant type and/or by different values in 
   * the record, QUAL, DP, MQ and AF1.
   * @param vcfView
   */
  public VCFFilter(final VCFview vcfView)
  {
    super("Variant Filter");

    final JPanel mainPanel = (JPanel)getContentPane();
    final JPanel btmPanel = new JPanel(new GridBagLayout());
    final JTabbedPane tabPane = new JTabbedPane();

    final List<HeaderLine> filterHdrLines = vcfView.getVcfReaders()[0].getFILTER();
    
    filterPanel = new FilteredPanel(filterHdrLines);
    JScrollPane ftrScroll = new JScrollPane(filterPanel);
    ftrScroll.setPreferredSize( new Dimension(ftrScroll.getPreferredSize().width, 150) );
    
    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
    tabPane.setPreferredSize(new Dimension(screen.width*9/22, screen.height*2/5));
    mainPanel.add(tabPane, BorderLayout.NORTH);
    mainPanel.add(ftrScroll, BorderLayout.CENTER);
    mainPanel.add(btmPanel, BorderLayout.SOUTH);

    ////
    ////
    //// Filter by type
    final Box typeBox = Box.createVerticalBox();
    tabPane.addTab("Type", typeBox);
    JLabel typeLabel = new JLabel("VARIANT TYPE");
    typeLabel.setFont(typeLabel.getFont().deriveFont(Font.BOLD));
    typeBox.add(typeLabel);

    final JCheckBox showSyn = new JCheckBox("Synonymous", vcfView.showSynonymous);
    typeBox.add(showSyn);
    showSyn.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        setFlagFilter(HeaderLine.FILTER_SYN, "SYN", showSyn.getText(), showSyn.isSelected());
        vcfView.showSynonymous = showSyn.isSelected();
        vcfView.repaint();
      }
    });

    final JCheckBox showNonSyn = new JCheckBox("Non-synonymous", vcfView.showNonSynonymous);
    typeBox.add(showNonSyn);
    showNonSyn.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        setFlagFilter(HeaderLine.FILTER_NONSYN, "NONSYN", showNonSyn.getText(), showNonSyn.isSelected());
        vcfView.showNonSynonymous = showNonSyn.isSelected();
        vcfView.repaint();
      }
    });

    final JCheckBox showDeletionsMenu = new JCheckBox("Deletions", vcfView.showDeletions);
    typeBox.add(showDeletionsMenu);
    showDeletionsMenu.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        setFlagFilter(HeaderLine.FILTER_DEL, "DEL", showDeletionsMenu.getText(), showDeletionsMenu.isSelected());
        vcfView.showDeletions = showDeletionsMenu.isSelected();
        vcfView.repaint();
      }
    });

    final JCheckBox showInsertionsMenu = new JCheckBox("Insertions", vcfView.showInsertions);
    typeBox.add(showInsertionsMenu);
    showInsertionsMenu.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        setFlagFilter(HeaderLine.FILTER_INS, "INS", showInsertionsMenu.getText(), showInsertionsMenu.isSelected());
        vcfView.showInsertions = showInsertionsMenu.isSelected();
        vcfView.repaint();
      }
    });

    final JCheckBox showMultiAllelesMenu = new JCheckBox("Multiple alleles", vcfView.showMultiAlleles);
    typeBox.add(showMultiAllelesMenu);
    showMultiAllelesMenu.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        setFlagFilter(HeaderLine.FILTER_MULTALL_FLAG, "MULTI_ALLELLES", "Multiple alleles", showMultiAllelesMenu.isSelected());
        vcfView.showMultiAlleles = showMultiAllelesMenu.isSelected();
        vcfView.repaint();
      }
    });
    
    final JCheckBox showNonOverlappingsMenu = new JCheckBox("Variants not overlapping CDS", vcfView.showNonOverlappings);
    typeBox.add(showNonOverlappingsMenu);
    showNonOverlappingsMenu.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        setFlagFilter(HeaderLine.FILTER_OVERLAP_FLAG, "NO-OVERLAP", showNonOverlappingsMenu.getText(), showNonOverlappingsMenu.isSelected());
        vcfView.showNonOverlappings = showNonOverlappingsMenu.isSelected();
        vcfView.repaint();
      }
    });
    
    final JCheckBox showNonVariantMenu = new JCheckBox("Non-Variants", vcfView.showNonVariants);
    typeBox.add(showNonVariantMenu);
    showNonVariantMenu.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        setFlagFilter(HeaderLine.FILTER_NV_FLAG, "NV", showNonVariantMenu.getText(), showNonVariantMenu.isSelected());
        vcfView.showNonVariants = showNonVariantMenu.isSelected();
        vcfView.repaint();
      }
    });
    setFlagFilter(HeaderLine.FILTER_NV_FLAG, "NV", showNonVariantMenu.getText(), showNonVariantMenu.isSelected());
    
    typeBox.add(Box.createVerticalStrut(20));
    final JCheckBox manualMenu = new JCheckBox("Manual Annotation", manualFilter);
    typeBox.add(manualMenu);
    manualMenu.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        manualFilter = manualMenu.isSelected();
        setFlagFilter(HeaderLine.FILTER_MANUAL, "MANUAL", manualMenu.getText(), manualMenu.isSelected());
        vcfView.repaint();
      }
    });
    setFlagFilter(HeaderLine.FILTER_MANUAL, "MANUAL", manualMenu.getText(), manualMenu.isSelected());

    typeBox.add(Box.createVerticalGlue());
    
    if(vcfView.getEntryGroup() == null || vcfView.getEntryGroup().getAllFeaturesCount() == 0)
    {
      showSyn.setEnabled(false);
      showNonSyn.setEnabled(false);
      showNonOverlappingsMenu.setEnabled(false);
    }
    
    ////
    ////
    //// Filter by property
    final GridBagConstraints c = new GridBagConstraints();
    c.gridx = 0;
    c.gridy = 0;
    c.anchor = GridBagConstraints.NORTHWEST;
    
    final List<HeaderLine> info = vcfView.getVcfReaders()[0].getINFO();

    final JPanel propPanel = new JPanel(new GridBagLayout());
    final JScrollPane jspProp = new JScrollPane(propPanel);
    propPanel.setBorder(BorderFactory.createLineBorder(Color.gray));

    tabPane.addTab("Info", jspProp);
    c.gridy += 1;
    propPanel.add(new JLabel(" "), c);
    c.gridy += 1;
    JLabel propLabel = new JLabel("VARIANT PROPERTY:");
    propLabel.setFont(propLabel.getFont().deriveFont(Font.BOLD));
    propPanel.add(propLabel, c);
    
    // min quality
    c.gridy += 1;
    c.gridx = 1;
    propPanel.add(new JLabel("MIN"), c);
    c.gridx += 1;
    propPanel.add(new JLabel("MAX"), c);;
    
    c.gridy += 1;
    c.gridx = 0;
    propPanel.add(new JLabel("Quality score (QUAL):"), c);
    final TextFieldFloat minQuality = new TextFieldFloat();
    minQuality.setColumns(8);
    c.gridx = 1;
    propPanel.add(minQuality, c); 
    final TextFieldFloat maxQuality = new TextFieldFloat();
    maxQuality.setColumns(8);
    c.gridx += 1;
    propPanel.add(maxQuality, c);
    
    String hdrLine = "##FILTER=<ID=QUAL,Type=Float,Number=1,Description=\"QUAL "
        + (minQuality.getText().equals("") ? "" : " < " + minQuality.getValue() + " ")
        + (maxQuality.getText().equals("") ? "" : " > " + maxQuality.getValue() + " ") + "\">";
    HeaderLine hLine = new HeaderLine(hdrLine, "FILTER_QUAL",
        AbstractVCFReader.getLineHash("FILTER", hdrLine));
    
    maxQuality.addKeyListener(new FilterListener(hLine, false, 0, 1));
    minQuality.addKeyListener(new FilterListener(hLine, true, 0, 1));
    

    ////
    ////
/*    final JTextField minDP = new JTextField(Integer.toString(MIN_DP), 8);
    final JTextField minMQ = new JTextField(Float.toString(MIN_MQ),8);
    final JTextField minAF1 = new JTextField(Float.toString(MIN_AF1),8);
    final JTextField maxCI95 = new JTextField(Float.toString(MAX_CI95),8);*/
    
    if(info.size() == 0)
    {
      //useHeader = false;
/*      // min DP
      c.gridy = c.gridy+1;
      c.gridx = 0;
      propPanel.add(new JLabel("Minimum combined depth across samples (DP):"), c);
      
      c.gridx = 1;
      propPanel.add(minDP, c);
      
      // min MQ
      c.gridy = c.gridy+1;
      c.gridx = 0;
      propPanel.add(new JLabel("Minimum RMS mapping quality (MQ):"), c);
      
      c.gridx = 1;
      propPanel.add(minMQ, c);
      
      // min AF1
      c.gridy = c.gridy+1;
      c.gridx = 0;
      propPanel.add(new JLabel("Minimum site frequency of strongest non-reference allele (AF1):"), c);
      
      c.gridx = 1;
      propPanel.add(minAF1, c);
      
      // max CI95
      c.gridy = c.gridy+1;
      c.gridx = 0;
      propPanel.add(new JLabel("Maximum 95% confidence interval variation from AF (CI95):"), c);

      c.gridx = 1;
      propPanel.add(maxCI95, c);*/
    }
    else
    {
      //
      createPanel(vcfView, propPanel, c, info, "INFO FIELDS:");
    }

    ////
    ////
    //// FORMAT - Genotype Filter
    final List<HeaderLine> format = vcfView.getVcfReaders()[0].getFORMAT();

    final JPanel formatPanel = new JPanel(new GridBagLayout());
    formatPanel.setBorder(BorderFactory.createLineBorder(Color.gray));

    tabPane.addTab("Genotype", formatPanel);
    createPanel(vcfView, formatPanel, c, format, "GENOTYPE FIELDS:");
    
    //
    c.gridy = c.gridy+1;
    c.gridx = 0;
    c.anchor = GridBagConstraints.EAST;
    JButton apply = new JButton("APPLY");
    apply.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        try
        {
/*          if(info.size() == 0)
          {
            MIN_DP = Integer.parseInt(minDP.getText());
            MIN_MQ = Float.parseFloat(minMQ.getText());
            MIN_AF1 = Float.parseFloat(minAF1.getText());
            MAX_CI95 = Float.parseFloat(maxCI95.getText());
          }*/
          vcfView.repaint();
        }
        catch(NumberFormatException ex)
        {
          JOptionPane.showMessageDialog(null, 
              ex.getMessage(), 
              "Format Error", JOptionPane.ERROR_MESSAGE);
        }
      }
    });
    btmPanel.add(apply, c);

    c.gridx = 1;
    c.anchor = GridBagConstraints.WEST;
    JButton close = new JButton("OK");
    close.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        vcfView.repaint();
        setVisible(false);
      }
    });
    btmPanel.add(close, c);
    
    pack();
    Utilities.centreFrame(this);
  }

  private void setFlagFilter(final int hdrIndex, final String id, final String descr, final boolean isSelected)
  {
    final String ID = id+"_FLAG";
    if(!isSelected) 
    {
      String hdrLine = "##FILTER=<ID="+id+",Type=Flag,Number=0,Description=\""+descr+"\">";
      HeaderLine hLine = new HeaderLine(hdrLine, HeaderLine.filterFlagStr[hdrIndex],
          AbstractVCFReader.getLineHash("FILTER", hdrLine));
      filterPanel.addFilter(ID, hLine, 0);
    }
    else
      filterPanel.removeFilter(ID);
    filterPanel.updateFilters();  
  }
  
  private void createPanel(final VCFview vcfView,
                           final JPanel panel, 
                           final GridBagConstraints c, 
                           final List<HeaderLine> headerLineList, 
                           final String name)
  {
    c.gridx = 0;
    c.gridy += 1;
    c.weightx = 1.d;
    c.weighty = 1.d;

    panel.add(new JLabel(" "), c);
    c.gridy += 1;
    JLabel nameLabel = new JLabel(name);
    nameLabel.setFont(panel.getFont().deriveFont(Font.BOLD));
    panel.add(nameLabel, c);
    
    c.gridy += 1;
    panel.add(new JLabel(" "), c);
    c.gridy += 1;

    for (int i = 0; i < headerLineList.size(); i++)
    {
      final HeaderLine hLine = headerLineList.get(i);
      final String type = hLine.getType();

      if (type != null && type.equals("String"))
        continue;

      int num = hLine.getNumber();

      c.gridy += 1;
      c.gridx = 0;
      for (int j = 0; j < num; j++)
      {
        c.gridx += 1;
        panel.add(new JLabel("MIN"), c);
        c.gridx += 1;
        panel.add(new JLabel("MAX"), c);
        c.gridx += 1;
      }
      
      c.gridy += 1;
      c.gridx = 0;
      
      JLabel lab = new JLabel(hLine.getID() + " (" + type + ")");
      lab.setToolTipText(hLine.getDescription());
      panel.add(lab, c);

      if (num == 0 && type.equals("Flag"))
      {
        final JCheckBox flag = new JCheckBox();
        flag.addActionListener(new ActionListener()
        {
          public void actionPerformed(ActionEvent e)
          {
            final String ID = hLine.getHeaderTypeStr()+":"+hLine.getID();
            if(flag.isSelected()) 
              filterPanel.addFilter(ID, hLine, 0);
            else
              filterPanel.removeFilter(ID);
            filterPanel.updateFilters();
          }
        });
        c.gridx += 1;
        panel.add(flag, c);
      }

      for (int j = 0; j < num; j++)
      {
        if (type != null && type.equals("Integer"))
        {
          final TextFieldInt min = new TextFieldInt();
          min.setColumns(8);
          final TextFieldInt max = new TextFieldInt();
          max.setColumns(8);
          c.gridx += 1;
          panel.add(min, c);
          c.gridx += 1;
          panel.add(max, c);
          
          max.addKeyListener(new FilterListener(hLine, false, j, num));
          min.addKeyListener(new FilterListener(hLine, true, j, num));
        }
        else if (type != null && type.equals("Float"))
        {
          TextFieldFloat min = new TextFieldFloat();
          min.setColumns(8);
          TextFieldFloat max = new TextFieldFloat();
          max.setColumns(8);
          c.gridx += 1;
          panel.add(min, c);
          c.gridx += 1;
          panel.add(max, c);
          
          max.addKeyListener(new FilterListener(hLine, false, j, num));
          min.addKeyListener(new FilterListener(hLine, true, j, num));
        }
        
        if(j<num-1)
        {
          c.gridx += 1;
          panel.add(new JLabel(":"), c);
        }
      }
    }

    // Filter out homozygous reference samples, GT - 0/0
    final JCheckBox showHomozygousMenu = new JCheckBox("Homozygous sample (GT)", vcfView.showHomozygous);
    c.gridx = 0;
    c.gridy += 1;
    c.gridwidth = 2;
    panel.add(showHomozygousMenu, c);
    showHomozygousMenu.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        setFlagFilter(HeaderLine.FILTER_HOMOZYG, "FILTER_HOMOZYG", "Homozygous", showHomozygousMenu.isSelected());
        vcfView.showHomozygous = showHomozygousMenu.isSelected();
        vcfView.repaint();
      }
    });
    c.gridwidth = 1;
    
    c.gridx = 100;
    c.weightx = 200.d;
    c.fill = GridBagConstraints.HORIZONTAL;
    panel.add(new JLabel(""), c);
    c.gridy += 1;
    c.weighty = 200.d;
    c.fill = GridBagConstraints.VERTICAL;
    panel.add(new JLabel(""), c);
    c.fill = GridBagConstraints.NONE;
  }
  
  /**
   * Check a VCF record against the hash of manually annotated variants
   * @param manualHash - hash containing manually annotated variants
   * @param record     - a VCF record
   * @param vcfIndex   - the file index for the VCF record
   * @return 0 - if no manual annotation
   *         1 - if manually passed
   *         2 - if manually failed
   */
  protected static int checkManualHash(final Map<String, Boolean> manualHash, 
                                       final VCFRecord record,
                                       final int vcfIndex, 
                                       final boolean applyManualFilter)
  {
    if(applyManualFilter || manualHash.size() == 0)
      return 0;
    
    final StringBuilder keyName = new StringBuilder();
    keyName.append(record.getPos());
    keyName.append(":");
    keyName.append(record.getChrom());
    keyName.append(":");
    keyName.append(vcfIndex);
    final String keyStr = keyName.toString();

    if(manualHash.containsKey(keyStr))
    {
      if(manualHash.get(keyStr))
        return 1;
      else
        return 2;
    }
    return 0;
  }
  
  protected static int checkManualHash(final Map<String, Boolean> manualHash, final VCFRecord record, final int vcfIndex)
  {
    return checkManualHash(manualHash, record, vcfIndex, manualFilter);
  }
  
  /**
   * Test for a given VCF record to see if it passes the filters.
   * @param manualHash
   * @param record
   * @param vcfReader
   * @param features
   * @param basePosition
   * @param sampleIndex
   * @return
   */
  protected static boolean passFilter(final Map<String, Boolean> manualHash, 
                                      final VCFRecord record, 
                                      final AbstractVCFReader vcfReader, 
                                      final FeatureVector features, 
                                      final int basePosition, 
                                      final int sampleIndex,
                                      final int vcfIndex)
  {
    final int manual = checkManualHash(manualHash, record, vcfIndex);
    if(manual == 1)
      return true;
    else if(manual == 2)
      return false;
    
    try
    {
      if(sampleIndex > -1)      // look at a specific sample and ignore
      {                         // those without values
        final String sample = record.getFormatValueForSample(sampleIndex);
        if(sample == null || sample.equals("."))
          return false;
      }
      
      if(record.getFilter().equals(".") || record.getFilter().equals("PASS"))
      {
      }
      else
      {
        // check if filter still being applied
        final String filterStr[] = SEMICOLON_PATTERN.split(record.getFilter());
        for(String s: filterStr)
        {
          if(FilteredPanel.getHeaderLineFiltersIDs().contains(s))
            return false;
        }
      }

      final Hashtable<String, RecordFilter> filters = FilteredPanel.getFilters();
      final Enumeration<String> enumFilter = filters.keys();
      while(enumFilter.hasMoreElements())
      {
        final String recid = enumFilter.nextElement();
        final RecordFilter recFilter = filters.get(recid);
        final String id = recFilter.getHeaderLine().getID();

        switch (recFilter.getHeaderLine().getHeaderType()) 
        {
          case HeaderLine.FILTER_NV_FLAG:  // FILTER non-variant
            if( record.getAlt().isNonVariant() )
              return false;
            break;
          case HeaderLine.FILTER_INS:
            if( record.getAlt().isInsertion(vcfReader.isVcf_v4()) )
              return false;
            break;
          case HeaderLine.FILTER_DEL:
            if( record.getAlt().isDeletion(vcfReader.isVcf_v4()) )
              return false;
            break;
          case HeaderLine.FILTER_OVERLAP_FLAG:  // FILTER no overlap
            if( !VCFview.isOverlappingFeature(features, basePosition) )
              return false;
            break;
          case HeaderLine.INFO_LINE:  // INFO line
            if (recFilter.getHeaderLine().isFlag())
            {
              if(record.containsInfoFlag(id))
                return false;
              continue;
            }
            else
            {
              String inf = record.getInfoValue(id);
              if(inf == null || !recFilter.pass(record, inf.split(","), vcfReader))
                return false;
            }
            break;
          case HeaderLine.FORMAT_LINE:  // FORMAT Genotype line
            final String samples[] = record.getFormatValues(id);
            if(samples == null)
              return false;
             if (recFilter.getHeaderLine().isFlag())
              return true;
              
            if(sampleIndex > -1) // look at a specific sample
            {
              if(samples[sampleIndex] == null || !recFilter.pass(record, samples[sampleIndex].split(","), vcfReader))
                return false;
            }
            else                 // look at all samples
            {
              for(int i=0; i<samples.length; i++)
              {
                if( samples[i] == null || !recFilter.pass(record, samples[i].split(","), vcfReader))
                  return false;
              }
            }
            break;
          case HeaderLine.FILTER_LINE:  // FILTER
            break;
          case HeaderLine.FILTER_QUAL:  // FILTER by quality score
            if( !recFilter.pass(record, new String[] { Float.toString(record.getQuality()) }, vcfReader))
              return false;
            break;
            
          case HeaderLine.FILTER_MULTALL_FLAG:
            if( record.getAlt().isMultiAllele(sampleIndex) )
              return false;
            break;
          case HeaderLine.FILTER_HOMOZYG:
            final String smpls[] = record.getFormatValues("GT");
            // look at a specific sample
            if(smpls != null && 
               sampleIndex > -1 && smpls[sampleIndex] != null && 
               HOMOZYGOUS_PATTERN.matcher(smpls[sampleIndex]).matches())
              return false;
            break;
          case HeaderLine.FILTER_NONSYN:
            if( !record.getAlt().isDeletion(vcfReader.isVcf_v4()) && 
                !record.getAlt().isInsertion(vcfReader.isVcf_v4()) && 
                record.getAlt().length() == 1 && 
                record.getRef().length() == 1)
           {
             short isSyn = record.getSynFlag(features, basePosition);
             if( (isSyn == 0 || isSyn == 2) )
               return false;
           }
           break;
          case HeaderLine.FILTER_SYN:
            if( !record.getAlt().isDeletion(vcfReader.isVcf_v4()) && 
                !record.getAlt().isInsertion(vcfReader.isVcf_v4()) && 
                record.getAlt().length() == 1 && 
                record.getRef().length() == 1)
           {
             short isSyn = record.getSynFlag(features, basePosition);
             if(isSyn == 1)
               return false;
           }
          default:
            break;
        }
      }
      return true;

/*      try
      {
        if(VCFFilter.MIN_DP > 0 && Integer.parseInt(record.getInfoValue("DP")) < VCFFilter.MIN_DP)
          return false;
      }
      catch(NullPointerException npe){}
    
      try
      {
        if(VCFFilter.MIN_MQ > 0 && Float.parseFloat(record.getInfoValue("MQ")) < VCFFilter.MIN_MQ)
          return false;
      }
      catch(NullPointerException npe){}

      try
      {
        // AF1 filtered - except for non-variant sites
        if((!record.getAlt().isNonVariant()) && VCFFilter.MIN_AF1 > 0 && Float.parseFloat(record.getInfoValue("AF1")) < VCFFilter.MIN_AF1)
          return false;
      }
      catch(NullPointerException npe){}

      try
      {
        String vals[] = COMMA_PATTERN.split(record.getInfoValue("CI95"));
        for(int i=0; i<vals.length; i++)
        {
          if(VCFFilter.MAX_CI95 < 10 && Float.parseFloat(vals[i]) > VCFFilter.MAX_CI95)
            return false;
        }
      }
      catch(NullPointerException npe){}*/
    }
    catch(NumberFormatException e)
    {
      System.err.println(e.getMessage()); 
    }

    return true;
  }
  
  /**
   * Test for a given VCF record to see if it passes the filters.
   * @param manualHash
   * @param record
   * @param vcfView
   * @param basePosition
   * @param features
   * @param vcfReader
   */
  protected static void setFilterString(final Map<String, Boolean> manualHash, 
                                        final VCFRecord record, 
                                        final VCFview vcfView, 
                                        final int basePosition, 
                                        final FeatureVector features, 
                                        final AbstractVCFReader vcfReader,
                                        final int vcfIndex)
  {
    
    
    final int manual = checkManualHash(manualHash, record, vcfIndex);
    if(manual == 1)
    {
      record.setFilter("PASS");
      return;
    }

    if(record.getFilter().equals(".") || record.getFilter().equals("PASS"))
    {
    }
    else
    {
      record.setFilter("");
      // check if filter still being applied
      final List<String> ids = FilteredPanel.getHeaderLineFiltersIDs();
      for(int i=0; i<ids.size(); i++)
      {
        record.appendFilter(ids.get(i));
        if(i<ids.size()-1)
          record.appendFilter(";");
      }
    }

    try
    {  
      // INFO, FORMAT
      if(FilteredPanel.getFilters().size() > 0)
      {
        Hashtable<String, RecordFilter> filters = FilteredPanel.getFilters();
        Enumeration<String> enumFilter = filters.keys();
        while(enumFilter.hasMoreElements())
        {
          String recid = enumFilter.nextElement();
          RecordFilter recFilter = filters.get(recid);
          String id = recFilter.getHeaderLine().getID();

          switch (recFilter.getHeaderLine().getHeaderType()) 
          {
            case HeaderLine.FILTER_NV_FLAG:  // FILTER non-variant
              if( record.getAlt().isNonVariant() )
                record.appendFilter(id);
              break;
            case HeaderLine.FILTER_INS:
              if( record.getAlt().isInsertion(vcfReader.isVcf_v4()) )
                record.appendFilter(id);
              break;
            case HeaderLine.FILTER_DEL:
              if( record.getAlt().isDeletion(vcfReader.isVcf_v4()) )
                record.appendFilter(id);
              break;
            case HeaderLine.FILTER_OVERLAP_FLAG:  // FILTER no overlap
              if( !VCFview.isOverlappingFeature(features, basePosition) )
                record.appendFilter(id);
            break;
            case HeaderLine.INFO_LINE:  // INFO line
              if (recFilter.getHeaderLine().isFlag())
              {
                if (recFilter.getHeaderLine().isFlag())
                {
                  if (record.containsInfoFlag(id))
                    record.appendFilter(id);
                }
              }
              else if ( record.getInfoValue(id) == null || 
                       !recFilter.pass(record, record.getInfoValue(id).split(","), vcfReader))
                record.appendFilter(id);
              break;
              
            case HeaderLine.FORMAT_LINE:  // FORMAT Genotype line
              final String samples[] = record.getFormatValues(id);
              
              id = "sample" + recFilter.getHeaderLine().getID();
              if(samples == null)
              {
                record.appendFilter(id);
                break;
              }

              //if (recFilter.getHeaderLine().isFlag())
              //{
              //  return true;
              //}
              
              for(int i=0; i<samples.length; i++)
              {
                if( samples[i] == null || !recFilter.pass(record, samples[i].split(","), vcfReader))
                {
                  record.appendFilter(id);
                  break;
                }
              }
              
              break;
            case HeaderLine.FILTER_LINE:  // FILTER
              break;
            case HeaderLine.FILTER_QUAL:  // FILTER by quality score
              if( !recFilter.pass(record, new String[] { Float.toString(record.getQuality()) }, vcfReader))
                record.appendFilter(id);
              break;
            case HeaderLine.FILTER_MULTALL_FLAG:  // FILTER by quality score
              if( record.getAlt().isMultiAllele(-1) )
                record.appendFilter(id);
              break;
            case HeaderLine.FILTER_NONSYN:
              if( !record.getAlt().isDeletion(vcfReader.isVcf_v4()) && 
                  !record.getAlt().isInsertion(vcfReader.isVcf_v4()) && 
                  record.getAlt().length() == 1 && 
                  record.getRef().length() == 1)
             {
               short isSyn = record.getSynFlag(features, basePosition);
               if( (isSyn == 0 || isSyn == 2) )
                 record.appendFilter(id);
             }
              break;
            case HeaderLine.FILTER_SYN:
              if( !record.getAlt().isDeletion(vcfReader.isVcf_v4()) && 
                  !record.getAlt().isInsertion(vcfReader.isVcf_v4()) && 
                  record.getAlt().length() == 1 && 
                  record.getRef().length() == 1)
             {
               short isSyn = record.getSynFlag(features, basePosition);
               if(isSyn == 1)
                 record.appendFilter(id);
             }
            case HeaderLine.FILTER_MANUAL:
              if(manual == 2)
                record.appendFilter(id);
            default:
              break;
          }
        }

        if(  record.getFilter().length() == 0 ||
            (record.getFilter().length() == 1 && record.getFilter().equals(".")) )
          record.setFilter("PASS");
        return;
      }

      // OTHERS
/*      try
      {
        if(VCFFilter.MIN_DP > 0 && Integer.parseInt(record.getInfoValue("DP")) < VCFFilter.MIN_DP)
          record.appendFilter("DP");
      }
      catch(NullPointerException npe){}
    
      try
      {
        if(VCFFilter.MIN_MQ > 0 && Float.parseFloat(record.getInfoValue("MQ")) < VCFFilter.MIN_MQ)
          record.appendFilter("MQ");
      }
      catch(NullPointerException npe){}

      try
      {
        // AF1 filtered - except for non-variant sites
        if((!record.getAlt().isNonVariant()) && VCFFilter.MIN_AF1 > 0 && Float.parseFloat(record.getInfoValue("AF1")) < VCFFilter.MIN_AF1)
          record.appendFilter("AF1");
      }
      catch(NullPointerException npe){}

      try
      {
        String vals[] = COMMA_PATTERN.split(record.getInfoValue("CI95"));
        for(int i=0; i<vals.length; i++)
        {
          if(VCFFilter.MAX_CI95 < 10 && Float.parseFloat(vals[i]) > VCFFilter.MAX_CI95)
            record.appendFilter("CI95");
        }
      }
      catch(NullPointerException npe){}*/
    }
    catch(NumberFormatException e)
    {
      System.err.println(e.getMessage()); 
    }
    
    if(record.getFilter().length() == 0)
      record.setFilter("PASS");
  }
  
  

  
  /**
   * Filter listener for integers and floats
   */
  class FilterListener extends KeyAdapter
  {
    private HeaderLine hLine;
    private boolean isMin;
    private int index;
    private int NUMBER;
    
    public FilterListener(HeaderLine hLine, boolean isMin, int index, int NUMBER)
    {
      this.hLine = hLine;
      this.isMin = isMin;
      this.index = index;
      this.NUMBER = NUMBER;
    }

    public void keyReleased(KeyEvent e) 
    {
      String type = hLine.getType();
      if (type.equals("Integer"))
        setFilterForInt(e);
      else if (type.equals("Float"))
        setFilterForFloat(e);
      filterPanel.updateFilters();
    }

    private void setFilterForInt(KeyEvent e)
    {  
      String ID = hLine.getHeaderTypeStr()+":"+hLine.getID();
      TextFieldInt field = (TextFieldInt)e.getComponent();
      RecordFilter filter;
      
      Hashtable<String, RecordFilter> filters = FilteredPanel.getFilters();
      if(filters.containsKey(ID))
        filter = filters.get(ID);
      else
        filter = new RecordFilter(hLine, NUMBER);

      if(field.getText().trim().equals(""))
      {
        if(isMin)
          filter.minIVal[index] = Integer.MIN_VALUE;
        else
          filter.maxIVal[index] = Integer.MAX_VALUE;
      }
      else if(isMin)
        filter.minIVal[index] = field.getValue();
      else
        filter.maxIVal[index] = field.getValue();
      
      if( filter.minIVal[index] == Integer.MIN_VALUE && 
          filter.maxIVal[index] == Integer.MAX_VALUE)
        filters.remove(ID);
      else
        filters.put(ID, filter);
    }
    
    private void setFilterForFloat(KeyEvent e)
    {
      String ID = hLine.getHeaderTypeStr()+":"+hLine.getID();
      TextFieldFloat field = (TextFieldFloat)e.getComponent();
      RecordFilter filter;
      Hashtable<String, RecordFilter> filters = FilteredPanel.getFilters();
      if(filters.containsKey(ID))
        filter = filters.get(ID);
      else
        filter = new RecordFilter(hLine, NUMBER);

      if(field.getText().trim().equals(""))
      {
        if(isMin)
          filter.minFVal[index] = Float.MIN_VALUE;
        else
          filter.maxFVal[index] = Float.MAX_VALUE;
      }
      else if(isMin)
        filter.minFVal[index] = (float) field.getValue();
      else
        filter.maxFVal[index] = (float) field.getValue();
      
      if( filter.minFVal[index] == Float.MIN_VALUE && 
          filter.maxFVal[index] == Float.MAX_VALUE)
        filters.remove(ID);
      else
        filters.put(ID, filter);
    }
  }

}