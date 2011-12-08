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
import java.awt.Insets;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.List;
import java.util.regex.Pattern;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextField;

import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.circular.TextFieldFloat;
import uk.ac.sanger.artemis.circular.TextFieldInt;
import uk.ac.sanger.artemis.components.Utilities;

public class VCFFilter extends JFrame
{
  private static final long serialVersionUID = 1L;

  private static int MIN_DP = 0;
  private static float MIN_MQ = 0;
  private static float MIN_AF1 = 0;
  private static float MAX_CI95 = 10;
  private static Pattern COMMA_PATTERN = Pattern.compile(",");
  private FilteredPanel filterPanel;

  /**
   * Filter VCF records by the variant type and/or by different values in 
   * the record, QUAL, DP, MQ and AF1.
   * @param vcfView
   */
  public VCFFilter(final VCFview vcfView)
  {
    super("Variant Filter");
    GridBagConstraints c = new GridBagConstraints();
    GridBagConstraints cMain = new GridBagConstraints();
    JPanel mainPanel = (JPanel)getContentPane();
    JPanel btmPanel = new JPanel(new GridBagLayout());
    JPanel yBox = new JPanel(new GridBagLayout());
    JScrollPane jsp = new JScrollPane(yBox);
    
    final List<HeaderLine> filterHdrLines = vcfView.getVcfReaders()[0].getFILTER();
    
    filterPanel = new FilteredPanel(filterHdrLines);
    JScrollPane ftrScroll = new JScrollPane(filterPanel);
    ftrScroll.setPreferredSize( new Dimension(ftrScroll.getPreferredSize().width, 150) );
    
    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
    jsp.setPreferredSize(new Dimension(screen.width*9/20, screen.height/2));
    mainPanel.add(jsp, BorderLayout.NORTH);
    mainPanel.add(ftrScroll, BorderLayout.CENTER);
    mainPanel.add(btmPanel, BorderLayout.SOUTH);
    
    c.gridx = 0;
    c.gridy = 0;
    c.anchor = GridBagConstraints.WEST;

    ////
    ////
    //// Filter by type
    JPanel typePanel = new JPanel(new GridBagLayout());
    typePanel.setBorder(BorderFactory.createLineBorder(Color.gray));
    cMain.gridx = 0;
    cMain.gridy = 0;
    cMain.insets = new Insets(2, 2, 2, 2);
    cMain.anchor = GridBagConstraints.WEST;
    yBox.add(typePanel, cMain);
    JLabel typeLabel = new JLabel("TYPE:");
    typeLabel.setFont(typeLabel.getFont().deriveFont(Font.BOLD));
    typePanel.add(typeLabel, c);
    
    c.gridy = c.gridy + 1;
    final JCheckBox showSyn = new JCheckBox("Synonymous", vcfView.showSynonymous);
    typePanel.add(showSyn, c);
    showSyn.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        vcfView.showSynonymous = showSyn.isSelected();
        vcfView.repaint();
      }
    });

    c.gridy = c.gridy + 1;
    final JCheckBox showNonSyn = new JCheckBox("Non-synonymous", vcfView.showNonSynonymous);
    typePanel.add(showNonSyn, c);
    showNonSyn.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        vcfView.showNonSynonymous = showNonSyn.isSelected();
        vcfView.repaint();
      }
    });
    
    c.gridy = c.gridy + 1;
    final JCheckBox showDeletionsMenu = new JCheckBox("Deletions", vcfView.showDeletions);
    typePanel.add(showDeletionsMenu, c);
    showDeletionsMenu.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        vcfView.showDeletions = showDeletionsMenu.isSelected();
        vcfView.repaint();
      }
    });

    c.gridy = c.gridy + 1;
    final JCheckBox showInsertionsMenu = new JCheckBox("Insertions", vcfView.showInsertions);
    typePanel.add(showInsertionsMenu, c);
    showInsertionsMenu.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        vcfView.showInsertions = showInsertionsMenu.isSelected();
        vcfView.repaint();
      }
    });
    
    c.gridy = c.gridy + 1;
    final JCheckBox showMultiAllelesMenu = new JCheckBox("Multiple alleles", vcfView.showMultiAlleles);
    typePanel.add(showMultiAllelesMenu, c);
    showMultiAllelesMenu.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        vcfView.showMultiAlleles = showMultiAllelesMenu.isSelected();
        vcfView.repaint();
      }
    });
    
    c.gridy = c.gridy + 1;
    final JCheckBox showNonOverlappingsMenu = new JCheckBox("Variants not overlapping CDS", vcfView.showNonOverlappings);
    typePanel.add(showNonOverlappingsMenu, c);
    showNonOverlappingsMenu.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        vcfView.showNonOverlappings = showNonOverlappingsMenu.isSelected();
        vcfView.repaint();
      }
    });
    
    c.gridy = c.gridy + 1;
    final JCheckBox showNonVariantMenu = new JCheckBox("Non-Variants", vcfView.showNonVariants);
    typePanel.add(showNonVariantMenu, c);
    showNonVariantMenu.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        vcfView.showNonVariants = showNonVariantMenu.isSelected();
        vcfView.repaint();
      }
    });
    
    if(vcfView.getEntryGroup() == null || vcfView.getEntryGroup().getAllFeaturesCount() == 0)
    {
      showSyn.setEnabled(false);
      showNonSyn.setEnabled(false);
      showNonOverlappingsMenu.setEnabled(false);
    }
    
    ////
    ////
    //// Filter by property
    final List<HeaderLine> info = vcfView.getVcfReaders()[0].getINFO();
    
    JPanel propPanel = new JPanel(new GridBagLayout());
    propPanel.setBorder(BorderFactory.createLineBorder(Color.gray));
    cMain.gridy = cMain.gridy+1;
    yBox.add(propPanel, cMain);
    c.gridy = c.gridy+1;
    propPanel.add(new JLabel(" "), c);
    c.gridy = c.gridy+1;
    JLabel propLabel = new JLabel("PROPERTY:");
    propLabel.setFont(propLabel.getFont().deriveFont(Font.BOLD));
    propPanel.add(propLabel, c);
    
    // min quality
    c.gridy = c.gridy + 1;
    c.gridx = 1;
    propPanel.add(new JLabel("MIN"), c);
    c.gridx = c.gridx + 1;
    propPanel.add(new JLabel("MAX"), c);;
    
    c.gridy = c.gridy+1;
    c.gridx = 0;
    c.anchor = GridBagConstraints.WEST;
    propPanel.add(new JLabel("Quality score (QUAL):"), c);
    final TextFieldFloat minQuality = new TextFieldFloat();
    minQuality.setColumns(8);
    c.gridx = 1;
    propPanel.add(minQuality, c); 
    final TextFieldFloat maxQuality = new TextFieldFloat();
    maxQuality.setColumns(8);
    c.gridx += 1;
    propPanel.add(maxQuality, c); 

    final JTextField minDP = new JTextField(Integer.toString(MIN_DP), 8);
    final JTextField minMQ = new JTextField(Float.toString(MIN_MQ),8);
    final JTextField minAF1 = new JTextField(Float.toString(MIN_AF1),8);
    final JTextField maxCI95 = new JTextField(Float.toString(MAX_CI95),8);
    
    if(info.size() == 0)
    {
      // min DP
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
      propPanel.add(maxCI95, c);
    }
    else
    {
      //
      createPanel(propPanel, c, info, "INFO FIELDS:");
    }

    ////
    ////
    //// FORMAT - Genotype Filter
    final List<HeaderLine> format = vcfView.getVcfReaders()[0].getFORMAT();
    
    JPanel formatPanel = new JPanel(new GridBagLayout());
    formatPanel.setBorder(BorderFactory.createLineBorder(Color.gray));
    cMain.gridy = cMain.gridy+1;
    yBox.add(formatPanel, cMain);
    
    createPanel(formatPanel, c, format, "GENOTYPE FIELDS:");

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
          checkQualityFiltering(minQuality, maxQuality);

          if(info.size() == 0)
          {
            MIN_DP = Integer.parseInt(minDP.getText());
            MIN_MQ = Float.parseFloat(minMQ.getText());
            MIN_AF1 = Float.parseFloat(minAF1.getText());
            MAX_CI95 = Float.parseFloat(maxCI95.getText());
          }
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
        checkQualityFiltering(minQuality, maxQuality);
        vcfView.repaint();
        setVisible(false);
      }
    });
    btmPanel.add(close, c);
    
    pack();
    Utilities.centreFrame(this);
    setVisible(true);
  }
  
  /**
   * Check the values in the quality (QUAL) filtering fields and update the filter panel.
   * @param minQuality
   * @param maxQuality
   */
  private void checkQualityFiltering(final TextFieldFloat minQuality, final TextFieldFloat maxQuality)
  {
    Hashtable<String, RecordFilter> filters = FilteredPanel.getFilters();
    String ID = "QUAL:QUAL";
    if (minQuality.getText().equals("") && maxQuality.getText().equals(""))
      filters.remove(ID);
    else
    {
      String hdrLine = "##Filter=<ID=QUAL,Type=Float,Number=1,Description=\"QUAL "
          + (minQuality.getText().equals("") ? "" : " < " + minQuality.getValue() + " ")
          + (maxQuality.getText().equals("") ? "" : " > " + maxQuality.getValue() + " ") + "\">";
      HeaderLine hLine = new HeaderLine(hdrLine, "QUAL",
          AbstractVCFReader.getLineHash("FILTER", hdrLine));

      RecordFilter filter = new RecordFilter(hLine, 1);
      if (minQuality.getText().equals(""))
        filter.minFVal[0] = Float.MIN_VALUE;
      else
        filter.minFVal[0] = (float) minQuality.getValue();

      if (maxQuality.getText().equals(""))
        filter.maxFVal[0] = Float.MAX_VALUE;
      else
        filter.maxFVal[0] = (float) maxQuality.getValue();
      filters.put(hLine.getHeaderTypeStr() + ":" + hLine.getID(), filter);
    }
    filterPanel.updateFilters();
  }
  
  private void createPanel(JPanel panel, GridBagConstraints c, List<HeaderLine> headerLineList, String name)
  {
    c.gridx = 0;
    c.gridy = c.gridy+1;
    panel.add(new JLabel(" "), c);
    c.gridy = c.gridy+1;
    JLabel nameLabel = new JLabel(name);
    nameLabel.setFont(panel.getFont().deriveFont(Font.BOLD));
    panel.add(nameLabel, c);
    
    c.gridy = c.gridy+1;
    panel.add(new JLabel(" "), c);
    c.gridy = c.gridy+1;

    
    for (int i = 0; i < headerLineList.size(); i++)
    {
      final HeaderLine hLine = headerLineList.get(i);
      final String type = hLine.getType();
      if (type.equals("String"))
        continue;

      int num = hLine.getNumber();

      c.gridy = c.gridy + 1;
      c.gridx = 0;
      for (int j = 0; j < num; j++)
      {
        c.gridx = c.gridx + 1;
        panel.add(new JLabel("MIN"), c);
        c.gridx = c.gridx + 1;
        panel.add(new JLabel("MAX"), c);
        c.gridx = c.gridx + 1;
      }
      
      c.gridy = c.gridy + 1;
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
        c.gridx = c.gridx + 1;
        panel.add(flag, c);
      }

      for (int j = 0; j < num; j++)
      {
        if (type.equals("Integer"))
        {
          final TextFieldInt min = new TextFieldInt();
          min.setColumns(8);
          final TextFieldInt max = new TextFieldInt();
          max.setColumns(8);
          c.gridx = c.gridx + 1;
          panel.add(min, c);
          c.gridx = c.gridx + 1;
          panel.add(max, c);
          
          max.addKeyListener(new FilterListener(hLine, false, j, num));
          min.addKeyListener(new FilterListener(hLine, true, j, num));
        }
        else if (type.equals("Float"))
        {
          TextFieldFloat min = new TextFieldFloat();
          min.setColumns(8);
          TextFieldFloat max = new TextFieldFloat();
          max.setColumns(8);
          c.gridx = c.gridx + 1;
          panel.add(min, c);
          c.gridx = c.gridx + 1;
          panel.add(max, c);
          
          max.addKeyListener(new FilterListener(hLine, false, j, num));
          min.addKeyListener(new FilterListener(hLine, true, j, num));
        }
        
        if(j<num-1)
        {
          c.gridx = c.gridx + 1;
          panel.add(new JLabel(":"), c);
        }
      }
    }
  }
  
  /**
   * Test for a given VCF record to see if it passes the filters.
   * @param record
   * @return
   */
  protected static boolean passFilter(final VCFRecord record, final AbstractVCFReader vcfReader)
  {
    try
    {
      if(record.getFilter().equals(".") || record.getFilter().equals("PASS"))
      {
      }
      else
        return false;

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
            case HeaderLine.INFO_LINE:  // INFO line
              if (recFilter.getHeaderLine().isFlag())
              {
                if(record.containsInfoFlag(id))
                  return false;
                continue;
              }
              else if( record.getInfoValue(id) == null || 
                       !recFilter.pass(record, record.getInfoValue(id).split(","), vcfReader))
                return false;
              break;
            case HeaderLine.FORMAT_LINE:  // FORMAT Genotype line
              String samples[] = record.getFormatValues(id);
              if(samples == null)
                return false;

              if (recFilter.getHeaderLine().isFlag())
                return true;
              for(int i=0; i<samples.length; i++)
              {
                if( !recFilter.pass(record, samples[i].split(","), vcfReader))
                  return false;
              }

              break;
            case HeaderLine.FILTER_LINE:  // FILTER
              break;
            case HeaderLine.FILTER_QUAL:  // FILTER by quality score
              if( !recFilter.pass(record, new String[] { Float.toString(record.getQuality()) }, vcfReader))
                return false;
              break;
            default:
              break;
          }
        }
        return true;
      }
      
      try
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
      catch(NullPointerException npe){}
    }
    catch(NumberFormatException e)
    {
      System.err.println(e.getMessage()); 
    }

    return true;
  }
  
  /**
   * Test for a given VCF record to see if it passes the filters.
   * @param record
   * @param vcfView
   * @param basePosition
   * @param features
   * @param vcf_v4
   */
  protected static void setFilterString(final VCFRecord record, final VCFview vcfView, 
      final int basePosition, final FeatureVector features, final AbstractVCFReader vcfReader)
  {
    boolean vcf_v4 = vcfReader.isVcf_v4();
    record.setFilter("");

    try
    {
      // TYPE
/*      if(record.getQuality() < VCFFilter.MIN_QUALITY &&
         record.getQuality() > VCFFilter.MAX_QUALITY)
        record.appendFilter("QUAL");*/

      if(!vcfView.showDeletions && record.getAlt().isDeletion(vcf_v4))
        record.appendFilter("DEL");
      if(!vcfView.showInsertions && record.getAlt().isInsertion(vcf_v4))
        record.appendFilter("IN");
      if(!vcfView.showNonOverlappings && !vcfView.isOverlappingFeature(features, basePosition))
        record.appendFilter("OVERLAP");
      if(!vcfView.showNonVariants && record.getAlt().isNonVariant())
        record.appendFilter("NV");
      if(!vcfView.showMultiAlleles && record.getAlt().isMultiAllele())
        record.appendFilter("MA");

      if( (!vcfView.showSynonymous || !vcfView.showNonSynonymous) &&
           !record.getAlt().isDeletion(vcf_v4) && 
           !record.getAlt().isInsertion(vcf_v4) && 
           record.getAlt().length() == 1 && 
           record.getRef().length() == 1)
      {
        short isSyn = record.getSynFlag(features, basePosition);

        if(!vcfView.showSynonymous && isSyn == 1) 
          record.appendFilter("SYN");
        if(!vcfView.showNonSynonymous && (isSyn == 0 || isSyn == 2))
          record.appendFilter("NONSYN");
      }
      
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
            case HeaderLine.FORMAT_LINE:  // FORMAT
              break;
            case HeaderLine.FILTER_LINE:  // FILTER
              break;
            case HeaderLine.FILTER_QUAL:  // FILTER by quality score
              if( !recFilter.pass(record, new String[] { Float.toString(record.getQuality()) }, vcfReader))
                record.appendFilter(id);
              break;
            default:
              break;
          }
        }

        if(record.getFilter().length() == 0)
          record.setFilter("PASS");
        return;
      }

      // OTHERS
      try
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
      catch(NullPointerException npe){}
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