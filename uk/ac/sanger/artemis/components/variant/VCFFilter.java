package uk.ac.sanger.artemis.components.variant;

import java.awt.BorderLayout;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.List;
import java.util.regex.Pattern;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextField;

import uk.ac.sanger.artemis.circular.TextFieldFloat;
import uk.ac.sanger.artemis.circular.TextFieldInt;
import uk.ac.sanger.artemis.components.Utilities;

public class VCFFilter extends JFrame
{
  private static final long serialVersionUID = 1L;
  private static float MIN_QUALITY = 0;
  private static int MIN_DP = 0;
  private static float MIN_MQ = 0;
  private static float MIN_AF1 = 0;
  private static float MAX_CI95 = 10;
  
  private static Pattern COMMA_PATTERN = Pattern.compile(",");
  
  private static Hashtable<String, InfoFilter> filters = null;
  
  /**
   * Filter VCF records by the variant type and/or by different values in 
   * the record, QUAL, DP, MQ and AF1.
   * @param vcfView
   */
  public VCFFilter(final VCFview vcfView)
  {
    super("Variant Filter");
    GridBagConstraints c = new GridBagConstraints();
    
    JPanel mainPanel = (JPanel)getContentPane();
    JPanel btmPanel = new JPanel(new GridBagLayout());
    JPanel panel = new JPanel(new GridBagLayout());
    JScrollPane jsp = new JScrollPane(panel);
    mainPanel.add(jsp, BorderLayout.CENTER);
    mainPanel.add(btmPanel, BorderLayout.SOUTH);
    
    // Filter by type
    c.gridx = 0;
    c.gridy = 0;
    c.anchor = GridBagConstraints.WEST;
    JLabel typeLabel = new JLabel("TYPE:");
    typeLabel.setFont(typeLabel.getFont().deriveFont(Font.BOLD));
    panel.add(typeLabel, c);
    
    c.gridy = c.gridy + 1;
    final JCheckBox showSyn = new JCheckBox("Synonymous", vcfView.showSynonymous);
    panel.add(showSyn, c);
    showSyn.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        vcfView.showSynonymous = showSyn.isSelected();
        vcfView.repaint();
      }
    });

    c.gridy = c.gridy + 1;
    final JCheckBox showNonSyn = new JCheckBox("Non-synonymous", vcfView.showNonSynonymous);
    panel.add(showNonSyn, c);
    showNonSyn.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        vcfView.showNonSynonymous = showNonSyn.isSelected();
        vcfView.repaint();
      }
    });
    
    c.gridy = c.gridy + 1;
    final JCheckBox showDeletionsMenu = new JCheckBox("Deletions", vcfView.showDeletions);
    panel.add(showDeletionsMenu, c);
    showDeletionsMenu.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        vcfView.showDeletions = showDeletionsMenu.isSelected();
        vcfView.repaint();
      }
    });

    c.gridy = c.gridy + 1;
    final JCheckBox showInsertionsMenu = new JCheckBox("Insertions", vcfView.showInsertions);
    panel.add(showInsertionsMenu, c);
    showInsertionsMenu.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        vcfView.showInsertions = showInsertionsMenu.isSelected();
        vcfView.repaint();
      }
    });
    
    c.gridy = c.gridy + 1;
    final JCheckBox showMultiAllelesMenu = new JCheckBox("Multiple alleles", vcfView.showMultiAlleles);
    panel.add(showMultiAllelesMenu, c);
    showMultiAllelesMenu.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        vcfView.showMultiAlleles = showMultiAllelesMenu.isSelected();
        vcfView.repaint();
      }
    });
    
    c.gridy = c.gridy + 1;
    final JCheckBox showNonOverlappingsMenu = new JCheckBox("Variants not overlapping CDS", vcfView.showNonOverlappings);
    panel.add(showNonOverlappingsMenu, c);
    showNonOverlappingsMenu.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        vcfView.showNonOverlappings = showNonOverlappingsMenu.isSelected();
        vcfView.repaint();
      }
    });
    
    c.gridy = c.gridy + 1;
    final JCheckBox showNonVariantMenu = new JCheckBox("Non-Variants", vcfView.showNonVariants);
    panel.add(showNonVariantMenu, c);
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
    
    // Filter by property
    c.gridy = c.gridy+1;
    panel.add(new JLabel(" "), c);
    c.gridy = c.gridy+1;
    JLabel propLabel = new JLabel("PROPERTY:");
    propLabel.setFont(propLabel.getFont().deriveFont(Font.BOLD));
    panel.add(propLabel, c);
    
    // min quality
    c.gridy = c.gridy+1;
    c.gridx = 0;
    c.anchor = GridBagConstraints.WEST;
    panel.add(new JLabel("Minimum quality score (QUAL):"), c);
    final JTextField minQuality = new JTextField(Float.toString(MIN_QUALITY), 8);
    c.gridx = 1;
    panel.add(minQuality, c);    

    final JTextField minDP = new JTextField(Integer.toString(MIN_DP), 8);
    final JTextField minMQ = new JTextField(Float.toString(MIN_MQ),8);
    final JTextField minAF1 = new JTextField(Float.toString(MIN_AF1),8);
    final JTextField maxCI95 = new JTextField(Float.toString(MAX_CI95),8);
    
    final List<Hashtable<String, String>> info = vcfView.getVcfReaders()[0].getINFO();
    if(info.size() == 0)
    {
      // min DP
      c.gridy = c.gridy+1;
      c.gridx = 0;
      panel.add(new JLabel("Minimum combined depth across samples (DP):"), c);
      
      c.gridx = 1;
      panel.add(minDP, c);
      
      // min MQ
      c.gridy = c.gridy+1;
      c.gridx = 0;
      panel.add(new JLabel("Minimum RMS mapping quality (MQ):"), c);
      
      c.gridx = 1;
      panel.add(minMQ, c);
      
      // min AF1
      c.gridy = c.gridy+1;
      c.gridx = 0;
      panel.add(new JLabel("Minimum site frequency of strongest non-reference allele (AF1):"), c);
      
      c.gridx = 1;
      panel.add(minAF1, c);
      
      // max CI95
      c.gridy = c.gridy+1;
      c.gridx = 0;
      panel.add(new JLabel("Maximum 95% confidence interval variation from AF (CI95):"), c);

      c.gridx = 1;
      panel.add(maxCI95, c);
    }
    else
    {
      //
      c.gridy = c.gridy+1;
      panel.add(new JLabel(" "), c);
      c.gridy = c.gridy+1;
      c.gridx = 0;
      JLabel filterLabel = new JLabel("INFO Column:");
      filterLabel.setFont(propLabel.getFont().deriveFont(Font.BOLD));
      panel.add(filterLabel, c);
      
      filters = new Hashtable<String, InfoFilter>();
      
      for (int i = 0; i < info.size(); i++)
      {
        final Hashtable<String, String> hinfo = info.get(i);
        final String type = hinfo.get("Type");
        if (type.equals("String"))
          continue;

        int num = 1;
        String numStr = hinfo.get("Number");
        if (numStr != null)
        {
          try
          {
            num = Integer.parseInt(numStr);
          }
          catch(NumberFormatException e)
          {
            num = 1;
          }
        }

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
        
        JLabel lab = new JLabel(hinfo.get("ID") + " (" + type + ")");
        lab.setToolTipText(hinfo.get("Description"));
        panel.add(lab, c);

        if (num == 0 && type.equals("Flag"))
        {
          final JCheckBox flag = new JCheckBox();
          flag.addActionListener(new ActionListener()
          {
            public void actionPerformed(ActionEvent e)
            {
              if(flag.isSelected()) 
                filters.put(hinfo.get("ID"), new InfoFilter(type, 0));
              else
                filters.remove(hinfo.get("ID"));
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
            
            max.addKeyListener(new FilterListener(hinfo.get("ID"), type, false, j, num));
            min.addKeyListener(new FilterListener(hinfo.get("ID"), type, true, j, num));
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
            
            max.addKeyListener(new FilterListener(hinfo.get("ID"), type, false, j, num));
            min.addKeyListener(new FilterListener(hinfo.get("ID"), type, true, j, num));
          }
          
          if(j<num-1)
          {
            c.gridx = c.gridx + 1;
            panel.add(new JLabel(":"), c);
          }
        }
      }
    }
    
    //
    c.gridy = c.gridy+1;
    c.gridx = 0;
    c.anchor = GridBagConstraints.EAST;
    JButton apply = new JButton("Apply");
    apply.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        try
        {
          MIN_QUALITY = Float.parseFloat(minQuality.getText());
          
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
    JButton close = new JButton("Close");
    close.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        setVisible(false);
      }
    });
    btmPanel.add(close, c);
    
    pack();
    Utilities.centreFrame(this);
    setVisible(true);
  }
  
  
  /**
   * Test for a given VCF record to see if it passes the filters.
   * @param record
   * @return
   */
  protected static boolean passFilter(VCFRecord record)
  {
    try
    {
      if(record.getQuality() < VCFFilter.MIN_QUALITY)
        return false;

      if(filters != null)
      {
        Enumeration<String> enumFilter = filters.keys();
        while(enumFilter.hasMoreElements())
        {
          String id = enumFilter.nextElement();
          InfoFilter filter = filters.get(id);
          
          if( filter.type.equals("Flag"))
          {
            if(record.containsInfoFlag(id))
              return false;
            continue;
          }
          else if(record.getInfoValue(id) == null || !filter.pass(record.getInfoValue(id).split(",")))
            return false;
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
   * Filter listener for integers and floats
   */
  class FilterListener extends KeyAdapter
  {
    private String ID;
    private String type;
    private boolean isMin;
    private int index;
    private int NUMBER;
    
    public FilterListener(String ID, String type, boolean isMin, int index, int NUMBER)
    {
      this.ID = ID;
      this.type = type;
      this.isMin = isMin;
      this.index = index;
      this.NUMBER = NUMBER;
    }

    public void keyReleased(KeyEvent e) 
    {
      if (type.equals("Integer"))
        setFilterForInt(e);
      else if (type.equals("Float"))
        setFilterForFloat(e);
    }

    private void setFilterForInt(KeyEvent e)
    {  
      TextFieldInt field = (TextFieldInt)e.getComponent();
      InfoFilter filter;
      if(filters.containsKey(ID))
        filter = filters.get(ID);
      else
        filter = new InfoFilter(type, NUMBER);

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
      
/*      if( filter.minIVal == Integer.MIN_VALUE && 
          filter.maxIVal == Integer.MAX_VALUE)
        filters.remove(ID);
      else*/
        filters.put(ID, filter);
    }
    
    private void setFilterForFloat(KeyEvent e)
    {  
      TextFieldFloat field = (TextFieldFloat)e.getComponent();
      InfoFilter filter;
      if(filters.containsKey(ID))
        filter = filters.get(ID);
      else
        filter = new InfoFilter(type, NUMBER);

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
      
/*      if( filter.minFVal == Float.MIN_VALUE && 
          filter.maxFVal == Float.MAX_VALUE)
        filters.remove(ID);
      else*/
        filters.put(ID, filter);
    }
  }
  
  class InfoFilter
  {
    private String type;
    private int NUMBER;
    private int minIVal[];
    private int maxIVal[];
    
    private float minFVal[];
    private float maxFVal[];
    
    private InfoFilter(String type, int NUMBER)
    {
      this.type = type;
      this.NUMBER = NUMBER;
      
      if (type.equals("Integer"))
      {
        minIVal = new int[NUMBER];
        maxIVal = new int[NUMBER];
        
        for(int i=0; i<NUMBER; i++)
        {
          minIVal[i] = Integer.MIN_VALUE;
          maxIVal[i] = Integer.MAX_VALUE;
        }
      }
      else if(type.equals("Float"))
      {
        minFVal = new float[NUMBER];
        maxFVal = new float[NUMBER];
        
        for(int i=0; i<NUMBER; i++)
        {
          minFVal[i] = Float.MIN_VALUE;
          maxFVal[i] = Float.MAX_VALUE;
        }
      }
    }
    
    private boolean pass(String valStr[])
    {
      for (int i = 0; i < NUMBER; i++)
      {
        if (type.equals("Integer"))
        {
          int val = Integer.parseInt(valStr[i]);
          if (val < minIVal[i] || val > maxIVal[i])
            return false;
        }
        else if (type.equals("Float"))
        {
          float val = Float.parseFloat(valStr[i]);
          if (val < minFVal[i] || val > maxFVal[i])
            return false;
        }
      }
      return true;
    }
  }
}