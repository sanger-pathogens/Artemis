/* GffPanel.java
 * This file is part of Artemis
 *
 * Copyright (C) 2007  Genome Research Limited
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
 * $Header: //tmp/pathsoft/artemis/uk/ac/sanger/artemis/components/genebuilder/gff/PropertiesPanel.java,v 1.11 2009-08-17 12:50:42 tjc Exp $
 */

package uk.ac.sanger.artemis.components.genebuilder.gff;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Iterator;
import java.util.Set;
import java.util.Vector;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.ButtonGroup;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import javax.swing.border.Border;

import org.gmod.schema.cv.CvTerm;

import uk.ac.sanger.artemis.FeatureChangeEvent;
import uk.ac.sanger.artemis.FeatureChangeListener;
import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.chado.ChadoTransactionManager;
import uk.ac.sanger.artemis.components.genebuilder.GeneUtils;
import uk.ac.sanger.artemis.components.genebuilder.JExtendedComboBox;
import uk.ac.sanger.artemis.editor.MultiLineToolTipUI;
import uk.ac.sanger.artemis.io.ChadoCanonicalGene;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.StringVector;

public class PropertiesPanel extends JPanel
                      implements FeatureChangeListener
{
  private static final long serialVersionUID = 1L;
  private QualifierVector gffQualifiers;
  private JTextField uniquenameTextField;
  private JTextField primaryNameTextField;
  private JCheckBox obsoleteField;
  private JCheckBox partialField5prime;
  private JCheckBox partialField3prime;
  
  private Feature feature;
  /** controls if this panel is automatically closed or open */
  private boolean empty = true;
  /** track if feature isObsolete flag has changed */
  protected boolean obsoleteChanged = false;
  private boolean partialChanged = false;

  private ButtonGroup phaseButtonGroup = null;
  
  private boolean showNames   = true;
  private boolean showParent  = true;
  private boolean showOptions = true;
  private boolean showTimeLastModified = true;
  
  static
  {
    MultiLineToolTipUI.initialize();
  }
  
  public PropertiesPanel(final Feature feature)
  {
    this(feature, true, true, true, true);
  }
  
  public PropertiesPanel(final Feature feature,
                         boolean showNames,
                         boolean showParent,
                         boolean showOptions, 
                         boolean showTimeLastModified)
  {
    super(new FlowLayout(FlowLayout.LEFT));
    this.showNames   = showNames;
    this.showOptions = showOptions;
    this.showParent  = showParent;
    this.showTimeLastModified = showTimeLastModified;
    
    setBackground(Color.WHITE);
    updateFromFeature(feature);
  }
  
  /**
   * Return true if this is a CV qualifier
   * @param qualifier
   * @return
   */
  public static boolean isPropertiesTag(final Qualifier qualifier, final Feature feature)
  {
    if(qualifier.getName().equals("ID") ||
       qualifier.getName().equals("Name") ||
       qualifier.getName().equals("feature_id") ||
       qualifier.getName().equals("Parent") ||
       qualifier.getName().equals("Derives_from") ||
       qualifier.getName().equals("feature_relationship_rank") ||
       qualifier.getName().equals("timelastmodified") ||
       qualifier.getName().equals("isObsolete") ||
       qualifier.getName().equals("Start_range") ||
       qualifier.getName().equals("End_range") ||
       qualifier.getName().equals("codon_start") ||
       ChadoTransactionManager.isSynonymTag(qualifier.getName(), 
           (GFFStreamFeature)feature.getEmblFeature()))
      return true;
    return false;
  }
  
  private Component createGffQualifiersComponent()
  {
    empty = true;
    GridBagConstraints c = new GridBagConstraints();
    JPanel gridPanel = new JPanel(new GridBagLayout());
    gridPanel.setBackground(Color.white);

    if(showNames)
    {
      addNames(c, gridPanel);
      addSynonyms(c, gridPanel);
    }
    
    if(showParent)
    {
      addParent(c, gridPanel, "Parent");
      addParent(c, gridPanel, "Derives_from");
    }
    
    // phase of translation wrt / codon_start
    if(feature.getEntry().getEntryInformation().isValidQualifier(
       feature.getKey(), "codon_start"))
      addPhaseComponent(c, gridPanel);
    else
      phaseButtonGroup = null;

    // partial/obsolete options
    if(showOptions)
      addOptions(c, gridPanel);

    // add buttons and timelastmodified
    if(showTimeLastModified)
      addTimeLastModified(c, gridPanel);
    return gridPanel;
  }
  
  /**
   * Add uniquename and name to the panel.
   * @param c
   * @param gridPanel
   */
  private void addNames(GridBagConstraints c, JPanel gridPanel)
  {
    Qualifier idQualifier   = gffQualifiers.getQualifierByName("ID");
    Qualifier nameQualifier = gffQualifiers.getQualifierByName("Name");
   
    final String uniquename = idQualifier.getValues().get(0);
    uniquenameTextField = new JTextField(uniquename);
    uniquenameTextField.setPreferredSize(calcPreferredMaxTextFieldWidth());
    uniquenameTextField.setCaretPosition(0);
    
    JLabel idField = new JLabel("ID ");
    idField.setFont(getFont().deriveFont(Font.BOLD));
    idField.setHorizontalAlignment(SwingConstants.RIGHT);
    idField.setPreferredSize(calcPreferredLabelWidth());

    c.gridx = 0;
    c.gridy++;
    c.ipadx = 5;
    c.anchor = GridBagConstraints.EAST;
    gridPanel.add(idField, c);
    c.gridx = 1;
    c.ipadx = 0;
    c.anchor = GridBagConstraints.WEST;
    gridPanel.add(uniquenameTextField, c);

    Qualifier featIdQualifier = gffQualifiers.getQualifierByName("feature_id");
    if (featIdQualifier != null)
    {
      Qualifier timeQualifier = gffQualifiers.getQualifierByName("timelastmodified");
      String time = null;
      if (timeQualifier != null)
        time = timeQualifier.getValues().get(0);
      
      String parent = getParentString();
      String tt = "feature_id=" +
        featIdQualifier.getValues().get(0) +
        (parent == null ? "" : "\n"+parent)+
        (time == null ? "" : "\n"+time);

      idField.setToolTipText(tt);
      uniquenameTextField.setToolTipText(tt);
    }
    
    if(feature.getKey().getKeyString().equals(DatabaseDocument.EXONMODEL))
      uniquenameTextField.setEditable(false);

    if (!feature.getKey().getKeyString().equals(DatabaseDocument.EXONMODEL))
    {
      primaryNameTextField = new JTextField();
      primaryNameTextField.setPreferredSize(calcPreferredMaxTextFieldWidth());
      
      if (nameQualifier != null)
      {
        primaryNameTextField.setText((String) nameQualifier.getValues().get(0));
        empty = false;
      }

      c.gridx = 2;
      c.ipadx = 5;
      c.anchor = GridBagConstraints.EAST;
      JLabel lab = new JLabel("  Name");
      lab.setFont(getFont().deriveFont(Font.BOLD));
      gridPanel.add(lab, c);
      c.gridx = 3;
      c.ipadx = 0;
      c.anchor = GridBagConstraints.WEST;
      gridPanel.add(primaryNameTextField, c);
    }
    
    ActionListener addAction = new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        addSynonym();
      }  
    };
    
    AddButton addSynonymButton = new AddButton(addAction, "Add Synonym");
    c.gridx = 4;
    gridPanel.add(addSynonymButton, c);
  }

  /**
   * Get the parent feature name.
   * @return
   */
  private String getParentString()
  {
    Qualifier parentQual = gffQualifiers.getQualifierByName("Parent");
    if(parentQual != null && 
       parentQual.getValues().size() == 1)
      return "Parent: "+parentQual.getValues().get(0);

    Qualifier derivesFromQual = gffQualifiers.getQualifierByName("Derives_from");
    if(derivesFromQual != null && 
       derivesFromQual.getValues().size() == 1)
      return "Derives from: "+derivesFromQual.getValues().get(0);

    return null;
  }
  
  /**
   * Add Parent or Derives_from.
   * @param c
   * @param gridPanel
   */
  private void addParent(GridBagConstraints c,
                         JPanel gridPanel,
                         String parentName)
  {
    Qualifier parentQualifier = gffQualifiers.getQualifierByName(parentName);
    if(parentQualifier != null && 
       parentQualifier.getValues().size() == 1)
    {
      JLabel parentField = new JLabel(parentName);
      parentField.setFont(getFont().deriveFont(Font.BOLD));
      c.gridx = 0;
      c.gridy++;
      c.ipadx = 5;
      c.fill = GridBagConstraints.NONE;
      c.anchor = GridBagConstraints.EAST;
      gridPanel.add(parentField, c);
      
      JTextField parent = new JTextField(" "+ parentQualifier.getValues().get(0));
      parent.setPreferredSize(calcPreferredMaxTextFieldWidth());
      parent.setBorder(BorderFactory.createEmptyBorder());
      c.gridx = 1;
      c.anchor = GridBagConstraints.WEST;
      gridPanel.add(parent, c);
      return;
    }
  }
  
  /**
   * Add synonyms to the panel.
   * @param c
   * @param gridPanel
   * @param nrows
   */
  private void addSynonyms(GridBagConstraints c, JPanel gridPanel)
  {
    for(Qualifier qualifier: gffQualifiers)
    {
      if( ChadoTransactionManager.isSynonymTag(qualifier.getName(), 
          (GFFStreamFeature)feature.getEmblFeature()) &&
          isSystematicId(qualifier.getName()))
      {
        addSynonymComponent(qualifier, c, gridPanel);  
      }
    }
    
    for(Qualifier qualifier: gffQualifiers)
    {
      if( ChadoTransactionManager.isSynonymTag(qualifier.getName(), 
          (GFFStreamFeature)feature.getEmblFeature()) &&
          !isSystematicId(qualifier.getName()))
      {
        addSynonymComponent(qualifier, c, gridPanel);  
      }
    }
  }
  
  
  /**
   * Add partial and obsolete options to the panel.
   * @param c
   * @param gridPanel
   */
  private void addOptions(GridBagConstraints c, JPanel gridPanel)
  {
    final Qualifier isPartialQualfier5;
    final Qualifier isPartialQualfier3;
    if(feature.isForwardFeature())
    {
      isPartialQualfier5 = gffQualifiers.getQualifierByName("Start_range");
      isPartialQualfier3 = gffQualifiers.getQualifierByName("End_range");
    }
    else
    {
      isPartialQualfier3 = gffQualifiers.getQualifierByName("Start_range");
      isPartialQualfier5 = gffQualifiers.getQualifierByName("End_range");
    }

    Box optionsBox = Box.createHorizontalBox();
    partialField5prime = new JCheckBox("partial 5'", 
        ( isPartialQualfier5 != null ) ? true : false);
    Dimension d = calcPreferred(partialField5prime.getPreferredSize().width);
    partialField5prime.setPreferredSize(d);
    partialField5prime.setOpaque(false);
    partialField5prime.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        checkPartial();
      }
    });
    optionsBox.add(partialField5prime);

    partialField3prime = new JCheckBox("partial 3'", 
        ( isPartialQualfier3 != null ) ? true : false);
    partialField3prime.setPreferredSize(d);
    partialField3prime.setOpaque(false);
    partialField3prime.addActionListener(new ActionListener(){
      public void actionPerformed(ActionEvent e)
      {
        checkPartial();
      }
    });
    optionsBox.add(partialField3prime);

    Qualifier obsoleteQual = gffQualifiers.getQualifierByName("isObsolete");
    obsoleteField = new JCheckBox("obsolete", (obsoleteQual == null ? false : 
      Boolean.parseBoolean(obsoleteQual.getValues().get(0))));
    obsoleteField.setPreferredSize(calcPreferred(obsoleteField.getPreferredSize().width));
    obsoleteField.setOpaque(false);
    obsoleteField.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        int result = JOptionPane.showConfirmDialog(
            PropertiesPanel.this, "Change this feature to be "+
            (obsoleteField.isSelected() ? "obsolete!" : "not obsolete!"), 
            "Change obsolete option", JOptionPane.OK_CANCEL_OPTION);
        if(result == JOptionPane.CANCEL_OPTION)
          obsoleteField.setSelected(!obsoleteField.isSelected());
      }
    });
    optionsBox.add(obsoleteField);

    c.gridx = 0;
    c.anchor = GridBagConstraints.WEST;
    c.gridwidth = GridBagConstraints.REMAINDER;
    c.gridy++;
    gridPanel.add(optionsBox, c);
    c.gridwidth = 1;
  }
  
  private void addTimeLastModified(GridBagConstraints c, JPanel gridPanel)
  {
    Qualifier timeQualifier = gffQualifiers.getQualifierByName("timelastmodified");
    if (timeQualifier != null)
    {
      JLabel timeLabel = new JLabel(timeQualifier.getValues().get(0));
      timeLabel.setEnabled(false);
      timeLabel.setToolTipText("time last modified");
      
      c.gridy++;
      c.gridx = 4;
      c.gridwidth = GridBagConstraints.REMAINDER;
      c.anchor = GridBagConstraints.EAST;
      gridPanel.add(timeLabel,c);
    }   
  }
  
  public void updateFromFeature(final Feature feature)
  {
    this.feature = feature;
    removeAll();
    if(gffQualifiers != null)
      feature.removeFeatureChangeListener(this);
    gffQualifiers = feature.getQualifiers().copy();
    
    gffQualifiers = new QualifierVector();
    final QualifierVector qualifiers = feature.getQualifiers();  
    for(Qualifier qualifier: qualifiers) 
    {
      if(isPropertiesTag(qualifier, feature))
        gffQualifiers.addElement(qualifier.copy());
    }
   
    feature.addFeatureChangeListener(this);  
    add(createGffQualifiersComponent());
    revalidate();
    repaint();
  }

  /**
   * Get the latest (edited) property qualifiers
   * @return
   */
  public QualifierVector getGffQualifiers(final Feature feature)
  {
    // check editable components for changes
    Qualifier idQualifier = gffQualifiers.getQualifierByName("ID");
    if(showNames &&
        !idQualifier.getValues().get(0).equals(uniquenameTextField.getText()))
    {
      if(!uniquenameTextField.getText().equals(""))
      {
        final String newName = uniquenameTextField.getText().trim();
        final String oldName = ((String) (idQualifier.getValues().get(0))).trim();
        gffQualifiers.remove(idQualifier);
        idQualifier = new Qualifier("ID", newName);
        gffQualifiers.addElement(idQualifier);
        
        GFFStreamFeature gffFeature = (GFFStreamFeature)feature.getEmblFeature();
        if(gffFeature.getChadoGene() != null)
        {
          Set<uk.ac.sanger.artemis.io.Feature> children = 
        	  gffFeature.getChadoGene().getChildren(gffFeature);
          gffFeature.getChadoGene().updateUniqueName(oldName, newName, children);
        }
      }
    }

    if(!feature.getKey().getKeyString().equals(DatabaseDocument.EXONMODEL))
    {
      Qualifier nameQualifier = gffQualifiers.getQualifierByName("Name");
      if( (nameQualifier != null &&
       !nameQualifier.getValues().get(0).equals(primaryNameTextField.getText())) ||
       (primaryNameTextField != null && !primaryNameTextField.getText().equals("")))
      {
        gffQualifiers.remove(nameQualifier);   
        final String newName = primaryNameTextField.getText().trim();
        nameQualifier = new Qualifier("Name", newName);
        gffQualifiers.addElement(nameQualifier);
      }
    }
    
    if(phaseButtonGroup != null)
    {
      String selectionCmd = phaseButtonGroup.getSelection().getActionCommand();
      Qualifier phaseQualifier = gffQualifiers.getQualifierByName("codon_start");
      if(phaseQualifier == null)
      {
        if(!selectionCmd.equals(""))
        {
          phaseQualifier = new Qualifier("codon_start", selectionCmd);
          gffQualifiers.addElement(phaseQualifier);
        }
      }
      else
      {
        String oldPhase = phaseQualifier.getValues().get(0);
        if(!oldPhase.equals(phaseButtonGroup))
        {
          gffQualifiers.remove(phaseQualifier);
          if(!selectionCmd.equals(""))
          {
            phaseQualifier = new Qualifier("codon_start", selectionCmd);
            gffQualifiers.addElement(phaseQualifier);
          }
        }
      }
    }
    
    Qualifier isObsoleteQualifier = gffQualifiers.getQualifierByName("isObsolete");
    if(isObsoleteQualifier != null)
    {
      if(showOptions)
      {
        String isObsoleteOld = isObsoleteQualifier.getValues().get(0);
        String isObsoleteNew = Boolean.toString(obsoleteField.isSelected());

        if (!isObsoleteNew.equals(isObsoleteOld))
        {
          gffQualifiers.remove(isObsoleteQualifier);
          isObsoleteQualifier = new Qualifier("isObsolete", isObsoleteNew);
          gffQualifiers.addElement(isObsoleteQualifier);
          obsoleteChanged = true;
        }
      } 
    }
    
    Qualifier isPartial5primeQualifier; 
    if(feature.isForwardFeature())
      isPartial5primeQualifier = gffQualifiers.getQualifierByName("Start_range");
    else
      isPartial5primeQualifier = gffQualifiers.getQualifierByName("End_range");
    if(isPartial5primeQualifier != null)
    {
      if(showOptions && !partialField5prime.isSelected())
      {
        gffQualifiers.remove(isPartial5primeQualifier);
        partialChanged = true;
      }
    }
    else if(showOptions && partialField5prime.isSelected())
    {
      if(feature.isForwardFeature())
        gffQualifiers.addElement(new Qualifier("Start_range",".,."));
      else
        gffQualifiers.addElement(new Qualifier("End_range",".,."));
      partialChanged = true;
    }
    
    Qualifier isPartial3primeQualifier;
    if(feature.isForwardFeature())
      isPartial3primeQualifier = gffQualifiers.getQualifierByName("End_range");
    else
      isPartial3primeQualifier = gffQualifiers.getQualifierByName("Start_range");
    if(isPartial3primeQualifier != null)
    {
      if(showOptions && !partialField3prime.isSelected())
      {
        gffQualifiers.remove(isPartial3primeQualifier);
        partialChanged = true;
      }
    }
    else if(showOptions && partialField3prime.isSelected())
    {
      if(feature.isForwardFeature())
        gffQualifiers.addElement(new Qualifier("End_range",".,."));
      else
        gffQualifiers.addElement(new Qualifier("Start_range",".,."));
      partialChanged = true;
    }

    return gffQualifiers;
  }
  
  /**
   * If the partial/isObsolete qualifier for this feature has been changed this 
   * method updates the partial/isObsolete qualifier of the children feature.
   */
  public void updateSettings()
  {
    if(partialChanged)
      updatePartialSettings((GFFStreamFeature)feature.getEmblFeature());
    if(!obsoleteChanged)
      return;
    
    obsoleteChanged = false;
    updateObsoleteSettings((GFFStreamFeature)feature.getEmblFeature());
  }
  
  public static void updateObsoleteSettings(GFFStreamFeature gffFeature)
  {
    Qualifier isObsoleteQualifier = gffFeature.getQualifierByName("isObsolete");
    String isObsoleteNew = isObsoleteQualifier.getValues().get(0);
    gffFeature.setVisible(!isObsoleteNew.equals("true"));
    
    if(gffFeature.getChadoGene() == null)
      return;
    Set<uk.ac.sanger.artemis.io.Feature> children =
    	gffFeature.getChadoGene().getChildren(gffFeature);
 
    if(children.size() > 0)
    {
      Qualifier idQualifier = gffFeature.getQualifierByName("ID");
      int select = JOptionPane.showConfirmDialog(null, 
          "Make children of "+idQualifier.getValues().get(0)+"\n"+
          (isObsoleteNew.equals("true") ? "obsolete?" : "not obsolete?"), 
          "Update Children", 
          JOptionPane.YES_NO_OPTION);
      
      if(select == JOptionPane.YES_OPTION)
      {
        try
        {
          Iterator<uk.ac.sanger.artemis.io.Feature> it = children.iterator();
          while(it.hasNext())
          {
            GFFStreamFeature gffChildFeature = (GFFStreamFeature)it.next();
            Feature f = (Feature)gffChildFeature.getUserData();
            f.setQualifier(new Qualifier("isObsolete", isObsoleteNew));
            gffChildFeature.setVisible(true);    
            
            if(isObsoleteNew.equals("true") ||
               GeneUtils.isHiddenFeature( gffChildFeature.getKey().getKeyString() ))
              gffChildFeature.setVisible(false);
          }
        }
        catch(Exception e)
        {
          e.printStackTrace();
        }
      }    
    }  
  }

  /**
   * Change partial settings of child featuers.
   * @param gffFeature
   */
  private void updatePartialSettings(GFFStreamFeature gffFeature)
  {
    if(!partialChanged)
      return;
    partialChanged = false;
    
    if(gffFeature.getChadoGene() == null)
      return;
    Qualifier fminQualifier = gffFeature.getQualifierByName("Start_range");
    Qualifier fmaxQualifier = gffFeature.getQualifierByName("End_range");
    ChadoCanonicalGene chadoGene = gffFeature.getChadoGene();
    Set<uk.ac.sanger.artemis.io.Feature> children = chadoGene.getChildren(gffFeature);

    if(children.size() > 0)
    {
      int select = JOptionPane.showConfirmDialog(null, 
          "Update partial setting on child features?", 
          "Update Children", JOptionPane.YES_NO_OPTION);
      if(select != JOptionPane.YES_OPTION)
        return;
      try
      {
        Iterator<uk.ac.sanger.artemis.io.Feature> it = children.iterator();
        while(it.hasNext())
        {
          final GFFStreamFeature gffChildFeature = (GFFStreamFeature)it.next();
          final Feature f = (Feature)gffChildFeature.getUserData();
          final String keyStr = f.getKey().getKeyString();
          
          if(keyStr.equals("five_prime_UTR") || keyStr.equals("three_prime_UTR"))
          {
            String fName = chadoGene.getQualifier(gffChildFeature, "ID");
            boolean isFwd = !f.getLocation().isComplement();
            if(fName != null && !chadoGene.isFirstUtr(fName, isFwd))
              continue;
            
            if( (keyStr.equals("five_prime_UTR")  &&  isFwd) ||
                (keyStr.equals("three_prime_UTR") && !isFwd) )
            {
              if(fminQualifier != null)
                f.setQualifier(new Qualifier("Start_range",".,."));
              else
                f.removeQualifierByName("Start_range");
            }
            else
            {
              if(fmaxQualifier != null)
                f.setQualifier(new Qualifier("End_range",".,."));
              else
                f.removeQualifierByName("End_range");
            }
          }
          else
          {
            if(fminQualifier != null)
              f.setQualifier(new Qualifier("Start_range",".,."));
            else
              f.removeQualifierByName("Start_range");
            if(fmaxQualifier != null)
              f.setQualifier(new Qualifier("End_range",".,."));
            else
              f.removeQualifierByName("End_range");
          }
        }
      }
      catch(Exception e)
      {
        e.printStackTrace();
      }  
    }  
  }
  
  /**
   * Partial settings are made on the gene or transcript features.
   * Provide a warning if this is not the case.
   */
  private void checkPartial()
  {
    String keyStr = feature.getKey().getKeyString();
    if(keyStr.equals("polypeptide") || 
       keyStr.equals("CDS") || 
       keyStr.equals("pseudogenic_exon"))
    {
      JOptionPane.showMessageDialog(null, 
          "Partial settings should be updated on the transcript\n"+
          "or gene feature not a "+keyStr+". Please make this change\n"+
          "on the transcript or gene feature now.", 
          "Error", JOptionPane.WARNING_MESSAGE);
    }
  }
  
  private boolean isSystematicId(final String synonymType)
  {
    return (synonymType.indexOf("systematic_id") > -1);
  }
  
  private void removeSynonym(String synonymName, String qualifierValue)
  {   
    int select = JOptionPane.showConfirmDialog(null, 
        "Delete "+qualifierValue+"?",
        "Select synonym type", JOptionPane.OK_CANCEL_OPTION);
    
    if(select != JOptionPane.OK_OPTION)
      return;
    
    StringVector values =
      gffQualifiers.getQualifierByName(synonymName).getValues();

    if(values.size()==1)
      gffQualifiers.removeQualifierByName(synonymName);
    else
    {
      int index = gffQualifiers.indexOfQualifierWithName(synonymName);
      values.remove(qualifierValue);
      gffQualifiers.remove(index);
      gffQualifiers.add(index, new Qualifier(synonymName, values));
    }
    
    removeAll();
    add(createGffQualifiersComponent());
    revalidate();
    repaint();
  }
  
  private void addSynonym()
  {
    final Vector<CvTerm> synonyms = DatabaseDocument.getCvterms("", 
        ChadoTransactionManager.SYNONYM_TAG_CVNAME, false);
    final JExtendedComboBox list = new JExtendedComboBox(synonyms);
    final String options[] = { "CANCEL", "NEXT>"};   
    
    int select = JOptionPane.showOptionDialog(null, list,
        "Select synonym type",
         JOptionPane.YES_NO_CANCEL_OPTION,
         JOptionPane.QUESTION_MESSAGE,
         null, options, options[1]);
    
    if(select == 0)
      return;
    
    Box xBox = Box.createHorizontalBox();
    final String synonymName = ((CvTerm)list.getSelectedItem()).getName();
    final JLabel name = new JLabel( synonymName );
    xBox.add(name);
    
    final JTextField newSynonym = new JTextField(15);
    xBox.add(newSynonym);
    
    final JCheckBox current = new JCheckBox("make current", true);
    xBox.add(current);
    
    select = JOptionPane.showConfirmDialog(null, xBox, 
        "Input name", JOptionPane.OK_CANCEL_OPTION);
    
    if(select == JOptionPane.CANCEL_OPTION || newSynonym.getText().equals(""))
      return;
    
    Qualifier synonymQualifier = gffQualifiers.getQualifierByName(synonymName);
    
    String newSynonymValue = newSynonym.getText();
    if(!current.isSelected())
      newSynonymValue = newSynonymValue + ";current=false";
    
    if(synonymQualifier == null)
    {
      synonymQualifier = new Qualifier(synonymName, newSynonymValue);
      gffQualifiers.add(synonymQualifier);
    }
    else
      synonymQualifier.addValue(newSynonymValue);
    
    final StringVector newValues = synonymQualifier.getValues();
    int index = gffQualifiers.indexOfQualifierWithName(synonymName);
    if(index == -1)
      gffQualifiers.setQualifier(new Qualifier(synonymName,newValues));
    else
    {
      gffQualifiers.remove(index);
      gffQualifiers.add(index, new Qualifier(synonymName,newValues));
    }
    removeAll();
    add(createGffQualifiersComponent());
    revalidate();
  }
  
  /**
   * Add codon_start component to the properties panel
   * @param c
   * @param gridPanel
   */
  private void addPhaseComponent(final GridBagConstraints c, final JPanel gridPanel)
  {
    phaseButtonGroup = new ButtonGroup();

    JRadioButton phaseNone = new JRadioButton("Default", true);
    phaseNone.setOpaque(false);
    phaseNone.setActionCommand("");
    phaseButtonGroup.add(phaseNone);
    int codon_start = feature.getCodonStart();
       
    Box xBox = Box.createHorizontalBox();
    c.gridx = 0;
    c.gridy++;
    c.anchor = GridBagConstraints.EAST;
    c.fill = GridBagConstraints.NONE;
    c.ipadx = 5;
    
    JLabel lab = new JLabel("Codon Start");
    lab.setFont(getFont().deriveFont(Font.BOLD));
    gridPanel.add(lab, c);

    Qualifier qualifierCodonStart = gffQualifiers.getQualifierByName("codon_start");
    for(int i=1; i<4; i++)
    {
      String s = Integer.toString(i);
      JRadioButton phase = new JRadioButton(s);
      phase.setOpaque(false);
      phase.setActionCommand(s);
      phaseButtonGroup.add(phase);
      if(qualifierCodonStart != null && i == codon_start)
      {
        empty = false;
        phase.setSelected(true);
      }
      xBox.add(phase);
    }
    xBox.add(phaseNone);
    xBox.add(Box.createHorizontalGlue());
    c.gridx = 1;
    c.anchor = GridBagConstraints.WEST;
    c.ipadx = 0;
    c.gridwidth = 2;
    gridPanel.add(xBox, c);
    c.gridwidth = 1;
  }
  
  /**
   * Add the synonym components
   * @param qualifier
   * @param c
   * @param gridPanel
   */
  private void addSynonymComponent(final Qualifier qualifier, 
                                   final GridBagConstraints c, 
                                   final JPanel gridPanel)
  {
    empty = false;
    final StringVector values = qualifier.getValues();

    String name = qualifier.getName();
    if(name.equals("previous_systematic_id"))
      name = "prev_sys_id";
    final JLabel sysidField = new JLabel(name+" ");
    sysidField.setFont(getFont().deriveFont(Font.BOLD));
    sysidField.setHorizontalAlignment(SwingConstants.RIGHT);
    sysidField.setPreferredSize(calcPreferredLabelWidth());

    c.gridx = 0;
    c.gridy++;
    c.ipadx = 5;
    c.anchor = GridBagConstraints.EAST;
    gridPanel.add(sysidField, c);

    c.gridx = 1;
    c.anchor = GridBagConstraints.WEST;

    Box synBox = Box.createHorizontalBox();
    for (final String val: values)
    {
      String strs[] = val.split(";");
      JLabel syn = new JLabel(" "+ strs[0] + ";");
      syn.setPreferredSize(calcPreferred(syn.getPreferredSize().width));
      
      if (strs.length > 1 && strs[1].indexOf("current=false") > -1)
        syn.setEnabled(false);
      synBox.add(syn);

      ActionListener removeAction = new ActionListener()
      {
        public void actionPerformed(ActionEvent e)
        {
          removeSynonym(qualifier.getName(), val);
        }  
      };
      synBox.add(new RemoveButton(removeAction));
    }
    c.gridwidth = GridBagConstraints.REMAINDER;

    gridPanel.add(synBox, c);
    c.gridwidth = 1;
  }
  
  public void featureChanged(FeatureChangeEvent event)
  {
    updateFromFeature(event.getFeature());
  }

  public boolean isEmpty()
  {
    return empty;
  }
  
  public void setObsoleteChanged(boolean obsoleteChanged)
  {
    obsoleteField.setSelected(obsoleteChanged);
  }
  
  private Dimension calcPreferredLabelWidth()
  {
    int maxLabelWidth = new JLabel("prev_sys_id ").getPreferredSize().width;
    return calcPreferred(maxLabelWidth);
  }
  
  private Dimension calcPreferredMaxTextFieldWidth()
  {
    int maxLabelWidth = new JLabel("previous_systematic_id       ").getPreferredSize().width;
    return calcPreferred(maxLabelWidth);
  }
  
  private Dimension calcPreferred(int w)
  {
    FontMetrics fm = getFontMetrics(getFont());
    int preferredHeight = fm.getHeight()+fm.getDescent()+4;
    Dimension d = super.getPreferredSize();
    d.height = preferredHeight;
    d.width  = w;
    return d;
  }
  
  protected void makeBorder()
  {
    Border grayline = BorderFactory.createLineBorder(Color.gray);
    setBorder(BorderFactory.createTitledBorder(grayline, 
        feature.getKey().getKeyString()));
  }

}