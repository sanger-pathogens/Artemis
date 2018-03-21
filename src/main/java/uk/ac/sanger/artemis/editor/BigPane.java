/*
 *
 * created: Wed Aug 3 2004
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2000  Genome Research Limited
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

package uk.ac.sanger.artemis.editor;

import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Insets;
import java.awt.Toolkit;
import java.util.Hashtable;

import javax.swing.Box;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JDesktopPane;
import javax.swing.JFrame;
import javax.swing.JInternalFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JRadioButtonMenuItem;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JTabbedPane;
import javax.swing.JToolBar;

import uk.ac.sanger.artemis.Feature;
import uk.ac.sanger.artemis.FeatureVector;
import uk.ac.sanger.artemis.chado.ChadoTransactionManager;
import uk.ac.sanger.artemis.components.QualifierTextArea;
import uk.ac.sanger.artemis.components.genebuilder.cv.CVPanel;
import uk.ac.sanger.artemis.components.genebuilder.ortholog.MatchPanel;
import uk.ac.sanger.artemis.io.ChadoCanonicalGene;
import uk.ac.sanger.artemis.io.EntryInformation;
import uk.ac.sanger.artemis.io.GFFStreamFeature;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierInfo;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.io.StreamQualifier;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.StringVector;


public class BigPane extends JFrame
{
  /** */
  private static final long serialVersionUID = 1L;
  protected static Font font    = new Font("Monospaced",Font.PLAIN,11);
  protected static Font font_sm = new Font("Monospaced",Font.PLAIN,10);
  protected static JCheckBoxMenuItem srsBrowser;
  protected static JCheckBoxMenuItem srsTabPane;
  protected static JCheckBoxMenuItem srsWin;
  protected static JInternalFrame srsFrame;
  protected static JCheckBox addNote = new JCheckBox("Add Note");
 
  public static int CACHE_SIZE = 300;
  public static int MAX_CACHE_SIZE = 1000;
  private QualifierTextArea qualifierTextArea;
  private DataViewInternalFrame dataView;
  //private FeatureVector overlapFeature;
  private Feature edit_feature;
  private JDesktopPane desktop = null;
  private MatchPanel matchForm;
  private CVPanel cvForm;
  private EntryInformation entryInformation;
  
  public BigPane(EntryInformation entryInformation)
  {
    super("Object Editor");
    addWindowListener(new winExit());
    setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
    MultiLineToolTipUI.initialize();
    this.entryInformation = entryInformation;
  }

  public void set(Hashtable dataFile, QualifierTextArea qualifierTextArea,
             FeatureVector overlapFeature, 
             final Feature edit_feature,
             final MatchPanel matchForm,
             final CVPanel cvForm) 
  {
    this.matchForm = matchForm; 
    this.cvForm    = cvForm;
    set(dataFile,qualifierTextArea.getText(),overlapFeature,edit_feature);
    
    this.qualifierTextArea = qualifierTextArea;
  }

  public void set(Hashtable dataFile, String qualifier_txt,
             FeatureVector overlapFeature,
             final Feature edit_feature)
  {
    if(matchForm != null || cvForm != null)
    {
      QualifierVector matchQualifiers = matchForm.getMatchQualifiers();
      
      if(matchQualifiers != null)
        qualifier_txt = getQualifierString(matchQualifiers);
      Qualifier productQualifier = cvForm.getCvQualifiers().getQualifierByName("product");
      
      if(productQualifier != null)
      {
        QualifierVector qv = new QualifierVector();
        qv.addQualifierValues(productQualifier);
        qualifier_txt = qualifier_txt + getQualifierString(qv) + "\n";
      }
    }
    
    
    //this.overlapFeature = overlapFeature;
    this.edit_feature   = edit_feature;
    addNote.setSelected(false);
    
    setFont(font);

    if(desktop == null)
    {
      desktop = new JDesktopPane();
      desktop.setDragMode(JDesktopPane.LIVE_DRAG_MODE);
      getContentPane().add(desktop);
    }

    //Make the big window be indented 80 pixels from each edge
    //of the screen.
    final int inset = 80;
    final Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
    
    if(screenSize.width > 1280)
      screenSize.width = 1280;
    
    setBounds(inset, inset,
              screenSize.width  - inset*2,
              screenSize.height - inset*2);

    final JScrollPane scrollEvidence = new JScrollPane();
    // data set
    final int hgt = getHeight()-85;
    final int wid = getWidth()/2-10;
    
    dataView = new DataViewInternalFrame(dataFile,desktop, scrollEvidence,
                                         wid,hgt,qualifier_txt,edit_feature);
    dataView.setLocation(5,0);
    dataView.setSize(wid,hgt);
    dataView.setVisible(true);
    desktop.add(dataView);

    // evidence
    final JInternalFrame evidence = new JInternalFrame("Evidence", true,
                                                      true, true, true);

    Box evidenceBox = dataView.getEvidenceBox();
    if(overlapFeature != null)
      evidenceBox.add(getOverlapFeatures(overlapFeature,desktop),0);
    evidenceBox.add(Box.createVerticalGlue());

    ScrollPanel scroller = new ScrollPanel();
    scroller.add(evidenceBox);
    scrollEvidence.setViewportView(scroller);
    evidence.getContentPane().add(scrollEvidence);
    evidence.setLocation(wid+10,0);
    evidence.setSize(wid,hgt);
    evidence.setVisible(true);
    desktop.add(evidence);   

    final JMenuBar menuBar = createMenuBar(desktop);
    setJMenuBar(menuBar);

    // toolbar
    final JToolBar toolBar = createToolbar();
    getContentPane().add(toolBar,BorderLayout.NORTH);

    setVisible(true);
  }

  /**
   *  Return a string containing one qualifier per line.  These are the
   *  original qualifiers, not the qualifiers from the qualifier_text_area.
   **/
  private String getQualifierString(QualifierVector qualifiers) 
  {
    final StringBuffer buffer = new StringBuffer();
    
    for(int qualifier_index = 0; qualifier_index < qualifiers.size();
        ++qualifier_index) 
    {
      final Qualifier this_qualifier = (Qualifier)qualifiers.elementAt(qualifier_index);


      final QualifierInfo qualifier_info =
          entryInformation.getQualifierInfo(this_qualifier.getName());

      final StringVector qualifier_strings =
                       StreamQualifier.toStringVector(qualifier_info, this_qualifier);

      for(int value_index = 0; value_index < qualifier_strings.size();
          ++value_index)
      {
        final String qualifier_string = (String)qualifier_strings.elementAt(value_index);
        buffer.append(qualifier_string + "\n");
      }
    }

    return buffer.toString();
  }

  /**
  *
  * Display for overlapping Pfam features.
  *
  */
  private Box getOverlapFeatures(FeatureVector overlapFeature,
                                 JDesktopPane desktop)
  {
    Box bdown = Box.createVerticalBox();
    bdown.add(new EvidenceViewer(edit_feature,overlapFeature,desktop));
    return bdown;
  }

  /**
  *
  * Create a toolbar
  * @return toolbar.
  *
  */
  private JToolBar createToolbar()
  {
    JToolBar toolBar = new JToolBar();
    
    JButton applyButt = new JButton("APPLY");
    applyButt.setToolTipText("Apply annotation changed to feature editor");
    applyButt.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        transferAnnotation(false);
      }
    });
    applyButt.setBackground(new Color(0,0,81));
    applyButt.setForeground(Color.red);
    applyButt.setBorderPainted(false);
    applyButt.setMargin(new Insets(0,0,0,0));
    applyButt.setFont(font);
    toolBar.add(applyButt);
    
    JButton closeButt = new JButton("CLOSE");
    closeButt.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        onClose();
      }
    });
    closeButt.setBackground(new Color(0,0,81));
    closeButt.setBorderPainted(false);
    closeButt.setMargin(new Insets(0,0,0,0));
    closeButt.setFont(font);
    toolBar.add(closeButt);

    addNote.addItemListener(new ItemListener()
    {
      public void itemStateChanged(ItemEvent e)
      {
        if(addNote.isSelected())
          dataView.updateNote();
        else
          dataView.deleteNote();
      }
    });
    toolBar.add(addNote);
    
    return toolBar;
  }

  /**
   * Transfer annotation data back to the feature editor
   */
  private void transferAnnotation(final boolean askToUpdate)
  {
    final String dataViewStringOriginal = dataView.getFeatureText().trim();
    
    ////////
    ////////  standard text area feature editor
    
    if(matchForm == null && cvForm == null)
    {
      final String oldTxt = qualifierTextArea.getText().trim();

      // changes have been made to feature annotation
      if(!oldTxt.equals(dataViewStringOriginal))
      {
        if(askToUpdate)
        {
          if(askToUpdate())
            qualifierTextArea.setText(dataViewStringOriginal);
        }
        else
          qualifierTextArea.setText(dataViewStringOriginal);
      }
      return;
    }
    
    ////////
    ////////  tabbed feature editor
    
    if(askToUpdate && !askToUpdate())
      return;

    QualifierVector matchQualifiers= new QualifierVector();
    Qualifier productQualifier = null;
    String otherString = null;
    
    ChadoCanonicalGene chadoGene =
      ((GFFStreamFeature)edit_feature.getEmblFeature()).getChadoGene();
    
    
    final StringVector v = StringVector.getStrings(dataViewStringOriginal, "\n");
    for(int i=0; i<v.size(); i++)
    {
      String value = (String)v.get(i);
      if(MatchPanel.isMatchTag(value) ||
         value.startsWith("/product"))
      {
        int index = value.indexOf('=');
        String key;
        if(index > -1)
        {
          key   = value.substring(1, index);
          value = value.substring(index+1);
        }
        else
        {
          key   = value;
          value = null;
        }

        if(value.startsWith("\""))
          value = value.substring(1);
        if(value.endsWith("\""))
          value = value.substring(0, value.length()-1);

        if(MatchPanel.isMatchTag(key))
        {
          Qualifier qualifier = matchQualifiers.getQualifierByName(key);
          if(qualifier == null)
            qualifier = new Qualifier(key);
          qualifier.addValue(value);
          matchQualifiers.setQualifier(qualifier);
        }
        else
        {
          if(productQualifier == null)
            productQualifier = new Qualifier("product");
          
          if( productQualifier.getValues() == null ||
              !containsProductValue(productQualifier.getValues(),value))
          {
            if(value.startsWith("term="))
              productQualifier.addValue(value);
            else
            {
              org.gmod.schema.cv.CvTerm cvTerm = 
              DatabaseDocument.getCvTermByCvPartAndCvTerm(value, ChadoTransactionManager.PRODUCT_CV);
            
              if(cvTerm == null)
              {
                int val = JOptionPane.showConfirmDialog(this, 
                  "This term is missing from the database:\n"+
                  value+"\nAdd this to the database?", 
                  "Add Product",JOptionPane.OK_CANCEL_OPTION);
                if(val == JOptionPane.CANCEL_OPTION)
                  return;
              }
              productQualifier.addValue("term="+value+";");
            }
          }
        }
      }
      else if(value.startsWith("/gene=") && chadoGene != null)
      {
        // Added as the Name qualifier on the gene feature, also
        // added to the peptide with first letter made uppercase
        String nameGene = value.substring(6);
        if(nameGene.startsWith("\""))
          nameGene = nameGene.substring(1);
        if(nameGene.endsWith("\""))
          nameGene = nameGene.substring(0, nameGene.length()-1);
        
        String namePep =
          Character.toUpperCase(nameGene.charAt(0)) + nameGene.substring(1);

        Qualifier qualifierGene = new Qualifier("Name", nameGene);
        Qualifier qualifierPep = new Qualifier("Name", namePep);
        try
        {
          ((Feature)chadoGene.getGene().getUserData()).setQualifier(qualifierGene);
          edit_feature.setQualifier(qualifierPep);
        }
        catch (Exception e)
        {
          e.printStackTrace();
        }
      }
      else
      {
        if(otherString == null)
          otherString = new String();
        otherString = otherString+value+"\n";
      }
    }
    matchForm.updateFromQualifiers(matchQualifiers,edit_feature);
    
    if(productQualifier != null)
    {
      QualifierVector cvQualifiers = cvForm.getCvQualifiers().copy();
      int productIndex = cvQualifiers.indexOfQualifierWithName("product");
      if(productIndex > -1)
        cvQualifiers.remove(productIndex);
      else
        productIndex = 0;

      cvQualifiers.add(productIndex, productQualifier);
      cvForm.updateFromQualifiers(cvQualifiers);
    }
    
    if(otherString != null)
    {
      otherString = qualifierTextArea.getText() + 
                    (qualifierTextArea.getText().endsWith("\n") ? "" : "\n") +
                    otherString;
      qualifierTextArea.setText(otherString);
    } 
  }
  
  /**
   * Check the values of a StringVector to see if it already contains
   * a new value.
   * @param values
   * @param newValue
   * @return
   */
  private boolean containsProductValue(StringVector values, String newValue)
  {
    for(int i=0; i<values.size(); i++)
    {
      String thisValue = (String)values.get(i);
      if(thisValue.startsWith("term="))
      {
        thisValue = thisValue.substring(5);
        if(thisValue.endsWith(";"))
          thisValue = thisValue.substring(0, thisValue.length()-1);
      }
      if(thisValue.equals(newValue))
        return true;
    }
    return false;
  }
  
  private boolean askToUpdate()
  {
    int ok = JOptionPane.showConfirmDialog(BigPane.this,
                        "Apply changes now?",
                        "Apply Changes",
                        JOptionPane.YES_NO_CANCEL_OPTION,
                        JOptionPane.QUESTION_MESSAGE);
    if(ok == JOptionPane.OK_OPTION)
      return true;
    return false;
  }
  
  /**
  *
  * Create a menu bar. 
  * @param desktop pane.
  * @return menu bar.
  *
  */
  private JMenuBar createMenuBar(final JDesktopPane desktop)
  {
    JMenuBar menuBar = new JMenuBar();
    JMenu fileMenu = new JMenu("File");
    menuBar.add(fileMenu);

    final JMenuItem reReadMenu = new JMenuItem("Re-read selected results");
    reReadMenu.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        dataView.reReadSelectedResults();
      }
    });
    fileMenu.add(reReadMenu);
    fileMenu.add(new JSeparator());
 
    JMenuItem applyMenu = new JMenuItem("Apply to Feature Editor");
    applyMenu.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        transferAnnotation(false);
      }
    });
    fileMenu.add(applyMenu);
    fileMenu.add(new JSeparator());
//
    final JMenuItem exitMenu = new JMenuItem("Close");
    exitMenu.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        onClose();
      }
    });
        
    fileMenu.add(exitMenu);

   //srs menu items
    final JMenu optionMenu = new JMenu("Options");
    menuBar.add(optionMenu);

    final JMenu srsMenu = new JMenu("Show SRS in");
    optionMenu.add(srsMenu);
    
    srsBrowser = new JCheckBoxMenuItem("Browser",false);
    srsTabPane = new JCheckBoxMenuItem("Tab Pane",true);
    srsWin     = new JCheckBoxMenuItem("New Window",false);

    srsMenu.add(srsBrowser);
    srsMenu.add(srsTabPane);
    srsMenu.add(srsWin);

   //drag mode
    final JMenu dragMenu = new JMenu("Drag Mode");
    optionMenu.add(dragMenu);

    JRadioButtonMenuItem liveDrag = new JRadioButtonMenuItem("Live",true);
    liveDrag.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        desktop.setDragMode(JDesktopPane.LIVE_DRAG_MODE);
      }
    });

    dragMenu.add(liveDrag);

    final JRadioButtonMenuItem outlineDrag = new JRadioButtonMenuItem("Outline",false);
    outlineDrag.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        desktop.setDragMode(JDesktopPane.OUTLINE_DRAG_MODE);
      }
    });
    dragMenu.add(outlineDrag);
    ButtonGroup buttGroup = new ButtonGroup();
    buttGroup.add(liveDrag);
    buttGroup.add(outlineDrag);

    // cache size
    /*final JMenuItem cacheMenu = new JMenuItem("Cache Size");
    cacheMenu.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        JTextField cacheSize = new JTextField(Integer.toString(CACHE_SIZE));
        int select = JOptionPane.showConfirmDialog(BigPane.this, 
                                 cacheSize, "Cache Size", 
                                 JOptionPane.OK_CANCEL_OPTION);
        if(select == JOptionPane.CANCEL_OPTION)
          return;
        
        CACHE_SIZE = Integer.parseInt(cacheSize.getText());
        
        HitInfo[] cacheHits = FastaTextPane.cacheHits;
        
        FastaTextPane.cacheHits = new HitInfo[CACHE_SIZE];
        for(int i=0; i<cacheHits.length; i++)
        {
          if(i >= FastaTextPane.cacheHits.length)
            break;
          FastaTextPane.cacheHits[i] = cacheHits[i];
        }
      }
    });
    optionMenu.add(cacheMenu);*/
    return menuBar;
  }

  /**
  *
  * Routine to call when the editor is closed.
  *
  */
  private void onClose()
  {
    // remember the splitpane divider locations
    dataView.setDataDividerLocation();
    dataView.setAnnotationDividerLocation();
    
    transferAnnotation(true);
    // stop getz processes
    setVisible(false);
    dataView.stopGetz();
    dataView.dispose();
    BigPane.srsFrame = null;
    dispose();
  }

  /**
  *
  * Set up the tabbed SRS frame
  *
  */ 
  protected static void setUpSRSFrame(int hgt, JDesktopPane desktop)
  {
    BigPane.srsFrame = new JInternalFrame("SRS",
                                           true, //resizable
                                           true, //closable
                                           true, //maximizable
                                           true);//iconifiable
    BigPane.srsFrame.setDefaultCloseOperation(JInternalFrame.HIDE_ON_CLOSE);
    
    BigPane.srsFrame.setLocation(0,0);
    BigPane.srsFrame.setSize(500,hgt);
    final JTabbedPane jtab = new JTabbedPane();
    BigPane.srsFrame.getContentPane().add(jtab);

    final JMenuBar menuBar = new JMenuBar();
    final CommonMenu cmen = new CommonMenu(BigPane.srsFrame);
    menuBar.add(cmen);
    final JMenuItem closeTabMenu = new JMenuItem("Close tab");
    cmen.add(closeTabMenu);
    closeTabMenu.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        int select = jtab.getSelectedIndex(); 
        if(select > -1)
          jtab.removeTabAt(select);
      }
    });

    BigPane.srsFrame.setJMenuBar(menuBar);
    desktop.add(BigPane.srsFrame);
  }

  /**
  *
  * Extends WindowAdapter to close window.
  *
  */
  class winExit extends WindowAdapter
  {
     public void windowClosing(WindowEvent we)
     {
       onClose();
     }
  }

}
