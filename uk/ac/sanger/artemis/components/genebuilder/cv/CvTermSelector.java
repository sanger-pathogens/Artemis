package uk.ac.sanger.artemis.components.genebuilder.cv;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Point;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.util.Collections;
import java.util.Vector;
import java.util.regex.Pattern;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTextField;

import org.gmod.schema.cv.CvTerm;

import uk.ac.sanger.artemis.chado.ChadoTransactionManager;
import uk.ac.sanger.artemis.components.database.DatabaseEntrySource;
import uk.ac.sanger.artemis.components.genebuilder.JExtendedComboBox;
import uk.ac.sanger.artemis.util.DatabaseDocument;


class CvTermSelector extends JPanel
                     implements ActionListener
{
  private static final long serialVersionUID = 1L;
  private JExtendedComboBox cvCombo;
  private JExtendedComboBox editableCvCombo;
  private JComboBox goAspect = new JComboBox(new String[]{ "F", "P", "C" });
  private JTextField keyWord = new JTextField(20);
  private JCheckBox ignoreCase = new JCheckBox("Ignore case",false);
  private JExtendedComboBox termList = new JExtendedComboBox(true);
  private JExtendedComboBox evidenceList = new JExtendedComboBox(GoBox.evidenceCodes[1]);
  private GridBagConstraints c = new GridBagConstraints();
  private int termRow = 0;
  private Dimension d = new Dimension(500,termList.getPreferredSize().height);
  private CvTerm cvTerm;
  private Vector<CvTerm> terms;
  private String cv;
  private String evidenceCode;
  
  // components used when adding a new term to the database
  private JTextField newTerm;
  private JTextField definition;
  private JCheckBox addToAnnotation;
  
  public CvTermSelector()
  {
    super(new GridBagLayout());

    setBackground(Color.white);
    createComponentToAddCvTerm();
  }

  /**
   * Add components to the panel for adding CvTerm's to
   * the annotation.
   */
  private void createComponentToAddCvTerm()
  {
    createCvComboBox();
    int row = 0;
    c.gridx = 0;
    c.gridy = row;
    c.anchor = GridBagConstraints.EAST;
    add(new JLabel("CV: "), c);
    
    c.gridx = 1;
    c.anchor = GridBagConstraints.WEST;
    add(cvCombo, c);
    
    c.gridy = ++row;
    c.gridx = 1;
    add(goAspect, c);

    // keyword
    keyWord.setSelectionStart(0);
    keyWord.setSelectionEnd(keyWord.getText().length());
    keyWord.setSelectedTextColor(Color.blue);
    keyWord.setMinimumSize(d);
    keyWord.setPreferredSize(d);
    keyWord.addActionListener(new ActionListener(){
      // carry out search when enter key is pressed
      public void actionPerformed(ActionEvent arg0)
      {
        searchCvTerms();
      }
    });
    c.gridy = ++row;
    c.gridx = 0;
    c.anchor = GridBagConstraints.EAST;
    add(new JLabel("Keywords (or GOid): "),c);
    c.gridx = 1;
    c.anchor = GridBagConstraints.WEST;
    add(keyWord,c);

    c.gridy = ++row;
    c.gridwidth = 1;
    add(ignoreCase,c);
    
    // search button
    c.gridx = 0;
    c.gridy = ++row;
    c.gridwidth = 1;
    JButton search = new JButton("Search");
    search.setActionCommand("SEARCH");
    search.addActionListener(this);
    add(search,c);
    
    // term selection
    termList.setPreferredSize(d);
    termList.setMaximumSize(d);
    c.gridy = ++row;
    c.gridwidth = 2;
    termRow = row;
    add(termList,c);
    
    // evidence
    evidenceList.setSelectedItem("NR \t:: Not Recorded");
    c.gridx = 0;
    c.gridy = ++row;
    c.gridwidth = 2;
    add(evidenceList,c);
  }

  /**
   * Add components to the panel for adding a CvTerm to
   * the database.
   */
  private void createComponentToAddCvTermToDb()
  {
    removeAll();

    final java.util.List<String> cvNames = 
      DatabaseDocument.getCvControledCurationNames();
    
    editableCvCombo = 
      new JExtendedComboBox(new Vector<String>(cvNames));
    editableCvCombo.addItem(JExtendedComboBox.SEPARATOR);
    editableCvCombo.addItem(ChadoTransactionManager.PRODUCT_CV);
    editableCvCombo.setPreferredSize(d);
    
    c.anchor = GridBagConstraints.WEST;
    c.gridx = 0;
    c.gridy = 0;
    c.gridwidth = 1;
    
    add(new JLabel("Controlled vocabulary:"), c);
    
    c.gridx = 1;
    add(new JLabel("Term to add:"), c);
    
    c.gridy = 1;
    c.gridx = 0;
    add(editableCvCombo, c);

    c.gridx = 1;
    newTerm = new JTextField(20);
    newTerm.setMinimumSize(newTerm.getPreferredSize());
    newTerm.setPreferredSize(d);
    add(newTerm, c);
    
    c.gridy = 2;
    c.gridx = 0;
    add(new JLabel("Definition (optional):"), c);
    
    c.gridy = 3;
    c.gridwidth = 2;
    c.fill = GridBagConstraints.HORIZONTAL;
    definition = new JTextField();
    add(definition, c);
    
    c.gridy = 4;
    c.gridwidth = 1;
    addToAnnotation = new JCheckBox("Add to annotation", false);
    add(addToAnnotation, c);
    
    repaint();
    revalidate();
  }

  /**
   * Create the JComboBox of the CV's.
   */
  private void createCvComboBox()
  {
    cvCombo = 
      new JExtendedComboBox(ChadoTransactionManager.CV_NAME);
 
    cvCombo.setActionCommand("CV");
    cvCombo.addActionListener(this);
    final java.util.List<String> cvNames = 
      DatabaseDocument.getCvControledCurationNames();
    cvCombo.addItem(JExtendedComboBox.SEPARATOR);
    for(int i=0; i<cvNames.size(); i++)
      cvCombo.addItem(cvNames.get(i));
    
    cvCombo.addItem(JExtendedComboBox.SEPARATOR);
    cvCombo.addItem(JExtendedComboBox.SEPARATOR);
    
    final String ADD_TERM = "Add term ...";
    cvCombo.addItem(ADD_TERM);
    
    Dimension d = cvCombo.getPreferredSize();
    d = new Dimension(180,(int)d.getHeight());
    cvCombo.setPreferredSize(d);
  }

  public void actionPerformed(ActionEvent e)
  {
    if(e.getActionCommand().equals("CV"))
    {
      if( ((String)cvCombo.getSelectedItem()).equals("Add term ...") )
      {
        createComponentToAddCvTermToDb();
        return;
      }
      
      if( !((String)cvCombo.getSelectedItem()).equals("GO") )
      {
        goAspect.setVisible(false);
        evidenceList.setVisible(false);
      }
      else
      {
        goAspect.setVisible(true);
        evidenceList.setVisible(true);
      }
      termList.removeAllItems();
    }
    else if(e.getActionCommand().equals("SEARCH"))
      searchCvTerms();
  }

  /**
   * Utility to return the name of the CV in the database.
   * @return
   */
  protected String getCvName()
  {
    String cvName = (String) cvCombo.getSelectedItem();
    if (cvName.equals("GO"))
    {
      if (((String) goAspect.getSelectedItem()).equals("F"))
        cvName = "molecular_function";
      else if (((String) goAspect.getSelectedItem()).equals("P"))
        cvName = "biological_process";
      else if (((String) goAspect.getSelectedItem()).equals("C"))
        cvName = "cellular_component";
    }
    else if (cvName.equals("product"))
      cvName = DatabaseDocument.PRODUCTS_TAG_CVNAME;
    else if (cvName.equals("controlled_curation"))
      cvName = DatabaseDocument.CONTROLLED_CURATION_TAG_CVNAME;
    else if (cvName.equals("class"))
      cvName = DatabaseDocument.RILEY_TAG_CVNAME;
    return cvName;
  }
  
  /**
   * Search the CvTerms for the selected CV and keywords.
   */
  private void searchCvTerms()
  {
    String key = keyWord.getText().trim();
    
    if (((String) cvCombo.getSelectedItem()).equals("GO") &&
        Pattern.matches("^\\d{7}$", key))
    {
      terms = getTermsForGoId();
    }
    else
    {
      String cvName = getCvName();
      terms = 
        DatabaseDocument.getCvterms(key, 
                                    cvName, ignoreCase.isSelected());
      Collections.sort(terms, new CvTermsComparator());
    }
    
    remove(termList);
    termList = new JExtendedComboBox(terms, true);
    termList.setPreferredSize(d);
    termList.setMaximumSize(d);
    
    if(((CvTerm)termList.getSelectedItem()).getName().equals("") &&
        termList.getItemCount() > 2)
      termList.setSelectedIndex(1);
    c.gridy = termRow;
    c.gridwidth = 2;
    c.gridx = 0;
    add(termList,c);

    repaint();
    revalidate();
  }
  
  private Vector<CvTerm> getTermsForGoId()
  {
    CvTerm cvTerm = 
      DatabaseDocument.getCvtermFromGoId(keyWord.getText().trim());

    terms = new Vector<CvTerm>(); 
    if(cvTerm != null)
      terms.add(cvTerm);
    return terms;
  }

  /**
   * Show the CvTerm selector in an dialog.
   * @param frame
   * @param doc
   * @return
   */
  protected boolean showOptions(final JFrame frame, final DatabaseDocument doc)
  {
    final JOptionPane optionPane = new JOptionPane(
        this,
        JOptionPane.QUESTION_MESSAGE, JOptionPane.OK_CANCEL_OPTION);

    final JDialog dialog = new JDialog(frame, "CV Term Selector", true);
    dialog.setContentPane(optionPane);
    dialog.setDefaultCloseOperation(JDialog.DO_NOTHING_ON_CLOSE);
    
    optionPane.addPropertyChangeListener(new PropertyChangeListener()
    {
      public void propertyChange(PropertyChangeEvent e)
      {
        String prop = e.getPropertyName();

        if (dialog.isVisible() && (e.getSource() == optionPane)
            && (prop.equals(JOptionPane.VALUE_PROPERTY)))
        {
          int value = ((Integer) optionPane.getValue()).intValue();
          if (value == JOptionPane.OK_OPTION)
          {
            if (editableCvCombo != null)
            {
              cvTerm = addCvTermToDb(doc);
              if (addToAnnotation.isSelected())
              {
                cv = (String) editableCvCombo.getSelectedItem();
                if (DatabaseDocument.PRODUCTS_TAG_CVNAME.equals(cv))
                  cv = "product";
              }
            }
            else
            {
              if (termList == null)
                cvTerm = null;
              else
                cvTerm = getCvTermFromSelectedItem();

              cv = (String) cvCombo.getSelectedItem();
              if (evidenceList != null && evidenceList.getSelectedIndex() > -1)
                evidenceCode = GoBox.evidenceCodes[2][evidenceList.getSelectedIndex()];
            }
          }
          dialog.setVisible(false);
          dialog.dispose();
        }
      }
    });
    dialog.pack();
    centerDialog(dialog);
    dialog.setMinimumSize(getPreferredSize());
    dialog.setVisible(true);

    int value = ((Integer) optionPane.getValue()).intValue();  
    if(value == JOptionPane.OK_OPTION && cvTerm != null)
      return true;
    else
      return false;
  }

  /**
   * Center the given dialog on the screen.
   * @param dialog
   */
  private void centerDialog(JDialog dialog)
  {
    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
    final int x_position =(screen.width - dialog.getSize().width) / 2;
    int y_position =(screen.height - dialog.getSize().height) / 2;
    if(y_position < 10) 
      y_position = 10;

    dialog.setLocation(new Point(x_position, y_position));
  }

  /**
   * Return the CvTerm for a selected item in the list of terms.
   * @return
   */
  private CvTerm getCvTermFromSelectedItem()
  {
    if (termList.getSelectedItem() instanceof String)
    {
      final String selectedStr = (String) termList.getSelectedItem();
      for (int i = 0; i < terms.size(); i++)
      {
        CvTerm thisCvTerm = (CvTerm) terms.get(i);
        if (thisCvTerm.getName().equals(selectedStr))
          return thisCvTerm;
      }
    }
    else
      return (CvTerm) termList.getSelectedItem();
    return null;
  }

  /**
   * Add the Selected term to the database.
   * @param doc
   * @return
   */
  private CvTerm addCvTermToDb(final DatabaseDocument doc)
  {
    final java.util.List<String> cvNames = 
      DatabaseDocument.getCvControledCurationNames();
    final String db;
    if(cvNames.contains((String)editableCvCombo.getSelectedItem()))
      db = ChadoTransactionManager.CONTROLLED_CURATION_DB;
    else
      db = ChadoTransactionManager.PRODUCT_DB;
    
    cvTerm = ChadoTransactionManager.getCvTerm(
        newTerm.getText().trim(), 
        (String)editableCvCombo.getSelectedItem(), 
        definition.getText().trim(), db);
    try
    {
      doc.insertCvTerm(cvTerm);
      return cvTerm;
    }
    catch(RuntimeException re)
    {
      JOptionPane.showMessageDialog(null, re.getMessage(),
          "Problems Writing to Database ",
          JOptionPane.ERROR_MESSAGE);
    }
    return null;
  }

  // getters
  protected CvTerm getCvTerm()
  {
    return cvTerm;
  }

  protected String getCv()
  {
    return cv;
  }
  
  protected String getEvidenceCode()
  {
    return evidenceCode;
  }
  
  public static void main(final String args[])
  {
    DatabaseEntrySource src = new DatabaseEntrySource();
    src.setLocation(true);

    DatabaseDocument doc = new DatabaseDocument(src.getLocation(), src.getPfield());
    doc.getCvTermsByCvName("molecular_function");
    
    CvTermSelector cvTermSelector = new CvTermSelector();
    if(cvTermSelector.showOptions(null, doc))
    {
      System.out.println(cvTermSelector.getCvTerm().getName());
    }

    System.exit(0);
  }
}