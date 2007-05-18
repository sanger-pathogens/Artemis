/* Similarityox.java
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
 */

package uk.ac.sanger.artemis.components.genebuilder.ortholog;

import java.awt.Color;
import java.awt.Component;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Vector;

import javax.swing.DefaultCellEditor;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.SizeRequirements;
import javax.swing.UIManager;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableColumn; 
import javax.swing.table.TableModel;

import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.util.StringVector;

public class SimilarityTable extends AbstractMatchTable
{ 
  private int NUMBER_COLUMNS = 12;
  private Vector rowData   = new Vector();
  private Vector tableData = new Vector(NUMBER_COLUMNS);
  private JTable similarityTable;
  private JButton infoLevelButton = new JButton("Verbose");
  private Qualifier origQualifier;
  private boolean isChanged = false;
  
  //
  // column headings
  final static String ORGANISM_COL = "Organism";
  final static String HIT_COL = "Hit";
  final static String DESCRIPTION_COL = "Description";
  final static String EVALUE_COL = "E-value";
  final static String LENGTH_COL = "Length";
  final static String ID_COL = "ID";
  final static String QUERY_COL = "Query";
  final static String SUBJECT_COL = "Subject";
  final static String SCORE_COL = "Score";
  final static String OVERLAP_COL = "Overlap";
  final static String METHOD_COL = "Method";
  final static String REMOVE_BUTTON_COL = "";
  
  
  /**
   * Contruct a component for a similarity line
   * @param similarity
   * @param similarityString
   */
  protected SimilarityTable(final Qualifier simQualifier)
  {
    this.origQualifier = simQualifier;
    
    infoLevelButton.setOpaque(false);
    tableData.setSize(NUMBER_COLUMNS);
    
    tableData.setElementAt(ORGANISM_COL,0);
    tableData.setElementAt(HIT_COL,1);
    tableData.setElementAt(DESCRIPTION_COL,2);
    tableData.setElementAt(EVALUE_COL,3);
    tableData.setElementAt(LENGTH_COL,4);
    tableData.setElementAt(ID_COL,5);
    tableData.setElementAt(QUERY_COL,6);
    tableData.setElementAt(SUBJECT_COL,7);
    tableData.setElementAt(SCORE_COL,8);
    tableData.setElementAt(OVERLAP_COL,9);
    tableData.setElementAt(METHOD_COL,10);
    tableData.setElementAt(REMOVE_BUTTON_COL,11);
    
    
    // add rows of similarity
    StringVector sims = simQualifier.getValues();
    for(int i=0; i<sims.size(); i++)
      rowData.add(getRowData((String)sims.get(i), tableData));
    
    
    similarityTable = new JTable(rowData, tableData);
    
    TableColumn col = similarityTable.getColumn(LENGTH_COL);
    col.setPreferredWidth(30);
    col = similarityTable.getColumn(EVALUE_COL);
    col.setPreferredWidth(40);
    col = similarityTable.getColumn(ID_COL);
    col.setPreferredWidth(25);
    
    final TableColumn[] hideColumns = new TableColumn[5];
    hideColumns[0] = similarityTable.getColumn(QUERY_COL);
    hideColumns[1] = similarityTable.getColumn(SUBJECT_COL);
    hideColumns[2] = similarityTable.getColumn(SCORE_COL);
    hideColumns[3] = similarityTable.getColumn(OVERLAP_COL);
    hideColumns[4] = similarityTable.getColumn(METHOD_COL);
    
    for(int i=0; i<hideColumns.length; i++)
    {
      hideColumns[i].setMinWidth(0);
      hideColumns[i].setMaxWidth(0);
      hideColumns[i].setPreferredWidth(0);
    }
    
    infoLevelButton.addActionListener(new ActionListener()
    {
      private SizeRequirements size = new SizeRequirements(
          20, 40, 180, 0.f);
      public void actionPerformed(ActionEvent e)
      {
        // get the current size of the column
        SizeRequirements temp = new SizeRequirements();
        
        temp.preferred = hideColumns[0].getPreferredWidth();
        temp.minimum = hideColumns[0].getMinWidth();
        temp.maximum = hideColumns[0].getMaxWidth();
        
        // change the column size to the old size 
        for(int i=0; i<hideColumns.length; i++)
        {
          hideColumns[i].setMinWidth(size.minimum);
          hideColumns[i].setMaxWidth(size.maximum);
          hideColumns[i].setPreferredWidth(size.preferred);
        }
        // save the old size
        size = temp;
        
        if(infoLevelButton.getText().equals("Verbose"))
          infoLevelButton.setText("Brief");
        else
          infoLevelButton.setText("Verbose");
      } 
    });
    
    TableModel tableModel = similarityTable.getModel();
    // remove button column
    col = similarityTable.getColumn(REMOVE_BUTTON_COL);
    col.setMinWidth(8);
    col.setMaxWidth(40);
    col.setPreferredWidth(40);

    final SimilarityRenderer renderer = new SimilarityRenderer();
    TableColumn tc;

    for(int columnIndex = 0; columnIndex <tableModel.getColumnCount();
        columnIndex++) 
    {
      tc = similarityTable.getColumnModel().getColumn(columnIndex);
      tc.setCellRenderer(renderer);
      tc.setCellEditor(new CellEditing(new JTextField()));
    }
    
    // remove JButton column
    tc = similarityTable.getColumn(REMOVE_BUTTON_COL);
    tc.setCellEditor(new ButtonEditor(new JCheckBox(),
        (DefaultTableModel)similarityTable.getModel()));
  }
  
  /**
   * Build a vector of the row data
   * @param similarityString
   * @return
   */
  private Vector getRowData(final String similarityString,
                            final Vector tableData)
  {
    Vector row = new Vector(NUMBER_COLUMNS);
    row.setSize(NUMBER_COLUMNS);
    System.out.println(similarityString);
    StringVector sim = StringVector.getStrings(similarityString, ";");
    
    // organism
    if(sim.size() >=3)
    {
      int columnIndex = tableData.indexOf(ORGANISM_COL);
      row.setElementAt(((String)sim.get(2)).trim(), columnIndex);
    }
    
    // hit
    if(sim.size() >=2)
    {
      int columnIndex = tableData.indexOf(HIT_COL);
      row.setElementAt(((String)sim.get(1)).trim(), columnIndex);
    }
    
    // description
    if(sim.size() >=4)
    {
      int columnIndex = tableData.indexOf(DESCRIPTION_COL);
      row.setElementAt(((String)sim.get(3)).trim(), columnIndex);
    }
    
    // e-value
    String evalueString;
    if( !(evalueString=getField("E()=", similarityString)).equals("") )
    {
      int columnIndex = tableData.indexOf(EVALUE_COL);
      row.setElementAt(evalueString, columnIndex);
    }
    
    // length
    String lenString;
    if( !(lenString=getField("length", similarityString).trim()).equals("") )
    {
      int columnIndex = tableData.indexOf(LENGTH_COL);
      row.setElementAt(lenString, columnIndex);
    }
    
    String ungappedId;
    if( !(ungappedId=getField("ungapped id=", similarityString)).equals("") )
    {
      int columnIndex = tableData.indexOf(ID_COL);
      row.setElementAt(ungappedId, columnIndex);
    }
    
    String query;
    if( !(query=getField("query", similarityString).trim()).equals("") )
    {
      int columnIndex = tableData.indexOf(QUERY_COL);
      row.setElementAt(query, columnIndex);
    }
    
    String subject;
    if( !(subject=getField("subject", similarityString).trim()).equals("") )
    {
      int columnIndex = tableData.indexOf(SUBJECT_COL);
      row.setElementAt(subject, columnIndex);
    }
    
    String score;
    if( !(score=getField("score=", similarityString)).equals("") )
    {
      int columnIndex = tableData.indexOf(SCORE_COL);
      row.setElementAt(score, columnIndex);
    }
    
    String overlap;
    if( !(overlap=getField("overlap=", similarityString)).equals("") )
    {
      int columnIndex = tableData.indexOf(OVERLAP_COL);
      row.setElementAt(overlap, columnIndex);
    }
    
    int columnIndex = tableData.indexOf(METHOD_COL);
    row.setElementAt(((String)sim.get(0)).trim(), columnIndex);
    return row;
  }
  
  /**
   * Get the column index from the column name
   * @param columnName
   * @return
   */
  private int getColumnIndex(final String columnName)
  {
    int modelColumnIndex = getSimilarityTable().getColumn(columnName).getModelIndex();
    return getSimilarityTable().convertColumnIndexToView(modelColumnIndex);
  }


  /**
   * Get the JTable of similarity data
   * @return
   */
  protected JTable getSimilarityTable()
  {
    return similarityTable;
  }
  
  /**
   * Button to show/hide columns
   * @return
   */
  public JButton getInfoLevelButton()
  {
    return infoLevelButton;
  }
 
  protected boolean isQualifierChanged()
  {
    return isChanged;
  }

  protected void updateQualifier(QualifierVector qv)
  {
    StringVector values = origQualifier.getValues();
    values.removeAllElements();
    
    if(similarityTable.getRowCount() < 1)
      return;
    
    System.out.println("\nHERE:\n");
    for(int i=0; i<similarityTable.getRowCount(); i++)
    {
      String updatedQualifierString = updateQualifierString(i);
      values.add(updatedQualifierString);
      System.out.println(updatedQualifierString);
    }
    System.out.println("\n\n");
    
    int index = qv.indexOfQualifierWithName(origQualifier.getName());
    origQualifier = new Qualifier(origQualifier.getName(), values);
    qv.remove(index);
    qv.add(index, origQualifier);
  }
  
  private String updateQualifierString(final int row)
  {
    StringBuffer similarityStr = new StringBuffer(
             (String)similarityTable.getValueAt(row, getColumnIndex(METHOD_COL)) );   // method
    similarityStr.append(";");
    similarityStr.append(
             (String)similarityTable.getValueAt(row, getColumnIndex(HIT_COL)) );      // hit
    similarityStr.append(";");
    similarityStr.append(
             (String)similarityTable.getValueAt(row, getColumnIndex(ORGANISM_COL)) ); // organism
    similarityStr.append(";");
    similarityStr.append(
             (String)similarityTable.getValueAt(row, getColumnIndex(DESCRIPTION_COL)) ); // description
    similarityStr.append(";");
    similarityStr.append("length "+
             (String)similarityTable.getValueAt(row, getColumnIndex(LENGTH_COL)) ); // length
    similarityStr.append(";");
    similarityStr.append("E()="+
             (String)similarityTable.getValueAt(row, getColumnIndex(EVALUE_COL)) ); // evalue
    similarityStr.append(";");
    similarityStr.append("score="+
             (String)similarityTable.getValueAt(row, getColumnIndex(SCORE_COL)) );  // score
    similarityStr.append(";");
    similarityStr.append("query "+
             (String)similarityTable.getValueAt(row, getColumnIndex(QUERY_COL)) );  // query
    similarityStr.append(";");
    similarityStr.append("subject "+
             (String)similarityTable.getValueAt(row, getColumnIndex(SUBJECT_COL)) ); // subject
    similarityStr.append(";");
    similarityStr.append("ungapped id="+
             (String)similarityTable.getValueAt(row, getColumnIndex(ID_COL)) ); // ungapped id
    similarityStr.append(";");
    similarityStr.append("overlap="+
             (String)similarityTable.getValueAt(row, getColumnIndex(OVERLAP_COL)) ); // overlap
    
    return similarityStr.toString();
  }
  
  /**
   * Renderer for the Similarity cells
   */
  public class SimilarityRenderer extends DefaultTableCellRenderer
  {  
    /** */
    private static final long serialVersionUID = 1L;
    private int minHeight = -1;
    
    private final JTextArea hitTextArea = new JTextArea();
    private final JLabel evalue = new JLabel();
    private final JLabel length = new JLabel();
    private final JTextArea organismTextArea = new JTextArea();
    private final JTextArea descriptionTextArea = new JTextArea();
    private final JLabel ungappedId = new JLabel();
    private final JLabel queryCoord = new JLabel();
    private final JLabel subjCoord  = new JLabel();
    private final JLabel score      = new JLabel();
    private final JLabel overlap    = new JLabel();
    private final JLabel method     = new JLabel();
    private final JButton buttRemove = new JButton("X");
    private Color fgColor = new Color(139,35,35);
    
    public SimilarityRenderer() 
    {
      organismTextArea.setLineWrap(true);
      organismTextArea.setWrapStyleWord(true);
      
      hitTextArea.setLineWrap(true);
      hitTextArea.setWrapStyleWord(true);
      
      descriptionTextArea.setLineWrap(true);
      descriptionTextArea.setWrapStyleWord(true);
      
      buttRemove.setOpaque(false);
      buttRemove.setText("X");
      
      Font font = getFont().deriveFont(Font.BOLD);
      buttRemove.setFont(font);
      buttRemove.setToolTipText("REMOVE");
      
      buttRemove.setPreferredSize(new Dimension(
            20,20));
      buttRemove.setMaximumSize(buttRemove.getPreferredSize());
    }

    public Component getTableCellRendererComponent(
        JTable table,
        Object value,
        boolean isSelected,
        boolean hasFocus,
        final int row,
        final int column) 
    {
      Component c = null;
      String text = null;
      if(value != null)
        text = (String)value;
      Dimension dim;

      TableColumn tableCol;
      if(column == getColumnIndex(ORGANISM_COL))
      {
        organismTextArea.setText(text);

        tableCol = table.getColumnModel().getColumn(column);
        organismTextArea.setSize(tableCol.getWidth(), table.getRowHeight(row));

        dim = organismTextArea.getPreferredSize();
        minHeight = Math.max(minHeight, dim.height);
        
        c = organismTextArea;
      }
      else if(column == getColumnIndex(HIT_COL))
      {
        hitTextArea.setText(text);
        
        tableCol = table.getColumnModel().getColumn(column);
        hitTextArea.setSize(tableCol.getWidth(), table.getRowHeight(row));

        dim = hitTextArea.getPreferredSize();
        minHeight = Math.max(minHeight, dim.height);
        
        c = hitTextArea;
      }
      else if(column == getColumnIndex(DESCRIPTION_COL))
      {
        descriptionTextArea.setText(text);

        tableCol = table.getColumnModel().getColumn(column);
        descriptionTextArea.setSize(tableCol.getWidth(), table
            .getRowHeight(row));

        dim = descriptionTextArea.getPreferredSize();
        minHeight = Math.max(minHeight, dim.height);
        c = descriptionTextArea;
      }
      else if(column == getColumnIndex(EVALUE_COL))
      {
        evalue.setText(text);
        c = evalue;
      }
      else if(column == getColumnIndex(LENGTH_COL))
      {
        length.setText(text);
        c = length;
      }
      else if(column == getColumnIndex(ID_COL))
      {
        ungappedId.setText(text);
        c = ungappedId;
      }
      else if(column == getColumnIndex(QUERY_COL))
      {
        queryCoord.setText(text);
        c = queryCoord;
      }
      else if(column == getColumnIndex(SUBJECT_COL))
      {
        subjCoord.setText(text);
        c = subjCoord;
      }
      else if(column == getColumnIndex(SCORE_COL))
      {
        score.setText(text);
        c = score;
      }
      else if(column == getColumnIndex(OVERLAP_COL))
      {
        overlap.setText(text);
        c = overlap;
      }
      else if(column ==getColumnIndex(METHOD_COL))
      {
        method.setText(text);
        c = method;
      }
      else if(column == getColumnIndex(REMOVE_BUTTON_COL))
      {
        if(isSelected) 
        {
          buttRemove.setForeground(fgColor);
          buttRemove.setBackground(table.getSelectionBackground());
          
        } 
        else
        {
          buttRemove.setForeground(fgColor);
          buttRemove.setBackground(UIManager.getColor("Button.background"));
        }
        c = buttRemove;
      }
      else
      {
        throw new RuntimeException("invalid column! " + column +
                                   " " + text);
      }

      // adjust row height for columns with multiple lines
      if(column < 3)
      {
        if(table.getRowHeight(row) < minHeight)
          table.setRowHeight(row, minHeight);

        minHeight = -1;
      }
      
      return c;
    }
  }

 
  public class CellEditing extends DefaultCellEditor
  {
    /** */
    private static final long serialVersionUID = 1L;
    public CellEditing(JTextField textField)
    {
      super(textField);
    }    

    public boolean stopCellEditing()
    {
      isChanged = true;
      return super.stopCellEditing();
    }
  }
  
  /**
   *
   */
  public class ButtonEditor extends DefaultCellEditor 
  {
    /** */
    private static final long serialVersionUID = 1L;
    protected JButton buttRemove;
    private boolean   isPushed;
    private int selectedRow;
    private Color fgColor = new Color(139,35,35);
    private DefaultTableModel tableModel;
    
    public ButtonEditor(JCheckBox checkBox, final DefaultTableModel tableModel) 
    {
      super(checkBox);
      this.tableModel = tableModel;
      
      buttRemove = new JButton("X");

      buttRemove.setOpaque(false);
      Font font = buttRemove.getFont().deriveFont(Font.BOLD);
      buttRemove.setFont(font);
      buttRemove.setToolTipText("REMOVE");
      buttRemove.setForeground(fgColor);
      buttRemove.setPreferredSize(new Dimension(
          20,20));
      buttRemove.setMaximumSize(buttRemove.getPreferredSize());
      
      buttRemove.addActionListener(new ActionListener()
      {
        public void actionPerformed(ActionEvent e) 
        {
          fireEditingStopped();
          
        }
      });
    }
   
    public Component getTableCellEditorComponent(JTable table, Object value,
                     boolean isSelected, int row, int column)
    {
      if (isSelected) 
      {
        buttRemove.setForeground(fgColor);
        buttRemove.setBackground(table.getSelectionBackground()); 
      }
      else
      {
        buttRemove.setForeground(fgColor);
        buttRemove.setBackground(table.getBackground());
      }
      
      selectedRow = row;
      isPushed = true;
      return buttRemove;
    }
   
    public Object getCellEditorValue() 
    {
      if(isPushed)  
      {
        tableModel.removeRow(selectedRow);
        isChanged = true;
        return null;
      }
      isPushed = false;
      return new String("X") ;
    }
     
    public boolean stopCellEditing() 
    {
      isPushed = false;
      return super.stopCellEditing();
    }
   
    protected void fireEditingStopped() 
    {
      try
      {
        super.fireEditingStopped();
      }
      catch(ArrayIndexOutOfBoundsException e){}
    }
  }
}