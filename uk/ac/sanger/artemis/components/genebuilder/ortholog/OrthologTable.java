/* OrthologTable.java
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
import java.awt.FlowLayout;
import java.awt.Font;
import java.util.Vector;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import javax.swing.UIManager;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableColumn;
import javax.swing.table.TableModel;

import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;

public class OrthologTable extends AbstractMatchTable
{
  private int NUMBER_COLUMNS = 4;
  private Vector rowData   = new Vector();
  private Vector tableData = new Vector(NUMBER_COLUMNS);
  private JTable orthologTable;
  private JButton infoLevelButton = new JButton("Details");
  private Qualifier origQualifier;
  private boolean isChanged = false;
  
  //
  // column headings
  final static String ORGANISM_COL = "Organism";
  final static String ORTHO_COL = "Ortholog";
  final static String DESCRIPTION_COL = "Description";
  final static String REMOVE_BUTTON_COL = "";
  
  /**
   * Contruct a component for a similarity line
   * @param similarity
   * @param similarityString
   */
  protected OrthologTable(final Qualifier origQualifier)
  {
    this.origQualifier = origQualifier;
    
    infoLevelButton.setOpaque(false);
    tableData.setSize(NUMBER_COLUMNS);
    
    tableData.setElementAt(ORGANISM_COL,0);
    tableData.setElementAt(ORTHO_COL,1);
    tableData.setElementAt(DESCRIPTION_COL,2);
    tableData.setElementAt(REMOVE_BUTTON_COL,3);
    
    // add row data
    Vector thisRowData = new Vector();
    thisRowData.add("organism");
    thisRowData.add("blah");
    thisRowData.add("blah blah");
    rowData.add(thisRowData);
    Vector thisRowData2 = new Vector();
    thisRowData2.add("organism2");
    thisRowData2.add("blah2");
    thisRowData2.add("blah blah2");
    rowData.add(thisRowData2);
    
    orthologTable = new JTable(rowData, tableData);
    setTable(orthologTable);
    orthologTable.setColumnSelectionAllowed(false);
    orthologTable.setRowSelectionAllowed(true);
    
    TableModel tableModel = orthologTable.getModel();
    // remove button column
    TableColumn col = orthologTable.getColumn(REMOVE_BUTTON_COL);
    col.setMinWidth(35);
    col.setMaxWidth(40);
    col.setPreferredWidth(40);

    final OrthologRenderer renderer = new OrthologRenderer();

    for(int columnIndex = 0; columnIndex <tableModel.getColumnCount();
        columnIndex++) 
    {
      col = orthologTable.getColumnModel().getColumn(columnIndex);
      col.setCellRenderer(renderer);
      col.setCellEditor(new CellEditing(new JTextField()));
    }
    
    // remove JButton column
    col = orthologTable.getColumn(REMOVE_BUTTON_COL);
    col.setCellEditor(new ButtonEditor(new JCheckBox(),
        (DefaultTableModel)orthologTable.getModel()));
  }
  
  protected boolean isQualifierChanged()
  {
    // TODO Auto-generated method stub
    return false;
  }

  protected void updateQualifier(QualifierVector qv)
  {
    // TODO Auto-generated method stub
  } 
  
  

  /**
   * Renderer for the Ortholog cells
   */
  private class OrthologRenderer extends DefaultTableCellRenderer
  {  
    /** */
    private static final long serialVersionUID = 1L;
    private int minHeight = -1;
    
    private final JTextArea orthologTextArea = new JTextArea();
    private final JTextArea organismTextArea = new JTextArea();
    private final JTextArea descriptionTextArea = new JTextArea();
    private final JLabel buttRemove = new JLabel("X");
    private Color fgColor = new Color(139,35,35);
    
    public OrthologRenderer() 
    {
      
      organismTextArea.setLineWrap(true);
      organismTextArea.setWrapStyleWord(true);
      
      orthologTextArea.setLineWrap(true);
      orthologTextArea.setWrapStyleWord(true);
      
      descriptionTextArea.setLineWrap(true);
      descriptionTextArea.setWrapStyleWord(true);

      buttRemove.setOpaque(true);
      buttRemove.setText("X");
      
      Font font = getFont().deriveFont(Font.BOLD);
      buttRemove.setFont(font);
      buttRemove.setToolTipText("REMOVE");
      buttRemove.setHorizontalAlignment(SwingConstants.CENTER);
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
      else if(column == getColumnIndex(ORTHO_COL))
      {
        orthologTextArea.setText(text);
        
        tableCol = table.getColumnModel().getColumn(column);
        orthologTextArea.setSize(tableCol.getWidth(), table.getRowHeight(row));

        dim = orthologTextArea.getPreferredSize();
        minHeight = Math.max(minHeight, dim.height);
        
        c = orthologTextArea;
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
      
      // highlight on selection
      if(isSelected) 
        c.setBackground(table.getSelectionBackground()); 
      else
        c.setBackground(Color.white);
      
      return c;
    }
  }

 
}