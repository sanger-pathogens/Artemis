/* AbstractMatchTable.java
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
import java.awt.FontMetrics;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.DefaultCellEditor;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;
import javax.swing.table.DefaultTableColumnModel;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumn;

import uk.ac.sanger.artemis.io.QualifierVector;

abstract class AbstractMatchTable
{
  protected boolean isChanged = false;
  protected JTable table;
  
  
  protected abstract boolean isQualifierChanged();
  protected abstract void updateQualifier(final QualifierVector qv);
  
  public void setTable(JTable table)
  {
    this.table = table;
  }
  
  protected JTable getTable()
  {
    return table;
  }
  
  /**
   * Get the column index from the column name
   * @param columnName
   * @return
   */
  protected int getColumnIndex(final String columnName)
  {
    int modelColumnIndex = getTable().getColumn(columnName).getModelIndex();
    return getTable().convertColumnIndexToView(modelColumnIndex);
  }
  
  
  /**
   * Strip out the value of a field of interest from a qualifier string
   * 
   * @param fieldName
   * @param qualifierString
   * @return
   */
  protected String getField(final String fieldName, final String qualifierString)
  {
    String field = "";
    
    int ind1 = qualifierString.toLowerCase().indexOf(fieldName.toLowerCase());
    int ind2 = qualifierString.indexOf(";", ind1);
    
    int len = fieldName.length();

    if(ind2 > ind1 && ind1 > -1)
      field = qualifierString.substring(ind1+len,ind2);
    else if(ind1 > -1)
      field = qualifierString.substring(ind1+len);
    
    return field;
  }
  
  /**
   * Sets the preferred, min & max width of the column specified by columnIndex. 
   * The column will be just wide enough to show the column head and the widest 
   * cell in the column. margin pixels are added to the left and right
   * @param table
   * @param columnIndex
   * @param margin
   */
  protected void packColumn(JTable table, int columnIndex, int margin) 
  {
    DefaultTableColumnModel colModel = (DefaultTableColumnModel)table.getColumnModel();
    TableColumn col = colModel.getColumn(columnIndex);
    int width = 0;
    int maxWidth;
    
    // Get width of column header
    TableCellRenderer renderer = col.getHeaderRenderer();
    if(renderer == null)
      renderer = table.getTableHeader().getDefaultRenderer();
      
    Component comp = renderer.getTableCellRendererComponent(
          table, col.getHeaderValue(), false, false, 0, 0);
    //width = comp.getPreferredSize().width;
  
    String text = ((JLabel)comp).getText();
    Font font = comp.getFont();
    FontMetrics fontMetrics = comp.getFontMetrics ( font );

    width = SwingUtilities.computeStringWidth ( fontMetrics, text );
    
    // Get maximum width of column data
    for(int r=0; r<table.getRowCount(); r++) 
    {
      renderer = table.getCellRenderer(r, columnIndex);
      comp = renderer.getTableCellRendererComponent(
              table, table.getValueAt(r, columnIndex), false, false, r, columnIndex);
      
      text = ((JLabel)comp).getText();
      font = comp.getFont();
      fontMetrics = comp.getFontMetrics ( font );

      maxWidth = SwingUtilities.computeStringWidth ( fontMetrics, text );
      //  maxWidth = comp.getPreferredSize().width;
      width = Math.max(width, maxWidth);
    }
  
    // Add margin
    width += 2*margin;
  
    // Set the width
    col.setPreferredWidth(width);
    col.setMaxWidth(width);
    col.setMinWidth(width);
  }
  
  protected class CellEditing extends DefaultCellEditor
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
  protected class ButtonEditor extends DefaultCellEditor 
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
     buttRemove.setBorderPainted(false);
     buttRemove.setOpaque(false);
     Font font = buttRemove.getFont().deriveFont(Font.BOLD);
     buttRemove.setFont(font);
     buttRemove.setToolTipText("REMOVE");
     buttRemove.setForeground(fgColor);
     
     Dimension size = new Dimension(20,20);
     buttRemove.setPreferredSize(size);
     buttRemove.setMaximumSize(size);
     
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