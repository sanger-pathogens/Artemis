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
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Insets;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.StringSelection;
import java.awt.datatransfer.Transferable;
import java.awt.datatransfer.UnsupportedFlavorException;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;

import javax.swing.DefaultCellEditor;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.TransferHandler;
import javax.swing.table.DefaultTableColumnModel;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumn;

import uk.ac.sanger.artemis.components.genebuilder.GeneEdit;
import uk.ac.sanger.artemis.editor.BrowserControl;
import uk.ac.sanger.artemis.editor.DataCollectionPane;
import uk.ac.sanger.artemis.io.DatabaseDocumentEntry;
import uk.ac.sanger.artemis.io.EntryInformationException;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.util.DatabaseDocument;

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
  protected static String getField(final String fieldName, final String qualifierString)
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

      if(text == null)
        continue;
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
  
  
  /**
  *
  */
  protected class LinkEditor extends DefaultCellEditor 
  {
   /** */
   private static final long serialVersionUID = 1L;
   protected JButton linkButton = new JButton();
   private boolean   isPushed;
   private Color fgLinkColor = Color.BLUE;
   private DatabaseDocument doc;
   
   public LinkEditor(JCheckBox checkBox, final DefaultTableModel tableModel,
                     DatabaseDocument doc) 
   {
     super(checkBox);
     this.doc = doc;
     
     linkButton.setBorderPainted(false);
     linkButton.setOpaque(true);

     linkButton.setHorizontalAlignment(SwingConstants.LEFT);
     linkButton.setVerticalAlignment(SwingConstants.TOP);
     linkButton.setMargin(new Insets(0,1,0,1));
     
     linkButton.addActionListener(new ActionListener()
     {
       public void actionPerformed(ActionEvent e) 
       {
         fireEditingStopped(); 
       }
     });
   }
   
   /**
    * 
    * @param schema
    * @param uniquename
    * @return
    */
   private DatabaseDocumentEntry makeEntry(final String schema, 
       final String uniquename)
   {
     DatabaseDocumentEntry db_entry = null;
     DatabaseDocument newdoc = new DatabaseDocument(doc, 
             uniquename, schema, true);
     
     try
     {
       db_entry = new DatabaseDocumentEntry(newdoc, null);
     }
     catch(EntryInformationException e)
     {
       e.printStackTrace();
     }
     catch(IOException e)
     {
       e.printStackTrace();
     }
     catch(NullPointerException npe)
     {
       JOptionPane.showMessageDialog(null, schema+":"+uniquename+
           " not found!", "Warning", JOptionPane.WARNING_MESSAGE);
     }

     return db_entry;
   }
   
   public Component getTableCellEditorComponent(JTable table, Object value,
                    boolean isSelected, int row, int column)
   {
     linkButton.setText((String)value);
     if (isSelected) 
     {
       linkButton.setForeground(fgLinkColor);
       linkButton.setBackground(table.getSelectionBackground()); 
     }
     else
     {
       linkButton.setForeground(fgLinkColor);
       linkButton.setBackground(table.getBackground());
     }
     
     isPushed = true;
     return linkButton;
   }
  
   public Object getCellEditorValue() 
   {
     String link = linkButton.getText();
     if(isPushed)  
     {
       if(doc == null)
       {
         // open in default browser
         String srscmd = DataCollectionPane.getSrsSite()+"/wgetz?-e+["+link+"]";

         // link to uniprot accession
         int ind;
         if( (ind = srscmd.indexOf("UniProt:")) > -1)
           srscmd = srscmd.substring(0,ind+7)+"-acc:"+
                    srscmd.substring(ind+8);

         if(srscmd.indexOf("ebi.ac.uk") > -1)
           srscmd = srscmd + "+-vn+2";
         
         
         System.out.println(srscmd);
         BrowserControl.displayURL(srscmd);
       }
       else
       {  
          // open gene editor for this gene link
          linkButton.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
          String reference[] = link.split(":");
          DatabaseDocumentEntry entry = makeEntry(reference[0], reference[1]);

          if(entry != null)
          {
            entry.setReadOnly(true);
            GeneEdit.showGeneEditor(reference[0], reference[1], entry);
          }
        }
       isChanged = true;
       linkButton.setCursor(Cursor.getDefaultCursor());
       return link;
     }
     isPushed = false;
     return link;
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

  //
  //http://java.sun.com/docs/books/tutorial/uiswing/dnd/intro.html
  //
  protected abstract class StringTransferHandler extends TransferHandler
  {
    protected abstract String exportString(JComponent c);
    protected abstract void importString(JComponent c, String str);
    protected abstract void cleanup(JComponent c, boolean remove);

    protected Transferable createTransferable(JComponent c)
    {
      return new StringSelection(exportString(c));
    }

    public int getSourceActions(JComponent c)
    {
      return COPY_OR_MOVE;
    }

    public boolean importData(JComponent c, Transferable t)
    {
      if(canImport(c, t.getTransferDataFlavors()))
      {
        try
        {
          String str = (String) t.getTransferData(DataFlavor.stringFlavor);
          importString(c, str);
          return true;
        }
        catch(UnsupportedFlavorException ufe){}
        catch(IOException ioe){}
      }
      return false;
    }

    protected void exportDone(JComponent c, Transferable data, int action)
    {
      cleanup(c, action == MOVE);
    }

    public boolean canImport(JComponent c, DataFlavor[] flavors)
    {
      for(int i = 0; i < flavors.length; i++)
      {
        if(DataFlavor.stringFlavor.equals(flavors[i]))
          return true;
      }
      return false;
    }
  }

  protected class TableTransferHandler extends StringTransferHandler
  {
    /** */
    private static final long serialVersionUID = 1L;
    private int[] rows = null;
    private int addIndex = -1; //Location where items were added
    private int addCount = 0; //Number of items added.

    protected String exportString(JComponent c)
    {
      JTable table = (JTable) c;
      rows = table.getSelectedRows();
      int colCount = table.getColumnCount();

      StringBuffer buff = new StringBuffer();

      for(int i = 0; i < rows.length; i++)
      {
        for(int j = 0; j < colCount; j++)
        {
          Object val = table.getValueAt(rows[i], j);
          buff.append(val == null ? "" : val.toString());
          if(j != colCount - 1)
            buff.append(",");
        }
        if(i != rows.length - 1)
          buff.append("\n");
      }
      return buff.toString();
    }

    protected void importString(JComponent c, String str)
    {
      JTable target = (JTable) c;
      DefaultTableModel model = (DefaultTableModel) target.getModel();
      int index = target.getSelectedRow();

      //Prevent the user from dropping data back on itself.
      //For example, if the user is moving rows #4,#5,#6 and #7 and
      //attempts to insert the rows after row #5, this would
      //be problematic when removing the original rows.
      //So this is not allowed.
      if(rows != null && index >= rows[0] - 1 && index <= rows[rows.length - 1])
      {
        rows = null;
        return;
      }

      int max = model.getRowCount();
      if(index < 0)
        index = max;
      else
      {
        index++;
        if(index > max)
          index = max;
      }
      addIndex = index;
      String[] values = str.split("\n");
      addCount = values.length;
      int colCount = target.getColumnCount();
      for(int i = 0; i < values.length && i < colCount; i++)
        model.insertRow(index++, values[i].split(","));
    }

    protected void cleanup(JComponent c, boolean remove)
    {
      JTable source = (JTable) c;
      if(remove && rows != null)
      {
        DefaultTableModel model = (DefaultTableModel) source.getModel();

        //If we are moving items around in the same table, we
        //need to adjust the rows accordingly, since those
        //after the insertion point have moved.
        if(addCount > 0)
        {
          for(int i = 0; i < rows.length; i++)
          {
            if(rows[i] > addIndex)
              rows[i] += addCount;
          }
        }
        for(int i = rows.length - 1; i >= 0; i--)
          model.removeRow(rows[i]);
      }
      rows = null;
      addCount = 0;
      addIndex = -1;
    }
  }


}