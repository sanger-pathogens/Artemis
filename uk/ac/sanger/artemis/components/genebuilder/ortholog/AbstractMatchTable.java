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
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.io.IOException;
import java.util.Collection;
import java.util.Iterator;

import javax.swing.DefaultCellEditor;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.SwingConstants;
import javax.swing.SwingUtilities;
import javax.swing.TransferHandler;
import javax.swing.table.DefaultTableColumnModel;
import javax.swing.table.DefaultTableModel;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumn;
import javax.swing.text.JTextComponent;

import org.gmod.schema.sequence.FeatureLoc;

import uk.ac.sanger.artemis.Entry;
import uk.ac.sanger.artemis.components.ArtemisMain;
import uk.ac.sanger.artemis.components.EntryEdit;
import uk.ac.sanger.artemis.components.SwingWorker;
import uk.ac.sanger.artemis.components.Utilities;
import uk.ac.sanger.artemis.components.genebuilder.GeneEdit;
import uk.ac.sanger.artemis.editor.BrowserControl;
import uk.ac.sanger.artemis.editor.DataCollectionPane;
import uk.ac.sanger.artemis.io.DatabaseDocumentEntry;
import uk.ac.sanger.artemis.io.EntryInformationException;
import uk.ac.sanger.artemis.io.Qualifier;
import uk.ac.sanger.artemis.io.QualifierVector;
import uk.ac.sanger.artemis.io.Range;
import uk.ac.sanger.artemis.sequence.Bases;
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.InputStreamProgressEvent;
import uk.ac.sanger.artemis.util.InputStreamProgressListener;
import uk.ac.sanger.artemis.util.OutOfRangeException;
import uk.ac.sanger.artemis.util.StringVector;

abstract class AbstractMatchTable
{
  protected boolean isChanged = false;
  protected JTable table;
  protected QualifierVector origQualifiers;
  private JLabel status_line;
  private static org.apache.log4j.Logger logger4j = 
    org.apache.log4j.Logger.getLogger(AbstractMatchTable.class);
  protected abstract String updateQualifierString(int row);
  
  /**
   * Determine if qualifiers have changed
   */
  protected boolean isQualifierChanged()
  {
    return isChanged; 
  }
  
  /**
   * Override this if there are multiple qualifier types in the table
   * e.g. for ortholog, paralog
   * @param qualifierName
   * @param row
   * @return
   */
  protected boolean isRowOfType(final String qualifierName, final int row)
  {
    return true;
  }
  
  /**
   * Override this method if there are multiple qualifier types in the table
   * e.g. for ortholog, paralog. See OrthoParalogTable.getOtherValues().
   * @param origQualifier
   * @return
   */
  protected StringVector getOtherValues(final Qualifier origQualifier)
  {
    StringVector values = origQualifier.getValues();
    values.removeAllElements();
    return values;
  }
  
  /**
   * Update the qualifiers from the entries in the table
   * @param qv
   */
  protected void updateQualifier(final QualifierVector qv)
  {
    for(int i=0; i<origQualifiers.size(); i++)
    {
      Qualifier origQualifier = (Qualifier) origQualifiers.elementAt(i);
      StringVector values = getOtherValues(origQualifier);

      logger4j.debug("AbstractMatchTable.updateQualifier() new value:\n");
      for(int j = 0; j < getTable().getRowCount(); j++)
      {
        if(isRowOfType(origQualifier.getName(), j))
        {
          String updatedQualifierString = updateQualifierString(j);
          values.add(updatedQualifierString);
          logger4j.debug(updatedQualifierString);
        }
      }
      logger4j.debug("\n");

      if(values.size() < 1)
        qv.remove(origQualifier);
      else
      {
        int index = qv.indexOfQualifierWithName(origQualifier.getName());
        origQualifier = new Qualifier(origQualifier.getName(), values);
        qv.remove(index);
        qv.add(index, origQualifier);
      }
    }
  }
  
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
    try
    {
      int modelColumnIndex = getTable().getColumn(columnName).getModelIndex();
      return getTable().convertColumnIndexToView(modelColumnIndex);
    }
    catch(IllegalArgumentException e)
    {
      return -1;
    }
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
    
    if(field.startsWith("="))
      field = field.substring(1);
    return field;
  }
  
  protected void openArtemis(final DatabaseDocument doc, final int selectedRow)
  {
    int columnIndex = getColumnIndex(OrthoParalogTable.GENE_COL);
    
    String geneRef = (String)getTable().getValueAt(selectedRow, columnIndex);
    final String gene[] = geneRef.split(":");
    
    final org.gmod.schema.sequence.Feature geneFeature = doc.getFeatureByUniquename(gene[1]);

    Collection featureLocs = geneFeature.getFeatureLocsForFeatureId();
    Iterator it = featureLocs.iterator();
    final FeatureLoc featureLoc = (FeatureLoc)it.next();

    final JFrame progressFrame = progressReading();

    SwingWorker readWorker = new SwingWorker()
    {
      public Object construct()
      {
        try
        {
          int start = featureLoc.getFmin().intValue()-10000;
          if(start <= 0)
            start = 1;
          Range range = new Range(start,featureLoc.getFmax().intValue()+10000);
          final DatabaseDocument newDoc = new DatabaseDocument(
              doc, gene[0], geneFeature, range,
              stream_progress_listener);
          newDoc.setLazyFeatureLoad(false);
          
          DatabaseDocumentEntry db_entry = new DatabaseDocumentEntry(newDoc, null);
          Bases bases = new Bases(db_entry.getSequence());
          Entry entry = new Entry(bases, db_entry);

          final EntryEdit new_entry_edit = ArtemisMain.makeEntryEdit(entry);
          new_entry_edit.getGotoEventSource().gotoBase(featureLoc.getFmin().intValue());
          new_entry_edit.setVisible(true);
        }
        catch(EntryInformationException e)
        {
          e.printStackTrace();
        }
        catch(IOException e)
        {
          e.printStackTrace();
        }
        catch(OutOfRangeException e)
        {
          e.printStackTrace();
        }
        return null;
      }
      
      public void finished()
      {
        progressFrame.dispose();
      }
    };
    readWorker.start();
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

      if(comp instanceof JLabel)
        text = ((JLabel)comp).getText();
      else
        text = ((JTextComponent)comp).getText();
      
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
   private String text;
   private DatabaseDocument doc;
   
   public ButtonEditor(JCheckBox checkBox, final DefaultTableModel tableModel,
       final String text, final DatabaseDocument doc) 
   {
     super(checkBox);
     this.tableModel = tableModel;
     this.text = text;
     this.doc = doc;
     
     buttRemove = new JButton(text);
     buttRemove.setBorderPainted(false);
     buttRemove.setOpaque(false);
     Font font = buttRemove.getFont().deriveFont(Font.BOLD);
     buttRemove.setFont(font);
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
       if(text.equals("X"))
       {
         tableModel.removeRow(selectedRow);
         isChanged = true;
       }
       else
         openArtemis(doc, selectedRow);
       //return null;
     }
     isPushed = false;
     return text;
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
   * Let the user know the progress of reading from the database
   */
  private JFrame progressReading()
  {
    final JFrame fread = new JFrame();
    fread.setUndecorated(true);
    status_line = new JLabel("Loading ...                         ");
    status_line.setBackground(Color.white);
    status_line.setOpaque(true);
    fread.getContentPane().add(status_line);
    
    fread.pack();
    Utilities.centreFrame(fread);
    fread.setVisible(true);
    return fread;
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
     
     linkButton.addMouseListener(new MouseListener()
     {

      public void mouseClicked(MouseEvent e)
      {
        if(e.getClickCount() == 2 &&
           !e.isShiftDown() && !e.isPopupTrigger())
          fireEditingStopped(); 
      }

      public void mouseEntered(MouseEvent e)
      { 
      }

      public void mouseExited(MouseEvent e)
      {
      }

      public void mousePressed(MouseEvent e)
      {
      }

      public void mouseReleased(MouseEvent e)
      {
      }     
     });
   }
   
   public Component getTableCellEditorComponent(JTable table, Object value,
                    boolean isSelected, int row, int column)
   {
     linkButton.setActionCommand((String)value);
     String gene[] = ((String)value).split(":");
     if(doc == null)
       linkButton.setText((String)value);
     else
       linkButton.setText(gene[1]);
     
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
         
         BrowserControl.displayURL(srscmd);
       }
       else
       {  
          // open gene editor for this gene link
          linkButton.setCursor(Cursor.getPredefinedCursor(Cursor.WAIT_CURSOR));
          final String reference[] = linkButton.getActionCommand().split(":");
          DatabaseDocumentEntry entry = GeneEdit.makeGeneEntry(
          		reference[0], reference[1], doc, stream_progress_listener);

          if(entry != null)
          {
            //entry.setReadOnly(true);
            GeneEdit.showGeneEditor(reference[0], reference[1], entry);
          }
       }
       //isChanged = true;
       linkButton.setCursor(Cursor.getDefaultCursor());
       return linkButton.getActionCommand();
     }
     isPushed = false;
     
     if(doc == null)
       return link;
     else
       return linkButton.getActionCommand();
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
      return MOVE;
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

  /**
   *  An InputStreamProgressListener used to update the error label with the
   *  current number of chars read.
   **/
  protected final InputStreamProgressListener stream_progress_listener =
    new InputStreamProgressListener() 
  {
    public void progressMade(final InputStreamProgressEvent event) 
    {
      final int char_count = event.getCharCount();
      if(char_count == -1) 
        status_line.setText("");
      else 
        status_line.setText("chars read so far: " + char_count);
    }

    public void progressMade(String progress)
    {
      status_line.setText(progress);
    }
  };

}