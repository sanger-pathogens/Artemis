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
import java.awt.Cursor;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseMotionAdapter;
import java.util.Vector;

import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JLabel;
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
import uk.ac.sanger.artemis.util.DatabaseDocument;
import uk.ac.sanger.artemis.util.StringVector;

public class SimilarityTable extends AbstractMatchTable
{ 
  private int NUMBER_COLUMNS = 13;
  private Vector rowData   = new Vector();
  private Vector tableData = new Vector(NUMBER_COLUMNS);
  
  private JButton infoLevelButton = new JButton("Details");
  
  
  //
  // column headings
  final static String ORGANISM_COL = "Organism";
  final static String HIT_COL = "Hit";
  final static String HIT_DBXREF_COL = "DbXRef";
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
  protected SimilarityTable(final Qualifier simQualifier,
       final DatabaseDocument doc)
  {
    this.origQualifiers = new QualifierVector();
    this.origQualifiers.add(simQualifier);
    
    infoLevelButton.setOpaque(false);
    infoLevelButton.setHorizontalAlignment(SwingConstants.LEFT);
    tableData.setSize(NUMBER_COLUMNS);
    
    tableData.setElementAt(ORGANISM_COL,0);
    tableData.setElementAt(HIT_COL,1);
    tableData.setElementAt(HIT_DBXREF_COL,2);
    tableData.setElementAt(DESCRIPTION_COL,3);
    tableData.setElementAt(EVALUE_COL,4);
    tableData.setElementAt(LENGTH_COL,5);
    tableData.setElementAt(ID_COL,6);
    tableData.setElementAt(QUERY_COL,7);
    tableData.setElementAt(SUBJECT_COL,8);
    tableData.setElementAt(SCORE_COL,9);
    tableData.setElementAt(OVERLAP_COL,10);
    tableData.setElementAt(METHOD_COL,11);
    tableData.setElementAt(REMOVE_BUTTON_COL,12);
    
    // add rows of similarity
    StringVector sims = simQualifier.getValues();
    for(int i=0; i<sims.size(); i++)
      rowData.add(getRowData((String)sims.get(i), tableData));
    
    JTable similarityTable = new JTable(rowData, tableData);
    setTable(similarityTable);
    
    // set hand cursor
    similarityTable.addMouseMotionListener( new MouseMotionAdapter() 
    {
      private Cursor handCursor = Cursor.getPredefinedCursor(Cursor.HAND_CURSOR);
      public void mouseMoved(MouseEvent e) 
      {
        int col = table.columnAtPoint(e.getPoint());
        
        String colName = table.getColumnName(col);
     
        if(colName.equals(HIT_COL) || colName.equals(HIT_DBXREF_COL) ||
           colName.equals(REMOVE_BUTTON_COL)) 
          table.setCursor(handCursor);
        else 
          table.setCursor(Cursor.getDefaultCursor());  
      }
    });
    
    
    similarityTable.setColumnSelectionAllowed(false);
    similarityTable.setRowSelectionAllowed(true);
    
    packColumn(similarityTable, getColumnIndex(LENGTH_COL), 4);
    packColumn(similarityTable, getColumnIndex(EVALUE_COL), 4);
    packColumn(similarityTable, getColumnIndex(ID_COL), 4);
    packColumn(similarityTable, getColumnIndex(HIT_COL), 6);
    packColumn(similarityTable, getColumnIndex(HIT_DBXREF_COL), 6);

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
    }
    
    infoLevelButton.addActionListener(new ActionListener()
    {
      private boolean show = true;
      public void actionPerformed(ActionEvent e)
      {
        // change the column size 
        for(int i=0; i<hideColumns.length; i++)
        {
          if(show)
            packColumn(getTable(), getColumnIndex(
               (String) hideColumns[i].getHeaderValue()), 2);
          else
          {
            hideColumns[i].setMinWidth(0);
            hideColumns[i].setMaxWidth(0);
          }
        }
        show = !show;
        
        if(infoLevelButton.getText().equals("Details"))
          infoLevelButton.setText("Hide Details");
        else
          infoLevelButton.setText("Details");
      } 
    });
    
    TableModel tableModel = getTable().getModel();
    // remove button column
    TableColumn col = getTable().getColumn(REMOVE_BUTTON_COL);
    col.setMinWidth(35);
    col.setMaxWidth(40);
    col.setPreferredWidth(40);

    final SimilarityRenderer renderer = new SimilarityRenderer();

    for(int columnIndex = 0; columnIndex <tableModel.getColumnCount();
        columnIndex++) 
    {
      col = getTable().getColumnModel().getColumn(columnIndex);
      col.setCellRenderer(renderer);
      col.setCellEditor(new CellEditing(new JTextField()));
    }
    
    
    col = getTable().getColumn(HIT_COL);
    col.setCellEditor(new LinkEditor(new JCheckBox(),
        (DefaultTableModel)getTable().getModel(), null));
    
    col = getTable().getColumn(HIT_DBXREF_COL);
    col.setCellEditor(new LinkEditor(new JCheckBox(),
        (DefaultTableModel)getTable().getModel(), null));
    
    // remove JButton column
    col = getTable().getColumn(REMOVE_BUTTON_COL);
    col.setCellEditor(new ButtonEditor(new JCheckBox(),
        (DefaultTableModel)getTable().getModel(), "X", doc));
  }
  
  /**
   * Build a vector of the row data
   * @param similarityString
   * @return
   */
  private Vector getRowData(String similarityString,
                            final Vector tableData)
  {
    Vector row = new Vector(NUMBER_COLUMNS);
    row.setSize(NUMBER_COLUMNS);
    
    if(similarityString.startsWith("\""))
      similarityString = similarityString.substring(1);
    if(similarityString.endsWith("\""))
      similarityString = similarityString.substring(0,similarityString.length()-1);
      
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
      
      String hit = ((String)sim.get(1)).trim();
      
      if(hit.startsWith("with="))
        hit = hit.substring(5);
      
      final String hits[] = hit.split(" ");
      
      row.setElementAt(hits[0], columnIndex);
      
      if(hits.length > 1)
      {
        // dbxref 
        columnIndex = tableData.indexOf(HIT_DBXREF_COL);
        
        if(hits[1].startsWith("(") && hits[1].endsWith(")"))
          hits[1] = hits[1].substring(1,hits[1].length()-1);
        
        row.setElementAt(hits[1], columnIndex);
      }
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
    if( !(lenString=getField("length=", similarityString).trim()).equals("") )
    {
      int columnIndex = tableData.indexOf(LENGTH_COL);
      row.setElementAt(lenString, columnIndex);
    }
    else if( !(lenString=getField("length", similarityString).trim()).equals("") )
    {
      int columnIndex = tableData.indexOf(LENGTH_COL);
      row.setElementAt(lenString, columnIndex);
    }
    
    String ungappedId;
    if( !(ungappedId=getField("ungapped id", similarityString)).equals("") )
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
    else if(similarityString.indexOf("overlap;") > -1)
    {
      overlap = null;
      for(int i=0;i<sim.size(); i++)
      {
        String val = (String)sim.get(i); 
        if( val.endsWith("overlap") )
        {
          overlap = val;
          break;
        }
      }
      if(overlap != null)
      {
        int columnIndex = tableData.indexOf(OVERLAP_COL);
        row.setElementAt(overlap, columnIndex);
      }
    }
    
    int columnIndex = tableData.indexOf(METHOD_COL);
    row.setElementAt(((String)sim.get(0)).trim(), columnIndex);
    return row;
  }


  /**
   * Button to show/hide columns
   * @return
   */
  public JButton getInfoLevelButton()
  {
    return infoLevelButton;
  }

  
  /**
   * Called by AbstractMatchTable.updateQualifier()
   */
  protected String updateQualifierString(final int row)
  {
    StringBuffer similarityStr = new StringBuffer(
             (String)getTable().getValueAt(row, getColumnIndex(METHOD_COL)) );   // method
    similarityStr.append(";");
    similarityStr.append(
             (String)getTable().getValueAt(row, getColumnIndex(HIT_COL)) );      // hit
    
                                                                                 // dbxref
    String dbxref = (String)getTable().getValueAt(row, getColumnIndex(HIT_DBXREF_COL));
    if(!dbxref.equals(""))
      similarityStr.append(" ("+dbxref+")");
    
    similarityStr.append(";");
    similarityStr.append(
             (String)getTable().getValueAt(row, getColumnIndex(ORGANISM_COL)) ); // organism
    similarityStr.append(";");
    similarityStr.append(
             (String)getTable().getValueAt(row, getColumnIndex(DESCRIPTION_COL)) ); // description
    similarityStr.append(";");
    similarityStr.append("length "+
             (String)getTable().getValueAt(row, getColumnIndex(LENGTH_COL)) ); // length
    similarityStr.append(";");
    similarityStr.append("E()="+
             (String)getTable().getValueAt(row, getColumnIndex(EVALUE_COL)) ); // evalue
    similarityStr.append(";");
    
    if(getTable().getValueAt(row, getColumnIndex(SCORE_COL)) != null)
    {
      similarityStr.append("score="+
             (String)getTable().getValueAt(row, getColumnIndex(SCORE_COL)) );  // score
      similarityStr.append(";");
    }
    
    similarityStr.append("query "+
             (String)getTable().getValueAt(row, getColumnIndex(QUERY_COL)) );  // query
    similarityStr.append(";");
    similarityStr.append("subject "+
             (String)getTable().getValueAt(row, getColumnIndex(SUBJECT_COL)) ); // subject
    similarityStr.append(";");
    
    if(getTable().getValueAt(row, getColumnIndex(ID_COL)) != null)
    {
      similarityStr.append("ungapped id="+
             (String)getTable().getValueAt(row, getColumnIndex(ID_COL)) ); // ungapped id
      similarityStr.append(";");
    }
    
    if(getTable().getValueAt(row, getColumnIndex(OVERLAP_COL)) != null)
      similarityStr.append("overlap="+
             (String)getTable().getValueAt(row, getColumnIndex(OVERLAP_COL)) ); // overlap
    
    return similarityStr.toString();
  }
  
  /**
   * Check whether s qualifier string exists in a StringVector for that qualifier.
   * If the StringVector contains the hit, organism, description & e-value then
   * return true.
   * @param qualStr
   * @param qualStringVector
   * @return
   */
  public static boolean containsStringInStringVector(final String qualStr, 
                                                     final StringVector qualStringVector)
  {
    StringVector sim1 = StringVector.getStrings(qualStr, ";");
    for(int i=0; i<qualStringVector.size(); i++)
    {
      String thisStr = (String)qualStringVector.get(i);
      
      StringVector sim2 = StringVector.getStrings(thisStr, ";");
      
      // hit
      if( !((String)sim1.get(1)).equals((String)sim2.get(1)) )
        continue;      
      
      // organism
      if( !((String)sim1.get(2)).equals((String)sim2.get(2)) )
        continue;

      // description
      if( !((String)sim1.get(3)).equals((String)sim2.get(3)) )
        continue;
      
      // e-value
      final String evalueString1 = getField("E()=", qualStr);
      final String evalueString2 = getField("E()=", thisStr);
      if( !(evalueString1.equals(evalueString2)) )
        continue;
      
      return true; 
    }
    return false;
  }
  
  /**
   * Renderer for the Similarity cells
   */
  private class SimilarityRenderer extends DefaultTableCellRenderer
  {  
    /** */
    private static final long serialVersionUID = 1L;
    private int minHeight = -1;
    
    private final JLabel hit = new JLabel();
    private final JLabel hit_dbxref = new JLabel();
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
    private final JLabel buttRemove = new JLabel("X");
    private Color fgColor = new Color(139,35,35);
    private Color fgLinkColor = Color.BLUE;
    
    public SimilarityRenderer() 
    {
      evalue.setHorizontalAlignment(SwingConstants.RIGHT);
      evalue.setVerticalAlignment(SwingConstants.TOP);
      evalue.setOpaque(true);
      length.setHorizontalAlignment(SwingConstants.RIGHT);
      length.setVerticalAlignment(SwingConstants.TOP);
      length.setOpaque(true);
      ungappedId.setHorizontalAlignment(SwingConstants.RIGHT);
      ungappedId.setVerticalAlignment(SwingConstants.TOP);
      ungappedId.setOpaque(true);
      queryCoord.setHorizontalAlignment(SwingConstants.RIGHT);
      queryCoord.setVerticalAlignment(SwingConstants.TOP);
      queryCoord.setOpaque(true);
      subjCoord.setHorizontalAlignment(SwingConstants.RIGHT);
      subjCoord.setVerticalAlignment(SwingConstants.TOP);
      subjCoord.setOpaque(true);
      score.setHorizontalAlignment(SwingConstants.RIGHT);
      score.setVerticalAlignment(SwingConstants.TOP);
      score.setOpaque(true);
      overlap.setHorizontalAlignment(SwingConstants.RIGHT);
      overlap.setVerticalAlignment(SwingConstants.TOP);
      overlap.setOpaque(true);
      method.setHorizontalAlignment(SwingConstants.RIGHT);
      method.setVerticalAlignment(SwingConstants.TOP);
      method.setOpaque(true);
      
      organismTextArea.setLineWrap(true);
      organismTextArea.setWrapStyleWord(true);
      
      hit.setHorizontalAlignment(SwingConstants.CENTER);
      hit.setVerticalAlignment(SwingConstants.TOP);
      hit.setOpaque(true);
      hit_dbxref.setHorizontalAlignment(SwingConstants.CENTER);
      hit_dbxref.setVerticalAlignment(SwingConstants.TOP);
      hit_dbxref.setOpaque(true);
      
      descriptionTextArea.setLineWrap(true);
      descriptionTextArea.setWrapStyleWord(true);
      
      buttRemove.setOpaque(true);
      buttRemove.setText("X");
      
      Font font = getFont().deriveFont(Font.BOLD);
      buttRemove.setFont(font);
      buttRemove.setToolTipText("REMOVE");
      buttRemove.setHorizontalAlignment(SwingConstants.CENTER);
      buttRemove.setVerticalAlignment(SwingConstants.TOP);
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
        hit.setText(text);
        hit.setForeground(fgLinkColor);

        c = hit;
      }
      else if(column == getColumnIndex(HIT_DBXREF_COL))
      {
        hit_dbxref.setText(text);
        hit_dbxref.setForeground(fgLinkColor);
        
        c = hit_dbxref;
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
      if(column < 4)
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