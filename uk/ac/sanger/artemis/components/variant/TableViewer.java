/* TableViewer
 *
 * created: 2011
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2011 Genome Research Limited
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
package uk.ac.sanger.artemis.components.variant;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Comparator;
import java.util.Vector;

import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.table.TableModel;
import javax.swing.table.TableRowSorter;

import uk.ac.sanger.artemis.components.StickyFileChooser;


class TableViewer extends JFrame
{
  private static final long serialVersionUID = 1L;
  private Comparator<Integer> comparator = new Comparator<Integer>() 
  {
    public int compare(Integer i1, Integer i2) 
    {
        return i1.compareTo(i2);
    }
  };

  private TableRowSorter<TableModel> sorter;
  
  public TableViewer(final Vector<Vector<Object>> rowData, final Vector<String> columnData, String title)
  {
    super(title);
    final JTable variantData = new JTable(rowData, columnData);
    sorter = new TableRowSorter<TableModel>(variantData.getModel());
    variantData.setRowSorter(sorter);
    
    final JFrame f = new JFrame("Variant Overview");
    final JScrollPane jsp = new JScrollPane(variantData);
    final Dimension screen = Toolkit.getDefaultToolkit().getScreenSize();
    jsp.setPreferredSize(new Dimension(screen.width/3, screen.height/3));
    final JPanel overviewPane = new JPanel(new BorderLayout());
    overviewPane.add(jsp, BorderLayout.CENTER);
    
    final Box xBox = Box.createHorizontalBox();
    xBox.add(Box.createHorizontalGlue());
    
    final JButton save = new JButton("SAVE");
    save.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        save(variantData, columnData);
      }
    });
    xBox.add(save);
    
    final JButton close = new JButton("CLOSE");
    close.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        f.dispose();
      }
    });
    xBox.add(close);
    overviewPane.add(xBox, BorderLayout.SOUTH);
    f.getContentPane().add(overviewPane);
    f.pack();
    f.setVisible(true);
  }
  
  public void setIntegerRowSorter(int colIndex)
  {
    sorter.setComparator(colIndex, comparator);
  }
  
  private void save(final JTable table, final Vector<String> columnData)
  {
    StickyFileChooser fc = new StickyFileChooser();
    int status = fc.showSaveDialog(TableViewer.this);
    if(status != StickyFileChooser.APPROVE_OPTION)
      return;
    try
    {
      File f = fc.getSelectedFile();
      if(f.exists())
      {
        status = JOptionPane.showConfirmDialog(TableViewer.this, 
            f.getName()+" exists overwrite?", "Overwrite", JOptionPane.YES_NO_OPTION);
        if(status != JOptionPane.YES_OPTION)
          return;
      }
      
      FileWriter writer = new FileWriter(fc.getSelectedFile());
      for(String col: columnData)
        writer.write(col+"\t");
      writer.write("\n");
      
      int nrows = table.getRowCount();
      int ncols = table.getColumnCount();
      for(int i=0; i<nrows; i++)
      {
        for(int j=0; j<ncols; j++)
        {
          // use the sorted row index
          int rowIndex = table.convertRowIndexToModel(i); 
          writer.write(table.getModel().getValueAt(rowIndex, j)+"\t");
        }
        writer.write("\n");
      }

      writer.close();
    }
    catch (IOException e1)
    {
      JOptionPane.showMessageDialog(TableViewer.this, e1.getMessage(), 
          "Problem Writing", JOptionPane.ERROR_MESSAGE);
    } 
  }
}