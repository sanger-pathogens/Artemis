/*
 * Copyright (C) 2008  Genome Research Limited
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
 *  @author: Tim Carver
 */

package uk.ac.sanger.artemis.circular;


import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.DefaultCellEditor;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JColorChooser;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSeparator;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.KeyStroke;
import javax.swing.event.*;
import javax.swing.table.AbstractTableModel;
import javax.swing.table.TableCellRenderer;
import javax.swing.border.Border;
import java.awt.event.*;
import java.awt.*;
import java.util.Vector;
import java.awt.Dimension;


public class RestrictionEnzyme extends JPanel
                           implements TableModelListener
{

  private DNADraw dna;

  public RestrictionEnzyme(final DNADraw dna, final Vector restrictionEnzyme)
  {
    super();
    this.dna = dna;
 
    Box bdown = Box.createVerticalBox();
    Dimension d = new Dimension(100,25);
    bdown.add(Box.createVerticalStrut(4));

    Vector columnNames = new Vector();
    columnNames.add("Label");
    columnNames.add("Position");
    columnNames.add("Colour");

    final RETableModel reModel = 
               new RETableModel(restrictionEnzyme,columnNames);
    final JTable reTable = new JTable(reModel);
    reTable.getModel().addTableModelListener(this);
    setUpColorRenderer(reTable);
    setUpColorEditor(reTable);

    JScrollPane scrollRE = new JScrollPane(reTable);
    scrollRE.setPreferredSize(new Dimension(300,150));
    scrollRE.getViewport().setBackground(Color.white);
    bdown.add(scrollRE);

    bdown.add(Box.createVerticalStrut(4));
    bdown.add(new JSeparator());

// Data entry
    bdown.add(Box.createVerticalStrut(4));
    Box bacross = Box.createHorizontalBox();
    
    final JTextField reLabel = new JTextField("");
    reLabel.setMaximumSize(d);
    reLabel.setPreferredSize(d);
    bacross.add(reLabel);
    bacross.add(new JLabel("Label"));
    bacross.add(Box.createHorizontalStrut(4));

    final TextFieldInt pos = new TextFieldInt();
    pos.setMaximumSize(d);
    pos.setPreferredSize(d);
    bacross.add(pos);
    bacross.add(new JLabel("Position"));
    bacross.add(Box.createHorizontalStrut(4));
 
    final JButton reColour = setUpColorButton(Color.blue);
    bacross.add(reColour);
    bacross.add(Box.createHorizontalStrut(4));
    bacross.add(Box.createHorizontalGlue());
    bdown.add(bacross);

    bacross = Box.createHorizontalBox();
    JButton addRE = new JButton("Add RE");
    addRE.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        Vector re = new Vector();
        re.add(reLabel.getText());
        re.add(new Integer(pos.getValue()));
        re.add(reColour.getBackground());
        restrictionEnzyme.add(re);
        reModel.addRow();
      }
    });
    bacross.add(addRE);
    bacross.add(Box.createHorizontalGlue());
    bdown.add(bacross);

    add(bdown);
  }

  public void tableChanged(TableModelEvent e) 
  {
    if(dna != null)
      dna.repaint();
  }
  
  protected JMenuBar createMenuBar(final JFrame f)
  {
    JMenuBar menuBar = new JMenuBar();

    JMenu fileMenu = new JMenu("File");
    fileMenu.setMnemonic(KeyEvent.VK_F);
    menuBar.add(fileMenu);

    JMenuItem closeMenu = new JMenuItem("Close");
    closeMenu.setAccelerator(KeyStroke.getKeyStroke(
              KeyEvent.VK_E, ActionEvent.CTRL_MASK));

    closeMenu.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        f.dispose();
      }
    });
    fileMenu.add(closeMenu);

    return menuBar;
  }


  class ColorRenderer extends JLabel
                       implements TableCellRenderer 
  {
    Border unselectedBorder = null;
    Border selectedBorder = null;
    boolean isBordered = true;

    public ColorRenderer(boolean isBordered) 
    {
      super();
      this.isBordered = isBordered;
      setOpaque(true); //MUST do this for background to show up.
    }

    public Component getTableCellRendererComponent(
                       JTable table, Object color, 
                       boolean isSelected, boolean hasFocus,
                       int row, int column) 
    {
      setBackground((Color)color);
      if(isBordered) 
      {
        if(isSelected) 
        {
          if(selectedBorder == null) 
            selectedBorder = BorderFactory.createMatteBorder(2,5,2,5,
                                     table.getSelectionBackground());
          setBorder(selectedBorder);
        }
        else 
        {
          if(unselectedBorder == null) 
            unselectedBorder = BorderFactory.createMatteBorder(2,5,2,5,
                                                table.getBackground());
          setBorder(unselectedBorder);
        }
      }
      return this;
    }
  }


  private void setUpColorRenderer(JTable table)
  {
    table.setDefaultRenderer(Color.class,
                             new ColorRenderer(true));
  }


  private void setUpColorEditor(JTable table)
  {
    //First, set up the button that brings up the dialog.
    final JButton button = new JButton("") 
    {
      public void setText(String s) {
           //Button never shows text -- only color.
      }
    };
    button.setBackground(Color.white);
    button.setBorderPainted(false);
    button.setMargin(new Insets(0,0,0,0));

    //Now create an editor to encapsulate the button, and
    //set it up as the editor for all Color cells.
    final ColorEditor colorEditor = new ColorEditor(button);
    table.setDefaultEditor(Color.class, colorEditor);

   //Set up the dialog that the button brings up.
    final JColorChooser colorChooser = new JColorChooser();
    ActionListener okListener = new ActionListener() 
    {
      public void actionPerformed(ActionEvent e) 
      {
        colorEditor.currentColor = colorChooser.getColor();
      }
    };
    final JDialog dialog = JColorChooser.createDialog(button,
                                        "Pick a Color",
                                        true,
                                        colorChooser,
                                        okListener,
                                        null); //XXXDoublecheck this is OK

    //Here's the code that brings up the dialog.
    button.addActionListener(new ActionListener() 
    {
      public void actionPerformed(ActionEvent e) 
      {
        button.setBackground(colorEditor.currentColor);
        colorChooser.setColor(colorEditor.currentColor);
        dialog.show();
      }
    });
  }


  private JButton setUpColorButton(Color col)
  {
    
    //First, set up the button that brings up the dialog.
    final JButton button = new JButton("")
    {
      public void setText(String s) {
           //Button never shows text -- only color.
      }
    };
    button.setBackground(col);
    button.setBorderPainted(false);
    button.setMargin(new Insets(0,0,0,0));
    Dimension d = new Dimension(25,25);
    button.setPreferredSize(d);
    button.setMinimumSize(d);
    button.setMaximumSize(d);

   //Set up the dialog that the button brings up.
    final JColorChooser colorChooser = new JColorChooser();
    ActionListener okListener = new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
        button.setBackground(colorChooser.getColor());
      }
    };
    final JDialog dialog = JColorChooser.createDialog(button,
                                        "Pick a Color",
                                        true,
                                        colorChooser,
                                        okListener,
                                        null); //XXXDoublecheck this is OK

    //Here's the code that brings up the dialog.
    button.addActionListener(new ActionListener()
    {
      public void actionPerformed(ActionEvent e)
      {
//      button.setBackground(colorChooser.getColor());
        dialog.show();
      }
    });
    return button;
  }

  /*
  * The editor button that brings up the dialog.
  * We extend DefaultCellEditor for convenience,
  * even though it mean we have to create a dummy
  * check box.  Another approach would be to copy
  * the implementation of TableCellEditor methods
  * from the source code for DefaultCellEditor.
  */
  class ColorEditor extends DefaultCellEditor 
  {
    Color currentColor = null;

    public ColorEditor(JButton b) 
    {
      super(new JCheckBox()); //Unfortunately, the constructor
                              //expects a check box, combo box,
                              //or text field.
      editorComponent = b;
      setClickCountToStart(1); //This is usually 1 or 2.

      //Must do this so that editing stops when appropriate.
      b.addActionListener(new ActionListener() 
      {
        public void actionPerformed(ActionEvent e) 
        {
          fireEditingStopped();
        }
      });
     }

     protected void fireEditingStopped() 
     {
       super.fireEditingStopped();
     }

     public Object getCellEditorValue() 
     {
       return currentColor;
     }

     public Component getTableCellEditorComponent(JTable table, 
                                                  Object value,
                                                  boolean isSelected,
                                                  int row,
                                                  int column) 
     {
       ((JButton)editorComponent).setText(value.toString());
       currentColor = (Color)value;
       return editorComponent;
     }
  }

  class RETableModel extends AbstractTableModel 
  {

    private Vector rowData;
    private Vector columnNames;

    public RETableModel(Vector rowData,Vector columnNames)
    {
      super();
      this.columnNames = columnNames;
      this.rowData     = rowData;
    }

    public int getColumnCount() 
    {
      return columnNames.size();
    }
        
    public int getRowCount() 
    {
      return rowData.size();
    }

    public String getColumnName(int col) 
    {
      return (String)columnNames.elementAt(col);
    }

    public Object getValueAt(int row, int col) 
    {
      return ((Vector)rowData.elementAt(row)).elementAt(col);
    }

    /*
    * JTable uses this method to determine the default renderer/
    * editor for each cell.  If we didn't implement this method,
    * then the last column would contain text ("true"/"false"),
    * rather than a check box.
    */
    public Class getColumnClass(int c) 
    {
      return getValueAt(0,c).getClass();
    }

    /*
    * Don't need to implement this method unless your table's
    * editable.
    */
    public boolean isCellEditable(int row, int col) 
    {
      return true;
    }

    /*
    * Don't need to implement this method unless your table's
    * data can change.
    */
    public void setValueAt(Object value, int row, int col) 
    {
      Vector vrow = (Vector)rowData.elementAt(row);
      vrow.setElementAt(value,col);
      rowData.setElementAt(vrow,row);
      fireTableCellUpdated(row, col);
    }


    public void addRow()
    {
      fireTableRowsInserted(rowData.size(),rowData.size());
    }

    /**
    *
    * Delete a row from the table
    * @param row  number to delete
    * @return     true if deleted
    *
    */
    public boolean deleteRow(int row)
    {
      if(row < 0 || row>=rowData.size())
        return false;
      rowData.remove(row);
      fireTableRowsDeleted(row,row); 

      return true;
    }


  }

}

