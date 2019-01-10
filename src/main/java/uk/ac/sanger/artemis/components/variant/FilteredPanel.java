/*
 * created: 2011
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2011  Genome Research Limited
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
import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.List;
import java.util.Vector;

import javax.swing.BorderFactory;
import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.border.Border;

  class FilteredPanel extends JPanel
  {
    private static final long serialVersionUID = 1L;
    private static Hashtable<String, RecordFilter> filters = new Hashtable<String, RecordFilter>();
    private static List<HeaderLine> hdrFltLines;
    private static List<String> hdrFltLinesIDs;
    
    private Box hdrFilterBox = Box.createVerticalBox();
    private Box userFilterBox = Box.createVerticalBox();
    
    FilteredPanel(final List<HeaderLine> hdrFltLines)
    {
      FilteredPanel.hdrFltLines = hdrFltLines;
      setLayout(new BorderLayout());
      Border raisedbevel = BorderFactory.createRaisedBevelBorder();
      Border loweredbevel = BorderFactory.createLoweredBevelBorder();    
      
      setBorder(BorderFactory.createCompoundBorder(
          raisedbevel, loweredbevel));
      setBackground(Color.white);

      Box titleBox = Box.createHorizontalBox();
      JLabel filterLabel = new JLabel("Filter Overview:");
      filterLabel.setFont(filterLabel.getFont().deriveFont(Font.BOLD));
      titleBox.add(filterLabel);
      titleBox.add(Box.createHorizontalGlue());
      add(titleBox, BorderLayout.NORTH);
      
      Box ftYBox = Box.createVerticalBox();
      add(ftYBox, BorderLayout.CENTER);
      ftYBox.add(userFilterBox, BorderLayout.CENTER);
      ftYBox.add(hdrFilterBox, BorderLayout.SOUTH);
      ftYBox.add(Box.createVerticalGlue());
      
      addHeaderFilters();
    }

    private void addHeaderFilters()
    {
      if(hdrFltLines.size() > 0)
      {
        JLabel lab = new JLabel("");
        Dimension d = new Dimension(150, lab.getPreferredSize().height);
        
        // SHOW CURRENT FILTERS
        for (int i = 0; i < hdrFltLines.size(); i++)
        {
          final Box hdrFilterLine = Box.createHorizontalBox();
          final HeaderLine line = hdrFltLines.get(i);
          lab = new JLabel(line.getID());
          lab.setPreferredSize(d);
          hdrFilterLine.add(lab);

          JLabel des = new JLabel(line.getDescription());
          hdrFilterLine.add(des);

          hdrFilterLine.add(Box.createHorizontalStrut(10));
          JButton remove = new JButton("X");
          remove.setBackground(Color.red);
          remove.addActionListener(new ActionListener()
          {
            public void actionPerformed(ActionEvent arg0)
            {
              hdrFilterBox.remove(hdrFilterLine);
              hdrFltLines.remove(line);
              hdrFltLinesIDs = null;
              repaint();
            }
          });
          hdrFilterLine.add(remove);
          
          hdrFilterLine.add(Box.createHorizontalGlue());
          hdrFilterBox.add(hdrFilterLine);
        }
        
        hdrFilterBox.add(Box.createVerticalGlue());
      }
    }

    protected void updateFilters()
    {
      userFilterBox.removeAll();
      final Enumeration<String> enumFilter = filters.keys();
      JLabel lab = new JLabel("");
      Dimension d = new Dimension(150, lab.getPreferredSize().height);
      while(enumFilter.hasMoreElements())
      {
        String id = enumFilter.nextElement();
        RecordFilter filter = filters.get(id);

        Box filterLine = Box.createHorizontalBox();
        lab = new JLabel(id+" ");
        lab.setPreferredSize(d);
        filterLine.add(lab);

        filterLine.add(new JLabel(filter.toString()));
        filterLine.add(Box.createHorizontalGlue());
        userFilterBox.add(filterLine);
      }
      revalidate();
      repaint();
    }
    
    protected static List<HeaderLine> getHeaderLineFilters()
    {
      return hdrFltLines;
    }
    
    protected static List<String> getHeaderLineFiltersIDs()
    {
      if(hdrFltLinesIDs == null)
      {
        hdrFltLinesIDs = new Vector<String>();
        for(HeaderLine ln: hdrFltLines)
          hdrFltLinesIDs.add(ln.getID());
      }
      return hdrFltLinesIDs;
    }
    
    protected static String getHeader()
    {
      StringBuffer buff = new StringBuffer();
      for(HeaderLine ln: hdrFltLines)
        buff.append(ln.toString()+"\n");

      Enumeration<String> filterStr = filters.keys();
      while(filterStr.hasMoreElements())
      {
        RecordFilter recFilter = filters.get(filterStr.nextElement());
        buff.append("##FILTER=<ID=");
        if(recFilter.getHeaderLine().getHeaderType() == HeaderLine.FORMAT_LINE)
          buff.append("sample");
        buff.append(recFilter.getHeaderLine().getID());
        buff.append(",Description=\"");
        buff.append(recFilter.toString());
        buff.append("\">\n");
      }
      return buff.toString();
    }
    
    protected void addFilter(final String ID, final HeaderLine hLine, final int number)
    {
      filters.put(ID, new RecordFilter(hLine, number));
    }
    
    protected void removeFilter(final String ID)
    {
      filters.remove(ID);
    }
    
    protected static Hashtable<String, RecordFilter> getFilters()
    {
      return filters;
    }
  }