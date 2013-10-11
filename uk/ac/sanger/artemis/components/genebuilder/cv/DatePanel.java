/* DatePanel.java
 *
 * created: 2007
 *
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
 **/

package uk.ac.sanger.artemis.components.genebuilder.cv;

import java.awt.Color;
import java.awt.Dimension;
import java.util.Calendar;
import java.util.Date;

import javax.swing.JSpinner;
import javax.swing.SpinnerDateModel;
import javax.swing.SpinnerListModel;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;


/**
 * Date input panel
 */
public class DatePanel extends JSpinner
{
  private static final long serialVersionUID = 1L;

  public DatePanel(final String date, 
                   final int height)
  {
    this(date);
    setMaximumSize(
        new Dimension(getPreferredSize().width, height));
  }
  
  public DatePanel(final String date)
  {
    Date initDate = getDate(date);
    Calendar calendar = Calendar.getInstance();
    calendar.add(Calendar.YEAR, -100);
    Date earliestDate = calendar.getTime();
    calendar.add(Calendar.YEAR, 200);
    Date latestDate = calendar.getTime();
    
    if(initDate != null)
    {
      SpinnerDateModel model = new SpinnerDateModel(initDate,
                                 earliestDate,
                                 latestDate,
                                 Calendar.YEAR);
      setModel(model);
      setBackground(Color.white);
      setEditor(new JSpinner.DateEditor(this, "yyyy/MM/dd"));
    }
    else
    {
      SpinnerListModel model = new SpinnerListModel( new String[] { "", "----/--/--", "" } );
      setModel(model);
      setBackground(Color.white);
      setValue("----/--/--");
      model.addChangeListener(new ChangeListener()
      {
        public void stateChanged(ChangeEvent e)
        {
          if(getModel() instanceof SpinnerListModel)
          {
            setModel(new SpinnerDateModel());
            setEditor(new JSpinner.DateEditor(DatePanel.this, "yyyy/MM/dd"));
          }
        }  
      });
    }
  }
  
  /**
   * Get a Date object from a date string in the format 20061129
   * @param dateStr
   * @return
   */
  protected Date getDate(String dateStr)
  {
    Calendar cal = Calendar.getInstance();
    cal.clear();
    if(dateStr == null ||
       dateStr.equals("") ||
       dateStr.length() != 8)
      return null;
    
    int year  = Integer.parseInt(dateStr.substring(0,4));
    int month = Integer.parseInt(dateStr.substring(4, 6))-1;
    int day   = Integer.parseInt(dateStr.substring(6,8));
    
    cal.set(year,month,day);
    return cal.getTime();
  }
  
  protected String getText()
  { 
    if(getModel() instanceof SpinnerDateModel)
    {
      java.text.SimpleDateFormat sdf = new java.text.SimpleDateFormat(
      "yyyyMMdd");
      Date date = ((SpinnerDateModel)getModel()).getDate();
      return sdf.format(date);
    }
    return "";
  }
  
  public static String getDate()
  {
    java.text.SimpleDateFormat sdf = new java.text.SimpleDateFormat(
        "yyyyMMdd");
    Date date = new Date();
    return sdf.format(date);
  }
}
