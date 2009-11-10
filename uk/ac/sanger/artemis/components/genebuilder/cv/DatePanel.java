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
public class DatePanel
{
  private JSpinner spinner;
  
  public DatePanel(final String date, 
                   final int height)
  {
    this(date);
    spinner.setMaximumSize(
        new Dimension(spinner.getPreferredSize().width, height));
  }
  
  public DatePanel(final String date)
  {
    Date initDate = getDate(date);
    Calendar calendar = Calendar.getInstance();
    
    //if(initDate == null)
    //  initDate = calendar.getTime();
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
      spinner = new JSpinner(model);
      spinner.setBackground(Color.white);
      spinner.setEditor(new JSpinner.DateEditor(spinner, "yyyy/MM/dd"));
    }
    else
    {
      SpinnerListModel model = new SpinnerListModel( new String[] { "", "----/--/--", "" } );
      spinner = new JSpinner(model);
      spinner.setBackground(Color.white);
      spinner.setValue("----/--/--");
      model.addChangeListener(new ChangeListener()
      {
        public void stateChanged(ChangeEvent e)
        {
          if(spinner.getModel() instanceof SpinnerListModel)
          {
            spinner.setModel(new SpinnerDateModel());
            spinner.setEditor(new JSpinner.DateEditor(spinner, "yyyy/MM/dd"));
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
  
  /*class DateVerifier extends InputVerifier 
  {
    java.text.SimpleDateFormat sdf = new java.text.SimpleDateFormat(
            "yyyyMMdd");
    java.util.Calendar cal = java.util.Calendar.getInstance();
    public DateVerifier() 
    {
      sdf.setLenient(false);
    }

    public boolean verify(JComponent input) 
    {
      JFormattedTextField ftf = (JFormattedTextField) input;
      // allow null entry which will include slashes because of the
      // mask
      if(ftf.getText().trim().equals(""))
      {
        ftf.setValue( null );
        return true;
      }
        
      try 
      {
        cal.setTime(sdf.parse(ftf.getText()));
      }
      catch (Exception pe) 
      {
        return false;
      }
      return true;
    }
  }*/

  public JSpinner getDateSpinner()
  {
    return spinner;
  }
  
  public String getText()
  { 
    if(spinner.getModel() instanceof SpinnerDateModel)
    {
      java.text.SimpleDateFormat sdf = new java.text.SimpleDateFormat(
      "yyyyMMdd");
      Date date = ((SpinnerDateModel)spinner.getModel()).getDate();
      return sdf.format(date);
    }
    return "";
  }
  
  protected static String getDate()
  {
    java.text.SimpleDateFormat sdf = new java.text.SimpleDateFormat(
        "yyyyMMdd");
    Date date = new Date();
    return sdf.format(date);
  }
  
}