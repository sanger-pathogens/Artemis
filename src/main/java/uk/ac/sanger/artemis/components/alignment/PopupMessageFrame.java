/* PopupMessageFrame
 *
 * created: 2009
 *
 * This file is part of Artemis
 *
 * Copyright(C) 2009  Genome Research Limited
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
package uk.ac.sanger.artemis.components.alignment;

import java.awt.Dimension;
import java.awt.Point;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JTextArea;
import javax.swing.SwingWorker;

/**
 * Create undecorated popup frames.
 */
class PopupMessageFrame extends JFrame
{
  private static final long serialVersionUID = 1L;
  private JTextArea textArea = new JTextArea();
  private double thisHeight;
  
  PopupMessageFrame()
  {
    super();
    getRootPane().putClientProperty("Window.alpha", new Float(0.8f));
    textArea.setWrapStyleWord(true);
    textArea.setEditable(false);
    setUndecorated(true);
    getContentPane().add(textArea);
  }

  PopupMessageFrame(String msg)
  {
    this();
    textArea.setFont(textArea.getFont().deriveFont(14.f));
    textArea.setText(msg);
    setTitle(msg);
    pack();
    thisHeight = getSize().getHeight();
  }
  
  protected void showWaiting(final String msg, final JPanel mainPanel)
  {
    PopupThread hide = new PopupThread(true, msg, mainPanel);
    hide.execute();
  }

  protected void hideFrame()
  {
    PopupThread hide = new PopupThread(false, null, null);
    hide.execute();
  }

  protected void show(final String msg, final JPanel mainPanel,
      final long aliveTime)
  {
    PopupThread p = new PopupThread(msg, mainPanel, aliveTime) ;
    p.execute();
  }
  
  class PopupThread extends SwingWorker<String, Object> 
  {
    private boolean show;
    private String msg;
    private JPanel mainPanel;
    private long aliveTime = -1;
    
    PopupThread(boolean show, String msg, JPanel mainPanel) 
    {
      this.show = show;
      this.msg = msg;
      this.mainPanel = mainPanel;
    }
    
    PopupThread(String msg, JPanel mainPanel, long aliveTime) 
    {
      this(true, msg, mainPanel);
      this.aliveTime = aliveTime;
    }
    
    @Override
    protected String doInBackground() throws Exception
    {
      if(aliveTime > 0)
        showFrameForGivenTime();
      else if(show)
        showFrame();
      else
        hideFrame();
      return null;
    }
    
    private void hideFrame()
    {
      if(getSize().getHeight() < thisHeight)
        return;

      Dimension d = getSize();
      double hgt = thisHeight/100.d;
      try
      {
        for(int i=0; i<100; i++)
        {
          Thread.sleep(1);
          d.setSize(d.getWidth(), thisHeight-(hgt*i));
          setPreferredSize(d);
          pack();
        }
      }
      catch (InterruptedException e)
      {
        e.printStackTrace();
      }
      setVisible(false);
      
      d.setSize(d.getWidth(), thisHeight);
      setPreferredSize(d);
      pack();
    }
    
    private void showFrame()
    {
      textArea.setText(msg);
      pack();
      Point p = mainPanel.getLocationOnScreen();
      p.x += (mainPanel.getWidth() - PopupMessageFrame.this.getWidth()) / 2;
      p.y += mainPanel.getHeight() / 2;
      setLocation(p);
      setVisible(true);
    }
    
    private void showFrameForGivenTime()
    {
      textArea.setText(msg);
      pack();
      Point p = mainPanel.getLocationOnScreen();
      p.x += (mainPanel.getWidth() - PopupMessageFrame.this.getWidth()) / 2;
      p.y += mainPanel.getHeight() / 2;
      setLocation(p);
      
      setVisible(true);
      requestFocus();

      try
      {
        Thread.sleep(aliveTime);
      }
      catch (InterruptedException e)
      {
        e.printStackTrace();
      }
      setVisible(false);
    }
  }

}
