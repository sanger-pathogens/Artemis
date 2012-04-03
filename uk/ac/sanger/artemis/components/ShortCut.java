/* ShortCut.java
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
package uk.ac.sanger.artemis.components;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;

import javax.swing.JMenuItem;
import javax.swing.KeyStroke;

import uk.ac.sanger.artemis.Options;
import uk.ac.sanger.artemis.components.database.DatabaseJPanel;

/**
 * Used to cache shortcut keystroke between sessions
 */
class ShortCut implements Serializable
{
  private static final long serialVersionUID = 1L;
  private String shortCutName;
  private String menuName;
  private KeyStroke ks;
  private static org.apache.log4j.Logger logger4j = 
      org.apache.log4j.Logger.getLogger(DatabaseJPanel.class);
  private static String SHORTCUT_CACHE_PATH = Options.CACHE_PATH + File.separatorChar + "MenuShortCuts";
  private static String REGEXP_REPLACE = "[ /]";
    
  ShortCut(final String menuName, final String shortCutName, final KeyStroke ks)
  {
    this.menuName = menuName;
    this.shortCutName = shortCutName;
    this.ks = ks;
    writeCache();
  }
    
  private void writeCache()
  {
    try
    {
      File dir = new File(SHORTCUT_CACHE_PATH);
      if(!dir.exists())
        dir.mkdirs();
      FileOutputStream fos = new FileOutputStream(SHORTCUT_CACHE_PATH +
          File.separatorChar + menuName + "#" + shortCutName.replaceAll(REGEXP_REPLACE, "_"));
      ObjectOutputStream out = new ObjectOutputStream(fos);
      out.writeObject(ShortCut.this);
      out.close();
    }
    catch(Exception ex)
    {
      ex.getMessage();
    }
  }
  
  /**
   * Return true if shortcut_cache flag is set.
   * @return
   */
  protected static boolean usingCache()
  {
    return Options.getOptions().getPropertyTruthValue("shortcut_cache");
  }
  
  private static File[] getCachedShortCuts()
  {
    File cacheDir = new File(SHORTCUT_CACHE_PATH);
    if(cacheDir.exists())
      return cacheDir.listFiles();
    return null;
  }

  protected static void applyShortCutFromCache(final String menuName, final JMenuItem menuItem)
  {
    File[] shorts = getCachedShortCuts();
    if(shorts != null)
    {
      String shortCutName = menuItem.getText();
      for(int i=0; i<shorts.length; i++)
      {
        if(shorts[i].getName().equals(menuName + "#" + shortCutName.replaceAll(REGEXP_REPLACE, "_")) ||
           shorts[i].getName().equals(menuName + ":" + shortCutName.replaceAll(REGEXP_REPLACE, "_")))
        {
          try
          {
            FileInputStream fis = new FileInputStream(shorts[i]);
            ObjectInputStream in = new ObjectInputStream(fis);
            ShortCut sc = (ShortCut) in.readObject();
            in.close();
            menuItem.setAccelerator(sc.ks);
            
            logger4j.debug("USING CACHED SHORTCUT ("+menuItem.getText()+
                ") FROM: "+shorts[i].getAbsolutePath());
          }
          catch (FileNotFoundException e)
          {
            System.err.println(e.getMessage());
          }
          catch (IOException e)
          {
            System.err.println(e.getMessage());
          }
          catch (ClassNotFoundException e)
          {
            System.err.println(e.getMessage());
          }
        }
      }
    }
  }
}