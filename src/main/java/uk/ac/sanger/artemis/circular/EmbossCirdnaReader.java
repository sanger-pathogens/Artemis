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

import javax.swing.*;
import java.awt.Color;
import java.io.*;
import java.util.Vector;
import java.util.StringTokenizer;

/**
*
*
*/
public class EmbossCirdnaReader
{

  /** cirdna file */
  private File cirdnaFile;
  /** true if read ok */
  private boolean reading = false;
  /** */
  private Vector restrictionEnzyme = new Vector();
  /** genetic markers */
  private Vector block = new Vector();
  /** emboss colour scheme for cirdna */
  private Color[] embossColor = { 
                      Color.black,
                      Color.red,
                      Color.yellow,
                      Color.green,
                      Color.decode("#99CCFF"),   //AQUAMARINE
                      Color.decode("#FFCCCC"),   //pink
                      Color.decode("#FFFFCC"),   //wheat
                      Color.gray,
                      Color.decode("#993300"),   //brown
                      Color.blue,
                      Color.decode("#9933FF"),   //blueviolet
                      Color.cyan,
                      Color.decode("#33FFCC"),   //turqoise
                      Color.decode("#FF00FF"),   //magenta
                      Color.decode("#FF9966"),   //salmon
                      Color.white
                                 };
  int start = 0;
  int end = 0;


  public EmbossCirdnaReader()
  {
    SecurityManager sm = System.getSecurityManager();
    System.setSecurityManager(null);
    JFileChooser fc = new JFileChooser(System.getProperty("user.home"));
    System.setSecurityManager(sm);

    int returnVal = fc.showOpenDialog(fc);
  
    if(returnVal == JFileChooser.APPROVE_OPTION)
    {
      cirdnaFile = fc.getSelectedFile();
      readFile();
      reading = true;
    }
  } 

  /**
  *
  * @param cirdnaFile	cirdna file
  *
  */
  public EmbossCirdnaReader(File cirdnaFile)
  {
    this.cirdnaFile = cirdnaFile;
    readFile();
    reading = true;
  }


  /**
  *
  *
  */
  public boolean isReading()
  {
    return reading;
  }

  /**
  *
  * Read a cirdna file 
  *
  */
  public Vector readFile()
  {
    BufferedReader in = null;
    try
    {
      in = new BufferedReader(new FileReader(cirdnaFile));
      String line;
      while((line = in.readLine()) != null )
      {
        line = line.trim().toLowerCase();
        if(!line.equals(""))
        {
          StringTokenizer stok = new StringTokenizer(line," ");
          if(line.startsWith("start "))
          {
            stok.nextElement();
            start = Integer.parseInt((String)stok.nextElement());
          }
          else if(line.startsWith("end "))
          {  
            stok.nextElement();
            end = Integer.parseInt((String)stok.nextElement());
          }
          else if(line.startsWith("group"))
          {
            while((line = in.readLine()) != null )
            {
              line = line.trim().toLowerCase();
              if(line.startsWith("endgroup"))
                break;
              else if(line.startsWith("block "))
              {
                stok = new StringTokenizer(line," ");
                Vector marker = new Vector();
                stok.nextElement();
                Integer bstart = new Integer((String)stok.nextElement());
                Integer bend   = new Integer((String)stok.nextElement());
                Color col = Color.red;

                if(stok.hasMoreTokens())
                  col = embossColor[Integer.parseInt((String)stok.nextElement())];
                String name = in.readLine().trim();
                if(name.equals("endlabel"))
                  name = "";
                 
                marker.add(name);
                marker.add(bstart);
                marker.add(bend);
                marker.add(col);
                marker.add(new Float(10.f));
                marker.add(new Boolean(false));
                marker.add(new Boolean(false));
                block.add(marker);
              }
              else if(line.startsWith("tick"))
              {
                stok = new StringTokenizer(line," ");
 
                stok.nextElement();
                Integer pos = new Integer((String)stok.nextElement());
                Color col = Color.red;
                if(stok.hasMoreTokens())
                  col = embossColor[Integer.parseInt((String)stok.nextElement())];
                String name = in.readLine();     
                if(line.equals("endlabel"))
                  name = "";

                Vector re = new Vector();
                re.add(name);
                re.add(pos);
                re.add(col);
                restrictionEnzyme.add(re);
              }
            }
          }
          

        }
      }
    }
    catch (IOException e)
    {
      System.out.println("SequenceReader Error");
    }

    System.out.println("Start : "+start);
    System.out.println("End   : "+end);

    return null;
  }

 
  protected Vector getRestrictionEnzyme()
  { 
    return restrictionEnzyme;
  }
  

  protected Vector getBlock()
  {
    return block;
  }

  protected int getStart()
  {
    return start;
  }

  protected int getEnd()
  {
    return end;
  }

}

